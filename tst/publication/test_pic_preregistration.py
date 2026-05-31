#!/usr/bin/env python3
"""Regression tests for PIC Q-022/Q-023 preregistration sidecars."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import unittest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from tst.publication.pic_qualification_manifest import validate_schema  # noqa: E402


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
Q022_PATH = READINESS_DIR / "q022_independent_comparison_preregistration_2026-05-30.json"
Q023_PATH = READINESS_DIR / "q023_statistical_qualification_preregistration_2026-05-30.json"
Q018_LINKS_PATH = READINESS_DIR / "q018_claim_links_successor_2026-05-30.json"
Q022_ENTITY_MAP_PATH = (
    READINESS_DIR / "q022_xcmp_entity_micro_equation_map_2026-05-30.json"
)
Q022_DATASET_PROVENANCE_PATH = (
    READINESS_DIR / "q022_dataset_provenance_manifest_2026-05-30.json"
)
Q022_PRIVATE_INGEST_PATH = (
    READINESS_DIR / "q022_external_reference_private_ingest_2026-05-30.json"
)
Q039_PATH = READINESS_DIR / "q039_gizmo_rsol_decision_2026-05-30.json"
Q039_EXCLUSION_CANDIDATE_PATH = (
    READINESS_DIR / "q039_gizmo_documented_exclusion_candidate_2026-05-30.json"
)
Q022_DATASET_PROVENANCE_SCHEMA_PATH = (
    READINESS_DIR / "schemas" / "q022_dataset_provenance_manifest.schema.json"
)
Q022_EQUATION_MAP_SCHEMA_PATH = (
    READINESS_DIR / "schemas" / "q022_equation_normalization_map.schema.json"
)
Q022_TOLERANCE_TABLE_SCHEMA_PATH = (
    READINESS_DIR / "schemas" / "q022_tolerance_table.schema.json"
)
Q023_DRAFTS_PATH = READINESS_DIR / "q023_campaign_drafts_2026-05-30.json"
Q023_PROFILES_PATH = (
    READINESS_DIR / "q023_local_preregistration_profiles_2026-05-30.json"
)
Q023_SCHEMA_PATH = (
    READINESS_DIR / "schemas" / "q023_campaign_preregistration.schema.json"
)

Q022_REFERENCE_IDS = {
    "bai_2015_arxiv_1412.1087",
    "riquelme_spitkovsky_2009_arxiv_0810.4565",
    "gargate_2010_arxiv_1002.1701",
    "zacharegkas_2022_arxiv_2210.08072",
    "mignone_2018_arxiv_1804.01946",
    "van_marle_casse_marcowith_2018_doi_stx2509",
    "bai_2019_arxiv_1902.10219",
    "plotnikov_ostriker_bai_2021_arxiv_2102.11878",
    "sun_bai_zhao_2024_arxiv_2409.08592",
    "ji_hopkins_2022_arxiv_2111.14704",
    "entity_toolkit_snapshot_512998c4",
}

Q022_COMPARISON_IDS = {
    "XCMP-ATHENA-BAI-BELL-SHOCK",
    "XCMP-BELL-NONLINEAR-INDEPENDENT",
    "XCMP-BAI2019-CRSI",
    "XCMP-PLUTO-COUPLING",
    "XCMP-AMRVAC-SHOCK",
    "XCMP-GIZMO-RSOL-DECISION",
    "XCMP-ENTITY-MICRO",
    "XCMP-EXT-CRSI-IN-DAMPING",
    "XCMP-EXT-HALL-BELL",
    "XCMP-EXT-CRPAI-TRANSPORT",
}

Q022_AUTHORIZED_EXTRACTED_DATASET_ROOT = (
    "/lustre/orion/ast207/proj-shared/dfielding/PIC/reference_artifacts/"
    "q022_extracted_datasets_2026-05-30"
)


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _test_only_q022_ready_path_fixtures() -> tuple[
    dict[str, object], dict[str, object], dict[str, object]
]:
    """Return synthetic complete records for schema ready-path coverage only."""
    comparison_id = "XCMP-TEST-ONLY-READY-PATH"
    dataset_id = "Q022-DATASET-TEST-ONLY-READY-PATH"
    provenance = copy.deepcopy(_load(Q022_DATASET_PROVENANCE_PATH))
    provenance["authorized_extracted_dataset_root_status"] = (
        "present_extracted_dataset_ingest_requires_manifest_update"
    )
    provenance["dataset_candidates"] = [
        {
            "dataset_id": dataset_id,
            "comparison_id": comparison_id,
            "reference_ids": ["test_only_reference"],
            "source_records": ["test-only synthetic schema fixture"],
            "extraction_status": (
                "extracted_dataset_checksum_verified_pending_external_review"
            ),
            "redistribution_basis": "test-only synthetic redistribution basis",
            "reviewer_disposition": "pending external review",
            "blocked_input_ids": [],
            "extracted_dataset_locator": (
                f"{Q022_AUTHORIZED_EXTRACTED_DATASET_ROOT}/"
                "test-only-ready-path/dataset.csv"
            ),
            "extracted_dataset_sha256": "a" * 64,
            "extraction_method": "test-only synthetic extraction method",
            "extraction_script_sha256": "b" * 64,
            "extraction_uncertainty": "test-only synthetic uncertainty record",
        }
    ]
    provenance["manifest_status"] = "test_only_ready_path_schema_fixture"

    equation_map = {
        "schema_version": 2,
        "map_status": "frozen_before_comparison_run",
        "comparison_id": comparison_id,
        "dataset_provenance_id": dataset_id,
        "reference_ids": ["test_only_reference"],
        "athenak_physical_mode": "test_only_ready_path_mode",
        "matched_equations": ["test-only synthetic matched equation"],
        "intentional_mismatches": ["test-only synthetic excluded mismatch"],
        "unit_map": {"observable": "test-only synthetic unit map"},
        "normalization_map": {"observable": "test-only synthetic normalization"},
        "parameter_overlap": {"parameter": "test-only synthetic overlap"},
        "excluded_regimes": ["test-only synthetic excluded regime"],
        "source_records": ["test-only synthetic schema fixture"],
        "blocked_inputs": [],
        "reviewer_disposition": "pending external review",
    }
    tolerance_table = {
        "schema_version": 2,
        "comparison_id": comparison_id,
        "dataset_provenance_id": dataset_id,
        "freeze_status": "frozen_before_qualifying_run",
        "observable_families": ["test_only_observable_family"],
        "rows": [
            {
                "observable": "test_only_observable",
                "units": "test_only_units",
                "reference_value_or_extraction": (
                    "test-only synthetic extracted value with uncertainty"
                ),
                "expected_discretization_trend": (
                    "test-only synthetic convergent trend"
                ),
                "absolute_tolerance": 0.01,
                "relative_tolerance": 0.02,
                "tolerance_rationale": "test-only synthetic tolerance rationale",
                "cpu_mpi_gpu_reduction_variation": (
                    "test-only synthetic reduction variation bound"
                ),
                "failure_artifact_path": "test-only/failures/observable.json",
            }
        ],
        "blocked_inputs": [],
        "reviewer_disposition": "pending external review",
    }
    return provenance, equation_map, tolerance_table


class PicPreregistrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.q022 = _load(Q022_PATH)
        self.q023 = _load(Q023_PATH)

    def test_sidecars_remain_open_scaffolds(self) -> None:
        for doc, gate in ((self.q022, "Q-022"), (self.q023, "Q-023")):
            with self.subTest(gate=gate):
                self.assertEqual(doc["gate"], gate)
                self.assertEqual(doc["qualification_effect"],
                                 "preregistration_scaffold_only")
                self.assertEqual(doc["disposition"], "open")
                self.assertTrue(doc["open_items"])
                self.assertIn("open", doc[f"q{gate[-3:]}_status"])

    def test_q022_freezes_reference_and_comparison_inventory(self) -> None:
        references = self.q022["reference_families"]
        reference_ids = {item["reference_id"] for item in references}
        self.assertEqual(reference_ids, Q022_REFERENCE_IDS)
        self.assertEqual(len(reference_ids), len(references))
        comparisons = self.q022["comparison_matrix"]
        comparison_ids = {item["comparison_id"] for item in comparisons}
        self.assertEqual(comparison_ids, Q022_COMPARISON_IDS)
        self.assertEqual(len(comparison_ids), len(comparisons))
        for comparison in comparisons:
            with self.subTest(comparison=comparison["comparison_id"]):
                self.assertTrue(comparison["claim_ids"])
                self.assertTrue(comparison["reference_ids"])
                self.assertTrue(set(comparison["reference_ids"]) <= reference_ids)
                self.assertTrue(comparison["observable_families"])
                self.assertTrue(comparison["dataset_provenance_id"])
                self.assertTrue(comparison["equation_normalization_map"])
                self.assertTrue(comparison["tolerance_table"])
                if comparison["comparison_id"] == "XCMP-ENTITY-MICRO":
                    self.assertEqual(
                        comparison["equation_map_status"],
                        "bounded_map_recorded_external_review_open",
                    )
                    self.assertEqual(
                        comparison["normalization_map_status"],
                        "bounded_map_recorded_external_review_open",
                    )
                    self.assertEqual(
                        comparison["equation_normalization_map"],
                        "tst/publication/readiness/"
                        "q022_xcmp_entity_micro_equation_map_2026-05-30.json",
                    )
                    self.assertEqual(
                        comparison["tolerance_table_status"],
                        "blocked_pending_external_review",
                    )
                elif comparison["comparison_id"] == "XCMP-GIZMO-RSOL-DECISION":
                    self.assertEqual(
                        comparison["equation_map_status"],
                        "documented_exclusion_candidate_pending_external_review",
                    )
                    self.assertEqual(
                        comparison["normalization_map_status"],
                        "documented_exclusion_candidate_pending_external_review",
                    )
                    self.assertEqual(
                        comparison["tolerance_table_status"],
                        "not_applicable_documented_exclusion_candidate",
                    )
                else:
                    self.assertEqual(
                        comparison["equation_map_status"],
                        "blocked_pending_reference_specific_mapping_and_"
                        "external_review",
                    )
                    self.assertEqual(
                        comparison["normalization_map_status"],
                        "blocked_pending_reference_specific_mapping_and_"
                        "external_review",
                    )
                    self.assertEqual(
                        comparison["tolerance_table_status"],
                        "blocked_pending_reference_dataset_extraction_and_"
                        "external_review",
                    )

    def test_q022_freezes_tolerance_and_discrepancy_rules(self) -> None:
        policy = self.q022["tolerance_policy"]
        self.assertTrue(policy["reference_extraction_uncertainty_required"])
        self.assertTrue(policy["no_relaxation_after_inspection"])
        self.assertTrue(policy["threshold_revision_invalidates_candidate_dataset"])
        self.assertIn("open", policy["numeric_thresholds_status"])
        required_tolerance_fields = {
            "observable",
            "units",
            "reference_value_or_extraction",
            "expected_discretization_trend",
            "absolute_tolerance",
            "relative_tolerance",
            "tolerance_rationale",
            "cpu_mpi_gpu_reduction_variation",
            "failure_artifact_path",
        }
        self.assertEqual(set(policy["required_fields_before_comparison_run"]),
                         required_tolerance_fields)
        ledger = self.q022["discrepancy_ledger_schema"]
        self.assertIn("unresolved", ledger["classification_enum"])
        self.assertIn("rejected_claim", ledger["disposition_enum"])
        self.assertIn("reviewer", ledger["required_fields"])

    def test_q023_freezes_seed_interval_and_exclusion_policy(self) -> None:
        seeds = self.q023["seed_policy"]
        self.assertEqual(seeds["default_design"], "fixed_sample")
        self.assertIn("locally_frozen", seeds["qualifying_seed_list_status"])
        self.assertTrue(seeds["archive_every_attempted_seed"])
        self.assertTrue(seeds["archive_failed_and_outlier_seeds"])
        self.assertTrue(seeds["informal_optional_seed_addition_forbidden"])
        self.assertEqual(seeds["sequential_design"]["default_status"], "disabled")

        interval = self.q023["interval_policy"]
        self.assertEqual(interval["confidence_level"], 0.95)
        self.assertEqual(interval["bootstrap_resamples"], 10000)
        self.assertEqual(interval["bootstrap_rng_seed"], 20260530)

        exclusions = self.q023["exclusion_policy"]
        self.assertTrue(exclusions["archive_excluded_attempts"])
        self.assertTrue(exclusions["post_inspection_window_or_threshold_tuning_forbidden"])
        self.assertTrue(exclusions["candidate_dataset_invalidated_by_policy_revision"])

    def test_q023_requires_independent_recompute_and_campaign_freeze(self) -> None:
        template = self.q023["campaign_preregistration_template"]
        required = set(template["required_fields_before_qualifying_run"])
        self.assertTrue({
            "qualifying_seed_list",
            "estimators",
            "fit_windows",
            "exclusion_rules",
            "parameter_grid",
            "tolerances",
            "independent_recompute_plan",
        } <= required)
        recompute = self.q023["independent_recompute_policy"]
        self.assertEqual(recompute["status"], "frozen_required_before_claim_closure")
        self.assertIn("open", recompute["current_artifact_status"])
        self.assertIn("script_sha256", recompute["required_records"])

    def test_q023_campaign_claims_exist_and_remain_open(self) -> None:
        registry = _load(READINESS_DIR / "claims_registry.json")
        claims = {claim["claim_id"]: claim for claim in registry["claims"]}
        for claim_id in self.q023["campaigns_requiring_campaign_specific_freeze"]:
            with self.subTest(claim_id=claim_id):
                self.assertIn(claim_id, claims)
                self.assertEqual(claims[claim_id]["disposition"], "open")

    def test_q018_claim_link_scaffold_covers_registry(self) -> None:
        registry = _load(READINESS_DIR / "claims_registry.json")
        links = _load(Q018_LINKS_PATH)
        self.assertEqual(links["disposition"], "open")
        self.assertEqual(links["reviewer_assignment"], "pending external review")
        registry_ids = {claim["claim_id"] for claim in registry["claims"]}
        link_ids = {claim["claim_id"] for claim in links["claim_links"]}
        self.assertEqual(link_ids, registry_ids)
        self.assertEqual(len(link_ids), len(links["claim_links"]))
        for claim in links["claim_links"]:
            with self.subTest(claim_id=claim["claim_id"]):
                self.assertTrue(claim["applicability_envelope"])
                self.assertTrue(claim["limitations"])
                self.assertIsInstance(claim["evidence_links"], list)
                for evidence_link in claim["evidence_links"]:
                    with self.subTest(evidence_link=evidence_link):
                        self.assertTrue(
                            (REPO_ROOT / evidence_link).is_file(),
                            f"missing claim evidence link: {evidence_link}",
                        )

    def test_q018_registered_f1_bundle_tracks_pending_review_manifests(self) -> None:
        links = _load(Q018_LINKS_PATH)
        bundle = links["registered_f1_pending_external_review_bundle"]
        successor = _load(
            READINESS_DIR
            / "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        f2 = _load(
            READINESS_DIR
            / "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )["accepted_v2_execution"]
        self.assertEqual(
            bundle["status"],
            "frozen_pending_external_review_not_claim_qualified",
        )
        self.assertEqual(
            bundle["historical_bounded_f2_projection_policy_sha256"],
            f2["active_policy_sha256"],
        )
        self.assertEqual(
            bundle["historical_bounded_f2_projection_promotion_sha256"],
            f2["active_promotion_sha256"],
        )
        strict = _load(
            READINESS_DIR
            / "q027_manual_frontier_accounting_activation_2026-05-30.json"
        )["control_plane_transition"]
        self.assertEqual(
            bundle["live_strict_operational_policy_sha256"],
            strict["active_policy_sha256"],
        )
        self.assertEqual(
            bundle["live_strict_operational_promotion_sha256"],
            strict["active_promotion_sha256"],
        )
        expected = {
            execution["registered_science_authorization_id"]: (
                execution["qualification_manifest_path"],
                execution["qualification_manifest_sha256"],
            )
            for execution in [
                successor["accepted_gyro_v3_registered_execution"],
                successor["accepted_paper_coupling_v2_registered_execution"],
            ]
        }
        actual = {
            qualification["registered_science_authorization_id"]: (
                qualification["manifest_path"],
                qualification["manifest_sha256"],
            )
            for qualification in bundle["qualifications"]
        }
        self.assertEqual(actual, expected)
        self.assertTrue(
            all(
                qualification["review"] == "pending external review"
                for qualification in bundle["qualifications"]
            )
        )

    def test_q022_entity_map_matches_schema_and_stays_bounded(self) -> None:
        schema = _load(Q022_EQUATION_MAP_SCHEMA_PATH)
        entity_map = _load(Q022_ENTITY_MAP_PATH)
        validate_schema(entity_map, schema)
        self.assertEqual(entity_map["map_status"],
                         "bounded_map_recorded_external_review_open")
        self.assertEqual(entity_map["comparison_id"], "XCMP-ENTITY-MICRO")
        self.assertEqual(entity_map["reviewer_disposition"],
                         "pending external review")
        self.assertIn("MHD feedback", entity_map["excluded_regimes"])
        self.assertTrue(any("not claimed interchangeable" in mismatch
                            for mismatch in entity_map["intentional_mismatches"]))

    def test_q022_dataset_provenance_manifest_stays_fail_closed(self) -> None:
        schema = _load(Q022_DATASET_PROVENANCE_SCHEMA_PATH)
        provenance = _load(Q022_DATASET_PROVENANCE_PATH)
        validate_schema(provenance, schema)
        self.assertEqual(
            provenance["authorized_extracted_dataset_root"],
            Q022_AUTHORIZED_EXTRACTED_DATASET_ROOT,
        )
        self.assertEqual(
            provenance["authorized_extracted_dataset_root_status"],
            "absent_no_bulk_extracted_dataset_ingested",
        )
        policy = provenance["storage_policy"]
        self.assertEqual(policy["source_control_scope"],
                         "metadata_only_no_bulk_extracted_data")
        self.assertEqual(policy["allowed_bulk_storage_system"], "Orion")
        self.assertEqual(policy["forbidden_storage_systems"], ["Kronos"])
        self.assertTrue(policy["no_fabricated_extracted_measurements"])

        blocked_input_ids = set(provenance["blocked_inputs"])
        comparisons = {
            item["comparison_id"]: item for item in self.q022["comparison_matrix"]
        }
        candidates = provenance["dataset_candidates"]
        self.assertEqual(
            {item["comparison_id"] for item in candidates},
            set(comparisons),
        )
        self.assertEqual(
            len({item["dataset_id"] for item in candidates}),
            len(candidates),
        )
        for candidate in candidates:
            with self.subTest(dataset=candidate["dataset_id"]):
                comparison = comparisons[candidate["comparison_id"]]
                self.assertEqual(candidate["dataset_id"],
                                 comparison["dataset_provenance_id"])
                self.assertEqual(candidate["reference_ids"],
                                 comparison["reference_ids"])
                self.assertEqual(candidate["reviewer_disposition"],
                                 "pending external review")
                self.assertTrue(candidate["blocked_input_ids"])
                self.assertTrue(set(candidate["blocked_input_ids"])
                                <= blocked_input_ids)
                self.assertNotIn("extracted_dataset_locator", candidate)
                self.assertNotIn("extracted_dataset_sha256", candidate)

    def test_q022_ingest_metadata_records_absent_authorized_extraction_root(
        self,
    ) -> None:
        ingest = _load(Q022_PRIVATE_INGEST_PATH)
        provenance = ingest["dataset_provenance"]
        self.assertEqual(
            provenance["source_controlled_manifest"],
            "tst/publication/readiness/"
            "q022_dataset_provenance_manifest_2026-05-30.json",
        )
        self.assertEqual(
            provenance["authorized_orion_extracted_dataset_root"],
            Q022_AUTHORIZED_EXTRACTED_DATASET_ROOT,
        )
        self.assertEqual(
            provenance["inspection_status"],
            "absent_no_bulk_extracted_dataset_ingested",
        )
        self.assertEqual(provenance["bulk_extracted_dataset_files_recorded"], 0)
        self.assertEqual(provenance["allowed_bulk_storage_system"], "Orion")
        self.assertEqual(provenance["forbidden_storage_systems"], ["Kronos"])

    def test_q022_comparison_sidecars_validate_without_fabricated_rows(
        self,
    ) -> None:
        equation_schema = _load(Q022_EQUATION_MAP_SCHEMA_PATH)
        tolerance_schema = _load(Q022_TOLERANCE_TABLE_SCHEMA_PATH)
        for comparison in self.q022["comparison_matrix"]:
            with self.subTest(comparison=comparison["comparison_id"]):
                equation_map = _load(REPO_ROOT
                                     / comparison["equation_normalization_map"])
                tolerance_table = _load(REPO_ROOT
                                        / comparison["tolerance_table"])
                validate_schema(equation_map, equation_schema)
                validate_schema(tolerance_table, tolerance_schema)
                self.assertEqual(equation_map["comparison_id"],
                                 comparison["comparison_id"])
                self.assertEqual(equation_map["dataset_provenance_id"],
                                 comparison["dataset_provenance_id"])
                self.assertEqual(equation_map["reference_ids"],
                                 comparison["reference_ids"])
                self.assertEqual(equation_map["map_status"],
                                 comparison["equation_map_status"])
                self.assertEqual(tolerance_table["comparison_id"],
                                 comparison["comparison_id"])
                self.assertEqual(tolerance_table["dataset_provenance_id"],
                                 comparison["dataset_provenance_id"])
                self.assertEqual(tolerance_table["freeze_status"],
                                 comparison["tolerance_table_status"])
                self.assertEqual(tolerance_table["observable_families"],
                                 comparison["observable_families"])
                self.assertEqual(tolerance_table["rows"], [])
                self.assertTrue(tolerance_table["blocked_inputs"])
                self.assertEqual(tolerance_table["reviewer_disposition"],
                                 "pending external review")

    def test_q039_gizmo_documented_exclusion_candidate_stays_provisional(
        self,
    ) -> None:
        q039 = _load(Q039_PATH)
        candidate = _load(Q039_EXCLUSION_CANDIDATE_PATH)
        comparison = next(
            item for item in self.q022["comparison_matrix"]
            if item["comparison_id"] == "XCMP-GIZMO-RSOL-DECISION"
        )
        self.assertEqual(
            q039["linked_records"]["documented_exclusion_candidate"],
            "tst/publication/readiness/"
            "q039_gizmo_documented_exclusion_candidate_2026-05-30.json",
        )
        self.assertEqual(candidate["candidate_route"], "documented_exclusion")
        self.assertEqual(candidate["disposition"], "open")
        self.assertEqual(candidate["reviewer_disposition"],
                         "pending external review")
        self.assertEqual(
            comparison["tolerance_table_status"],
            "not_applicable_documented_exclusion_candidate",
        )
        self.assertTrue(any("No extracted measurement" in item
                            for item in candidate["non_claims"]))

    def test_q022_schemas_reject_unfilled_controls_promoted_to_ready(
        self,
    ) -> None:
        provenance_schema = _load(Q022_DATASET_PROVENANCE_SCHEMA_PATH)
        provenance = copy.deepcopy(_load(Q022_DATASET_PROVENANCE_PATH))
        provenance["dataset_candidates"][0]["extraction_status"] = (
            "extracted_dataset_checksum_verified_pending_external_review"
        )
        with self.assertRaises(ValueError):
            validate_schema(provenance, provenance_schema)

        comparison = self.q022["comparison_matrix"][0]
        equation_map = _load(REPO_ROOT / comparison["equation_normalization_map"])
        equation_map["map_status"] = "frozen_before_comparison_run"
        with self.assertRaises(ValueError):
            validate_schema(equation_map, _load(Q022_EQUATION_MAP_SCHEMA_PATH))

        tolerance_table = _load(REPO_ROOT / comparison["tolerance_table"])
        tolerance_table["freeze_status"] = "frozen_before_qualifying_run"
        with self.assertRaises(ValueError):
            validate_schema(tolerance_table,
                            _load(Q022_TOLERANCE_TABLE_SCHEMA_PATH))

    def test_q022_schemas_accept_complete_test_only_ready_path_fixtures(
        self,
    ) -> None:
        provenance, equation_map, tolerance_table = (
            _test_only_q022_ready_path_fixtures()
        )
        validate_schema(provenance, _load(Q022_DATASET_PROVENANCE_SCHEMA_PATH))
        validate_schema(equation_map, _load(Q022_EQUATION_MAP_SCHEMA_PATH))
        validate_schema(tolerance_table, _load(Q022_TOLERANCE_TABLE_SCHEMA_PATH))

        candidate = provenance["dataset_candidates"][0]
        self.assertEqual(candidate["comparison_id"], equation_map["comparison_id"])
        self.assertEqual(candidate["comparison_id"],
                         tolerance_table["comparison_id"])
        self.assertEqual(candidate["dataset_id"],
                         equation_map["dataset_provenance_id"])
        self.assertEqual(candidate["dataset_id"],
                         tolerance_table["dataset_provenance_id"])
        self.assertEqual(len(tolerance_table["rows"]), 1)

    def test_q023_drafts_cover_policy_campaigns_and_block_runs(self) -> None:
        drafts = _load(Q023_DRAFTS_PATH)
        profiles = _load(Q023_PROFILES_PATH)
        schema = _load(Q023_SCHEMA_PATH)
        self.assertEqual(drafts["disposition"], "open")
        self.assertEqual(
            drafts["profile_registry"],
            "tst/publication/readiness/"
            "q023_local_preregistration_profiles_2026-05-30.json",
        )
        campaign_claims = {
            claim_id
            for campaign in drafts["campaigns"]
            for claim_id in campaign["claim_ids"]
        }
        self.assertEqual(
            campaign_claims,
            set(self.q023["campaigns_requiring_campaign_specific_freeze"]),
        )
        self.assertEqual(len(drafts["campaigns"]), len(campaign_claims))
        self.assertIn("blocked", drafts["rule"])

        profile_fields = {
            "seed_profile": "seed_profiles",
            "sampling_rule_profile": "sampling_rule_profiles",
            "parameter_grid_profile": "grid_profiles",
            "fit_window_profile": "analysis_window_profiles",
            "tolerance_profile": "tolerance_profiles",
            "exclusion_profile": "exclusion_profiles",
            "independent_recompute_profile": "independent_recompute_profiles",
        }
        input_ids = {item["input_id"] for item in profiles["local_input_inventory"]}
        analyzer_ids = {
            item["analyzer_id"] for item in profiles["local_analyzer_inventory"]
        }
        for campaign in drafts["campaigns"]:
            with self.subTest(campaign=campaign["campaign_id"]):
                validate_schema(campaign, schema)
                self.assertEqual(
                    campaign["freeze_status"],
                    "locally_preregistered_run_blocked",
                )
                self.assertEqual(
                    campaign["candidate_bindings"]["status"],
                    "open_before_qualifying_run",
                )
                self.assertIn("open", campaign["failure_artifact_root"])
                self.assertTrue(campaign["open_dependencies"])
                self.assertTrue(set(campaign["local_input_ids"]) <= input_ids)
                self.assertTrue(set(campaign["local_analyzer_ids"]) <= analyzer_ids)
                for field, registry in profile_fields.items():
                    self.assertIn(campaign[field], profiles[registry])

    def test_q023_local_profiles_freeze_seed_and_recompute_plan(self) -> None:
        profiles = _load(Q023_PROFILES_PATH)
        seeds = profiles["seed_profiles"]["Q023-SEEDS-FIXED-8"]
        self.assertEqual(seeds["status"],
                         "frozen_before_qualifying_output_inspection")
        self.assertEqual(seeds["qualifying_seed_count"], 8)
        self.assertEqual(len(seeds["qualifying_seed_list"]), 8)
        self.assertEqual(len(set(seeds["qualifying_seed_list"])), 8)
        self.assertEqual(seeds["pilot_seed_count"], 2)
        self.assertEqual(len(seeds["pilot_seed_list"]), 2)
        self.assertFalse(
            set(seeds["qualifying_seed_list"]) & set(seeds["pilot_seed_list"])
        )

        deltaf = profiles["deltaf_validity_profiles"][
            "Q023-DELTAF-VALIDITY-BOUNDED"
        ]
        self.assertEqual(deltaf["max_absolute_deltaf_weight"], 0.5)
        self.assertTrue(deltaf["matched_reduced_nonlinear_fullf_required"])

        recompute = profiles["independent_recompute_profiles"][
            "Q023-RECOMPUTE-RAW-ARTIFACT-INDEPENDENT"
        ]
        self.assertIn("plan_frozen", recompute["status"])
        self.assertIn("must not import", recompute["implementation_rule"])
        self.assertIn("script_sha256", recompute["required_records"])
        self.assertIn("independent_metric",
                      recompute["required_metric_table_columns"])

    def test_q023_local_inventory_checksums_recompute(self) -> None:
        profiles = _load(Q023_PROFILES_PATH)
        for inventory_name in ("local_analyzer_inventory",
                               "local_input_inventory"):
            for item in profiles[inventory_name]:
                with self.subTest(inventory=inventory_name, path=item["path"]):
                    path = REPO_ROOT / item["path"]
                    self.assertTrue(path.is_file())
                    self.assertEqual(_sha256(path), item["sha256"])


if __name__ == "__main__":
    unittest.main()
