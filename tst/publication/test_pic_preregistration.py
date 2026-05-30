#!/usr/bin/env python3
"""Regression tests for PIC Q-022/Q-023 preregistration sidecars."""

from __future__ import annotations

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
    "XCMP-EXT-CRPAI-TRANSPORT",
}


def _load(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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
                else:
                    self.assertIn("open_placeholder",
                                  comparison["equation_map_status"])
                    self.assertIn("open_placeholder",
                                  comparison["normalization_map_status"])

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

    def test_q022_entity_map_matches_schema_and_stays_bounded(self) -> None:
        schema = _load(
            READINESS_DIR / "schemas" / "q022_equation_normalization_map.schema.json"
        )
        entity_map = _load(Q022_ENTITY_MAP_PATH)
        validate_schema(entity_map, schema)
        self.assertEqual(entity_map["comparison_id"], "XCMP-ENTITY-MICRO")
        self.assertEqual(entity_map["reviewer_disposition"],
                         "pending external review")
        self.assertIn("MHD feedback", entity_map["excluded_regimes"])
        self.assertTrue(any("not claimed interchangeable" in mismatch
                            for mismatch in entity_map["intentional_mismatches"]))

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
