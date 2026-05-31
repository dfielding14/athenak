#!/usr/bin/env python3
"""Regression for the bounded Q-032 Plotnikov applicability derivation."""

from __future__ import annotations

import json
from pathlib import Path
import unittest


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS = REPO_ROOT / "tst/publication/readiness"
DERIVATION = (
    READINESS / "q032_plotnikov_reduced_map_applicability_derivation_2026-05-30.md"
)
INGEST = READINESS / "q022_external_reference_private_ingest_2026-05-30.json"
PROVENANCE = READINESS / "q022_dataset_provenance_manifest_2026-05-30.json"
EQUATION_MAP = (
    READINESS / "q022_xcmp_ext_crsi_in_damping_equation_map_2026-05-30.json"
)
TOLERANCES = (
    READINESS / "q022_xcmp_ext_crsi_in_damping_tolerance_table_2026-05-30.json"
)
RUNTIME = READINESS / "q032_reduced_static_neutral_runtime_local_2026-05-30.json"
PARSER_HARDENING = (
    READINESS / "q032_plotnikov_damped_crsi_parser_hardening_successor_2026-05-30.json"
)

REFERENCE_ID = "plotnikov_ostriker_bai_2021_arxiv_2102.11878"
COMPARISON_ID = "XCMP-EXT-CRSI-IN-DAMPING"
DATASET_ID = "Q022-DATASET-XCMP-EXT-CRSI-IN-DAMPING"


def _load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


class Q032PlotnikovReducedMapApplicabilityDerivationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.derivation = DERIVATION.read_text(encoding="utf-8")

    def test_artifact_remains_explicitly_nonqualifying(self) -> None:
        for marker in (
            "artifact_role = source_local_applicability_derivation_only",
            "qualification_effect = none",
            "claim_closure = false",
            "plotnikov_matched_qualification = false",
            "No extracted reference values, matched equation mappings, or numeric",
            "tolerances are supplied by this artifact.",
        ):
            with self.subTest(marker=marker):
                self.assertIn(marker, self.derivation)

    def test_source_substep_is_not_wave_amplitude_or_energy_asymptote(self) -> None:
        for marker in (
            "D_src = exp(-nu_in dt)",
            "u_ion,perp' = exp(-nu_in dt) u_ion,perp",
            "K_ion,perp' = exp(-2 nu_in dt) K_ion,perp",
            "-Im(omega) = Gamma_d = nu_in / 2",
            "D_wave_amp(t) = exp(-nu_in t / 2)",
            "D_wave_energy(t) = exp(-nu_in t)",
            "D_wave_amp(dt) = sqrt(D_src(dt))",
            "D_wave_energy(dt) = D_src(dt)",
            "D_ion_perp_kinetic_energy(dt) = D_src(dt)^2",
            "conditional algebraic",
            "coincidence between different observables",
            "not permission to identify the",
            "source sink with magnetic wave energy.",
        ):
            with self.subTest(marker=marker):
                self.assertIn(marker, self.derivation)

    def test_phase_scrambling_requirement_stays_split_and_fail_closed(self) -> None:
        self.assertIn(
            "particle_phase_scrambling_for_source_local_sink_derivation = not_required",
            self.derivation,
        )
        self.assertIn(
            "particle_phase_scrambling_for_plotnikov_matched_campaign = "
            "unresolved_fail_closed",
            self.derivation,
        )
        self.assertIn(
            "either bind and validate the required phase-scrambling",
            self.derivation,
        )
        self.assertIn(
            "or obtain a reviewed, scoped justification for omitting it",
            self.derivation,
        )

    def test_archived_plotnikov_pdf_remains_a_source_reference_only(self) -> None:
        ingest = _load_json(INGEST)
        reference = next(
            item for item in ingest["artifacts"] if item["reference_id"] == REFERENCE_ID
        )
        self.assertEqual(reference["path"], "arxiv_2102.11878.pdf")
        self.assertEqual(reference["source_locator"], "arXiv:2102.11878")
        self.assertEqual(len(reference["sha256"]), 64)
        self.assertIn("source-reference artifacts", ingest["scope"])
        self.assertIn("not digitized comparison datasets", ingest["scope"])
        self.assertIn("source reference only", self.derivation)

    def test_q022_dataset_map_and_tolerances_remain_blocked_and_empty(self) -> None:
        provenance = _load_json(PROVENANCE)
        dataset = next(
            item
            for item in provenance["dataset_candidates"]
            if item["dataset_id"] == DATASET_ID
        )
        self.assertEqual(dataset["comparison_id"], COMPARISON_ID)
        self.assertEqual(dataset["reference_ids"], [REFERENCE_ID])
        self.assertEqual(
            dataset["extraction_status"], "blocked_extraction_input_unavailable"
        )
        self.assertEqual(dataset["reviewer_disposition"], "pending external review")

        equation_map = _load_json(EQUATION_MAP)
        self.assertEqual(equation_map["comparison_id"], COMPARISON_ID)
        self.assertEqual(
            equation_map["map_status"],
            "blocked_pending_reference_specific_mapping_and_external_review",
        )
        self.assertEqual(equation_map["matched_equations"], [])
        self.assertEqual(equation_map["intentional_mismatches"], [])
        self.assertEqual(equation_map["unit_map"], {})
        self.assertEqual(equation_map["normalization_map"], {})
        self.assertEqual(equation_map["parameter_overlap"], {})
        self.assertEqual(equation_map["excluded_regimes"], [])

        tolerances = _load_json(TOLERANCES)
        self.assertEqual(tolerances["comparison_id"], COMPARISON_ID)
        self.assertEqual(tolerances["dataset_provenance_id"], DATASET_ID)
        self.assertEqual(
            tolerances["freeze_status"],
            "blocked_pending_reference_dataset_extraction_and_external_review",
        )
        self.assertEqual(tolerances["rows"], [])

    def test_existing_q032_readiness_records_stay_nonqualifying(self) -> None:
        runtime = _load_json(RUNTIME)
        self.assertEqual(runtime["qualification_effect"], "none")
        self.assertFalse(runtime["claim_closure"])
        reduced_map = runtime["generator"]["reduced_map"]
        self.assertEqual(reduced_map["attenuation_factor"], "exp(-nu_in*dt)")
        self.assertEqual(
            reduced_map["physical_scope"],
            "reduced_static_neutral_high_frequency_transverse_friction_only",
        )
        result = runtime["source_local_runtime_probe"]["result"]
        self.assertEqual(result["qualification_effect"], "none")
        self.assertEqual(result["plotnikov_qualification"], "not_claimed")

        parser_hardening = _load_json(PARSER_HARDENING)
        self.assertFalse(parser_hardening["claim_closure"])
        boundary = parser_hardening["plotnikov_comparison_boundary"]
        self.assertEqual(boundary["dataset_status"],
                         "blocked_extraction_input_unavailable")
        self.assertEqual(
            boundary["equation_map_status"],
            "blocked_pending_reference_specific_mapping_and_external_review",
        )
        self.assertEqual(boundary["numeric_tolerance_rows"], [])
        launch_block = parser_hardening["launch_block"]
        self.assertFalse(launch_block["analyzer_always_emits_qualifying_evidence"])
        self.assertFalse(
            launch_block["analyzer_can_emit_matched_plotnikov_qualification"]
        )


if __name__ == "__main__":
    unittest.main()
