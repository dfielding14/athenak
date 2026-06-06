#!/usr/bin/env python3
"""Fail-closed Q-011 production-campaign contract successor.

This source-local contract binds one nine-output production deck through every
campaign stage, expands the runtime-critical source closure, and freezes the
remaining physical and engineering gates.  It does not authorize a launch,
qualifying-output inspection, scientific acceptance, or claim closure.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Mapping, Sequence


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"
SCHEMA_VERSION = 2
SUCCESSOR_ID = "q011_section54_production_campaign_contract_successor_v2"
CAMPAIGN_ID = "Q011-SECTION54-PRODUCTION-SCIENCE-SUCCESSOR-V2"
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
SELECTED_PRESSURE = {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0}

SUCCESSOR_DECK = (
    "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)
SUCCESSOR_DECK_LAUNCH_PATH = (
    "bindings/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)
SUCCESSOR_DECK_SHA256 = (
    "0e8589618fb8ffdd6da5d8c4aadba8308ac8ce74ac6d0d512c5f216ef7401d79"
)
PREREGISTRATION = (
    "tst/publication/readiness/"
    "q011_section54_qualifying_campaign_preregistration_successor_v4_2026-06-06.json"
)

GRID_VARIANTS = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
QUALIFYING_SEED_POOL = (
    23050101,
    23050102,
    23050103,
    23050104,
    23050105,
    23050106,
    23050107,
    23050108,
)
CORE_QUALIFYING_SEEDS = QUALIFYING_SEED_POOL[:3]
RESERVE_QUALIFYING_SEEDS = QUALIFYING_SEED_POOL[3:]
MIN_FROZEN_TRIAD_COUNT = 3
MAX_FROZEN_TRIAD_COUNT = 8
GLOBAL_NODE_HOUR_LEDGER_CAP = 10000.0
RESOURCE_CONTINGENCY_FRACTION = 0.15
RESOURCE_CONTINGENCY_COMPLETE_TRIAD_FLOOR = 1
PAIRED_INTERVAL_NOMINAL_COVERAGE = 0.80
ENGINEERING_SEEDS = (24060601, 24060602, 24060603, 24060604)
SEED_OVERRIDE_NAMES = (
    "particles/pic_random_seed",
    "problem/ps_inject_seed",
    "problem/ps_seed_noise_seed",
)
VARIANT_OVERRIDES = {
    "coarse_uniform_dx12": (
        "mesh_refinement/refinement=none",
        "mesh_refinement/num_levels=1",
        "problem/ps_enable_curvature_amr=false",
    ),
    "three_level_amr_root_dx12_finest_dx3": (
        "mesh_refinement/refinement=adaptive",
        "mesh_refinement/num_levels=3",
        "problem/ps_enable_curvature_amr=true",
    ),
    "fine_uniform_dx3": (
        "mesh/nx1=16000",
        "mesh/nx2=1040",
        "meshblock/nx1=20",
        "meshblock/nx2=20",
        "mesh_refinement/refinement=none",
        "mesh_refinement/num_levels=1",
        "problem/ps_enable_curvature_amr=false",
    ),
}

STAGE_ORDER = (
    "planner",
    "resource_pilot_materializer",
    "registered_launch_materializer",
    "campaign_execution",
    "admission",
    "analyzer",
)
STAGE_IMPLEMENTATIONS = {
    "planner": (
        "tst/publication/"
        "q011_section54_production_campaign_planner_successor_v2.py"
    ),
    "resource_pilot_materializer": (
        "tst/publication/"
        "q011_resource_scaling_preproduction_pilot_materializer_successor_v2.py"
    ),
    "registered_launch_materializer": (
        "tst/publication/q011_section54_registered_launch_materializer_successor_v2.py"
    ),
    "campaign_execution": (
        "tst/publication/q011_section54_qualifying_campaign_execution_successor_v2.py"
    ),
    "admission": (
        "tst/publication/"
        "q011_section54_production_science_admission_orchestration_successor_v2.py"
    ),
    "analyzer": (
        "tst/publication/q011_section54_production_science_analyzer_successor_v2.py"
    ),
}
STAGE_RECORD_TYPES = {
    "planner": (
        "none_source_local_entrypoint",
        "q011_section54_production_campaign_plan_successor_v2",
    ),
    "resource_pilot_materializer": (
        "q011_section54_production_campaign_plan_successor_v2",
        "q011_section54_resource_pilot_materialization_successor_v2",
    ),
    "registered_launch_materializer": (
        "q011_section54_resource_pilot_materialization_successor_v2",
        "q011_section54_registered_launch_review_candidates_successor_v2",
    ),
    "campaign_execution": (
        "q011_section54_registered_launch_review_candidates_successor_v2",
        "q011_section54_execution_handoffs_successor_v2",
    ),
    "admission": (
        "q011_section54_execution_handoffs_successor_v2",
        "q011_section54_source_local_admission_successor_v2",
    ),
    "analyzer": (
        "q011_section54_source_local_admission_successor_v2",
        "q011_section54_source_local_analysis_gate_packet_successor_v2",
    ),
}

EXPECTED_OUTPUTS: dict[str, dict[str, str]] = {
    "output1": {
        "file_type": "bin",
        "variable": "mhd_w_bcc",
        "id": "mhd_w_bcc",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output2": {
        "file_type": "bin",
        "variable": "prtcl_rho",
        "id": "prtcl_rho",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output3": {
        "file_type": "bin",
        "variable": "prtcl_jx",
        "id": "prtcl_jx",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output4": {
        "file_type": "bin",
        "variable": "prtcl_jy",
        "id": "prtcl_jy",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output5": {
        "file_type": "bin",
        "variable": "prtcl_jz",
        "id": "prtcl_jz",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output6": {
        "file_type": "bin",
        "variable": "mhd_j2",
        "id": "mhd_j2",
        "dt": "100.0",
        "ghost_zones": "false",
    },
    "output7": {
        "file_type": "pvtk",
        "variable": "prtcl_all",
        "id": "prtcl_all",
        "dt": "100.0",
    },
    "output8": {"file_type": "hst", "dt": "10.0"},
    "output9": {"file_type": "rst", "dt": "100.0"},
}

RUNTIME_CRITICAL_SOURCE_GROUPS = {
    "particle_push_deposition_tasks_and_boundaries": (
        "src/particles/particles.cpp",
        "src/particles/particles.hpp",
        "src/particles/particles_data_structs.hpp",
        "src/particles/particles_moments.cpp",
        "src/particles/particles_pushers.cpp",
        "src/particles/particles_tasks.cpp",
        "src/bvals/bvals_part.cpp",
        "src/bvals/bvals_tasks.cpp",
    ),
    "mhd_stages_eos_driver_and_task_scheduler": (
        "src/driver/driver.cpp",
        "src/driver/driver.hpp",
        "src/eos/eos.cpp",
        "src/eos/eos.hpp",
        "src/eos/ideal_mhd.cpp",
        "src/mhd/mhd.cpp",
        "src/mhd/mhd.hpp",
        "src/mhd/mhd_corner_e.cpp",
        "src/mhd/mhd_ct.cpp",
        "src/mhd/mhd_fluxes.cpp",
        "src/mhd/mhd_fofc.cpp",
        "src/mhd/mhd_newdt.cpp",
        "src/mhd/mhd_tasks.cpp",
        "src/mhd/mhd_update.cpp",
        "src/mhd/rsolvers/llf_mhd.hpp",
        "src/srcterms/srcterms.cpp",
        "src/srcterms/srcterms.hpp",
        "src/srcterms/srcterms_newdt.cpp",
        "src/tasklist/task_list.hpp",
    ),
    "amr_refinement_migration_and_load_balancing": (
        "src/bvals/bvals.cpp",
        "src/bvals/bvals.hpp",
        "src/bvals/bvals_cc.cpp",
        "src/bvals/bvals_fc.cpp",
        "src/bvals/bvals_mom.cpp",
        "src/bvals/flux_correct_cc.cpp",
        "src/bvals/flux_correct_fc.cpp",
        "src/bvals/prolongation.cpp",
        "src/bvals/prolong_prims.cpp",
        "src/mesh/build_tree.cpp",
        "src/mesh/load_balance.cpp",
        "src/mesh/mesh.cpp",
        "src/mesh/mesh.hpp",
        "src/mesh/mesh_refinement.cpp",
        "src/mesh/mesh_refinement.hpp",
        "src/mesh/meshblock.cpp",
        "src/mesh/meshblock.hpp",
        "src/mesh/meshblock_pack.cpp",
        "src/mesh/meshblock_pack.hpp",
        "src/mesh/meshblock_tree.cpp",
        "src/mesh/meshblock_tree.hpp",
        "src/mesh/prolongation.hpp",
        "src/mesh/restriction.hpp",
    ),
    "shock_problem_and_nine_output_path": (
        "src/outputs/basetype_output.cpp",
        "src/outputs/binary.cpp",
        "src/outputs/derived_variables.cpp",
        "src/outputs/history.cpp",
        "src/outputs/io_wrapper.cpp",
        "src/outputs/io_wrapper.hpp",
        "src/outputs/outputs.cpp",
        "src/outputs/outputs.hpp",
        "src/outputs/restart.cpp",
        "src/outputs/restart_utils.cpp",
        "src/outputs/restart_utils.hpp",
        "src/outputs/vtk_prtcl.cpp",
        "src/pgen/pgen.cpp",
        "src/pgen/pgen.hpp",
        "src/pgen/tests/pic_parallel_shock.cpp",
    ),
    "campaign_and_analysis_contracts": (
        SUCCESSOR_DECK,
        PREREGISTRATION,
        "tst/publication/analyze_q011_section54_outputs.py",
        "tst/publication/q011_section54_artifacts.py",
        "tst/publication/q011_section54_model.py",
        "tst/publication/q011_section54_morphology_diagnostic_successor_v1.py",
        "tst/publication/q011_section54_particles.py",
        "tst/publication/q011_section54_production_science_successor_v1.py",
        "tst/publication/q011_section54_production_campaign_contract_successor_v2.py",
        *tuple(STAGE_IMPLEMENTATIONS[stage] for stage in STAGE_ORDER),
    ),
}

PROTECTED_HISTORICAL_ARTIFACTS = {
    "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput": (
        "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
    ),
    (
        "tst/publication/readiness/"
        "q011_section54_qualifying_campaign_preregistration_successor_v3_2026-06-06.json"
    ): "1216fc0fcaa78fe423855b9c2ad2039597ef2ed722c46ff8605d184fcedcf1c0",
    "tst/publication/readiness/q011_resource_scaling_preproduction_plan_2026-06-06.json": (
        "f6b059183cd3dd1cfbca55ef2d02300d597bac4e9ac5851aecc3ae79abdebfb5"
    ),
    "tst/publication/q011_resource_scaling_preproduction_pilot_materializer.py": (
        "b8c0e3199b13245409fcef56af2e4c7fbc9e1da675d52f6b0f907a42015784bb"
    ),
    "tst/publication/q011_section54_registered_launch_materializer.py": (
        "8c0a159234238cf602aa5a5b8d9be79d8caad60e06514e92274319b73fbe9f13"
    ),
    "tst/publication/q011_section54_qualifying_campaign_execution.py": (
        "6236898c50e32d1ef5d5bcaf130d638ea7cd572f487095ef73c8d5fea654080e"
    ),
    (
        "tst/publication/"
        "q011_section54_production_science_admission_orchestration_successor_v1.py"
    ): "4768e935078bf7727bc104db17110d8956500da20c21f60a14fda9fcf3f4e6ef",
    "tst/publication/q011_section54_production_science_successor_v1.py": (
        "c18b61990b4e95858da983b50b31334861291c50137440c5fde2b370a4f22b73"
    ),
    "tst/publication/q011_section54_morphology_diagnostic_successor_v1.py": (
        "71cc6e620159f3124fda35327fe04dd4ace0609297d00e035d36b5be62fab314"
    ),
    "tst/publication/test_q011_resource_scaling_preproduction_pilot_materializer.py": (
        "bafc98625e5f739b01cc6ce85f2a32f1cd9e34df0c029f281108c28c5e18af58"
    ),
    "tst/publication/test_q011_section54_registered_launch_materializer.py": (
        "8c0b4d8110ef3353911b053d28f217ecccccfa62d4910ce71b777c2cc8f2665d"
    ),
    "tst/publication/test_q011_section54_qualifying_campaign_execution.py": (
        "27068601ad6685f58241fe000b07369f1a39db05054e5e6c8bb64102105a8e3b"
    ),
    (
        "tst/publication/"
        "test_q011_section54_production_science_admission_orchestration_successor_v1.py"
    ): "57d391f69babd915822cb8e9cad613849f5e3a0ff0400e6687aa002d36f06505",
    "tst/publication/test_q011_section54_production_science_successor_v1.py": (
        "41d0c9d3e8a7efcdd3050b7f4ada13ff2e8f37f985f4dd6cab7cbbd0672e2a5e"
    ),
    "tst/publication/test_q011_section54_morphology_diagnostic_successor_v1.py": (
        "e08c10739dffd52b26a7f83b13fc08496908dcd2f335a32f88393024d9ad96bd"
    ),
    "tst/publication/test_q011_section54_qualifying_campaign_preregistration.py": (
        "28a617b5deec956531752a455f591de44104fed299be5ae6d9b4e9beec0d1e6a"
    ),
}

AUTHORIZATION_BOUNDARY = {
    "launch_authorized": False,
    "scheduler_submission_authorized": False,
    "live_policy_mutation_authorized": False,
    "qualifying_output_inspection_authorized": False,
    "scientific_acceptance_authorized": False,
    "claim_closure_authorized": False,
}

THRESHOLD_SOURCE_CATEGORIES = frozenset(
    {
        "physical-limit",
        "analysis-quality",
        "conservative applicability",
        "engineering closure",
        "literature comparison",
    }
)

def _policy_rule(
    statement: str, source_category: str, rationale: str
) -> dict[str, str]:
    return {
        "rule": statement,
        "source_category": source_category,
        "rationale": rationale,
    }


def campaign_size_resource_design() -> dict[str, object]:
    """Return the preregistered paired-triad resource design."""
    return {
        "variants": list(GRID_VARIANTS),
        "ordered_qualifying_seed_pool": list(QUALIFYING_SEED_POOL),
        "mandatory_core_qualifying_seeds": list(CORE_QUALIFYING_SEEDS),
        "preregistered_reserve_qualifying_seeds": list(RESERVE_QUALIFYING_SEEDS),
        "paired_seed_rule": _policy_rule(
            "Every selected seed must retain a complete coarse-AMR-fine triad.",
            "analysis-quality",
            (
                "Complete paired triads preserve within-seed grid comparisons and "
                "prevent resource pressure from selectively dropping a variant."
            ),
        ),
        "ordered_prefix_rule": _policy_rule(
            "The frozen seed set must be an ordered prefix of the preregistered pool.",
            "analysis-quality",
            (
                "The ordered prefix prevents post-pilot or post-output selection of "
                "apparently favorable random realizations."
            ),
        ),
        "resource_only_freeze_rule": _policy_rule(
            (
                "After excluded exact-deck scaling and I/O pilots pass, and before "
                "any qualifying launch or output inspection, freeze N as the largest "
                "affordable complete paired-triad prefix."
            ),
            "engineering closure",
            (
                "The freeze uses resource evidence only, preserving the largest "
                "affordable preregistered paired design without science-output-driven "
                "sample-size selection."
            ),
        ),
        "post_qualifying_execution_rule": _policy_rule(
            (
                "After the first qualifying launch begins, expansion, contraction, "
                "seed dropping, and science-output-driven sample-size changes are "
                "forbidden; retries retain the frozen seed and complete triad."
            ),
            "analysis-quality",
            (
                "An immutable frozen matrix prevents optional stopping and selective "
                "removal while preserving paired comparisons after runtime failures."
            ),
        ),
        "core_fit_failure_rule": _policy_rule(
            (
                "If the mandatory core plus contingency does not fit the measured "
                "ledger, fail closed and launch no qualifying production attempt."
            ),
            "engineering closure",
            (
                "The mandatory core cannot be weakened to fit the budget, and the "
                "global ledger cap cannot be exceeded."
            ),
        ),
        "numeric_rules": [
            {
                "rule_id": "mandatory_core_triad_count",
                "value": MIN_FROZEN_TRIAD_COUNT,
                "source_category": "analysis-quality",
                "rationale": (
                    "Three paired seeds are the minimum preregistered repeatability "
                    "core; this is an analysis-design choice, not a literature sample "
                    "size."
                ),
            },
            {
                "rule_id": "mandatory_core_attempt_count",
                "value": MIN_FROZEN_TRIAD_COUNT * len(GRID_VARIANTS),
                "source_category": "analysis-quality",
                "rationale": (
                    "Nine attempts are the complete three-variant realization of the "
                    "mandatory three-seed paired core, not a literature-derived count."
                ),
            },
            {
                "rule_id": "maximum_preregistered_triad_count",
                "value": MAX_FROZEN_TRIAD_COUNT,
                "source_category": "engineering closure",
                "rationale": (
                    "Eight paired triads cap the preregistered reserve pool without "
                    "promising an unaffordable full matrix; this is a project resource "
                    "choice with no literature origin."
                ),
            },
            {
                "rule_id": "maximum_preregistered_attempt_count",
                "value": MAX_FROZEN_TRIAD_COUNT * len(GRID_VARIANTS),
                "source_category": "engineering closure",
                "rationale": (
                    "Twenty-four is only the maximum preregistered design envelope, "
                    "not a promised run count or a literature-derived sample size."
                ),
            },
            {
                "rule_id": "global_node_hour_ledger_cap",
                "value": GLOBAL_NODE_HOUR_LEDGER_CAP,
                "source_category": "engineering closure",
                "rationale": (
                    "The global ten-thousand node-hour cap is an AthenaK project "
                    "resource-control limit and has no scientific literature origin."
                ),
            },
            {
                "rule_id": "resource_contingency_fraction",
                "value": RESOURCE_CONTINGENCY_FRACTION,
                "source_category": "engineering closure",
                "rationale": (
                    "The fifteen-percent proportional reserve is an AthenaK "
                    "engineering recovery allowance selected before qualifying "
                    "execution and has no scientific literature origin."
                ),
            },
            {
                "rule_id": "resource_contingency_complete_triad_floor",
                "value": RESOURCE_CONTINGENCY_COMPLETE_TRIAD_FLOOR,
                "source_category": "engineering closure",
                "rationale": (
                    "At least one measured complete-triad cost is reserved so a "
                    "runtime retry does not require seed dropping; this is an "
                    "engineering allowance with no literature origin."
                ),
            },
        ],
    }


def resource_only_freeze_definition() -> dict[str, object]:
    """Return the exact resource-only freeze contract."""
    return {
        "record_type": "q011_section54_resource_only_campaign_size_freeze_successor_v2",
        "required_timing": _policy_rule(
            (
                "after excluded exact-deck scaling and I/O pilots pass and before "
                "any qualifying launch, qualifying output inspection, or qualifying "
                "production reservation"
            ),
            "engineering closure",
            (
                "This ordering isolates campaign-size selection from qualifying "
                "science outputs and prevents unbudgeted production reservations."
            ),
        ),
        "required_resource_evidence": _policy_rule(
            (
                "bind the measured global ledger snapshot, the complete reservation "
                "manifest for every other registered campaign, and the exact-deck "
                "pilot completion receipt"
            ),
            "engineering closure",
            (
                "The freeze must account for measured use and all prior obligations "
                "before assigning the remaining ledger to Q011."
            ),
        ),
        "selection_rule": _policy_rule(
            (
                "select the largest integer N in the preregistered range whose "
                "complete paired-triad production cost and stated contingency fit "
                "the global ledger after all other obligations"
            ),
            "engineering closure",
            (
                "Largest-affordable selection uses only pre-launch resource evidence "
                "and never observed scientific outcomes."
            ),
        ),
        "budget_equation": _policy_rule(
            (
                "global consumed + global currently reserved after all other "
                "campaign reservations + unreserved Q011 nonbaseline obligations + "
                "N times measured conservative complete-triad cost + max(one "
                "complete-triad cost, contingency fraction times N times measured "
                "complete-triad cost) <= global ledger cap"
            ),
            "engineering closure",
            (
                "The equation reserves all measured and registered obligations plus "
                "a retry and proportional recovery allowance before production."
            ),
        ),
        "prohibited_inputs": _policy_rule(
            (
                "qualifying measurements, gate results, physical observables, and "
                "science-output-derived costs or preferences"
            ),
            "analysis-quality",
            (
                "Excluding qualifying science information prevents optional stopping "
                "and outcome-dependent campaign-size selection."
            ),
        ),
        "post_freeze_rule": campaign_size_resource_design()[
            "post_qualifying_execution_rule"
        ],
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def reporting_uncertainty_policy() -> dict[str, object]:
    """Return the preregistered honest small-N reporting policy."""
    return {
        "all_frozen_N": _policy_rule(
            (
                "report every paired seed value, every preregistered paired contrast, "
                "and the observed min-max range without seed dropping"
            ),
            "analysis-quality",
            (
                "Transparent paired reporting exposes realization variability and "
                "prevents a summary statistic from hiding an unfavorable seed."
            ),
        ),
        "N_equals_3": _policy_rule(
            (
                "report paired values, paired contrasts, and ranges; make only "
                "narrow claims about the three preregistered realizations and do not "
                "make broad population inference"
            ),
            "analysis-quality",
            (
                "Three paired realizations do not support a useful distribution-free "
                "population interval, so the honest scope is descriptive and narrow."
            ),
        ),
        "N_greater_than_3": _policy_rule(
            (
                "before launch freeze the exact two-sided distribution-free "
                "order-statistic interval for the median paired contrast, report its "
                "achieved finite-sample coverage, document ties, and avoid "
                "asymptotic-normal or standard-error claims"
            ),
            "analysis-quality",
            (
                "An exact order-statistic interval remains auditable for the small "
                "frozen seed ensemble and does not assume Gaussian paired contrasts."
            ),
        ),
        "inference_scope": _policy_rule(
            (
                "the interval describes repeatability over the preregistered seed "
                "ensemble and does not establish universal astrophysical-parameter "
                "population inference"
            ),
            "conservative applicability",
            (
                "Seed replication quantifies numerical-realization sensitivity but "
                "does not sample the broader physical parameter population."
            ),
        ),
        "numeric_rules": [
            {
                "rule_id": "descriptive_only_triad_count",
                "value": MIN_FROZEN_TRIAD_COUNT,
                "source_category": "analysis-quality",
                "rationale": (
                    "At three paired seeds the exact median interval at the chosen "
                    "coverage is unavailable, so only paired values and ranges support "
                    "narrow claims; this is not a literature threshold."
                ),
            },
            {
                "rule_id": "interval_estimator_minimum_triad_count",
                "value": MIN_FROZEN_TRIAD_COUNT + 1,
                "source_category": "analysis-quality",
                "rationale": (
                    "Four is the smallest frozen N admitting the preregistered exact "
                    "two-sided median interval at the chosen coverage; this is a "
                    "finite-sample analysis property, not a literature sample size."
                ),
            },
            {
                "rule_id": "paired_interval_nominal_coverage",
                "value": PAIRED_INTERVAL_NOMINAL_COVERAGE,
                "source_category": "analysis-quality",
                "rationale": (
                    "Eighty-percent nominal coverage guarantees an exact interval for "
                    "every allowed N above three and permits an inner order-statistic "
                    "interval at N seven and eight; achieved exact coverage is always "
                    "reported and no literature origin is implied."
                ),
            },
        ],
    }


GATE_DEFINITIONS: dict[str, dict[str, object]] = {
    "nonrelativistic_accelerated_tail_validity": {
        "required_products": ["prtcl_all"],
        "metric": "macro_weighted_q999_particle_speed_over_artificial_C_and_maximum_speed_over_C",
        "acceptance_criteria": [
            {
                "criterion_id": "q999_speed_over_C",
                "operator": "<=",
                "value": 0.1,
                "source_category": "conservative applicability",
                "rationale": (
                    "AthenaK-selected applicability bound requiring the weighted "
                    "accelerated tail to remain well below the artificial particle "
                    "light speed; it is not a literature-derived tolerance."
                ),
            },
            {
                "criterion_id": "maximum_speed_over_C",
                "operator": "<=",
                "value": 0.3,
                "source_category": "conservative applicability",
                "rationale": (
                    "AthenaK-selected fail-safe bound on rare particles so a small "
                    "weighted tail cannot hide a relativistic outlier; it is not a "
                    "literature-derived tolerance."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [500.0, 1200.0],
            "source_category": "literature comparison",
            "rationale": (
                "These are the Sun and Bai Section 5.4 comparison epochs; they are "
                "sampling times, not numerical acceptance tolerances."
            ),
        },
    },
    "emax_acceleration_history": {
        "required_products": ["prtcl_all"],
        "metric": "macro_weighted_q999_chi_history_with_linear_fit",
        "acceptance_criteria": [
            {
                "criterion_id": "minimum_finite_snapshot_count",
                "operator": ">=",
                "value": 5,
                "source_category": "analysis-quality",
                "rationale": (
                    "Five snapshots are an AthenaK-selected minimum for a late-time "
                    "trend fit with residual diagnostics; this is not a literature "
                    "tolerance."
                ),
            },
            {
                "criterion_id": "q999_chi_t1200_minus_t500",
                "operator": ">",
                "value": 0.0,
                "source_category": "analysis-quality",
                "rationale": (
                    "A positive endpoint difference is the minimal directional "
                    "quality check for an acceleration claim, not a literature "
                    "prediction of the acceleration rate."
                ),
            },
            {
                "criterion_id": "q999_chi_linear_fit_slope",
                "operator": ">",
                "value": 0.0,
                "source_category": "analysis-quality",
                "rationale": (
                    "A positive fitted slope is the minimal analysis-quality "
                    "condition for calling the measured trend acceleration."
                ),
            },
            {
                "criterion_id": "q999_chi_linear_fit_R2",
                "operator": ">=",
                "value": 0.8,
                "source_category": "analysis-quality",
                "rationale": (
                    "The R2 bound is an AthenaK-selected trend-quality screen chosen "
                    "before output inspection; no cited paper supplies this tolerance."
                ),
            },
            {
                "criterion_id": "nondecreasing_interval_fraction",
                "operator": ">=",
                "value": 0.6,
                "source_category": "analysis-quality",
                "rationale": (
                    "The interval-fraction bound is an AthenaK-selected robustness "
                    "screen allowing turbulent variability; it is not literature "
                    "derived."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [500.0, 700.0, 900.0, 1100.0, 1200.0],
            "source_category": "analysis-quality",
            "rationale": (
                "The endpoints match literature-comparison epochs while the three "
                "intermediate samples are an AthenaK-selected minimum late-time fit "
                "cadence; the times are not pass tolerances."
            ),
        },
    },
    "acceleration_efficiency": {
        "required_products": ["mhd_w_bcc", "prtcl_all", "hst"],
        "metric": "shock_injected_CR_energy_gain_over_integrated_upstream_kinetic_energy_flux",
        "acceptance_criteria": [
            {
                "criterion_id": "acceleration_efficiency_range",
                "operator": "closed_interval",
                "value": [0.0, 1.0],
                "source_category": "physical-limit",
                "rationale": (
                    "For the preregistered efficiency definition, a finite fraction "
                    "outside the available upstream kinetic-energy budget is "
                    "physically inadmissible; this is a definition-level bound, not "
                    "a literature fit."
                ),
            },
            {
                "criterion_id": "admitted_CR_energy_gain",
                "operator": ">",
                "value": 0.0,
                "source_category": "physical-limit",
                "rationale": (
                    "Positive admitted CR energy gain is the physical sign condition "
                    "for an acceleration-efficiency measurement."
                ),
            },
            {
                "criterion_id": "history_snapshot_relative_closure",
                "operator": "<=",
                "value": 0.1,
                "source_category": "engineering closure",
                "rationale": (
                    "The ten-percent cross-reducer closure bound is an AthenaK "
                    "engineering consistency screen selected before output "
                    "inspection; it is not literature derived."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [500.0, 1200.0],
            "source_category": "literature comparison",
            "rationale": (
                "These are published Section 5.4 comparison epochs and do not imply "
                "a literature-derived efficiency tolerance."
            ),
        },
    },
    "compression_and_jump_conditions": {
        "required_products": ["mhd_w_bcc"],
        "metric": "detected_front_relative_compression_and_mass_momentum_energy_flux_residuals",
        "acceptance_criteria": [
            {
                "criterion_id": "median_late_compression_ratio",
                "operator": "closed_interval",
                "value": [3.5, 4.5],
                "source_category": "conservative applicability",
                "rationale": (
                    "The interval is an AthenaK-selected applicability screen around "
                    "the strong gamma=5/3 adiabatic-shock expectation of four; it is "
                    "not a tolerance quoted by Sun and Bai."
                ),
            },
            {
                "criterion_id": "absolute_normalized_mass_momentum_energy_flux_residuals",
                "operator": "<=",
                "value": 0.1,
                "source_category": "engineering closure",
                "rationale": (
                    "The residual bound is an AthenaK engineering closure screen "
                    "selected before output inspection; no cited source supplies it."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [900.0, 1000.0, 1100.0, 1200.0],
            "source_category": "analysis-quality",
            "rationale": (
                "The late-time window is an AthenaK-selected stationarity screen "
                "ending at the published t=1200 comparison epoch; it is not a "
                "literature tolerance."
            ),
        },
    },
    "bell_scale_power_current_field_correlation_and_precursor_morphology": {
        "required_products": [
            "mhd_w_bcc",
            "prtcl_rho",
            "prtcl_jx",
            "prtcl_jy",
            "prtcl_jz",
            "mhd_j2",
            "prtcl_all",
        ],
        "metric": (
            "upstream_transverse_magnetic_power_spectrum_characteristic_scale_"
            "growth_deposited_current_field_statistics_and_bound_t500_morphology"
        ),
        "acceptance_criteria": [
            {
                "criterion_id": "spectra_and_correlations_finite",
                "operator": "required",
                "value": True,
                "source_category": "analysis-quality",
                "rationale": "Finite diagnostics are a prerequisite for any Bell comparison.",
            },
            {
                "criterion_id": "upstream_transverse_magnetic_power_minus_t0_noise_floor",
                "operator": ">",
                "value": 0.0,
                "source_category": "analysis-quality",
                "rationale": (
                    "Positive excess over the measured initial noise floor is the "
                    "minimal analysis-quality condition for resolved amplification; "
                    "it is not a literature amplitude tolerance."
                ),
            },
            {
                "criterion_id": "characteristic_scale_t1200_minus_t500",
                "operator": ">",
                "value": 0.0,
                "source_category": "literature comparison",
                "rationale": (
                    "Growth of the characteristic scale is a qualitative Bell "
                    "literature comparison observable; the zero sign boundary is not "
                    "a claimed literature numerical tolerance."
                ),
            },
            {
                "criterion_id": "bound_precursor_morphology_inventory",
                "operator": "required",
                "value": True,
                "source_category": "literature comparison",
                "rationale": (
                    "Cavities, filaments, corrugation, current, magnetic field, and "
                    "refinement overlays are qualitative comparison products; no "
                    "numeric morphology tolerance is claimed."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [0.0, 500.0, 1200.0],
            "source_category": "literature comparison",
            "rationale": (
                "t=0 supplies the measured noise baseline and t=500/t=1200 are the "
                "published Section 5.4 comparison epochs; these are sampling times."
            ),
        },
    },
    "ideal_vs_detected_shock_sensitivity": {
        "required_products": ["mhd_w_bcc", "prtcl_all"],
        "metric": "paired_key_observables_under_ideal_surface_and_negative_gradient_detected_front",
        "acceptance_criteria": [
            {
                "criterion_id": "relative_q999_chi_difference",
                "operator": "<=",
                "value": 0.1,
                "source_category": "conservative applicability",
                "rationale": (
                    "The ten-percent bound is an AthenaK-selected classifier "
                    "sensitivity screen, not a literature-derived tolerance."
                ),
            },
            {
                "criterion_id": "relative_acceleration_efficiency_difference",
                "operator": "<=",
                "value": 0.1,
                "source_category": "conservative applicability",
                "rationale": (
                    "The ten-percent bound is an AthenaK-selected classifier "
                    "sensitivity screen, not a literature-derived tolerance."
                ),
            },
            {
                "criterion_id": "absolute_compression_ratio_difference",
                "operator": "<=",
                "value": 0.25,
                "source_category": "conservative applicability",
                "rationale": (
                    "The compression-difference bound is an AthenaK-selected "
                    "classifier applicability screen chosen before output inspection."
                ),
            },
            {
                "criterion_id": "classifier_and_particle_count_inventory",
                "operator": "required",
                "value": True,
                "source_category": "engineering closure",
                "rationale": (
                    "Archiving both classifier populations is required to make the "
                    "sensitivity comparison independently auditable."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [500.0, 1200.0],
            "source_category": "literature comparison",
            "rationale": (
                "These are published Section 5.4 comparison epochs; they are not "
                "classifier-sensitivity tolerances."
            ),
        },
    },
    "exact_successor_deck_scaling_and_io_pilots": {
        "required_products": [
            "scheduler_accounting",
            "runtime_telemetry",
            "active_cell_history",
            "meshblock_history",
            "particle_count_and_updates",
            "memory_high_water_mark",
            "output_timers",
            "artifact_byte_counts",
        ],
        "metric": "exact_nine_output_successor_deck_resource_and_io_calibration",
        "acceptance_criteria": [
            {
                "criterion_id": "complete_variant_count",
                "operator": "==",
                "value": 3,
                "source_category": "engineering closure",
                "rationale": (
                    "All coarse, AMR, and fine variants must be calibrated so no "
                    "production variant inherits an unmeasured resource model."
                ),
            },
            {
                "criterion_id": "memory_high_water_fraction",
                "operator": "<=",
                "value": 0.8,
                "source_category": "engineering closure",
                "rationale": (
                    "The twenty-percent memory reserve is an AthenaK engineering "
                    "failure margin, not a literature-derived threshold."
                ),
            },
            {
                "criterion_id": "projected_walltime_fraction",
                "operator": "<=",
                "value": 0.8,
                "source_category": "engineering closure",
                "rationale": (
                    "The twenty-percent scheduler-walltime reserve is an AthenaK "
                    "engineering recovery margin, not a literature threshold."
                ),
            },
            {
                "criterion_id": "complete_pilot_consumed_node_hours",
                "operator": "<=",
                "value": 500.0,
                "source_category": "engineering closure",
                "rationale": (
                    "The cap is the preregistered Q011 engineering-pilot budget "
                    "inside the project-wide 10000 node-hour ceiling; it is unrelated "
                    "to a scientific literature tolerance."
                ),
            },
        ],
        "snapshot_selection": {
            "values_omega0_inverse": [],
            "source_category": "engineering closure",
            "rationale": (
                "Resource pilots use phase-specific runtime endpoints rather than "
                "scientific snapshot-selection thresholds."
            ),
        },
    },
}

_SHA256 = re.compile(r"[0-9a-f]{64}")
_ASSIGNMENT = re.compile(r"^(\s*)([A-Za-z0-9_]+)(\s*=\s*)([^#\n]*?)(\s*(?:#.*)?)$")


class ProductionCampaignContractError(ValueError):
    """Reject drifted, incomplete, or authority-bearing successor contracts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProductionCampaignContractError(message)


def canonical_json_bytes(value: object) -> bytes:
    """Return deterministic finite JSON bytes."""
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ProductionCampaignContractError(
            "contract contains noncanonical JSON values"
        ) from error


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _strict_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(
            _strict_equal(left[key], right[key]) for key in left
        )
    if type(left) is list:
        return len(left) == len(right) and all(
            _strict_equal(lvalue, rvalue) for lvalue, rvalue in zip(left, right)
        )
    return left == right


def _stable_regular_bytes(path: Path, *, label: str) -> bytes:
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute() and lexical == path, f"{label}: path must be absolute")
    try:
        _require(path.resolve(strict=True) == path, f"{label}: path uses a symlink alias")
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ProductionCampaignContractError(f"{label}: file unavailable") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label}: expected regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(path, follow_symlinks=False)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        _require(
            all(getattr(before, name) == getattr(after, name) for name in stable)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size,
            f"{label}: file changed while reading",
        )
        return bytes(payload)
    finally:
        os.close(descriptor)


def _repo_payload(relative: str, *, label: str) -> bytes:
    path = PurePosixPath(relative)
    _require(
        not path.is_absolute()
        and path.as_posix() == relative
        and relative != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe repository-relative path",
    )
    return _stable_regular_bytes(REPO_ROOT / relative, label=label)


def _binding(relative: str, payload: bytes) -> dict[str, str]:
    return {"path": relative, "sha256": hashlib.sha256(payload).hexdigest()}


def _decode_json(payload: bytes, *, label: str) -> dict[str, Any]:
    def reject_constant(value: str) -> None:
        raise ProductionCampaignContractError(
            f"{label}: non-finite JSON constant {value}"
        )

    def reject_duplicates(pairs: Sequence[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ProductionCampaignContractError(f"{label}: invalid JSON") from error
    _require(type(value) is dict, f"{label}: expected JSON object")
    return value


def parse_deck(text: str) -> dict[str, dict[str, str]]:
    """Parse one AthenaK input deck with duplicate rejection."""
    blocks: dict[str, dict[str, str]] = {}
    current: str | None = None
    for line_number, raw in enumerate(text.splitlines(), 1):
        stripped = raw.split("#", 1)[0].strip()
        if not stripped:
            continue
        if stripped.startswith("<") and stripped.endswith(">"):
            current = stripped[1:-1].strip()
            _require(
                bool(current) and current not in blocks,
                f"deck line {line_number}: duplicate or empty block",
            )
            blocks[current] = {}
            continue
        _require(current is not None and "=" in stripped, f"deck line {line_number}: malformed")
        name, value = (item.strip() for item in stripped.split("=", 1))
        _require(
            bool(name) and bool(value) and name not in blocks[current],
            f"deck line {line_number}: duplicate or empty assignment",
        )
        blocks[current][name] = value
    return blocks


def validate_successor_deck(path: Path | None = None) -> dict[str, object]:
    """Validate exact bytes, physical mode, negative-gradient marker, and outputs."""
    selected = REPO_ROOT / SUCCESSOR_DECK if path is None else Path(path)
    payload = _stable_regular_bytes(selected, label="Q011 nine-output successor deck")
    digest = hashlib.sha256(payload).hexdigest()
    _require(digest == SUCCESSOR_DECK_SHA256, "Q011 nine-output successor deck SHA-256 drifted")
    try:
        blocks = parse_deck(payload.decode("utf-8"))
    except UnicodeDecodeError as error:
        raise ProductionCampaignContractError("Q011 successor deck is not UTF-8") from error
    outputs = {name: values for name, values in blocks.items() if name.startswith("output")}
    _require(outputs == EXPECTED_OUTPUTS, "Q011 successor deck nine-output contract drifted")
    marker = blocks.get("q011_section54_production_science_successor_v1", {})
    _require(
        marker.get("shock_front_gradient") == "negative"
        and marker.get("full_state_output") == "mhd_w_bcc"
        and marker.get("qualification_effect")
        == "none_no_launch_no_policy_authorization_no_claim_closure",
        "Q011 successor deck diagnostic marker drifted",
    )
    problem = blocks.get("problem", {})
    particles = blocks.get("particles", {})
    refinement = blocks.get("mesh_refinement", {})
    _require(
        problem.get("ps_p0") == "1.0"
        and problem.get("ps_u0") == "30.0"
        and problem.get("ps_b0") == "1.0"
        and problem.get("ps_eta") == "1.0e-3"
        and problem.get("ps_shock_speed_model") == "ideal_surface"
        and problem.get("ps_enable_injection") == "true"
        and problem.get("ps_enable_gas_subtraction") == "true",
        "Q011 successor deck physical baseline drifted",
    )
    _require(
        particles.get("pusher") == "boris_tsc"
        and particles.get("pic_physical_mode") == PHYSICAL_MODE
        and particles.get("pic_cr_light_speed") == "10000.0"
        and particles.get("deposit_moments") == "true"
        and particles.get("couple_moments_to_mhd") == "true"
        and particles.get("couple_moments_momentum_to_mhd") == "true"
        and particles.get("couple_moments_energy_to_mhd") == "true",
        "Q011 successor deck PIC coupling drifted",
    )
    _require(
        refinement.get("refinement") == "adaptive"
        and refinement.get("num_levels") == "3",
        "Q011 successor deck AMR baseline drifted",
    )
    return {
        "binding": {"path": SUCCESSOR_DECK, "sha256": digest},
        "launch_path": SUCCESSOR_DECK_LAUNCH_PATH,
        "output_count": len(outputs),
        "outputs": outputs,
        "physical_mode": PHYSICAL_MODE,
        "selected_pressure": dict(SELECTED_PRESSURE),
        "shock_front_gradient": "negative",
        "load_balance_cost_per_particle": float(
            particles["pic_load_balance_cost_per_particle"]
        ),
    }


def protected_historical_bindings() -> list[dict[str, str]]:
    """Verify that every protected predecessor remains byte-for-byte unchanged."""
    bindings: list[dict[str, str]] = []
    for relative in sorted(PROTECTED_HISTORICAL_ARTIFACTS):
        payload = _repo_payload(relative, label=f"protected historical artifact {relative}")
        digest = hashlib.sha256(payload).hexdigest()
        _require(
            digest == PROTECTED_HISTORICAL_ARTIFACTS[relative],
            f"protected historical artifact drifted: {relative}",
        )
        bindings.append({"path": relative, "sha256": digest})
    return bindings


def runtime_source_closure() -> dict[str, object]:
    """Bind the runtime-critical PIC, MHD, AMR, load-balance, and analysis path."""
    groups: dict[str, list[dict[str, str]]] = {}
    all_members: list[dict[str, str]] = []
    seen: set[str] = set()
    for group, paths in RUNTIME_CRITICAL_SOURCE_GROUPS.items():
        members: list[dict[str, str]] = []
        for relative in paths:
            _require(relative not in seen, f"runtime source closure duplicates {relative}")
            seen.add(relative)
            payload = _repo_payload(relative, label=f"runtime source closure {relative}")
            member = _binding(relative, payload)
            members.append(member)
            all_members.append(member)
        groups[group] = members
    all_members = sorted(all_members, key=lambda item: item["path"])
    return {
        "groups": groups,
        "member_count": len(all_members),
        "members": all_members,
        "sha256": canonical_sha256(all_members),
    }


def preregistration_binding() -> dict[str, object]:
    """Validate and bind the versioned negative-gradient preregistration."""
    payload = _repo_payload(PREREGISTRATION, label="Q011 successor preregistration")
    value = _decode_json(payload, label="Q011 successor preregistration")
    statistical_resource = validate_statistical_resource_design()
    _require(
        value.get("record_type")
        == "q011_section54_qualifying_campaign_preregistration_successor_v4"
        and value.get("schema_version") == 2
        and value.get("production_deck") == SUCCESSOR_DECK
        and value.get("source_local_gate_ids") == list(GATE_DEFINITIONS)
        and value.get("threshold_provenance_policy", {}).get("allowed_source_categories")
        == sorted(THRESHOLD_SOURCE_CATEGORIES)
        and value.get("campaign_size_resource_design")
        == statistical_resource["campaign_size_resource_design"]
        and value.get("resource_only_freeze_definition")
        == statistical_resource["resource_only_freeze_definition"]
        and value.get("reporting_uncertainty_policy")
        == statistical_resource["reporting_uncertainty_policy"]
        and value.get("shock_front", {}).get("detector")
        == "unique_strongest_negative_density_gradient"
        and value.get("authorization") == AUTHORIZATION_BOUNDARY,
        "Q011 successor preregistration contract drifted",
    )
    return {"binding": _binding(PREREGISTRATION, payload), "record": value}


def _require_finite_json_numbers(value: object, *, label: str) -> None:
    if type(value) is bool or value is None or type(value) is str:
        return
    if type(value) is int:
        return
    if type(value) is float:
        _require(math.isfinite(value), f"{label}: numeric value must be finite")
        return
    if type(value) is list:
        for index, item in enumerate(value):
            _require_finite_json_numbers(item, label=f"{label}[{index}]")
        return
    _require(False, f"{label}: unsupported threshold value type")


def _validate_policy_rule(value: object, *, label: str) -> None:
    _require(
        type(value) is dict
        and set(value) == {"rule", "source_category", "rationale"},
        f"{label}: policy-rule schema drifted",
    )
    _require(
        type(value["rule"]) is str and len(value["rule"].strip()) >= 24,
        f"{label}: policy rule is missing",
    )
    _require(
        value["source_category"] in THRESHOLD_SOURCE_CATEGORIES,
        f"{label}: source category is missing or unsupported",
    )
    _require(
        type(value["rationale"]) is str and len(value["rationale"].strip()) >= 24,
        f"{label}: rationale is missing",
    )


def _validate_numeric_rules(value: object, *, label: str) -> None:
    _require(type(value) is list and bool(value), f"{label}: numeric rules missing")
    seen: set[str] = set()
    for index, rule in enumerate(value):
        rule_label = f"{label}[{index}]"
        _require(
            type(rule) is dict
            and set(rule) == {"rule_id", "value", "source_category", "rationale"},
            f"{rule_label}: numeric-rule schema drifted",
        )
        rule_id = rule["rule_id"]
        _require(
            type(rule_id) is str and bool(rule_id) and rule_id not in seen,
            f"{rule_label}: numeric-rule ID is empty or duplicated",
        )
        seen.add(rule_id)
        _require_finite_json_numbers(rule["value"], label=f"{rule_label}/value")
        _require(
            rule["source_category"] in THRESHOLD_SOURCE_CATEGORIES,
            f"{rule_label}: source category is missing or unsupported",
        )
        _require(
            type(rule["rationale"]) is str and len(rule["rationale"].strip()) >= 24,
            f"{rule_label}: rationale is missing",
        )
        if rule["source_category"] == "literature comparison":
            _require(
                "not" in rule["rationale"].lower(),
                f"{rule_label}: literature-comparison rule implies unsupported origin",
            )


def validate_statistical_resource_design() -> dict[str, dict[str, object]]:
    """Validate resource-freeze and honest-reporting provenance."""
    design = campaign_size_resource_design()
    _require(
        len(GRID_VARIANTS) == 3
        and len(QUALIFYING_SEED_POOL) == MAX_FROZEN_TRIAD_COUNT
        and CORE_QUALIFYING_SEEDS == QUALIFYING_SEED_POOL[:MIN_FROZEN_TRIAD_COUNT]
        and RESERVE_QUALIFYING_SEEDS == QUALIFYING_SEED_POOL[MIN_FROZEN_TRIAD_COUNT:]
        and len(set(QUALIFYING_SEED_POOL)) == len(QUALIFYING_SEED_POOL)
        and set(QUALIFYING_SEED_POOL).isdisjoint(ENGINEERING_SEEDS),
        "Q011 qualifying paired-triad seed design drifted",
    )
    for key in (
        "paired_seed_rule",
        "ordered_prefix_rule",
        "resource_only_freeze_rule",
        "post_qualifying_execution_rule",
        "core_fit_failure_rule",
    ):
        _validate_policy_rule(design[key], label=f"campaign_size_resource_design/{key}")
    _validate_numeric_rules(
        design["numeric_rules"], label="campaign_size_resource_design/numeric_rules"
    )

    freeze = resource_only_freeze_definition()
    for key in (
        "required_timing",
        "required_resource_evidence",
        "selection_rule",
        "budget_equation",
        "prohibited_inputs",
        "post_freeze_rule",
    ):
        _validate_policy_rule(freeze[key], label=f"resource_only_freeze_definition/{key}")
    _require(
        freeze["authorization"] == AUTHORIZATION_BOUNDARY,
        "Q011 resource-only freeze definition acquired authority",
    )

    reporting = reporting_uncertainty_policy()
    for key in ("all_frozen_N", "N_equals_3", "N_greater_than_3", "inference_scope"):
        _validate_policy_rule(reporting[key], label=f"reporting_uncertainty_policy/{key}")
    _validate_numeric_rules(
        reporting["numeric_rules"], label="reporting_uncertainty_policy/numeric_rules"
    )
    return {
        "campaign_size_resource_design": design,
        "resource_only_freeze_definition": freeze,
        "reporting_uncertainty_policy": reporting,
    }


def reporting_rule_for_frozen_n(triad_count: int) -> dict[str, object]:
    """Return the exact preregistered reporting rule for one frozen N."""
    _require(
        type(triad_count) is int
        and MIN_FROZEN_TRIAD_COUNT <= triad_count <= MAX_FROZEN_TRIAD_COUNT,
        "Q011 frozen triad count is outside the preregistered range",
    )
    if triad_count == MIN_FROZEN_TRIAD_COUNT:
        return {
            "frozen_triad_count": triad_count,
            "estimator_id": "none_descriptive_paired_values_and_range_only",
            "interval_estimator_frozen": False,
            "broad_population_inference_authorized": False,
            "source_category": "analysis-quality",
            "rationale": reporting_uncertainty_policy()["N_equals_3"]["rationale"],
        }

    maximum_lower_index = (triad_count + 1) // 2
    selected_lower_index = 0
    selected_coverage = 0.0
    for lower_index in range(1, maximum_lower_index + 1):
        lower_tail = sum(
            math.comb(triad_count, index) for index in range(lower_index)
        ) / (2**triad_count)
        coverage = 1.0 - 2.0 * lower_tail
        if coverage >= PAIRED_INTERVAL_NOMINAL_COVERAGE:
            selected_lower_index = lower_index
            selected_coverage = coverage
    _require(
        selected_lower_index > 0,
        "Q011 frozen N does not admit the preregistered exact paired interval",
    )
    return {
        "frozen_triad_count": triad_count,
        "estimator_id": "exact_two_sided_order_statistic_interval_for_median_paired_contrast",
        "interval_estimator_frozen": True,
        "lower_order_statistic_index_1_based": selected_lower_index,
        "upper_order_statistic_index_1_based": triad_count - selected_lower_index + 1,
        "nominal_coverage": PAIRED_INTERVAL_NOMINAL_COVERAGE,
        "achieved_finite_sample_coverage": selected_coverage,
        "ties": "document_and_report_without_post_freeze_estimator_substitution",
        "broad_population_inference_authorized": False,
        "source_category": "analysis-quality",
        "rationale": reporting_uncertainty_policy()["N_greater_than_3"]["rationale"],
    }


RESOURCE_FREEZE_INPUT_KEYS = frozenset(
    {
        "global_ledger_snapshot_binding",
        "all_other_registered_campaign_reservations_manifest_binding",
        "exact_deck_scaling_io_pilot_completion_binding",
        "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots",
        "global_ledger_currently_reserved_node_hours_after_all_other_campaigns",
        "q011_unreserved_nonbaseline_obligations_node_hours",
        "measured_conservative_complete_paired_triad_node_hours",
        "exact_deck_scaling_io_pilots_passed",
        "exact_deck_scaling_io_pilots_excluded_from_qualifying_science",
        "all_other_registered_campaigns_reserved",
        "qualifying_production_reservations_present",
        "qualifying_launch_started",
        "qualifying_output_inspected",
    }
)


def _normalize_registered_binding(value: object, *, label: str) -> dict[str, str]:
    _require(
        type(value) is dict and set(value) == {"path", "sha256"},
        f"{label}: binding schema drifted",
    )
    relative = value["path"]
    digest = value["sha256"]
    path = PurePosixPath(relative) if type(relative) is str else PurePosixPath(".")
    _require(
        type(relative) is str
        and bool(relative)
        and not path.is_absolute()
        and path.as_posix() == relative
        and relative != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe registered binding path",
    )
    _require(
        type(digest) is str and _SHA256.fullmatch(digest) is not None,
        f"{label}: invalid registered binding SHA-256",
    )
    return {"path": relative, "sha256": digest}


def _normalize_nonnegative_node_hours(value: object, *, label: str) -> float:
    _require(
        type(value) in {int, float}
        and math.isfinite(float(value))
        and float(value) >= 0.0,
        f"{label}: expected finite nonnegative node hours",
    )
    return float(value)


def _normalize_resource_freeze_inputs(value: object) -> dict[str, object]:
    _require(
        type(value) is dict and set(value) == RESOURCE_FREEZE_INPUT_KEYS,
        "Q011 resource-only freeze inputs drifted or contain prohibited fields",
    )
    normalized: dict[str, object] = {}
    for key in (
        "global_ledger_snapshot_binding",
        "all_other_registered_campaign_reservations_manifest_binding",
        "exact_deck_scaling_io_pilot_completion_binding",
    ):
        normalized[key] = _normalize_registered_binding(value[key], label=key)
    for key in (
        "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots",
        "global_ledger_currently_reserved_node_hours_after_all_other_campaigns",
        "q011_unreserved_nonbaseline_obligations_node_hours",
        "measured_conservative_complete_paired_triad_node_hours",
    ):
        normalized[key] = _normalize_nonnegative_node_hours(value[key], label=key)
    _require(
        normalized["measured_conservative_complete_paired_triad_node_hours"] > 0.0,
        "Q011 measured conservative complete-triad cost must be positive",
    )
    required_true = (
        "exact_deck_scaling_io_pilots_passed",
        "exact_deck_scaling_io_pilots_excluded_from_qualifying_science",
        "all_other_registered_campaigns_reserved",
    )
    required_false = (
        "qualifying_production_reservations_present",
        "qualifying_launch_started",
        "qualifying_output_inspected",
    )
    for key in required_true:
        _require(value[key] is True, f"Q011 resource-only freeze requires {key}=true")
        normalized[key] = True
    for key in required_false:
        _require(value[key] is False, f"Q011 resource-only freeze requires {key}=false")
        normalized[key] = False
    return normalized


def build_resource_only_freeze_record(inputs: object) -> dict[str, object]:
    """Freeze the largest affordable paired-triad prefix from resource data only."""
    normalized = _normalize_resource_freeze_inputs(inputs)
    consumed = normalized[
        "global_ledger_cumulative_consumed_node_hours_after_excluded_pilots"
    ]
    reserved = normalized[
        "global_ledger_currently_reserved_node_hours_after_all_other_campaigns"
    ]
    nonbaseline = normalized["q011_unreserved_nonbaseline_obligations_node_hours"]
    triad_cost = normalized["measured_conservative_complete_paired_triad_node_hours"]
    fixed_obligations = consumed + reserved + nonbaseline
    evaluations: list[dict[str, object]] = []
    for triad_count in range(MIN_FROZEN_TRIAD_COUNT, MAX_FROZEN_TRIAD_COUNT + 1):
        production = triad_count * triad_cost
        contingency = max(
            RESOURCE_CONTINGENCY_COMPLETE_TRIAD_FLOOR * triad_cost,
            RESOURCE_CONTINGENCY_FRACTION * production,
        )
        projected = fixed_obligations + production + contingency
        evaluations.append(
            {
                "triad_count": triad_count,
                "complete_paired_triad_production_node_hours": production,
                "stated_contingency_node_hours": contingency,
                "projected_global_node_hours": projected,
                "fits_global_ledger": projected <= GLOBAL_NODE_HOUR_LEDGER_CAP,
            }
        )
    affordable = [
        evaluation for evaluation in evaluations if evaluation["fits_global_ledger"]
    ]
    _require(
        bool(affordable)
        and affordable[0]["triad_count"] == MIN_FROZEN_TRIAD_COUNT,
        "Q011 mandatory three-triad core plus contingency does not fit the ledger",
    )
    selected = affordable[-1]["triad_count"]
    _require(
        all(
            evaluation["fits_global_ledger"]
            == (evaluation["triad_count"] <= selected)
            for evaluation in evaluations
        ),
        "Q011 resource-only affordability is not a monotone paired-triad prefix",
    )
    return {
        "record_type": "q011_section54_resource_only_campaign_size_freeze_successor_v2",
        "schema_version": SCHEMA_VERSION,
        "campaign_id": CAMPAIGN_ID,
        "status": "resource_only_freeze_complete_execution_blocked",
        "qualification_effect": (
            "campaign_size_freeze_only_no_launch_no_qualifying_output_inspection_"
            "no_scientific_acceptance_no_claim_closure"
        ),
        "resource_only_inputs": normalized,
        "resource_only_freeze_definition_sha256": canonical_sha256(
            resource_only_freeze_definition()
        ),
        "campaign_size_resource_design_sha256": canonical_sha256(
            campaign_size_resource_design()
        ),
        "fixed_obligations_node_hours": fixed_obligations,
        "global_node_hour_cap": {
            "value": GLOBAL_NODE_HOUR_LEDGER_CAP,
            "source_category": "engineering closure",
            "rationale": (
                "The global cap is the registered AthenaK project resource limit and "
                "has no scientific literature origin."
            ),
        },
        "candidate_budget_evaluations": evaluations,
        "frozen_triad_count": selected,
        "selected_qualifying_seeds": list(QUALIFYING_SEED_POOL[:selected]),
        "unselected_preregistered_reserve_seeds": list(QUALIFYING_SEED_POOL[selected:]),
        "frozen_attempt_count": selected * len(GRID_VARIANTS),
        "reporting_rule": reporting_rule_for_frozen_n(selected),
        "post_qualifying_execution_rule": campaign_size_resource_design()[
            "post_qualifying_execution_rule"
        ],
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_resource_only_freeze_record(value: object) -> dict[str, object]:
    """Require a deterministic, resource-only, pre-launch campaign-size freeze."""
    _require(type(value) is dict, "Q011 resource-only freeze record must be an object")
    expected = build_resource_only_freeze_record(value.get("resource_only_inputs"))
    _require(
        _strict_equal(value, expected),
        "Q011 resource-only freeze record drifted, used science inputs, or acquired authority",
    )
    return expected


def validate_gate_definition_threshold_provenance() -> dict[str, dict[str, object]]:
    """Require an explicit allowed provenance category and rationale for every gate."""
    for gate_id, definition in GATE_DEFINITIONS.items():
        _require(
            set(definition)
            == {
                "required_products",
                "metric",
                "acceptance_criteria",
                "snapshot_selection",
            },
            f"{gate_id}: gate-definition keys drifted",
        )
        criteria = definition["acceptance_criteria"]
        _require(type(criteria) is list and bool(criteria), f"{gate_id}: criteria missing")
        seen: set[str] = set()
        for index, criterion in enumerate(criteria):
            label = f"{gate_id}/acceptance_criteria[{index}]"
            _require(
                type(criterion) is dict
                and set(criterion)
                == {
                    "criterion_id",
                    "operator",
                    "value",
                    "source_category",
                    "rationale",
                },
                f"{label}: criterion schema drifted",
            )
            criterion_id = criterion["criterion_id"]
            _require(
                type(criterion_id) is str
                and bool(criterion_id)
                and criterion_id not in seen,
                f"{label}: criterion ID is empty or duplicated",
            )
            seen.add(criterion_id)
            _require(
                criterion["source_category"] in THRESHOLD_SOURCE_CATEGORIES,
                f"{label}: threshold source category is missing or unsupported",
            )
            rationale = criterion["rationale"]
            _require(
                type(rationale) is str and len(rationale.strip()) >= 24,
                f"{label}: threshold rationale is missing",
            )
            _require_finite_json_numbers(criterion["value"], label=f"{label}/value")
            if criterion["source_category"] == "literature comparison":
                lowered = rationale.lower()
                _require(
                    "comparison" in lowered
                    and (
                        "not" in lowered
                        or "no numeric" in lowered
                        or "no numerical" in lowered
                    ),
                    f"{label}: literature-comparison criterion implies an unsupported tolerance",
                )
        selection = definition["snapshot_selection"]
        _require(
            type(selection) is dict
            and set(selection)
            == {"values_omega0_inverse", "source_category", "rationale"},
            f"{gate_id}: snapshot-selection schema drifted",
        )
        _require(
            selection["source_category"] in THRESHOLD_SOURCE_CATEGORIES,
            f"{gate_id}: snapshot-selection source category is missing",
        )
        _require(
            type(selection["rationale"]) is str
            and len(selection["rationale"].strip()) >= 24,
            f"{gate_id}: snapshot-selection rationale is missing",
        )
        _require_finite_json_numbers(
            selection["values_omega0_inverse"],
            label=f"{gate_id}/snapshot_selection/values_omega0_inverse",
        )
    return GATE_DEFINITIONS


def source_local_gate_records() -> list[dict[str, object]]:
    """Return frozen, unevaluated, fail-closed gate records."""
    records = []
    definitions = validate_gate_definition_threshold_provenance()
    for gate_id, definition in definitions.items():
        records.append(
            {
                "record_type": "q011_section54_source_local_gate_record_successor_v2",
                "schema_version": SCHEMA_VERSION,
                "gate_id": gate_id,
                "status": "blocked_pending_registered_runtime_evidence",
                "passed": False,
                "measurements_inspected": False,
                "qualifying_evidence": False,
                "definition": definition,
                "blockers": [
                    "registered_runtime_evidence_absent",
                    "raw_artifact_inventory_absent",
                    "exact_successor_deck_scaling_and_io_pilot_not_closed",
                ],
                "authorization": dict(AUTHORIZATION_BOUNDARY),
            }
        )
    return records


def validate_source_local_gate_records(value: object) -> list[dict[str, object]]:
    expected = source_local_gate_records()
    _require(
        _strict_equal(value, expected),
        "Q011 source-local gate records drifted or acquired authority/results",
    )
    return expected


def source_local_blockers(resource_only_freeze_record: object | None = None) -> list[str]:
    blockers = [
        "exact_combined_MHD_CR_conservation_not_integrated_for_dynamic_AMR_and_load_balancing",
        "registered_MPI_HIP_dynamic_AMR_migration_restart_and_load_balance_evidence_absent",
        "particle_load_balance_cost_per_particle_is_zero_and_requires_pilot_selection",
        "registered_production_attempts_and_admitted_raw_artifacts_absent",
        "all_source_local_scientific_gate_records_are_unevaluated",
    ]
    if resource_only_freeze_record is None:
        blockers[3:3] = [
            "exact_nine_output_successor_deck_scaling_and_io_pilots_not_executed",
            "resource_only_campaign_size_freeze_not_recorded_before_qualifying_launch",
        ]
    else:
        validate_resource_only_freeze_record(resource_only_freeze_record)
    return blockers


def _stage_record(stage_id: str, basis: Mapping[str, object]) -> dict[str, object]:
    _require(stage_id in STAGE_ORDER, f"unknown Q011 successor stage {stage_id!r}")
    index = STAGE_ORDER.index(stage_id)
    input_type, output_type = STAGE_RECORD_TYPES[stage_id]
    return {
        "record_type": f"q011_section54_{stage_id}_contract_successor_v2",
        "schema_version": SCHEMA_VERSION,
        "stage_id": stage_id,
        "stage_index": index,
        "status": "source_local_fail_closed",
        "implementation": basis["stage_implementations"][stage_id],
        "predecessor_stage": None if index == 0 else STAGE_ORDER[index - 1],
        "input_record_type": input_type,
        "emitted_record_type": output_type,
        "production_deck": basis["production_deck"],
        "production_deck_launch_path": SUCCESSOR_DECK_LAUNCH_PATH,
        "required_output_contract_sha256": basis["output_contract_sha256"],
        "runtime_source_closure_sha256": basis["runtime_source_closure_sha256"],
        "physical_preregistration": basis["physical_preregistration"],
        "source_local_gate_records_sha256": basis["source_local_gate_records_sha256"],
        "campaign_size_resource_design_sha256": basis[
            "campaign_size_resource_design_sha256"
        ],
        "resource_only_freeze_definition_sha256": basis[
            "resource_only_freeze_definition_sha256"
        ],
        "reporting_uncertainty_policy_sha256": basis[
            "reporting_uncertainty_policy_sha256"
        ],
        "blockers": list(basis["blockers"]),
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def _basis() -> dict[str, object]:
    deck = validate_successor_deck()
    closure = runtime_source_closure()
    preregistration = preregistration_binding()
    gates = source_local_gate_records()
    statistical_resource = validate_statistical_resource_design()
    implementations = {
        stage: _binding(
            STAGE_IMPLEMENTATIONS[stage],
            _repo_payload(
                STAGE_IMPLEMENTATIONS[stage],
                label=f"Q011 successor stage implementation {stage}",
            ),
        )
        for stage in STAGE_ORDER
    }
    return {
        "production_deck": deck["binding"],
        "output_contract_sha256": canonical_sha256(deck["outputs"]),
        "runtime_source_closure_sha256": closure["sha256"],
        "physical_preregistration": preregistration["binding"],
        "source_local_gate_records_sha256": canonical_sha256(gates),
        "campaign_size_resource_design_sha256": canonical_sha256(
            statistical_resource["campaign_size_resource_design"]
        ),
        "resource_only_freeze_definition_sha256": canonical_sha256(
            statistical_resource["resource_only_freeze_definition"]
        ),
        "reporting_uncertainty_policy_sha256": canonical_sha256(
            statistical_resource["reporting_uncertainty_policy"]
        ),
        "stage_implementations": implementations,
        "blockers": source_local_blockers(),
    }


def build_stage_contract(stage_id: str) -> dict[str, object]:
    """Build one exact source-local stage contract."""
    return _stage_record(stage_id, _basis())


def validate_stage_contract(stage_id: str, value: object) -> dict[str, object]:
    expected = build_stage_contract(stage_id)
    _require(_strict_equal(value, expected), f"Q011 {stage_id} stage contract drifted")
    return expected


def build_pipeline_contract() -> dict[str, object]:
    """Build the complete source-local end-to-end successor contract."""
    deck = validate_successor_deck()
    historical = protected_historical_bindings()
    closure = runtime_source_closure()
    preregistration = preregistration_binding()
    gates = source_local_gate_records()
    statistical_resource = validate_statistical_resource_design()
    basis = _basis()
    stages = [_stage_record(stage, basis) for stage in STAGE_ORDER]
    return {
        "record_type": "q011_section54_production_campaign_pipeline_contract_successor_v2",
        "schema_version": SCHEMA_VERSION,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": CAMPAIGN_ID,
        "status": "source_local_contract_complete_execution_and_science_blocked",
        "qualification_effect": (
            "contract_repair_only_no_launch_no_policy_no_qualifying_output_"
            "inspection_no_science_no_claim_closure"
        ),
        "physical_mode": PHYSICAL_MODE,
        "selected_pressure": dict(SELECTED_PRESSURE),
        "production_deck": deck,
        "campaign_size_resource_design": statistical_resource[
            "campaign_size_resource_design"
        ],
        "resource_only_freeze_definition": statistical_resource[
            "resource_only_freeze_definition"
        ],
        "reporting_uncertainty_policy": statistical_resource[
            "reporting_uncertainty_policy"
        ],
        "stage_order": list(STAGE_ORDER),
        "stages": stages,
        "runtime_source_closure": closure,
        "negative_gradient_preregistration": preregistration["binding"],
        "source_local_gate_records": gates,
        "threshold_provenance_policy": {
            "allowed_source_categories": sorted(THRESHOLD_SOURCE_CATEGORIES),
            "numeric_threshold_rule": (
                "every numeric science or admission criterion must carry one explicit "
                "source category and rationale; literature-comparison criteria must "
                "state that no unsupported literature tolerance is implied"
            ),
        },
        "protected_historical_artifacts": historical,
        "blockers": source_local_blockers(),
        "authorization": dict(AUTHORIZATION_BOUNDARY),
    }


def validate_pipeline_contract(value: object) -> dict[str, object]:
    expected = build_pipeline_contract()
    _require(
        _strict_equal(value, expected),
        "Q011 production-campaign pipeline contract drifted or acquired authority",
    )
    return expected


def expected_attempt_argv(attempt_id: str, variant: str, seed: int) -> list[str]:
    """Return the exact launch-prohibited argv for one preregistered attempt."""
    _require(
        variant in GRID_VARIANTS,
        "Q011 attempt variant is not in the preregistered paired-triad matrix",
    )
    _require(
        seed in QUALIFYING_SEED_POOL,
        "Q011 attempt seed is not in the preregistered paired-triad pool",
    )
    _require(
        re.fullmatch(r"[a-z0-9][a-z0-9._-]{0,127}", attempt_id) is not None,
        "Q011 attempt ID is unsafe",
    )
    return [
        "-i",
        SUCCESSOR_DECK_LAUNCH_PATH,
        "-d",
        f"<registered-orion-raw-root>/{attempt_id}",
        f"job/basename={attempt_id}",
        "problem/ps_p0=1.0",
        *(f"{name}={seed}" for name in SEED_OVERRIDE_NAMES),
        *VARIANT_OVERRIDES[variant],
    ]


def render_successor_deck(overrides: Mapping[tuple[str, str], str]) -> bytes:
    """Render an engineering-only deck while preserving all unspecified bytes."""
    payload = _repo_payload(SUCCESSOR_DECK, label="Q011 nine-output successor deck")
    text = payload.decode("utf-8")
    remaining = dict(overrides)
    current: str | None = None
    rendered: list[str] = []
    for raw in text.splitlines():
        stripped = raw.split("#", 1)[0].strip()
        if stripped.startswith("<") and stripped.endswith(">"):
            current = stripped[1:-1].strip()
            rendered.append(raw)
            continue
        match = _ASSIGNMENT.fullmatch(raw)
        if match is None or current is None:
            rendered.append(raw)
            continue
        prefix, name, separator, _, suffix = match.groups()
        key = (current, name)
        if key not in remaining:
            rendered.append(raw)
            continue
        rendered.append(f"{prefix}{name}{separator}{remaining.pop(key)}{suffix}")
    _require(not remaining, f"Q011 engineering deck overrides missing keys: {sorted(remaining)}")
    result = ("\n".join(rendered) + "\n").encode("utf-8")
    actual = parse_deck(result.decode("utf-8"))
    base = parse_deck(text)
    expected = {block: dict(values) for block, values in base.items()}
    for (block, name), value in overrides.items():
        expected[block][name] = value
    _require(actual == expected, "Q011 rendered engineering deck drifted outside overrides")
    return result


def main() -> None:
    print(canonical_json_bytes(build_pipeline_contract()).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
