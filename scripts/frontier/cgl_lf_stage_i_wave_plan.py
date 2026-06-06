#!/usr/bin/env python3
"""Emit an evidence-bound, deterministic, read-only CGL-LF Stage I wave plan.

The planner never launches work and never mutates campaign state.  It consumes
only one promoted schema-2 recost recommendation and its exact independent
review/publication audit.  That publication is the sole current authority for
the wave profiles, promoted controller/source identity, projection,
reconciliation, and any R17 readiness chain.  Current source identity must in
turn be authorized by the exact independently reviewed F-118 supersession
chain.  The planner independently re-authenticates its immutable F-116/F-115
predecessors, canonical stores, and
preserves the immutable historical F-115 R03 continuation checks without
treating F-115 as current controller identity.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation
import hashlib
import io
import json
import math
import os
from pathlib import Path
import pwd
import re
import stat
import struct
import subprocess
import sys
from typing import Iterable


EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
CAMPAIGN = "mks24_stage_i_cgl_lf_reproduction"

CANONICAL_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
FROZEN_SOURCE = Path("/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422")
CONTROLLER_HELPER = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-df/scripts/frontier/cgl_lf_stage_i.py"
)
SQUEUE = Path("/usr/bin/squeue")
SCONTROL = Path("/usr/bin/scontrol")
GIT = Path("/usr/lib/git/git")
GIT_EXEC_PATH = Path("/usr/lib/git")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
ACCOUNT = "AST207"
PARTITION = "batch"

SOURCE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
EXECUTABLE_SHA256 = "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
BUILD_FILE_SHA256 = {
    "athena.sha256": "8dedf6bbb74aaa1d62a0f1d9c46b050d254085d33843134d8b4278ee5777d4ab",
    "athena-config.txt": "45abe6e6a8952f300808f941425855808f47efdda43fd2233db74732370e84f1",
    "environment.txt": "d6ae3869da25a670290a8f4d023861c69640418fc4ad287e0750bf345cdebbe7",
    "CMakeCache.txt": "66cafe04d8f6cf5c74139d81fc7d3c5ba7cf875639653dfff33be3bde9f4b3f1",
}

MAX_LANES = 4
MAX_NODES = 10
RANKS_PER_NODE = 8
CPUS_PER_TASK = 7
MAX_SEGMENT_SECONDS = 2 * 60 * 60
SHUTDOWN_MARGIN_SECONDS = 10 * 60
STAGE_I_BUDGET_NODE_HOURS = Decimal("1400")
PROJECT_BUDGET_NODE_HOURS = Decimal("4000")
NODE_HOUR_PROJECTION_METHOD = "observed-stage-i-scoped-node-hour-rate-v2"
STATE_MAX_AGE = timedelta(minutes=20)
R17_READINESS_MAX_AGE = timedelta(hours=24)
FUTURE_SKEW = timedelta(minutes=5)
R17_REQUIRED_RETENTION_BYTES = 958271710272
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
RESTART_MESH_HEADER_SIZE = 252
RESTART_TIME_OFFSET = 232
RESTART_TIME_FORMAT = "<d"
ALLOWED_RESTART_MARKER_MODES = frozenset(("full_precision", "legacy_default_precision"))
CONTINUATION_MASS_TOLERANCE = 1.0e-12
CONTINUATION_MAX_NORMALIZED_CT_DIVB = 1.0e-12
CONTINUATION_ACTIVITY_ABSOLUTE_GT = 1.0e-6
CONTINUATION_ACTIVITY_NORMALIZED_GT = 1.0e-8
CONTINUATION_PLASMA_POLICY = "stage-i-clean-partial-continuation-v2"
FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY = (
    "stage-i-frozen-e03-exact-retained-continuation-migration-v1"
)
FROZEN_E03_CT_DIVERGENCE_REASON = (
    "the qualified historical E03 executable did not retain normalized CT "
    "divB in its exact user-history schema"
)
R12_HISTORICAL_PARTIAL_JOB_ID = "4766856"
R12_HISTORICAL_PARTIAL_SEGMENT = "s00_rankio_t0_t0p25"
R12_FRESH_RERUN_SEGMENT = "s01_rankio_t0_t0p12"
R12_FRESH_RERUN_TARGET = Decimal("0.12")
R12_FRESH_RERUN_NODES = 4
R12_FRESH_RERUN_RANKS = 32
R12_FRESH_RERUN_WALLTIME = "02:00:00"
R12_FRESH_RERUN_ATHENA_WALLTIME = "01:50:00"
R12_FRESH_RERUN_POLICY = "stage-i-r12-fresh-rerun-after-inventory-only-partial-v2"
R12_CONTINUATION_ALIGNMENT = Decimal("0.02")
STANDARD_CONTINUATION_ALIGNMENT = Decimal("0.25")
SCOPED_MEASUREMENT_BASIS_KEYS = (
    "job_id",
    "case_id",
    "segment",
    "actual_node_hours",
    "observed_cells",
    "observed_simulation_interval",
    "normalized_node_hours_per_cell_per_simulation_time",
)
SCOPED_CASE_BREAKDOWN_KEYS = (
    "matrix_full_case_node_hours_reference_only",
    "authenticated_progress_fraction",
    "remaining_simulation_time",
    "projected_cells",
    "projection_measurement_basis",
    "observed_rate_projected_remaining_node_hours",
    "authorized_profile_reserved_node_hours",
    "projected_remaining_node_hours",
)
SCOPED_BUDGET_KEYS = (
    "method",
    "measurement_basis",
    "actual_stage_i_node_hours",
    "authorized_wave_reserved_node_hours",
    "actual_plus_authorized_wave_node_hours",
    "computed_remaining_stage_i_node_hours",
    "computed_stage_i_total_node_hours",
    "promoted_stage_i_envelope_node_hours",
    "project_ceiling_node_hours",
    "computed_stage_i_margin_node_hours",
    "case_breakdown",
)
REQUEST_REVIEW_IDENTITY_ASSURANCE = "declared-process-independence-non-cryptographic"
REQUEST_REVIEW_IDENTITY_LIMITATION = (
    "Reviewer identity and process independence are declared evidence, "
    "not cryptographically proven."
)
CONTROLLER_CLEAN_PARTIAL_CHECKS = {
    "required_time_reached": False,
    "strict_lf_failure_counters_zero": True,
    "snapshots_retained": True,
    "terminal_snapshot_retained": True,
    "restart_retained": True,
    "terminal_restart_physical_time_matches_final": True,
    "plasma_continuation_policy": True,
}
FROZEN_E03_NO_MAX_NDIV_CHECKS = {
    key: value
    for key, value in CONTROLLER_CLEAN_PARTIAL_CHECKS.items()
    if key != "plasma_continuation_policy"
}
CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS = frozenset(
    (
        "schema_version", "execution_epoch", "inspected_utc", "manifest", "job_id",
        "case_id", "segment", "required_time", "final_time",
        "maximum_strict_failure_counts", "checks", "accepted",
        "clean_for_continuation", "mhd_history", "user_history",
        "plasma_continuation_policy", "plasma_continuation_evidence",
        "snapshots", "snapshot_times", "restarts", "restart_times",
        "restart_time_marker_modes", "terminal_restart", "terminal_restart_time",
        "restart_time_marker_bypass", "final_hardwall_projection_count",
    )
)
FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS = (
    CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS
    - {"plasma_continuation_policy", "plasma_continuation_evidence"}
)
FROZEN_E03_NO_MAX_NDIV_MIGRATIONS = {
    ("R03", "4762472"): {
        "segment": "s00_rankio_t0_t0p5",
        "ranks": 8,
        "manifest": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/manifest/prepared_run.json"
            ),
            "mode": 0o644,
            "size_bytes": 34030,
            "sha256": "ad3428b51945d65bed2b9ddeb337f4972bdc5b6ec07ac0cd6105ae6368055444",
        },
        "inspection": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/manifest/segment_inspection.json"
            ),
            "mode": 0o644,
            "size_bytes": 23551,
            "sha256": "24de23eaac2805eb7a0aa6e1c86ee6f1f133b573c9815548052fc618c5d3a643",
        },
        "mhd_history": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/output/"
                "E03_forcing_policy_paper_standard_active_alfvenic_beta100_"
                "s00_rankio_t0_t0p5.mhd.hst"
            ),
            "mode": 0o644,
            "size_bytes": 14145,
            "sha256": "7cb77ecf6d756b5a64c5bb02d917a5e75482df94b67f678cc34a68e05299f5af",
        },
        "user_history": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/output/"
                "E03_forcing_policy_paper_standard_active_alfvenic_beta100_"
                "s00_rankio_t0_t0p5.user.hst"
            ),
            "mode": 0o644,
            "size_bytes": 10621,
            "sha256": "766a79819ad14653042f578af6aa7f0ab8feba2abf5ca640fc4009a528001551",
        },
        "independent_validation": {
            "path": "accounting/4762472.stage_i.independent_validation.json",
            "mode": 0o644,
            "size_bytes": 24565,
            "sha256": "da15fc8a74013fd0f13276b3421ccc827ce0319c350cf265b4548297bfa1185e",
            "schema_version": 1,
            "record_type": "stage-i-binary-aware-clean-partial-validation",
        },
        "independent_review": None,
    },
    ("R12", "4766856"): {
        "segment": "s00_rankio_t0_t0p25",
        "ranks": 32,
        "manifest": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/manifest/prepared_run.json"
            ),
            "mode": 0o644,
            "size_bytes": 78976,
            "sha256": "ae1bc8256b3713dcf70ec99dde05a60d9d018650766e80ba29d49214c5e9bd19",
        },
        "inspection": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/manifest/segment_inspection.json"
            ),
            "mode": 0o644,
            "size_bytes": 66333,
            "sha256": "f901e51978a2e07d90728a10e8ca0f2c5abc393020a0c3ce8a4064ee5c914ab0",
        },
        "mhd_history": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/output/"
                "E03_forcing_policy_paper_heat_flux_beta10_strong_"
                "s00_rankio_t0_t0p25.mhd.hst"
            ),
            "mode": 0o644,
            "size_bytes": 6936,
            "sha256": "19f2fd19c003f644c9a53f0ba7bed5a754de46dd905d374eaf6cc65661e6d9b5",
        },
        "user_history": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/output/"
                "E03_forcing_policy_paper_heat_flux_beta10_strong_"
                "s00_rankio_t0_t0p25.user.hst"
            ),
            "mode": 0o644,
            "size_bytes": 5212,
            "sha256": "06507bd49e5b10f9db2c75f8858d70187b5ea7b0ce9a6e139db35aecd45d88aa",
        },
        "independent_validation": {
            "path": "accounting/4766856.stage_i.independent_validation.json",
            "mode": 0o444,
            "size_bytes": 278070,
            "sha256": "1d69d1b8f6cb35928f2334336399f3f7eee07aa511bf9ef36490670fe4da27e5",
            "schema_version": 3,
            "record_type": "stage-i-independent-segment-validation",
        },
        "independent_review": {
            "path": (
                "accounting/4766856.stage_i.independent_validation.json."
                "independent_review.json"
            ),
            "mode": 0o444,
            "size_bytes": 5628,
            "sha256": "45c20708e5619bcd46b3ec16f82f517c8acee415c1d2002ec0d43386c89d94f9",
        },
    },
}
FROZEN_E03_EXECUTABLE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
FROZEN_E03_EXECUTABLE_SHA256 = (
    "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
)
R17_MAX_NORMALIZED_CT_DIVB_TEXT = "1e-12"
R17_SCIENTIFIC_CHECKS = {
    "exact_endpoint": True,
    "complete_rank_inventory": True,
    "finite_synchronized_histories": True,
    "mass_conserved": True,
    "strict_lf_failure_counters_zero": True,
    "hard_volume_zero": True,
    "snapshot_hard_bounds_independently_verified": True,
    "normalized_ct_divb_below_threshold": True,
    "interval_cap_counts_valid": True,
    "nontrivial_forcing_and_pressure_work": True,
    "case_aware_policy_passed": True,
    "snapshot_cadence_complete": True,
    "terminal_snapshot_unique": True,
    "snapshots_structurally_complete": True,
    "restart_headers_authenticated": True,
    "terminal_restart_unique": True,
    "restart_load_smoke_passed": True,
}
INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION = (
    "Reviewer roles, agent identifiers, and process separation are retained "
    "declarations; exact artifact digests authenticate reviewed bytes but do not "
    "cryptographically authenticate a human or agent identity."
)
ACTIVE_SCHEDULER_STATES = frozenset(
    ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED")
)
STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
ACCOUNT_SACCT_FIELDS = (
    "JobIDRaw", "JobName", "State", "ExitCode", "NNodes", "ElapsedRaw",
    "Submit", "Start", "End", "Partition", "Account", "User",
)
ACCOUNT_SCHEDULER_HEADER = "|".join(ACCOUNT_SACCT_FIELDS)
ACCOUNT_MISSING_TIMESTAMPS = frozenset(("", "Unknown", "N/A", "None"))

R03_F115_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F115_source_bundle_recovery_supersession_evidence.json"
)
R03_F115_SHA256 = "cb50beb064678a9446ac33801a0023547d06c432d8c59b6bd3fbe34b11cf0391"
R03_F115_PUBLICATION_AUDIT_SHA256 = (
    "5923e3872b1d4a84a147bcd1d81fddc20ee79b72bfc410683d083039781f5b1f"
)
R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256 = (
    "6fcd19f9267f36332742f1f103968098216bd8e6b42fa9821964ab8df704bacd"
)
R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256 = (
    "a78357ed90e593809b1a82a641d2b40d651b16569ec943fc1781940e569b8440"
)
R03_F114_SEGMENT = "s01_rankio_t0p312823_t0p5"
R03_F115_SEGMENT = "s02_rankio_t0p312823_t0p5"
R03_F114_SOURCE_BUNDLE = (
    "source-archives/athenak-feature-cgl-through-c7e4fa30e.bundle"
)
R03_F114_SOURCE_BUNDLE_SHA256 = (
    "d94d559108470157f07981c9da8fec343a992128c0c4333df87181b81f0c505e"
)
R03_F114_CONTROLLER_REVISION = "c7e4fa30ea7162e4d5ce070a45dfec5b56ca2052"
R03_F115_CONTROLLER_REVISION = "5834a91e448a69ec0df5d011b7be3fe786666806"
R03_F115_CONTROLLER_SHA256 = (
    "0a6b30a60ba70faf32dae722c4472e96acb38031d8c1a4b52b816fc93e58e53a"
)
R03_F115_SOURCE_BUNDLE = "source-archives/athenak-feature-cgl-through-5834a91e4.bundle"
R03_F115_SOURCE_BUNDLE_SHA256 = (
    "a6aa40f8f3350d65022be6d898d5575a60185b05bba7b72f184c3f2ebd401ec8"
)
F115_SOURCE_BUNDLE_REQUIRED_REVISIONS = (
    SOURCE_REVISION,
    R03_F114_CONTROLLER_REVISION,
    "b0d3e8d526d8f3b4e5977333000db24c09d4eab1",
    "38aedd2a65c3c11855721858f5dbcce20bae11e4",
    R03_F115_CONTROLLER_REVISION,
)
F115_PROVENANCE_REVIEWER = {
    "agent_id": "019e9618-899e-7142-a748-ba57bf3efdd3",
    "identity": "OpenAI Codex independent read-only F115 provenance-security reviewer",
}
F115_PROVENANCE_REVIEW_SCOPE = (
    "Exact reviewed F115 candidate bytes and exact published F115 byte binding.",
    "All 19 declared F115 path, SHA-256, mode, and link bindings.",
    "Replacement source bundle identity, complete-history verification, strict object "
    "validation, active catalog integrity, and corrupt-C7 exclusion.",
    "Submitted-cancellation no-start evidence, immutable authorization packet, "
    "cancellation replay, retained cancelled s01 identity, and prohibition on job or "
    "segment reuse.",
    "Canonical reconciliation counts, zero-active-reservation state, zero-pending-"
    "transaction state, and qualification binding.",
    "F114-to-F115 source-binding and segment-identity supersession, sole-next-segment "
    "profile, prepare-only boundary, and non-broadening authorization.",
)
F115_PLASMA_REVIEWER = {
    "agent_id": "openai-codex-gpt-5",
    "identity": "OpenAI Codex, GPT-5 coding agent",
    "role": "independent read-only plasma/scientific continuation reviewer",
}
F115_PLASMA_AUTHORIZATION_LIMITATIONS = {
    "authorization_broadening": False,
    "limitations": [
        "Approval is limited to publication of the exact reviewed F-115 candidate bytes "
        "and preparation of only the sole authorized R03 s02 continuation after "
        "publication and clean reconciliation.",
        "Approval does not authorize reuse of cancelled job 4766485 or the retained s01 "
        "directory.",
        "Approval does not change the R03 parent, restart state, target, resources, "
        "executable, input, matrix, qualification, or F-113 bounded-concurrency and "
        "1400-node-hour policy.",
        "Submission must use the controller check-submit and submit paths; direct sbatch "
        "is not authorized.",
    ],
    "sole_authorized_case": "R03",
    "sole_authorized_segment": R03_F115_SEGMENT,
}
R17_READINESS_RELATIVE = Path(
    "accounting/mks24_stage_i_E03_forcing_policy_R17_readiness_evidence.json"
)
F116_CURRENT_SOURCE_AUTHORITY_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F116_current_source_authority_supersession_evidence.json"
)
F116_CANONICAL_SHA256 = {
    "evidence": "6cdbf9e4d10f1282744c6274aa3ef08afec4c510420296837fdbbdfcefe30a2a",
    "provenance_review": "e9731aab8305505e058c68ae8bb61c9ec5ff4885bde1bee3162719c41ab9bafd",
    "plasma_review": "bc4da6897263843f5d233f996539ff047d6ef465b775a9b5cdc8f864de96a6b8",
    "publication_audit": "3a6168e3039c02656b38ebdfcadffc07a2f151b1a474084307a80a341ba83096",
}
F118_CURRENT_SOURCE_AUTHORITY_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F118_current_source_authority_supersession_evidence.json"
)
F116_BRIDGE_REVISION = "36140ea825cb853b298714c27720440fdab60b9e"
F116_BRIDGE_SHA256 = "2c2f57a166877387244dd5bb6bdf87beb12492ea075a7431939b78e5df7307a0"
F116_BRIDGE_SOURCE_BUNDLE = "source-archives/athenak-feature-cgl-through-36140ea82.bundle"
F116_PRODUCTION_REQUIRED_REVISIONS = (
    SOURCE_REVISION,
    R03_F114_CONTROLLER_REVISION,
    "b0d3e8d526d8f3b4e5977333000db24c09d4eab1",
    "38aedd2a65c3c11855721858f5dbcce20bae11e4",
    R03_F115_CONTROLLER_REVISION,
    "469d38a841ef25d5c044713071adccd55ff90fef",
    "e1f4f4a0b62c3649b4d80a25a188b991f856ebbe",
    F116_BRIDGE_REVISION,
)
F116_REQUIRED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
F116_PUBLISHER_RELATIVE = "scripts/frontier/cgl_lf_stage_i_source_authority.py"
F116_AUTHORIZATION = {
    "current_source_selection_authorized": True,
    "source_authority_publication_authorized": True,
    "prepare_authorized": False,
    "submit_authorized": False,
    "direct_sbatch_authorized": False,
    "scheduler_mutation_authorized": False,
    "stage_i_execution_state_mutation_authorized": False,
    "scientific_configuration_change_authorized": False,
    "historical_manifest_rebinding_authorized": False,
}
F116_SCOPE_PRESERVES = [
    "The immutable F-115 evidence, reviews, publication audit, and historical R03 s02 authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F116_SCOPE_DOES_NOT_AUTHORIZE = [
    "prepare",
    "submit",
    "direct sbatch",
    "scheduler mutation",
    "Stage I execution-state mutation",
    "scientific configuration change",
    "historical manifest rebinding",
]
F116_PUBLICATION_REQUIREMENTS = {
    "published_evidence_mode": "0444",
    "published_review_mode": "0444",
    "published_audit_mode": "0444",
    "published_links": 1,
    "publication_audit_is_authority_commit_marker": True,
    "recovery_required_after_interruption": True,
}
F116_VALIDATION_CLAIMS = {
    "historical_f115_chain": "passed",
    "bridge_bundle_complete_history": "passed",
    "final_bundle_complete_history": "passed",
    "final_bundle_single_head_tip": "passed",
    "final_bundle_required_revisions": "passed",
    "committed_tool_bytes": "passed",
    "corrupt_c7_exclusion_preserved": True,
}
F116_PUBLICATION_METHOD = (
    "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
)
F118_SCOPE_PRESERVES = [
    "The immutable F-116 evidence, reviews, publication audit, selected source bundle, and nested F-115 historical authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F118_SCOPE_DOES_NOT_AUTHORIZE = F116_SCOPE_DOES_NOT_AUTHORIZE
F118_AUTHORIZATION = F116_AUTHORIZATION
F118_PUBLICATION_REQUIREMENTS = F116_PUBLICATION_REQUIREMENTS
F118_VALIDATION_CLAIMS = {
    "historical_f116_chain": "passed",
    "historical_f115_chain": "passed",
    "bridge_bundle_complete_history": "passed",
    "predecessor_current_source_bundle_complete_history": "passed",
    "final_bundle_complete_history": "passed",
    "final_bundle_single_head_tip": "passed",
    "final_bundle_required_revisions": "passed",
    "committed_tool_bytes": "passed",
    "corrupt_c7_exclusion_preserved": True,
}
F118_PUBLICATION_METHOD = F116_PUBLICATION_METHOD
F118_REQUIRED_TOOLS = F116_REQUIRED_TOOLS
F116_CORRUPT_C7_NAME = "athenak-feature-cgl-through-c7e4fa30e.bundle"

R03 = "R03"
R17 = "R17"
LOWER_CASES = tuple(f"R{index:02d}" for index in range(4, 17))
PROFILE_CASES = LOWER_CASES + (R17,)
ALL_CASES = tuple(f"R{index:02d}" for index in range(2, 18))
CONCURRENT_CASES = frozenset(f"R{index:02d}" for index in range(3, 17))
ACTIVE_CASES = frozenset(("R02", "R03", "R04", "R05", *tuple(f"R{i:02d}" for i in range(10, 18))))
PASSIVE_CASES = frozenset(("R06", "R07", "R08", "R09"))
FINITE_LIMITER_CASES = frozenset(("R14", "R15"))
HARDWALL_CASES = frozenset(ALL_CASES) - FINITE_LIMITER_CASES

INITIAL_TARGETS = {
    **{case_id: Decimal("0.25") for case_id in PROFILE_CASES},
    "R06": Decimal("0.50"),
    "R12": R12_FRESH_RERUN_TARGET,
    "R16": Decimal("1.50"),
}
MAX_REVIEWED_INCREMENTS = {
    **{case_id: Decimal("1.0") for case_id in LOWER_CASES},
    "R16": Decimal("3.0"),
    R17: Decimal("0.25"),
}
PRODUCTION_PRIORITY = (
    R03,
    "R04",
    "R12",
    "R16",
    "R06",
    "R10",
    "R11",
    "R13",
    "R05",
    "R07",
    "R14",
    "R15",
    "R08",
    "R09",
)
ROLLING_PRIORITY = PRODUCTION_PRIORITY[4:]
ALLOWED_NODES = {
    "R02": frozenset((1,)),
    R03: frozenset((1,)),
    **{case_id: frozenset((1, 2, 4)) for case_id in LOWER_CASES},
    "R16": frozenset((1, 2)),
    R17: frozenset((8,)),
}

EXPECTED_CASES = {
    "R02": ("paper_standard_active_alfvenic_beta10", "cgl_lf_paper_standard_active_alfvenic_beta10.athinput", "192x192x384", "c310509aa1638418bab7117427e8cca06e5a651c60e15212e4763baad5101206"),
    "R03": ("paper_standard_active_alfvenic_beta100", "cgl_lf_paper_standard_active_alfvenic_beta100.athinput", "192x192x384", "997f449abb3c2e4d1de50509efffa127c6230e3ecb2632612200010b8fa6c0b0"),
    "R04": ("paper_standard_active_random_beta10", "cgl_lf_paper_standard_active_random_beta10.athinput", "192x192x384", "571eea2ccec5d069ccb1b49d132ba4c5b8bdac7ce15ee8d8a04c92327c453137"),
    "R05": ("paper_standard_active_random_beta100", "cgl_lf_paper_standard_active_random_beta100.athinput", "192x192x384", "6527a2072105d515287701c904c3af2cebb07e10ab0b45c3fb41f29d3467622c"),
    "R06": ("paper_standard_passive_alfvenic_beta10", "cgl_lf_paper_standard_passive_alfvenic_beta10.athinput", "192x192x384", "c6e038c198b23bf20a7cf1e2a90fa82e83e5ec544492ddc64e1cf2bb535a38ae"),
    "R07": ("paper_standard_passive_alfvenic_beta100", "cgl_lf_paper_standard_passive_alfvenic_beta100.athinput", "192x192x384", "72ed6a3b38342d00855f34a7c7eb67102b44cd22b6c12ff2448ddb6f141e5145"),
    "R08": ("paper_standard_passive_random_beta10", "cgl_lf_paper_standard_passive_random_beta10.athinput", "192x192x384", "3535aca5e3d47dc383e353cff262d79c56cf3c083a3ac6024d04ddaa5d353780"),
    "R09": ("paper_standard_passive_random_beta100", "cgl_lf_paper_standard_passive_random_beta100.athinput", "192x192x384", "2ec05b247d917259c8306911cfd31d6be4e1c2772d6cdea544db883edfc35f2c"),
    "R10": ("paper_compressive_active_random_beta1", "cgl_lf_paper_compressive_active_random_beta1.athinput", "192x192x384", "93d7ed9b0846019f59bce34feeb9b22b098afffdfc2dbaa59aaf5913655cd705"),
    "R11": ("paper_compressive_active_random_beta100_sonic", "cgl_lf_paper_compressive_active_random_beta100_sonic.athinput", "192x192x384", "ab0dba23f80ea2d6c173ecbb724dbcb332b13d3b1751b8776ad05313f9ccaf6c"),
    "R12": ("paper_heat_flux_beta10_strong", "cgl_lf_paper_heat_flux_beta10_strong.athinput", "192x192x384", "98ddea4b4f7fec18cc40abdbf5f7c8ba5b583a91f84f4dae00e1411f23d42e7c"),
    "R13": ("paper_heat_flux_beta10_weak", "cgl_lf_paper_heat_flux_beta10_weak.athinput", "192x192x384", "a190129ee34c46a4fe83a392a19120a58b9a9a4064742ac58511441370c59c35"),
    "R14": ("paper_nulim_beta100_20", "cgl_lf_paper_nulim_beta100_20.athinput", "192x192x384", "2b8d5837f8a7f3070ca2ef56f8b44e53d048839918185a8eef155cb084736578"),
    "R15": ("paper_nulim_beta100_200", "cgl_lf_paper_nulim_beta100_200.athinput", "192x192x384", "9a698a60bef4c558ccee4635d69c3acbf3fea478401bd633d09bbc0526d943d8"),
    "R16": ("paper_scale_separation_active_alfvenic_beta10_nperp96", "cgl_lf_paper_scale_separation_beta10_nperp96.athinput", "96x96x192", "c0ac4b54248e8f8dfb0f5fd34c0cfb4414b5330529cbf2836961c5277af3f2d1"),
    "R17": ("paper_scale_separation_active_alfvenic_beta10_nperp384", "cgl_lf_paper_scale_separation_beta10_nperp384.athinput", "384x384x768", "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04"),
}

POLICY_TEXT = {
    "U": (
        "Require the exact endpoint; finite synchronized MHD and user histories; "
        "relative mass drift and MHD/user mass mismatch <= 1e-12; zero lf_dfloor, "
        "lf_pfloor, lf_nonfin, lf_nonpos, lf_hardbd, and hard_vol at every retained "
        "row; positive interval lf_nstage and lf_qface with bounded cap increments; "
        "finite LF heat and pressure-work ledgers; retain any emitted magnetic-"
        "divergence or CT diagnostics for independent scientific review without "
        "treating this planner as authentication of a CT metric or numerical threshold; "
        "complete ranked products; and one loadable exact terminal snapshot and restart "
        "group."
    ),
    "A": (
        "For active CGL, require segment and whole-case "
        "|Delta E - Delta force_work| / scale < 1e-8 and finite nonzero applied "
        "pressure work over the developed window."
    ),
    "P": (
        "For passive Delta, do not apply active-CGL total-energy closure; require "
        "lf_cpwrk == lf_cawrk == 0 throughout and retain finite forcing work."
    ),
    "H": (
        "For hardwall cases, require nonnegative monotonic lf_hwproj and review "
        "developed-turbulence activity; zero activity requires physics review rather "
        "than automatic acceptance."
    ),
    "F": (
        "For finite-limiter cases, require lf_hwproj == 0 and finite threshold "
        "occupancy and nu_eff that distinguish the configured limiter rate."
    ),
}

LEDGER_COLUMNS = (
    "execution_epoch", "job_id", "submitted_utc", "completed_utc", "case_id",
    "case_name", "segment", "state", "exit_code", "nodes", "requested_walltime",
    "elapsed_seconds", "reserved_node_hours", "actual_node_hours",
    "cumulative_stage_i_node_hours", "executable_revision", "executable_sha256",
    "input_revision", "input_file", "output_dir", "result", "notes",
)
RESERVATION_REQUIRED_COLUMNS = frozenset(
    (
        "execution_epoch", "manifest", "case_id", "case_name", "segment", "nodes",
        "requested_walltime", "reserved_node_hours", "state", "prepared_utc",
    )
)
RESERVATION_OPTIONAL_COLUMNS = frozenset(
    ("execution_intent_sha256", "job_id", "actual_node_hours", "result", "notes")
)
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
GIT_REVISION_RE = re.compile(r"^[0-9a-f]{40}$")
JOB_RE = re.compile(r"^[1-9][0-9]*$")
ACCOUNT_JOB_RE = re.compile(r"^[1-9][0-9]*(?:_[0-9]+|\+[0-9]+)?$")
DIGITS_RE = re.compile(r"^[0-9]+$")
RECOST_ARTIFACT_RE = re.compile(
    r"^mks24_stage_i_E03_forcing_policy_F(?P<number>[0-9]+)_recost_evidence\.json$"
)
SCHEDULER_NAIVE_RE = re.compile(r"^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}$")
SEGMENT_RE = re.compile(
    r"^s(?P<index>[0-9]+)_rankio_t(?P<start>[0-9]+(?:p[0-9]+)?)"
    r"_t(?P<target>[0-9]+(?:p[0-9]+)?)$"
)
WALLTIME_RE = re.compile(
    r"(?P<hours>[0-9]{2,}):(?P<minutes>[0-9]{2}):(?P<seconds>[0-9]{2})"
)
BATCH_SCRIPT_DIGEST_PLACEHOLDER = "0" * 64
BATCH_SCRIPT_DIGEST_PATTERN = re.compile(r"(?m)^BATCH_SCRIPT_SHA256=([0-9a-f]{64})$")
UTC_SUMMARY_RE = re.compile(r"^- Updated UTC: `([^`]+)`$", re.MULTILINE)
EPOCH_SUMMARY_RE = re.compile(r"^- Execution epoch: `([^`]+)`$", re.MULTILINE)
RECORDED_ROW_RE = re.compile(
    r"^\| `(?P<job>[0-9]+)` \| `(?P<case>R[0-9]{2})/(?P<segment>[^`/]+)` "
    r"\| (?P<state>[^|]+?) \| `(?P<hours>[0-9]+(?:\.[0-9]+)?)` "
    r"\| (?P<result>[a-z_]+) \|$"
)
ACTIVE_ROW_RE = re.compile(
    r"^\| `(?P<case>R[0-9]{2})/(?P<segment>[^`/]+)` "
    r"\| `(?P<nodes>[0-9]+)` \| `(?P<walltime>[0-9]{2,}:[0-9]{2}:[0-9]{2})` "
    r"\| `(?P<hours>[0-9]+(?:\.[0-9]+)?)` \| (?P<state>[a-z_]+) \|$"
)
SUMMARY_VALUE_PATTERNS = {
    "actual": re.compile(r"^- E03-forcing-policy Stage I actual use: `([^`]+)` node-hours$", re.MULTILINE),
    "reserved": re.compile(r"^- Active segment reservations: `([^`]+)` node-hours$", re.MULTILINE),
    "stage_remaining": re.compile(r"^- Unreserved E03 mapped-matrix remainder: `([^`]+)` node-hours$", re.MULTILINE),
    "project_remaining": re.compile(r"^- Incremental project remainder after active E03 Stage I use: `([^`]+)` node-hours$", re.MULTILINE),
}


def canonical_json(value: object) -> bytes:
    """Return stable bytes for one JSON-compatible value."""

    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("ascii")


def sha256_bytes(value: bytes) -> str:
    """Return a lowercase SHA-256 digest."""

    return hashlib.sha256(value).hexdigest()


def current_utc() -> datetime:
    """Return the current UTC time; isolated for deterministic tests."""

    return datetime.now(timezone.utc)


def duplicate_rejecting_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Build one JSON object while rejecting duplicate keys."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"JSON object contains duplicate key {key!r}")
        result[key] = value
    return result


def normalized_path(path: Path) -> Path:
    """Return an absolute lexical path without resolving symlinks."""

    return Path(os.path.abspath(os.fspath(path.expanduser())))


def component_identity(path: Path, label: str) -> tuple[tuple[str, int, int, int], ...]:
    """Return stable identities while rejecting every symlink path component."""

    current = Path(path.anchor)
    components = [current]
    for part in path.parts[1:]:
        current /= part
        components.append(current)
    result = []
    for component in components:
        try:
            profile = component.lstat()
        except OSError as error:
            raise ValueError(f"{label} path component is missing: {component}") from error
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symlink component: {component}")
        result.append((str(component), profile.st_dev, profile.st_ino, profile.st_mode))
    return tuple(result)


def read_stable_regular_file(path: Path, label: str) -> tuple[bytes, dict[str, object]]:
    """Read one stable, owner-controlled regular file without following symlinks."""

    path = normalized_path(path)
    components_before = component_identity(path, label)
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ValueError(f"{label} is not a readable non-symlink file: {path}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"{label} is not a regular file: {path}")
        if before.st_uid != os.geteuid() or before.st_nlink != 1:
            raise ValueError(f"{label} lacks the required owner/link profile: {path}")
        if stat.S_IMODE(before.st_mode) & 0o022:
            raise ValueError(f"{label} is group/world writable: {path}")
        blocks: list[bytes] = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        after_read = os.fstat(descriptor)
        try:
            final_path = path.lstat()
        except OSError as error:
            raise ValueError(f"{label} disappeared after it was read: {path}") from error
        if component_identity(path, label) != components_before:
            raise ValueError(f"{label} path components changed while it was read: {path}")
        final = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if any(
        getattr(snapshot, key) != getattr(before, key)
        for snapshot in (after_read, final_path, final)
        for key in stable_fields
    ):
        raise ValueError(f"{label} or its pathname changed while it was read: {path}")
    payload = b"".join(blocks)
    return payload, {
        "path": str(path),
        "sha256": sha256_bytes(payload),
        "size_bytes": len(payload),
        "mtime_ns": final.st_mtime_ns,
        "mode": f"{stat.S_IMODE(final.st_mode):04o}",
        "uid": final.st_uid,
        "nlink": final.st_nlink,
    }


def read_json_value(path: Path, label: str) -> tuple[object, dict[str, object]]:
    """Read stable unambiguous UTF-8 JSON."""

    payload, evidence = read_stable_regular_file(path, label)
    try:
        value = json.loads(payload.decode("utf-8"), object_pairs_hook=duplicate_rejecting_object)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"{label} is not unambiguous UTF-8 JSON: {path}") from error
    return value, evidence


def read_json_file(path: Path, label: str) -> tuple[dict[str, object], dict[str, object]]:
    """Read one stable JSON object."""

    value, evidence = read_json_value(path, label)
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain a JSON object")
    return value, evidence


def require_exact_keys(value: dict[str, object], keys: Iterable[str], label: str) -> None:
    """Require one exact object schema."""

    expected = set(keys)
    actual = set(value)
    if actual != expected:
        raise ValueError(
            f"{label} keys differ; missing={sorted(expected - actual)}, "
            f"extra={sorted(actual - expected)}"
        )


def require_nonempty_string(value: object, label: str) -> str:
    """Require one nonempty string."""

    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string")
    return value


def declared_process_independence_assurance(
    role_agents: dict[str, str], label: str
) -> dict[str, object]:
    """Represent strict declared process independence without identity overclaim."""

    if len(role_agents) < 2:
        raise ValueError(f"{label} must declare at least two distinct process roles")
    retained = {
        require_nonempty_string(role, f"{label} role"): require_nonempty_string(
            agent, f"{label} agent"
        )
        for role, agent in role_agents.items()
    }
    if len(retained) != len(role_agents) or len(set(retained.values())) != len(retained):
        raise ValueError(f"{label} roles or agents are not strictly distinct")
    return {
        "basis": "declared-process-independence",
        "strict_distinct_role_and_agent_declarations": True,
        "declared_role_agents": retained,
        "cryptographic_identity_verified": False,
        "non_cryptographic_limitation": INDEPENDENT_REVIEW_NON_CRYPTOGRAPHIC_LIMITATION,
    }


def require_git_revision(value: object, label: str) -> str:
    """Require one exact lowercase full Git object ID."""

    revision = require_nonempty_string(value, label)
    if GIT_REVISION_RE.fullmatch(revision) is None:
        raise ValueError(f"{label} must be a lowercase 40-hex Git revision")
    return revision


def require_sha256(value: object, label: str) -> str:
    """Require one exact lowercase SHA-256 digest."""

    digest = require_nonempty_string(value, label)
    if SHA256_RE.fullmatch(digest) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return digest


def require_utc(value: object, label: str) -> datetime:
    """Require one timezone-aware UTC timestamp."""

    text = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")
    return parsed


def require_scheduler_timestamp(value: object, label: str) -> datetime:
    """Parse controller-retained Slurm time, including its exact legacy naive form."""

    text = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(text.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        if SCHEDULER_NAIVE_RE.fullmatch(text) is None:
            raise ValueError(f"{label} has an ambiguous legacy scheduler timestamp")
        parsed = parsed.replace(tzinfo=datetime.now().astimezone().tzinfo)
    return parsed.astimezone(timezone.utc)


def require_fresh(
    evidence: dict[str, object],
    timestamp: object,
    label: str,
    max_age: timedelta = STATE_MAX_AGE,
) -> datetime:
    """Require timestamp and filesystem mtime to be current and mutually coherent."""

    observed = require_utc(timestamp, f"{label} timestamp")
    now = current_utc()
    if observed > now + FUTURE_SKEW or now - observed > max_age:
        raise ValueError(f"{label} timestamp is stale or implausibly future")
    mtime = datetime.fromtimestamp(int(evidence["mtime_ns"]) / 1_000_000_000, timezone.utc)
    if mtime > now + FUTURE_SKEW or now - mtime > max_age:
        raise ValueError(f"{label} filesystem evidence is stale or implausibly future")
    if abs(mtime - observed) > max_age:
        raise ValueError(f"{label} timestamp and filesystem mtime disagree")
    return observed


def decimal_value(value: object, label: str) -> Decimal:
    """Require one finite nonnegative exact decimal."""

    if isinstance(value, bool) or not isinstance(value, (str, int, float, Decimal)):
        raise ValueError(f"{label} must be an exact decimal")
    try:
        result = Decimal(str(value))
    except InvalidOperation as error:
        raise ValueError(f"{label} must be an exact decimal") from error
    if not result.is_finite() or result < 0:
        raise ValueError(f"{label} must be finite and nonnegative")
    return result


def recost_decimal_value(
    value: object, label: str, *, allow_negative: bool = False
) -> Decimal:
    """Require one finite non-exponent decimal string from a recost artifact."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be encoded as a decimal string")
    try:
        result = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"{label} must be encoded as a decimal string") from error
    if (
        not result.is_finite()
        or (result < 0 and not allow_negative)
        or format(result, "f") != value
    ):
        raise ValueError(f"{label} is not a valid recost decimal string")
    return result


def decimal_text(value: Decimal) -> str:
    """Return a plain canonical decimal string."""

    result = format(value.normalize(), "f")
    if "." in result:
        result = result.rstrip("0").rstrip(".")
    return result or "0"


def time_token(value: Decimal) -> str:
    """Return a Stage I segment-name time token."""

    return decimal_text(value).replace(".", "p")


def parse_time_token(value: str, label: str) -> Decimal:
    """Parse one Stage I segment-name time token."""

    return decimal_value(value.replace("p", "."), label)


def parse_segment(segment: str, label: str) -> tuple[int, Decimal, Decimal]:
    """Parse one exact rank-local Stage I segment identifier."""

    match = SEGMENT_RE.fullmatch(segment)
    if match is None:
        raise ValueError(f"{label} has unsupported segment identifier {segment!r}")
    start = parse_time_token(match.group("start"), f"{label} start")
    target = parse_time_token(match.group("target"), f"{label} target")
    if target <= start or target > Decimal("10"):
        raise ValueError(f"{label} has invalid exact interval")
    return int(match.group("index")), start, target


def resolution_cell_count(value: object, label: str) -> int:
    """Parse one exact three-dimensional matrix resolution."""

    retained = require_nonempty_string(value, label)
    match = re.fullmatch(r"([1-9][0-9]*)x([1-9][0-9]*)x([1-9][0-9]*)", retained)
    if match is None:
        raise ValueError(f"{label} must use NxNxN positive-integer dimensions")
    cells = math.prod(int(item) for item in match.groups())
    if cells <= 0:
        raise ValueError(f"{label} cell count is invalid")
    return cells


def walltime_seconds(value: object, label: str) -> int:
    """Require one positive normalized walltime."""

    text = require_nonempty_string(value, label)
    match = WALLTIME_RE.fullmatch(text)
    if match is None:
        raise ValueError(f"{label} must use HH:MM:SS")
    minutes = int(match.group("minutes"))
    seconds = int(match.group("seconds"))
    if minutes >= 60 or seconds >= 60:
        raise ValueError(f"{label} has invalid minutes or seconds")
    total = int(match.group("hours")) * 3600 + minutes * 60 + seconds
    if total <= 0:
        raise ValueError(f"{label} must be positive")
    return total


def expected_path(relative: str) -> Path:
    """Return one exact canonical-root path."""

    return normalized_path(CANONICAL_ROOT / relative)


def frozen_matrix_path() -> Path:
    """Return the exact frozen matrix path."""

    return normalized_path(FROZEN_SOURCE / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json")


def input_path(case_id: str) -> Path:
    """Return one exact frozen source-input path."""

    return normalized_path(FROZEN_SOURCE / "inputs/cgl_lf_paper" / EXPECTED_CASES[case_id][1])


def ledger_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv")


def reservations_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json")


def qualification_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_qualification_approval.json")


def transaction_store_paths() -> tuple[Path, Path]:
    """Return exact Stage I and recost transaction-store paths."""

    accounting = expected_path("accounting")
    return (
        accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions",
        accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions",
    )


def r03_f115_path() -> Path:
    """Return the exact published F-115 R03 authority path."""

    return expected_path(R03_F115_RELATIVE.as_posix())


def r03_f115_publication_audit_path() -> Path:
    """Return the exact published F-115 publication-audit path."""

    return expected_path(f"{R03_F115_RELATIVE.as_posix()}.publication_audit.json")


def r03_f115_provenance_security_review_path() -> Path:
    """Return the exact published F-115 provenance/security review path."""

    return expected_path(f"{R03_F115_RELATIVE.as_posix()}.provenance_security_review.json")


def r03_f115_plasma_scientific_review_path() -> Path:
    """Return the exact published F-115 plasma/scientific review path."""

    return expected_path(f"{R03_F115_RELATIVE.as_posix()}.plasma_scientific_review.json")


def f116_current_source_authority_path() -> Path:
    """Return the exact published F-116 current-source-authority path."""

    return expected_path(F116_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix())


def f116_provenance_security_review_path() -> Path:
    """Return the exact published F-116 provenance/security review path."""

    return expected_path(
        f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.provenance_security_review.json"
    )


def f116_plasma_scientific_review_path() -> Path:
    """Return the exact published F-116 plasma/scientific review path."""

    return expected_path(
        f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.plasma_scientific_review.json"
    )


def f116_publication_audit_path() -> Path:
    """Return the exact published F-116 publication-audit path."""

    return expected_path(
        f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.publication_audit.json"
    )


def f118_current_source_authority_path() -> Path:
    """Return the exact published F-118 current-source-authority path."""

    return expected_path(F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix())


def f118_provenance_security_review_path() -> Path:
    """Return the exact published F-118 provenance/security review path."""

    return expected_path(
        f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.provenance_security_review.json"
    )


def f118_plasma_scientific_review_path() -> Path:
    """Return the exact published F-118 plasma/scientific review path."""

    return expected_path(
        f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.plasma_scientific_review.json"
    )


def f118_publication_audit_path() -> Path:
    """Return the exact published F-118 publication-audit path."""

    return expected_path(
        f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix()}.publication_audit.json"
    )


def promoted_source_bundle_path(controller_revision: str) -> Path:
    """Return the exact canonical bundle path for one reviewed promoted revision."""

    revision = require_git_revision(controller_revision, "promoted controller revision")
    return expected_path(f"source-archives/athenak-feature-cgl-through-{revision[:9]}.bundle")


def f115_source_bundle_path() -> Path:
    """Return the exact immutable historical R03 F-115 source-bundle path."""

    return expected_path(R03_F115_SOURCE_BUNDLE)


def recost_independent_review_path(recost_path: Path) -> Path:
    """Return one schema-2 recost artifact's exact independent-review path."""

    return recost_path.with_name(f"{recost_path.name}.independent_review.json")


def recost_publication_audit_path(recost_path: Path) -> Path:
    """Return one schema-2 recost artifact's exact publication-audit path."""

    return recost_path.with_name(f"{recost_path.name}.publication_audit.json")


def r17_readiness_path() -> Path:
    """Return the sole recost-compatible R17 readiness artifact path."""

    return expected_path(R17_READINESS_RELATIVE.as_posix())


def controller_repository_path() -> Path:
    """Return the exact repository containing the promoted controller helper."""

    return normalized_path(CONTROLLER_HELPER.parents[2])


def require_trusted_system_executable(path: Path, label: str) -> str:
    """Require one absolute root-owned, single-link, non-writable executable."""

    path = normalized_path(path)
    component_identity(path, label)
    profile = path.stat()
    mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != 0
        or profile.st_nlink != 1
        or mode & 0o022
        or mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted system executable profile")
    return str(path)


def hardened_child_environment() -> dict[str, str]:
    """Return a scheduler child environment without caller execution controls."""

    return {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def hardened_git_environment() -> dict[str, str]:
    """Return Git's complete caller-independent execution environment."""

    return {
        "GIT_CONFIG_GLOBAL": os.devnull,
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_CONFIG_SYSTEM": os.devnull,
        "GIT_EXEC_PATH": str(GIT_EXEC_PATH),
        "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0",
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def git_read_only(
    repository: Path,
    arguments: list[str],
    *,
    descriptor: int | None = None,
) -> subprocess.CompletedProcess:
    """Run one bounded, environment-sanitized, non-shell Git query."""

    git = require_trusted_system_executable(GIT, "Git executable")
    try:
        return subprocess.run(
            [
                git, "--no-replace-objects", "-C", str(repository), *arguments
            ],
            check=False,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            env=hardened_git_environment(),
            pass_fds=(() if descriptor is None else (descriptor,)),
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError("read-only Git source-bundle validation failed to execute") from error


def validate_source_bundle_coverage(
    path: Path,
    *,
    bundle_sha256: str,
    controller_revision: str,
    controller_sha256: str,
    required_revisions: Iterable[str],
    expected_evidence: dict[str, object] | None = None,
    allow_branch_ref: bool,
    label: str,
) -> dict[str, object]:
    """Independently require one exact usable self-contained controller bundle."""

    path = normalized_path(path)
    bundle_sha256 = require_sha256(bundle_sha256, f"{label} SHA-256")
    controller_revision = require_git_revision(
        controller_revision, f"{label} controller revision"
    )
    controller_sha256 = require_sha256(controller_sha256, f"{label} controller SHA-256")
    revisions = tuple(
        require_git_revision(revision, f"{label} required revision")
        for revision in required_revisions
    )
    if not revisions or len(set(revisions)) != len(revisions):
        raise ValueError(f"{label} required revisions are empty or duplicated")

    payload, evidence = read_stable_regular_file(path, label)
    if evidence["sha256"] != bundle_sha256:
        raise ValueError(f"{label} digest changed")
    if expected_evidence is not None and evidence != expected_evidence:
        raise ValueError(f"{label} changed before Git validation")
    try:
        header = payload.split(b"\n\n", 1)[0].decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} header is not UTF-8") from error
    if not header or header[0] not in {"# v2 git bundle", "# v3 git bundle"}:
        raise ValueError(f"{label} is not a Git bundle")
    if any(line.startswith("-") for line in header[1:]):
        raise ValueError(f"{label} is not self-contained")
    advertised_header = []
    for line in header[1:]:
        if line.startswith("@"):
            continue
        match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
        if match is None:
            raise ValueError(f"{label} advertised reference is malformed")
        advertised_header.append((match.group(1), match.group(2)))
    if len(advertised_header) != 1:
        raise ValueError(f"{label} must advertise exactly one tip")
    advertised_revision, advertised_name = advertised_header[0]
    branch_ref = (
        advertised_name.startswith("refs/heads/")
        and len(advertised_name) > len("refs/heads/")
        and not any(character.isspace() for character in advertised_name)
    )
    if (
        advertised_revision != controller_revision
        or not (advertised_name == "HEAD" or (allow_branch_ref and branch_ref))
    ):
        raise ValueError(f"{label} does not advertise the exact promoted controller tip")

    components_before = component_identity(path, f"{label} Git validation")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        if sha256_bytes(b"".join(blocks)) != bundle_sha256:
            raise ValueError(f"{label} descriptor digest changed")
        os.lseek(descriptor, 0, os.SEEK_SET)
        descriptor_path = f"/proc/self/fd/{descriptor}"
        repository = controller_repository_path()
        verify = git_read_only(
            repository, ["bundle", "verify", descriptor_path], descriptor=descriptor
        )
        if verify.returncode:
            raise ValueError(f"{label} fails git bundle verify")
        heads = git_read_only(
            repository, ["bundle", "list-heads", descriptor_path], descriptor=descriptor
        )
        if heads.returncode:
            raise ValueError(f"{label} heads cannot be listed")
        try:
            advertised = [
                tuple(line.split(" ", 1))
                for line in heads.stdout.decode("utf-8").splitlines()
                if line
            ]
        except UnicodeDecodeError as error:
            raise ValueError(f"{label} heads are not UTF-8") from error
        if advertised != advertised_header:
            raise ValueError(f"{label} header and Git-advertised heads differ")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if any(getattr(before, field) != getattr(after, field) for field in stable_fields):
        raise ValueError(f"{label} changed during Git validation")
    if component_identity(path, f"{label} Git validation") != components_before:
        raise ValueError(f"{label} path changed during Git validation")
    _, retained = read_stable_regular_file(path, label)
    if retained != evidence:
        raise ValueError(f"{label} changed after Git validation")

    repository = controller_repository_path()
    for revision in revisions:
        if git_read_only(
            repository, ["cat-file", "-e", f"{revision}^{{commit}}"]
        ).returncode:
            raise ValueError(f"{label} required revision is unknown: {revision}")
        if git_read_only(
            repository,
            ["merge-base", "--is-ancestor", revision, controller_revision],
        ).returncode:
            raise ValueError(f"{label} required revision is not covered by its tip: {revision}")
    try:
        helper_relative = CONTROLLER_HELPER.relative_to(repository).as_posix()
    except ValueError as error:
        raise ValueError(f"{label} controller helper is outside its repository") from error
    committed_helper = git_read_only(
        repository, ["show", f"{controller_revision}:{helper_relative}"]
    )
    if committed_helper.returncode or sha256_bytes(committed_helper.stdout) != controller_sha256:
        raise ValueError(f"{label} tip does not bind the declared controller helper bytes")
    return {
        "evidence": evidence,
        "git_bundle_verify": "passed",
        "self_contained": True,
        "advertised_heads": [
            {"revision": revision, "name": name} for revision, name in advertised
        ],
        "required_revisions": list(revisions),
        "committed_controller_helper_sha256": controller_sha256,
    }


def validate_f115_source_bundle_coverage(
    expected_evidence: dict[str, object] | None = None,
) -> dict[str, object]:
    """Require the immutable historical R03 F-115 bundle and its exact HEAD."""

    return validate_source_bundle_coverage(
        f115_source_bundle_path(),
        bundle_sha256=R03_F115_SOURCE_BUNDLE_SHA256,
        controller_revision=R03_F115_CONTROLLER_REVISION,
        controller_sha256=R03_F115_CONTROLLER_SHA256,
        required_revisions=F115_SOURCE_BUNDLE_REQUIRED_REVISIONS,
        expected_evidence=expected_evidence,
        allow_branch_ref=False,
        label="R03 F-115 source bundle",
    )


def executable_path() -> Path:
    return expected_path("build/frontier-hip-9e07542281e4-cpe25.09-cce20-rocm6.4.2/src/athena")


def build_manifest_path() -> Path:
    return expected_path("runs/build-manifests/9e07542281e4-cpe25.09-cce20-rocm6.4.2")


def validate_binding(
    binding: object,
    evidence: dict[str, object],
    expected: Path,
    label: str,
) -> None:
    """Require an exact file binding at an exact trusted path."""

    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding must be an object")
    require_exact_keys(
        binding, ("path", "sha256", "size_bytes", "mtime_ns", "mode", "uid", "nlink"),
        f"{label} binding",
    )
    if normalized_path(Path(str(binding["path"]))) != normalized_path(expected):
        raise ValueError(f"{label} path is not the exact trusted path")
    if binding != evidence:
        raise ValueError(f"{label} evidence changed after review")


def directory_snapshot(path: Path, label: str) -> dict[str, object]:
    """Return one stable owner-controlled directory inventory."""

    path = normalized_path(path)
    components_before = component_identity(path, label)
    try:
        before = path.stat()
    except OSError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if not stat.S_ISDIR(before.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")
    if before.st_uid != os.geteuid() or stat.S_IMODE(before.st_mode) & 0o022:
        raise ValueError(f"{label} lacks the required owner/write profile: {path}")
    try:
        entries = sorted(item.name for item in os.scandir(path))
    except OSError as error:
        raise ValueError(f"{label} cannot be inventoried: {path}") from error
    after = path.stat()
    if component_identity(path, label) != components_before:
        raise ValueError(f"{label} path changed while inventoried: {path}")
    stable_fields = ("st_dev", "st_ino", "st_mode", "st_uid", "st_mtime_ns", "st_ctime_ns")
    if any(getattr(before, key) != getattr(after, key) for key in stable_fields):
        raise ValueError(f"{label} changed while inventoried: {path}")
    return {
        "path": str(path),
        "device": after.st_dev,
        "inode": after.st_ino,
        "mode": f"{stat.S_IMODE(after.st_mode):04o}",
        "uid": after.st_uid,
        "mtime_ns": after.st_mtime_ns,
        "ctime_ns": after.st_ctime_ns,
        "entries": entries,
    }


def discover_transaction_state() -> dict[str, dict[str, object]]:
    """Discover both transaction stores and reject any in-flight journal."""

    labels = ("stage_i", "recost")
    result = {
        label: directory_snapshot(path, f"{label} transaction store")
        for label, path in zip(labels, transaction_store_paths())
    }
    pending = [
        f"{label}:{entry}"
        for label, snapshot in result.items()
        for entry in snapshot["entries"]
    ]
    if pending:
        raise ValueError(f"canonical transaction stores are not empty: {', '.join(pending)}")
    return result


def collect_file_evidence_records(value: object) -> dict[str, dict[str, object]]:
    """Collect every exact retained-file evidence object from a nested value."""

    evidence_keys = {"path", "sha256", "size_bytes", "mtime_ns", "mode", "uid", "nlink"}
    result: dict[str, dict[str, object]] = {}

    def visit(item: object) -> None:
        if isinstance(item, dict):
            if set(item) == evidence_keys:
                path = str(item["path"])
                prior = result.get(path)
                if prior is not None and prior != item:
                    raise ValueError(f"global evidence contains conflicting bindings for {path}")
                result[path] = dict(item)
            for child in item.values():
                visit(child)
        elif isinstance(item, list):
            for child in item:
                visit(child)

    visit(value)
    return result


def revalidate_global_boundary(
    value: object,
    transaction_state: dict[str, dict[str, object]],
    submitted_scheduler: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Re-read every input and live scheduler fact before emitting a plan."""

    records = collect_file_evidence_records(value)
    for path, expected in sorted(records.items()):
        _, current = read_stable_regular_file(Path(path), "global TOCTOU evidence")
        if current != expected:
            raise ValueError(f"global evidence changed before plan completion: {path}")
    current_transactions = discover_transaction_state()
    if current_transactions != transaction_state:
        raise ValueError("transaction stores changed before plan completion")
    for job_id, expected in sorted(submitted_scheduler.items()):
        scheduler = expected["scheduler"]
        current = authenticate_live_scheduler_job(
            job_id, int(scheduler["nodes"]), str(scheduler["job_name"])
        )
        if current != scheduler:
            raise ValueError(f"live Slurm state changed before plan completion for job {job_id}")
        if sha256_bytes(collect_live_batch_script(job_id)) != expected["batch_script_sha256"]:
            raise ValueError(f"Slurm-stored batch script changed before plan completion for job {job_id}")
    return {
        "revalidated_file_count": len(records),
        "transaction_stores": current_transactions,
        "submitted_scheduler_jobs": submitted_scheduler,
    }


def policy_codes(case_id: str) -> tuple[str, ...]:
    """Return the exact case-aware physics policy codes."""

    if case_id in PASSIVE_CASES:
        return ("U", "P", "H")
    if case_id in FINITE_LIMITER_CASES:
        return ("U", "A", "F")
    if case_id in ACTIVE_CASES and case_id in HARDWALL_CASES:
        return ("U", "A", "H")
    raise ValueError(f"case {case_id} lacks an exact physics policy")


def model_role(case_id: str) -> str:
    """Return the scientific role of one mapped case."""

    if case_id in PASSIVE_CASES:
        return "isothermal_mhd_with_passive_cgl_delta_landau_fluid"
    if case_id in FINITE_LIMITER_CASES:
        return "active_cgl_mhd_landau_fluid_with_finite_limiter"
    return "active_cgl_mhd_landau_fluid_with_hardwall"


def validate_matrix(
    matrix: dict[str, object], evidence: dict[str, object]
) -> tuple[dict[str, dict[str, object]], dict[str, dict[str, object]]]:
    """Require the exact frozen matrix semantics and every exact input record."""

    if normalized_path(Path(str(evidence["path"]))) != frozen_matrix_path():
        raise ValueError("matrix path is not the exact frozen matrix path")
    if evidence["sha256"] != MATRIX_SHA256:
        raise ValueError("matrix digest differs from the known frozen matrix")
    if matrix.get("schema_version") != 1 or matrix.get("campaign") != CAMPAIGN:
        raise ValueError("matrix identity differs from the frozen Stage I matrix")
    cases = matrix.get("cases")
    if not isinstance(cases, list):
        raise ValueError("matrix cases must be a list")
    indexed: dict[str, dict[str, object]] = {}
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("matrix case must be an object")
        case_id = require_nonempty_string(case.get("id"), "matrix case id")
        if case_id in indexed:
            raise ValueError(f"matrix contains duplicate case {case_id}")
        indexed[case_id] = case
    if set(indexed) != set(ALL_CASES):
        raise ValueError("matrix case IDs differ from the exact R02-R17 mapping")
    input_evidence: dict[str, dict[str, object]] = {}
    for case_id, (name, filename, resolution, digest) in EXPECTED_CASES.items():
        case = indexed[case_id]
        expected_input = f"inputs/cgl_lf_paper/{filename}"
        estimate = decimal_value(case.get("estimated_node_hours"), f"matrix {case_id} estimate")
        if (
            case.get("name") != name
            or case.get("input") != expected_input
            or case.get("resolution") != resolution
            or estimate <= 0
        ):
            raise ValueError(f"matrix {case_id} semantics differ from the frozen mapping")
        case["_estimated_node_hours"] = estimate
        _, retained = read_stable_regular_file(input_path(case_id), f"{case_id} frozen input")
        if retained["sha256"] != digest:
            raise ValueError(f"{case_id} frozen input digest differs from the known record")
        input_evidence[case_id] = retained
    return indexed, input_evidence


def validate_static_provenance(
    recost: dict[str, object],
    matrix_evidence: dict[str, object],
    input_evidence: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Authenticate current promoted identity from the reviewed schema-2 recost."""

    provenance = recost.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError("recost provenance must be an object")
    required = {
        "generator_sha256",
        "generator_revision",
        "stage_i_helper_sha256",
        "stage_i_helper_revision",
        "matrix_sha256",
        "matrix_revision",
        "source_bundle_sha256",
        "source_bundle_verified_revisions",
        "source_authority",
        "qualification_approval_sha256",
    }
    if not required.issubset(provenance):
        raise ValueError("recost provenance omits current promoted identity")
    current_source = validate_f118_current_source_authority(provenance["source_authority"])
    controller_helper = current_source["controller_helper"]
    generator = current_source["generator"]
    matrix_identity = current_source["matrix"]
    source_bundle = current_source["source_bundle"]
    controller_revision = controller_helper["revision"]
    controller_sha256 = controller_helper["sha256"]
    generator_revision = generator["revision"]
    matrix_revision = matrix_identity["revision"]
    source_bundle_sha256 = source_bundle["sha256"]
    revisions = tuple(source_bundle["verified_revisions"])
    if provenance["matrix_sha256"] != MATRIX_SHA256 or matrix_evidence["sha256"] != MATRIX_SHA256:
        raise ValueError("recost matrix identity differs from the frozen Stage I matrix")
    if (
        provenance["stage_i_helper_revision"] != controller_revision
        or provenance["stage_i_helper_sha256"] != controller_sha256
        or provenance["generator_revision"] != generator_revision
        or provenance["generator_sha256"] != generator["sha256"]
        or provenance["matrix_revision"] != matrix_revision
        or provenance["matrix_sha256"] != matrix_identity["sha256"]
        or provenance["source_bundle_sha256"] != source_bundle_sha256
        or provenance["source_bundle_verified_revisions"] != list(revisions)
    ):
        raise ValueError("recost current-source identity differs from F118 authority")

    repository = controller_repository_path()
    matrix_relative = "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    committed_matrix = git_read_only(repository, ["show", f"{matrix_revision}:{matrix_relative}"])
    if committed_matrix.returncode or sha256_bytes(committed_matrix.stdout) != MATRIX_SHA256:
        raise ValueError("recost matrix revision does not bind the frozen matrix bytes")

    _, executable_evidence = read_stable_regular_file(executable_path(), "qualified executable")
    if executable_evidence["sha256"] != EXECUTABLE_SHA256:
        raise ValueError("qualified executable digest changed")
    build_evidence: dict[str, dict[str, object]] = {}
    build_payloads: dict[str, bytes] = {}
    for filename, digest in BUILD_FILE_SHA256.items():
        path = build_manifest_path() / filename
        payload, retained = read_stable_regular_file(path, f"build manifest {filename}")
        if retained["sha256"] != digest:
            raise ValueError(f"build manifest {filename} digest changed")
        build_evidence[filename] = retained
        build_payloads[filename] = payload
    if build_payloads["athena.sha256"].decode("utf-8").split()[0] != EXECUTABLE_SHA256:
        raise ValueError("build manifest does not bind the qualified executable")
    if re.search(
        rf"(?m)^git_revision={re.escape(SOURCE_REVISION)}\s*$",
        build_payloads["environment.txt"].decode("utf-8"),
    ) is None:
        raise ValueError("build manifest does not bind the qualified source revision")

    qualification, qualification_evidence = read_json_file(
        qualification_path(), "qualification approval"
    )
    if (
        qualification_evidence["sha256"] != provenance["qualification_approval_sha256"]
        or qualification.get("schema_version") != 1
        or qualification.get("execution_epoch") != EXECUTION_EPOCH
        or qualification.get("approved_executable") != str(executable_path())
        or qualification.get("approved_executable_revision") != SOURCE_REVISION
        or qualification.get("approved_executable_sha256") != EXECUTABLE_SHA256
        or qualification.get("build_manifest") != str(build_manifest_path())
    ):
        raise ValueError("qualification approval differs from the promoted build")

    return {
        "controller_helper": {
            "path": str(CONTROLLER_HELPER),
            "revision": controller_revision,
            "sha256": controller_sha256,
            "evidence": controller_helper["evidence"],
        },
        "source_revision": SOURCE_REVISION,
        "source_bundle": {
            "path": source_bundle["path"],
            "sha256": source_bundle_sha256,
            "verified_revisions": list(revisions),
            "evidence": source_bundle["evidence"],
            "independent_validation": source_bundle["independent_validation"],
        },
        "current_source_authority": provenance["source_authority"],
        "current_source_authority_evidence": current_source["authority_evidence"],
        "current_source_authority_review_assurance": current_source[
            "independent_review_assurance"
        ],
        "matrix": matrix_evidence,
        "executable": {
            "path": str(executable_path()),
            "revision": SOURCE_REVISION,
            "sha256": EXECUTABLE_SHA256,
            "evidence": executable_evidence,
        },
        "build_manifest": {
            "path": str(build_manifest_path()),
            "files": build_evidence,
        },
        "qualification": qualification_evidence,
        "inputs": {
            case_id: {
                "path": str(input_path(case_id)),
                "revision": SOURCE_REVISION,
                "sha256": EXPECTED_CASES[case_id][3],
                "evidence": input_evidence[case_id],
            }
            for case_id in ALL_CASES
        },
    }


def summary_section(text: str, heading: str, next_heading: str | None) -> str:
    """Return one unique Markdown section body."""

    marker = f"## {heading}\n"
    if text.count(marker) != 1:
        raise ValueError(f"controller summary must contain exactly one {heading} section")
    body = text.split(marker, 1)[1]
    if next_heading is not None:
        next_marker = f"## {next_heading}\n"
        if body.count(next_marker) != 1:
            raise ValueError(f"controller summary must contain exactly one {next_heading} section")
        body = body.split(next_marker, 1)[0]
    return body


def parse_table_rows(section: str, pattern: re.Pattern[str], label: str) -> list[dict[str, str]]:
    """Parse all non-header Markdown rows from one summary section."""

    rows = []
    for line in section.splitlines():
        if not line.startswith("|"):
            continue
        if line.startswith("| ---") or line.startswith("| Job ") or line.startswith("| Case/segment "):
            continue
        match = pattern.fullmatch(line)
        if match is None:
            raise ValueError(f"controller summary contains malformed {label} row: {line}")
        rows.append(match.groupdict())
    return rows


def unique_match(pattern: re.Pattern[str], text: str, label: str) -> str:
    """Return one unique regular-expression capture."""

    matches = pattern.findall(text)
    if len(matches) != 1:
        raise ValueError(f"controller summary must contain exactly one {label}")
    return matches[0]


def parse_summary(payload: bytes) -> dict[str, object]:
    """Parse the exact controller budget summary subset used by the planner."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("controller summary must be UTF-8") from error
    updated_utc = unique_match(UTC_SUMMARY_RE, text, "Updated UTC")
    if unique_match(EPOCH_SUMMARY_RE, text, "execution epoch") != EXECUTION_EPOCH:
        raise ValueError("controller summary execution epoch is invalid")
    if (
        "Up to 4 distinct R03-R16 cases may be active concurrently" not in text
        or "R17 remains exclusive and last" not in text
        or "at most one prepared submission packet" not in text
    ):
        raise ValueError("controller summary lacks the promoted wave policy")
    recorded = parse_table_rows(
        summary_section(text, "Recorded Segments", "Active Reservations"),
        RECORDED_ROW_RE, "recorded-segment",
    )
    active = parse_table_rows(
        summary_section(text, "Active Reservations", None), ACTIVE_ROW_RE, "active-reservation"
    )
    values = {
        key: decimal_value(unique_match(pattern, text, key), f"summary {key}")
        for key, pattern in SUMMARY_VALUE_PATTERNS.items()
    }
    return {"updated_utc": updated_utc, "recorded": recorded, "active": active, **values}


def read_ledger(path: Path) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Read the exact controller ledger schema."""

    payload, evidence = read_stable_regular_file(path, "canonical ledger")
    try:
        reader = csv.DictReader(io.StringIO(payload.decode("utf-8")))
        if tuple(reader.fieldnames or ()) != LEDGER_COLUMNS:
            raise ValueError("canonical ledger columns differ from the controller schema")
        rows = list(reader)
    except (UnicodeDecodeError, csv.Error) as error:
        raise ValueError("canonical ledger is not valid UTF-8 CSV") from error
    return rows, evidence


def validate_resource_values(
    case_id: str,
    nodes_value: object,
    walltime_value: object,
    athena_walltime_value: object,
    reserved_value: object,
    label: str,
) -> tuple[int, int, int, Decimal]:
    """Require exact controller node, walltime, margin, and reservation policy."""

    if isinstance(nodes_value, bool) or not isinstance(nodes_value, int):
        try:
            nodes = int(str(nodes_value))
        except ValueError as error:
            raise ValueError(f"{label} nodes are invalid") from error
    else:
        nodes = nodes_value
    if nodes not in ALLOWED_NODES[case_id]:
        raise ValueError(f"{label} nodes are outside the controller node policy")
    slurm_seconds = walltime_seconds(walltime_value, f"{label} Slurm walltime")
    athena_seconds = walltime_seconds(athena_walltime_value, f"{label} Athena walltime")
    if slurm_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError(f"{label} Slurm walltime exceeds the controller maximum")
    if athena_seconds > slurm_seconds - SHUTDOWN_MARGIN_SECONDS:
        raise ValueError(f"{label} Athena walltime lacks the controller shutdown margin")
    reserved = decimal_value(reserved_value, f"{label} reserved node-hours")
    expected_reserved = Decimal(nodes * slurm_seconds) / Decimal(3600)
    if abs(reserved - expected_reserved) > Decimal("0.000001"):
        raise ValueError(f"{label} reserved node-hours differ from the allocation")
    return nodes, slurm_seconds, athena_seconds, expected_reserved


def expected_job_name(case_id: str, segment: str) -> str:
    """Return the exact production Slurm job name."""

    value = f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_{case_id}_{segment}"
    value = re.sub(r"[^A-Za-z0-9_]+", "_", value)[:60]
    if re.fullmatch(r"[A-Za-z0-9_]+", value) is None:
        raise ValueError("generated Stage I job name is unsafe")
    return value


def scheduler_owner() -> str:
    """Return the effective local identity Slurm must report."""

    return pwd.getpwuid(os.geteuid()).pw_name


def collect_live_scheduler_job(job_id: str) -> dict[str, object]:
    """Read exactly one active top-level Slurm allocation."""

    if JOB_RE.fullmatch(job_id) is None:
        raise ValueError("live scheduler query has invalid job ID")
    squeue = require_trusted_system_executable(SQUEUE, "squeue executable")
    try:
        completed = subprocess.run(
            [
                squeue, "-j", job_id, "-h",
                "-t", ",".join(sorted(ACTIVE_SCHEDULER_STATES)),
                "-o", "%i|%j|%T|%D|%u|%a",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=hardened_child_environment(),
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"live squeue query failed for job {job_id}") from error
    lines = [line for line in completed.stdout.splitlines() if line]
    if len(lines) != 1:
        raise ValueError(f"live squeue must return exactly one active row for job {job_id}")
    fields = lines[0].split("|")
    if len(fields) != 6:
        raise ValueError("live squeue row has invalid columns")
    observed_job, job_name, state, nodes_text, owner, account = fields
    try:
        nodes = int(nodes_text)
    except ValueError as error:
        raise ValueError("live squeue node count is invalid") from error
    if (
        observed_job != job_id
        or state not in ACTIVE_SCHEDULER_STATES
        or nodes <= 0
        or owner != scheduler_owner()
        or account.casefold() != ACCOUNT.casefold()
    ):
        raise ValueError("live squeue row has invalid active-job identity")
    return {
        "job_id": observed_job,
        "job_name": job_name,
        "state": state,
        "nodes": nodes,
        "owner": owner,
        "account": account,
    }


def collect_live_batch_script(job_id: str) -> bytes:
    """Read exact batch-script bytes retained by the Slurm controller."""

    if JOB_RE.fullmatch(job_id) is None:
        raise ValueError("scheduler batch-script query has invalid job ID")
    scontrol = require_trusted_system_executable(SCONTROL, "scontrol executable")
    try:
        completed = subprocess.run(
            [scontrol, "write", "batch_script", job_id, "-"],
            check=True,
            capture_output=True,
            env=hardened_child_environment(),
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError(f"Slurm batch-script query failed for job {job_id}") from error
    payload = completed.stdout
    if isinstance(payload, str):
        payload = payload.encode("utf-8")
    if not isinstance(payload, bytes) or not payload:
        raise ValueError("Slurm returned an empty or invalid stored batch script")
    return payload


def normalized_batch_script_sha256(payload: bytes) -> str:
    """Return the controller-normalized self-digest of one batch script."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("prepared batch script is not UTF-8") from error
    matches = BATCH_SCRIPT_DIGEST_PATTERN.findall(text)
    if len(matches) != 1:
        raise ValueError("prepared batch script lacks one self-digest")
    normalized = BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}", text
    )
    return sha256_bytes(normalized.encode("utf-8"))


def authenticate_live_scheduler_job(job_id: str, nodes: int, job_name: str) -> dict[str, object]:
    """Require live Slurm to retain the exact submitted allocation."""

    record = collect_live_scheduler_job(job_id)
    if (
        record.get("job_id") != job_id
        or record.get("nodes") != nodes
        or record.get("job_name") != job_name
        or record.get("state") not in ACTIVE_SCHEDULER_STATES
        or record.get("owner") != scheduler_owner()
        or str(record.get("account", "")).casefold() != ACCOUNT.casefold()
    ):
        raise ValueError(f"live Slurm job {job_id} differs from its submitted manifest")
    return record


def authenticate_submitted_job(info: dict[str, object]) -> dict[str, object]:
    """Bind one submitted manifest to live Slurm and scheduler-stored script bytes."""

    manifest = info["manifest"]
    command = manifest["command"]
    paths = manifest.get("paths")
    if not isinstance(paths, dict):
        raise ValueError(f"{info['case_id']}/{info['segment']} submitted manifest lacks paths")
    batch_path = manifest_expected_path(info["case_id"], info["segment"]).parent / "cgl_lf_stage_i.sbatch"
    if normalized_path(Path(str(paths.get("batch_script", "")))) != batch_path:
        raise ValueError(f"{info['case_id']}/{info['segment']} batch-script path is not exact")
    payload, evidence = read_stable_regular_file(batch_path, "submitted batch script")
    digest = require_nonempty_string(command.get("batch_script_sha256"), "batch-script digest")
    if SHA256_RE.fullmatch(digest) is None or normalized_batch_script_sha256(payload) != digest:
        raise ValueError(f"{info['case_id']}/{info['segment']} batch-script digest is invalid")
    job_id = str(manifest["job_id"])
    live_payload = collect_live_batch_script(job_id)
    if live_payload != payload:
        raise ValueError(f"{info['case_id']}/{info['segment']} Slurm-stored batch script differs")
    scheduler = authenticate_live_scheduler_job(
        job_id, int(info["nodes"]), expected_job_name(info["case_id"], info["segment"])
    )
    result = {
        "scheduler": scheduler,
        "batch_script": evidence,
        "batch_script_sha256": sha256_bytes(payload),
        "normalized_batch_script_sha256": digest,
    }
    info["submitted_authentication"] = result
    return result


def manifest_expected_path(case_id: str, segment: str) -> Path:
    """Return one exact canonical manifest path."""

    return expected_path(
        f"runs/mks24-stage-i/{EXECUTION_EPOCH}/{case_id}/{segment}/manifest/prepared_run.json"
    )


def validate_manifest(
    manifest: dict[str, object],
    evidence: dict[str, object],
    cases: dict[str, dict[str, object]],
    promoted: dict[str, object],
) -> dict[str, object]:
    """Validate one canonical manifest sufficiently for read-only wave planning."""

    if manifest.get("schema_version") != 3 or manifest.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("manifest identity is invalid")
    if manifest.get("project_root") != str(CANONICAL_ROOT):
        raise ValueError("manifest project root differs from the canonical root")
    run = manifest.get("run")
    allocation = manifest.get("allocation")
    command = manifest.get("command")
    if not isinstance(run, dict) or not isinstance(allocation, dict) or not isinstance(command, dict):
        raise ValueError("manifest lacks run/allocation/command objects")
    case_id = require_nonempty_string(run.get("case_id"), "manifest case_id")
    segment = require_nonempty_string(run.get("segment"), "manifest segment")
    if case_id not in ALL_CASES:
        raise ValueError(f"manifest contains unsupported case {case_id}")
    if normalized_path(Path(str(evidence["path"]))) != manifest_expected_path(case_id, segment):
        raise ValueError(f"{case_id}/{segment} manifest path is not canonical")
    expected_case = EXPECTED_CASES[case_id]
    if run.get("case_name") != expected_case[0] or run.get("resolution") != expected_case[2]:
        raise ValueError(f"{case_id}/{segment} run identity differs from the frozen matrix")
    index, start, target = parse_segment(segment, f"{case_id} manifest")
    if (
        command.get("input_revision") != SOURCE_REVISION
        or command.get("input_sha256") != expected_case[3]
        or command.get("matrix_sha256") != MATRIX_SHA256
        or command.get("executable_revision") != SOURCE_REVISION
        or command.get("executable_sha256") != EXECUTABLE_SHA256
        or command.get("build_manifest") != str(build_manifest_path())
    ):
        raise ValueError(f"{case_id}/{segment} scientific/build provenance is invalid")
    overrides = command.get("overrides")
    if (
        not isinstance(overrides, list)
        or len(overrides) != 1
        or not isinstance(overrides[0], str)
        or "=" not in overrides[0]
        or overrides[0].split("=", 1)[0] != "time/tlim"
        or decimal_value(overrides[0].split("=", 1)[1], "manifest time/tlim override")
        != target
    ):
        raise ValueError(f"{case_id}/{segment} time/tlim override is invalid")
    if decimal_value(command.get("time_tlim_target"), "manifest time_tlim_target") != target:
        raise ValueError(f"{case_id}/{segment} command target differs from its segment")
    nodes, _, _, reserved = validate_resource_values(
        case_id,
        allocation.get("nodes"),
        allocation.get("requested_walltime"),
        command.get("athena_walltime"),
        allocation.get("reserved_node_hours"),
        f"{case_id}/{segment}",
    )
    if (
        allocation.get("ranks_per_node") != RANKS_PER_NODE
        or allocation.get("cpus_per_task") != CPUS_PER_TASK
        or allocation.get("requested_seconds")
        != walltime_seconds(allocation.get("requested_walltime"), "manifest requested walltime")
    ):
        raise ValueError(f"{case_id}/{segment} resource shape differs from controller policy")
    state = manifest.get("state")
    if state not in {"prepared", "submitted", "recorded", "cancelled"}:
        raise ValueError(f"{case_id}/{segment} manifest state is invalid")
    prepared_utc = require_utc(manifest.get("prepared_utc"), f"{case_id}/{segment} prepared_utc")
    if prepared_utc > current_utc() + FUTURE_SKEW:
        raise ValueError(f"{case_id}/{segment} prepared_utc is implausibly future")
    if state == "prepared" and manifest.get("job_id") is not None:
        raise ValueError(f"{case_id}/{segment} prepared manifest unexpectedly has a job ID")
    if state in {"submitted", "recorded"} and JOB_RE.fullmatch(str(manifest.get("job_id", ""))) is None:
        raise ValueError(f"{case_id}/{segment} submitted/recorded job ID is invalid")
    if state in {"prepared", "submitted"}:
        utility = command.get("production_utility")
        bundle = command.get("source_bundle")
        promoted_helper = promoted["controller_helper"]
        promoted_bundle = promoted["source_bundle"]
        if utility != {
            "committed": True,
            "path": str(CONTROLLER_HELPER),
            "revision": promoted_helper["revision"],
            "sha256": promoted_helper["sha256"],
        }:
            raise ValueError(f"{case_id}/{segment} active helper provenance is not promoted")
        if not isinstance(bundle, dict) or (
            bundle.get("path") != promoted_bundle["path"]
            or bundle.get("sha256") != promoted_bundle["sha256"]
            or bundle.get("verified_revisions") != promoted_bundle["verified_revisions"]
        ):
            raise ValueError(f"{case_id}/{segment} active source-bundle provenance is not promoted")
    return {
        "manifest": manifest,
        "evidence": evidence,
        "path": str(evidence["path"]),
        "case_id": case_id,
        "segment": segment,
        "index": index,
        "start": start,
        "target": target,
        "state": state,
        "nodes": nodes,
        "reserved": reserved,
    }


def restart_loadability(payload: bytes, path: Path) -> dict[str, object]:
    """Validate the qualified restart parameter boundary and binary time header."""

    marker = b"<par_end>\n"
    end = payload[:MAX_RESTART_PARAMETER_DUMP_BYTES].find(marker)
    if end < 0:
        raise ValueError(f"restart lacks a loadable <par_end> parameter boundary: {path}")
    try:
        text = payload[:end].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"restart parameter dump is not UTF-8: {path}") from error
    block = ""
    markers = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "restart_time":
                markers.append(value.strip())
    if len(markers) != 1:
        raise ValueError(f"restart must contain one time/restart_time marker: {path}")
    payload_offset = end + len(marker)
    binary_offset = payload_offset + RESTART_TIME_OFFSET
    if len(payload) < payload_offset + RESTART_MESH_HEADER_SIZE:
        raise ValueError(f"restart binary mesh header is truncated: {path}")
    binary_time = float(struct.unpack_from(RESTART_TIME_FORMAT, payload, binary_offset)[0])
    if not math.isfinite(binary_time):
        raise ValueError(f"restart binary physical time is not finite: {path}")
    marker_text = markers[0]
    try:
        marker_time = float(marker_text)
    except ValueError as error:
        raise ValueError(f"restart time marker is not numeric: {path}") from error
    if not math.isfinite(marker_time):
        raise ValueError(f"restart time marker is not finite: {path}")
    if marker_text == format(binary_time, ".17g"):
        marker_mode = "full_precision"
    elif marker_text == format(binary_time, ".6g"):
        marker_mode = "legacy_default_precision"
    else:
        raise ValueError(f"restart marker does not authenticate binary physical time: {path}")
    if marker_mode not in ALLOWED_RESTART_MARKER_MODES:
        raise ValueError(f"restart marker precision is not qualified: {path}")
    return {
        "parameter_dump_size_bytes": payload_offset,
        "binary_time_offset_bytes": binary_offset,
        "binary_time_format": RESTART_TIME_FORMAT,
        "marker_mode": marker_mode,
        "time": binary_time,
    }


def validate_restart_group(
    terminal: object,
    *,
    expected_root: Path,
    expected_count: int,
    expected_time: Decimal,
    label: str,
) -> dict[str, object]:
    """Authenticate one unique, complete, synchronized rank-local restart group."""

    if not isinstance(terminal, dict):
        raise ValueError(f"{label} terminal restart evidence is missing")
    path = require_nonempty_string(terminal.get("path"), f"{label} restart path")
    digest = require_nonempty_string(terminal.get("sha256"), f"{label} restart sha256")
    rank_files = terminal.get("rank_files")
    if (
        terminal.get("storage") != "per_rank"
        or SHA256_RE.fullmatch(digest) is None
        or not isinstance(rank_files, list)
        or len(rank_files) != expected_count
    ):
        raise ValueError(f"{label} terminal restart evidence is malformed")
    normalized_files = []
    seen: set[Path] = set()
    common_name: str | None = None
    binary_times = []
    for rank, item in enumerate(rank_files):
        if not isinstance(item, dict):
            raise ValueError(f"{label} restart rank-file record is malformed")
        item_path = normalized_path(Path(str(item.get("path", ""))))
        expected_parent = expected_root / f"rank_{rank:08d}"
        if item_path.parent != expected_parent:
            raise ValueError(f"{label} restart rank-file inventory is not exact")
        if item_path in seen:
            raise ValueError(f"{label} restart rank-file inventory contains duplicates")
        seen.add(item_path)
        if common_name is None:
            common_name = item_path.name
        elif item_path.name != common_name:
            raise ValueError(f"{label} restart rank-file names disagree")
        item_digest = require_nonempty_string(item.get("sha256"), f"{label} restart rank sha256")
        size = item.get("size_bytes")
        if (
            SHA256_RE.fullmatch(item_digest) is None
            or isinstance(size, bool)
            or not isinstance(size, int)
            or size <= 0
        ):
            raise ValueError(f"{label} restart rank-file metadata is invalid")
        payload, evidence = read_stable_regular_file(item_path, f"{label} restart rank {rank}")
        if evidence["sha256"] != item_digest or evidence["size_bytes"] != size:
            raise ValueError(f"{label} restart rank-file bytes differ from retained metadata")
        loadability = restart_loadability(payload, item_path)
        binary_times.append(Decimal(str(loadability["time"])))
        normalized_files.append(
            {
                "path": str(item_path),
                "sha256": item_digest,
                "size_bytes": size,
                "evidence": evidence,
                "loadability": loadability,
            }
        )
    primary = normalized_path(Path(path))
    terminal_size = terminal.get("size_bytes")
    if (
        primary != Path(normalized_files[0]["path"])
        or normalized_files[0]["sha256"] != digest
        or terminal_size != normalized_files[0]["size_bytes"]
    ):
        raise ValueError(f"{label} terminal restart primary record differs from rank evidence")
    if any(abs(value - expected_time) > Decimal("1e-12") for value in binary_times):
        raise ValueError(f"{label} restart binary times differ from the inspected endpoint")
    if any(abs(value - binary_times[0]) > Decimal("1e-12") for value in binary_times[1:]):
        raise ValueError(f"{label} restart rank-file binary times disagree")
    return {
        "restart_file": str(primary),
        "restart_sha256": digest,
        "restart_time": decimal_text(expected_time),
        "restart_files": normalized_files,
    }


def terminal_restart(info: dict[str, object]) -> dict[str, object]:
    """Return byte-authenticated restart provenance for one recorded terminal manifest."""

    cached = info.get("restart_validation")
    if isinstance(cached, dict):
        return cached
    manifest = info["manifest"]
    inspection = manifest.get("scientific_inspection")
    if not isinstance(inspection, dict):
        raise ValueError(f"{info['case_id']}/{info['segment']} lacks scientific inspection")
    final_time = decimal_value(inspection.get("final_time"), "inspection final_time")
    if decimal_value(inspection.get("terminal_restart_time"), "terminal restart time") != final_time:
        raise ValueError("terminal restart time differs from the inspected endpoint")
    validated = validate_restart_group(
        inspection.get("terminal_restart"),
        expected_root=expected_path(
            f"runs/mks24-stage-i/{EXECUTION_EPOCH}/{info['case_id']}/{info['segment']}/output/rst"
        ),
        expected_count=int(info["nodes"]) * RANKS_PER_NODE,
        expected_time=final_time,
        label=f"{info['case_id']}/{info['segment']}",
    )
    result = {"parent_manifest": info["evidence"], **validated}
    continuation_evidence = info.get("clean_partial_continuation_evidence")
    if isinstance(continuation_evidence, dict):
        result["clean_partial_continuation_evidence"] = continuation_evidence
    info["restart_validation"] = result
    return result


def validate_parent_link(current: dict[str, object], prior: dict[str, object]) -> None:
    """Require exact authenticated continuation-parent provenance."""

    command = current["manifest"]["command"]
    parent = command.get("parent_segment")
    restart = terminal_restart(prior)
    prior_manifest = prior["manifest"]
    prior_result = prior_manifest["accounting"]["result"]
    prior_final = decimal_value(prior_manifest["scientific_inspection"]["final_time"], "prior final_time")
    expected = {
        "case_id": prior["case_id"],
        "executable_sha256": EXECUTABLE_SHA256,
        "execution_epoch": EXECUTION_EPOCH,
        "final_time": float(prior_final),
        "input_sha256": EXPECTED_CASES[prior["case_id"]][3],
        "manifest": prior["path"],
        "restart_files": [item["path"] for item in restart["restart_files"]],
        "restart_sha256": restart["restart_sha256"],
        "restart_time": float(prior_final),
        "result": prior_result,
        "segment": prior["segment"],
    }
    if parent != expected:
        raise ValueError(f"{current['case_id']}/{current['segment']} parent provenance is invalid")
    if (
        command.get("source_restart_file") != restart["restart_file"]
        or command.get("restart_sha256") != restart["restart_sha256"]
    ):
        raise ValueError(f"{current['case_id']}/{current['segment']} restart provenance is invalid")
    if abs(current["start"] - prior_final) > Decimal("0.000001"):
        raise ValueError(f"{current['case_id']}/{current['segment']} start token differs from its parent")


def authenticate_controller_retained_file(
    value: object, label: str
) -> tuple[bytes, dict[str, object]]:
    """Authenticate one exact file record emitted by the production controller."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} record must be an object")
    require_exact_keys(value, ("path", "size_bytes", "sha256"), f"{label} record")
    path = normalized_path(Path(require_nonempty_string(value["path"], f"{label} path")))
    payload, evidence = read_stable_regular_file(path, label)
    if (
        value["sha256"] != evidence["sha256"]
        or value["size_bytes"] != evidence["size_bytes"]
    ):
        raise ValueError(f"{label} bytes differ from retained controller metadata")
    return payload, evidence


def authenticate_controller_product(
    value: object, root: Path, expected_ranks: int, label: str
) -> list[dict[str, object]]:
    """Authenticate one exact rank-local product record emitted by the controller."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} product must be an object")
    require_exact_keys(
        value,
        ("path", "size_bytes", "sha256", "storage", "rank_files"),
        f"{label} product",
    )
    rank_files = value["rank_files"]
    if value["storage"] != "per_rank" or not isinstance(rank_files, list):
        raise ValueError(f"{label} is not an exact rank-local controller product")
    if len(rank_files) != expected_ranks:
        raise ValueError(f"{label} rank-local inventory is incomplete")
    retained = []
    common_name = None
    for rank, item in enumerate(rank_files):
        _, evidence = authenticate_controller_retained_file(item, f"{label} rank {rank}")
        path = normalized_path(Path(str(evidence["path"])))
        if path.parent != root / f"rank_{rank:08d}":
            raise ValueError(f"{label} rank-local path inventory differs")
        if common_name is None:
            common_name = path.name
        elif path.name != common_name:
            raise ValueError(f"{label} rank-local file names disagree")
        retained.append({key: item[key] for key in ("path", "size_bytes", "sha256")})
    if any(value[key] != retained[0][key] for key in ("path", "size_bytes", "sha256")):
        raise ValueError(f"{label} representative differs from rank zero")
    return retained


def parse_controller_history(payload: bytes, label: str) -> dict[str, list[float]]:
    """Parse exact finite Athena history bytes using the controller's grammar."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8") from error
    labels: list[str] = []
    rows: list[list[float]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            found = re.findall(r"\[(\d+)\]=([^\s]+)", line)
            if found:
                indices = [int(index) for index, _ in found]
                candidate = [name for _, name in found]
                if (
                    indices != list(range(1, len(candidate) + 1))
                    or len(candidate) != len(set(candidate))
                    or (labels and candidate != labels)
                ):
                    raise ValueError(f"{label} labels are ambiguous or invalid")
                labels = candidate
            continue
        elif line.strip():
            try:
                rows.append([float(value) for value in line.split()])
            except ValueError as error:
                raise ValueError(f"{label} contains a nonnumeric row") from error
    if (
        not labels
        or len(labels) != len(set(labels))
        or not rows
        or any(len(row) != len(labels) for row in rows)
        or any(not math.isfinite(value) for row in rows for value in row)
    ):
        raise ValueError(f"{label} is incomplete or non-finite")
    return {
        name: [row[index] for row in rows] for index, name in enumerate(labels)
    }


def continuation_plasma_evidence(
    case_id: str,
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    *,
    allow_frozen_e03_no_max_ndiv: bool = False,
) -> dict[str, object]:
    """Reproduce the controller's complete clean-partial plasma policy evidence."""

    mhd_required = {
        "time", "mass", "tot-E", "lf_nstage", "lf_qface", "lf_qprcap",
        "lf_qpr10", "lf_qpecap", "lf_qpe10", "lf_qprwrk", "lf_qpewrk",
        "lf_cpwrk", "lf_cawrk", "lf_hwproj", *STRICT_LF_FAILURE_COLUMNS,
    }
    user_required = {"time", "mass", "hard_vol", "force_work"}
    missing_mhd = sorted(mhd_required - set(mhd))
    missing_user = sorted(user_required - set(user))
    if missing_mhd or missing_user:
        raise ValueError(
            "continuation histories lack required plasma columns: "
            f"MHD={missing_mhd}, user={missing_user}"
        )
    has_max_ndiv = "max_ndiv" in user
    if not has_max_ndiv and not allow_frozen_e03_no_max_ndiv:
        raise ValueError(
            "continuation normalized CT divB is irreducibly unavailable for the "
            "qualified historical E03 executable: user history lacks max_ndiv, "
            "retained mhd_w_bcc snapshots are cell-centered, and retained restart "
            "face fields have no authenticated CT-divergence interpretation"
        )
    lengths = {len(values) for values in [*mhd.values(), *user.values()]}
    if len(lengths) != 1 or not lengths or next(iter(lengths)) < 2:
        raise ValueError("continuation histories have inconsistent or insufficient rows")
    if mhd["time"] != user["time"] or any(
        right <= left for left, right in zip(mhd["time"], mhd["time"][1:])
    ):
        raise ValueError("continuation histories are not exactly synchronized and monotonic")
    floor = 1.0e-30

    def relative_drift(values: list[float]) -> float:
        return max(abs(value - values[0]) for value in values) / max(
            abs(values[0]), floor
        )

    mhd_mass_drift = relative_drift(mhd["mass"])
    user_mass_drift = relative_drift(user["mass"])
    mass_mismatch = max(
        abs(left - right) / max(abs(left), abs(right), floor)
        for left, right in zip(mhd["mass"], user["mass"])
    )
    maximum_divb = max(user["max_ndiv"]) if has_max_ndiv else None
    if (
        mhd_mass_drift > CONTINUATION_MASS_TOLERANCE
        or user_mass_drift > CONTINUATION_MASS_TOLERANCE
        or mass_mismatch > CONTINUATION_MASS_TOLERANCE
    ):
        raise ValueError("continuation mass conservation exceeds accepted policy")
    if maximum_divb is not None and (
        maximum_divb < 0.0 or maximum_divb >= CONTINUATION_MAX_NORMALIZED_CT_DIVB
    ):
        raise ValueError("continuation normalized CT divB exceeds accepted policy")
    if any(value != 0.0 for value in user["hard_vol"]):
        raise ValueError("continuation hard_vol is nonzero")
    strict_maxima = {name: max(mhd[name]) for name in STRICT_LF_FAILURE_COLUMNS}
    if any(value != 0.0 for name in STRICT_LF_FAILURE_COLUMNS for value in mhd[name]):
        raise ValueError("continuation strict LF failure counter is nonzero")
    count_columns = [
        "lf_nstage", "lf_qface", "lf_qprcap", "lf_qpr10", "lf_qpecap",
        "lf_qpe10", "lf_hwproj",
        *[name for name in ("lf_mirror", "lf_firehs") if name in mhd],
    ]
    for name in count_columns:
        values = mhd[name]
        if (
            any(value < 0.0 or not value.is_integer() for value in values)
            or any(right < left for left, right in zip(values, values[1:]))
        ):
            raise ValueError(f"continuation {name} violates count policy")
    for name in ("lf_qprcap", "lf_qpr10", "lf_qpecap", "lf_qpe10"):
        if any(value > qface for value, qface in zip(mhd[name], mhd["lf_qface"])):
            raise ValueError(f"continuation {name} exceeds lf_qface")
        if any(
            right - left > qright - qleft
            for left, right, qleft, qright in zip(
                mhd[name], mhd[name][1:], mhd["lf_qface"], mhd["lf_qface"][1:]
            )
        ):
            raise ValueError(f"continuation {name} increment exceeds lf_qface")
    if (
        mhd["lf_nstage"][-1] <= mhd["lf_nstage"][0]
        or mhd["lf_qface"][-1] <= mhd["lf_qface"][0]
    ):
        raise ValueError("continuation LF stage or qface diagnostics did not advance")
    state_scale = max(abs(mhd["tot-E"][0]), abs(mhd["tot-E"][-1]), floor)

    def activity(delta: float, label: str) -> dict[str, float]:
        absolute = abs(delta)
        normalized = absolute / state_scale
        if (
            absolute <= CONTINUATION_ACTIVITY_ABSOLUTE_GT
            or normalized <= CONTINUATION_ACTIVITY_NORMALIZED_GT
        ):
            raise ValueError(f"continuation {label} did not exceed accepted policy")
        return {
            "delta": delta,
            "absolute": absolute,
            "state_normalized": normalized,
        }

    forcing = activity(
        user["force_work"][-1] - user["force_work"][0], "forcing work"
    )
    passive = case_id in PASSIVE_CASES
    finite_limiter = case_id in FINITE_LIMITER_CASES
    if passive:
        if any(value != 0.0 for name in ("lf_cpwrk", "lf_cawrk") for value in mhd[name]):
            raise ValueError("passive continuation contains pressure-work feedback")
        pressure_activity = None
        closure = None
    else:
        pressure_delta = max(
            (
                mhd[name][-1] - mhd[name][0]
                for name in ("lf_cpwrk", "lf_cawrk")
            ),
            key=abs,
        )
        pressure_activity = activity(pressure_delta, "active pressure work")
        delta_energy = mhd["tot-E"][-1] - mhd["tot-E"][0]
        delta_force = user["force_work"][-1] - user["force_work"][0]
        closure = abs(delta_energy - delta_force) / max(
            abs(delta_energy), abs(delta_force), floor
        )
        if closure >= 1.0e-8:
            raise ValueError("continuation forcing-work closure exceeds accepted policy")
    if finite_limiter and any(value != 0.0 for value in mhd["lf_hwproj"]):
        raise ValueError("finite-limiter continuation contains hardwall projections")
    evidence = {
        "schema_version": 2,
        "policy": CONTINUATION_PLASMA_POLICY,
        "case_id": case_id,
        "continuation_authorized": True,
        "checks": {
            "finite_synchronized_histories": True,
            "mass_conservation": True,
            "normalized_ct_divb": True,
            "strict_lf_policy": True,
            "forcing_policy": True,
            "pressure_feedback_policy": True,
            "limiter_policy": True,
        },
        "normalized_ct_divb_evidence": {
            "source": "authenticated user history max_ndiv",
            "comparison": "strictly-less-than",
            "threshold": CONTINUATION_MAX_NORMALIZED_CT_DIVB,
            "maximum": maximum_divb,
            "passed": True,
        },
        "measurements": {
            "mhd_mass_relative_drift": mhd_mass_drift,
            "user_mass_relative_drift": user_mass_drift,
            "mhd_user_mass_relative_mismatch": mass_mismatch,
            "normalized_ct_divb_max": maximum_divb,
            "strict_lf_maxima": strict_maxima,
            "forcing_activity": forcing,
            "pressure_activity": pressure_activity,
            "forcing_work_closure_normalized_residual": closure,
            "final_hardwall_projection_count": mhd["lf_hwproj"][-1],
        },
    }
    if maximum_divb is None:
        evidence["checks"]["normalized_ct_divb"] = False
        evidence["normalized_ct_divb_evidence"] = None
        evidence["measurements"]["normalized_ct_divb_max"] = None
    return evidence


def read_frozen_e03_migration_binding(
    binding: dict[str, object], label: str
) -> tuple[Path, bytes, dict[str, object]]:
    """Read one exact owner-controlled frozen-E03 migration artifact."""

    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding is invalid")
    relative_text = require_nonempty_string(binding.get("path"), f"{label} path")
    relative = Path(relative_text)
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or relative.as_posix() != relative_text
    ):
        raise ValueError(f"{label} path is not normalized and root-relative")
    mode = binding.get("mode")
    size_bytes = binding.get("size_bytes")
    expected_sha256 = require_sha256(binding.get("sha256"), f"{label} SHA-256")
    if (
        isinstance(mode, bool)
        or not isinstance(mode, int)
        or isinstance(size_bytes, bool)
        or not isinstance(size_bytes, int)
        or size_bytes < 1
    ):
        raise ValueError(f"{label} binding is invalid")
    path = expected_path(relative_text)
    retained, profile = read_stable_regular_file(path, label)
    if (
        profile["mode"] != f"{mode:04o}"
        or profile["size_bytes"] != size_bytes
        or profile["sha256"] != expected_sha256
    ):
        raise ValueError(f"{label} differs from the controller migration binding")
    return path, retained, {
        "path": str(path),
        "sha256": expected_sha256,
        "size_bytes": size_bytes,
        "mode": f"{mode:04o}",
        "links": 1,
    }


def read_frozen_e03_migration_json(
    binding: dict[str, object], label: str
) -> tuple[dict[str, object], dict[str, object]]:
    """Read and decode one exact JSON frozen-E03 migration artifact."""

    _, payload, public_binding = read_frozen_e03_migration_binding(binding, label)
    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=duplicate_rejecting_object
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"{label} is not unambiguous UTF-8 JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain an object")
    return value, public_binding


def validate_frozen_e03_independent_validation(
    contract: dict[str, object],
    validation: dict[str, object],
    review: dict[str, object] | None,
    bindings: dict[str, dict[str, object] | None],
) -> None:
    """Require the controller/recost exact independent migration validation."""

    manifest = contract["manifest"]
    inspection = contract["inspection"]
    validation_contract = contract["independent_validation"]
    if not all(
        isinstance(item, dict)
        for item in (manifest, inspection, validation_contract)
    ):
        raise ValueError("frozen-E03 migration contract differs")
    case_id = str(contract["case_id"])
    job_id = str(contract["job_id"])
    segment = str(contract["segment"])
    expected_common = {
        "case_id": case_id,
        "job_id": job_id,
        "segment": segment,
        "schema_version": validation_contract["schema_version"],
        "record_type": validation_contract["record_type"],
        "validation_accepted": True,
        "inspection_sha256": inspection["sha256"],
    }
    if any(validation.get(key) != value for key, value in expected_common.items()):
        raise ValueError("frozen-E03 independent validation differs")
    product_sha256 = validation.get("product_sha256")
    product_sizes = validation.get("product_sizes")
    if not isinstance(product_sha256, dict) or not isinstance(product_sizes, dict):
        raise ValueError("frozen-E03 independent validation lacks product inventory")
    for key in ("mhd_history", "user_history"):
        history = bindings.get(key)
        if (
            not isinstance(history, dict)
            or product_sha256.get(history["path"]) != history["sha256"]
            or product_sizes.get(history["path"]) != history["size_bytes"]
        ):
            raise ValueError("frozen-E03 independent validation history binding differs")
    if validation.get("strict_failure_maxima") != {
        name: 0 for name in STRICT_LF_FAILURE_COLUMNS
    }:
        raise ValueError("frozen-E03 independent validation strict LF evidence differs")
    if case_id == "R03":
        if (
            validation.get("clean_for_continuation") is not True
            or validation.get("formal_inspection_accepted") is not False
            or validation.get("product_count") != validation.get("expected_product_count")
            or review is not None
        ):
            raise ValueError("R03 frozen-E03 migration validation differs")
        return
    executable = {
        "revision": FROZEN_E03_EXECUTABLE_REVISION,
        "sha256": FROZEN_E03_EXECUTABLE_SHA256,
    }
    expected_ct = {
        "authorizing": False,
        "executable_revision": executable["revision"],
        "executable_sha256": executable["sha256"],
        "reason": FROZEN_E03_CT_DIVERGENCE_REASON,
        "status": "historically_unavailable",
    }
    expected_face = {
        "authorizing": False,
        "reason": (
            "native restart values are fully finite-decoded, but face-centered "
            "magnetic fields are not independently interpreted for CT divergence"
        ),
        "status": "unavailable",
    }
    independent = validation.get("independent_authorization")
    if (
        validation.get("record_role") != "historical_read_only_non_authorizing"
        or validation.get("authorization_effect") != "none"
        or validation.get("normalized_ct_divb_validation") != expected_ct
        or validation.get("restart_face_field_validation") != expected_face
        or not isinstance(independent, dict)
        or independent.get("authorizing") is not False
        or independent.get("status") != "non_authorizing"
        or not isinstance(review, dict)
    ):
        raise ValueError("R12 frozen-E03 migration validation differs")
    reviewed = review.get("reviewed_validation")
    verified = review.get("verified_validation")
    validation_binding = bindings.get("independent_validation")
    if (
        review.get("schema_version") != 1
        or review.get("record_type")
        != "stage-i-independent-segment-validation-independent-review"
        or review.get("case_id") != case_id
        or review.get("job_id") != job_id
        or review.get("segment") != segment
        or review.get("validation_accepted") is not True
        or review.get("decision") != "approved-for-publication"
        or not isinstance(reviewed, dict)
        or not isinstance(validation_binding, dict)
        or reviewed.get("sha256") != validation_binding["sha256"]
        or not isinstance(verified, dict)
        or verified.get("case_id") != case_id
        or verified.get("job_id") != job_id
        or verified.get("segment") != segment
        or verified.get("inspection_sha256") != inspection["sha256"]
        or verified.get("validation_accepted") is not True
        or verified.get("non_authorizing") is not True
        or verified.get("authorization_effect") != "none"
    ):
        raise ValueError("R12 frozen-E03 independent review binding differs")


def require_frozen_e03_no_max_ndiv_migration(
    info: dict[str, object],
    inspection: dict[str, object],
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    expected_ranks: int,
) -> dict[str, object]:
    """Reproduce the controller/recost exact retained R03/R12 migration."""

    case_id = str(info["case_id"])
    segment = str(info["segment"])
    manifest = info["manifest"]
    job_id = str(manifest.get("job_id", ""))
    retained_contract = FROZEN_E03_NO_MAX_NDIV_MIGRATIONS.get((case_id, job_id))
    if retained_contract is None:
        raise ValueError(
            f"manifest {job_id} is not eligible for frozen-E03 no-max_ndiv migration"
        )
    contract = {**retained_contract, "case_id": case_id, "job_id": job_id}
    command = manifest.get("command")
    run = manifest.get("run")
    allocation = manifest.get("allocation")
    accounting = manifest.get("accounting")
    if (
        manifest.get("project_root") != str(CANONICAL_ROOT)
        or manifest.get("execution_epoch") != EXECUTION_EPOCH
        or manifest.get("state") != "recorded"
        or not isinstance(command, dict)
        or command.get("executable_revision") != FROZEN_E03_EXECUTABLE_REVISION
        or command.get("executable_sha256") != FROZEN_E03_EXECUTABLE_SHA256
        or not isinstance(run, dict)
        or run.get("case_id") != case_id
        or run.get("segment") != segment
        or segment != contract["segment"]
        or not isinstance(allocation, dict)
        or allocation.get("nodes") * RANKS_PER_NODE != expected_ranks
        or expected_ranks != contract["ranks"]
        or not isinstance(accounting, dict)
        or accounting.get("result") != "clean_partial"
        or accounting.get("job_id") != job_id
        or "max_ndiv" in user
    ):
        raise ValueError(
            f"manifest {job_id} frozen-E03 no-max_ndiv migration provenance differs"
        )
    exact_manifest, manifest_binding = read_frozen_e03_migration_json(
        contract["manifest"], "frozen-E03 migration manifest"
    )
    exact_inspection, inspection_binding = read_frozen_e03_migration_json(
        contract["inspection"], "frozen-E03 migration inspection"
    )
    retained_manifest = {
        key: value for key, value in manifest.items() if not key.startswith("_")
    }
    if exact_manifest != retained_manifest or exact_inspection != inspection:
        raise ValueError("frozen-E03 migration manifest or inspection bytes differ")
    if (
        exact_manifest.get("scientific_inspection") != exact_inspection
        or set(inspection) != FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS
        or inspection.get("checks") != FROZEN_E03_NO_MAX_NDIV_CHECKS
        or inspection.get("manifest") != manifest_binding["path"]
        or normalized_path(Path(str(info["path"]))) != Path(manifest_binding["path"])
    ):
        raise ValueError("frozen-E03 migration inspection identity differs")
    bindings: dict[str, dict[str, object] | None] = {
        "manifest": manifest_binding,
        "inspection": inspection_binding,
    }
    histories = {"mhd_history": mhd, "user_history": user}
    for key in histories:
        path, payload, public_binding = read_frozen_e03_migration_binding(
            contract[key], f"frozen-E03 migration {key}"
        )
        bindings[key] = public_binding
        if (
            inspection.get(key)
            != {
                "path": str(path),
                "size_bytes": public_binding["size_bytes"],
                "sha256": public_binding["sha256"],
            }
            or parse_controller_history(payload, f"frozen-E03 migration {key}")
            != histories[key]
        ):
            raise ValueError(f"frozen-E03 migration {key} binding differs")
    validation, validation_binding = read_frozen_e03_migration_json(
        contract["independent_validation"],
        "frozen-E03 migration independent validation",
    )
    bindings["independent_validation"] = validation_binding
    review_contract = contract["independent_review"]
    if review_contract is None:
        review = None
        bindings["independent_review"] = None
    else:
        if not isinstance(review_contract, dict):
            raise ValueError("frozen-E03 migration independent review differs")
        review, review_binding = read_frozen_e03_migration_json(
            review_contract, "frozen-E03 migration independent review"
        )
        bindings["independent_review"] = review_binding
    validate_frozen_e03_independent_validation(contract, validation, review, bindings)
    executable = {
        "revision": FROZEN_E03_EXECUTABLE_REVISION,
        "sha256": FROZEN_E03_EXECUTABLE_SHA256,
    }
    normalized_ct_divb_evidence = {
        "ct_divergence_claimed": False,
        "status": "historically_unavailable",
        "reason": FROZEN_E03_CT_DIVERGENCE_REASON,
        "authorizing": False,
        "source": "authenticated frozen-E03 migration contract",
        "qualified_executable": executable,
        "independent_validation": validation_binding,
        "independent_review": bindings["independent_review"],
    }
    evidence = continuation_plasma_evidence(
        case_id, mhd, user, allow_frozen_e03_no_max_ndiv=True
    )
    r12_historical_inventory = case_id == "R12"
    evidence["policy"] = FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
    evidence["continuation_authorized"] = False
    evidence["continuation_eligible"] = not r12_historical_inventory
    evidence["eligibility_only"] = True
    evidence["normalized_ct_divb_evidence"] = normalized_ct_divb_evidence
    evidence["ct_divergence_claimed"] = False
    evidence["ct_divergence_reason"] = FROZEN_E03_CT_DIVERGENCE_REASON
    evidence["authorization_basis"] = (
        (
            "none; historical frozen-E03 evidence is inventory-only and cannot "
            "authorize continuation"
        )
        if r12_historical_inventory
        else (
            "none; frozen-E03 migration evidence is historical inventory only and "
            "cannot authorize continuation"
        )
    )
    evidence["checks"]["frozen_e03_exact_migration"] = True
    evidence["migration_contract"] = {
        "schema_version": 1,
        "policy": FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY,
        "case_id": case_id,
        "job_id": job_id,
        "segment": segment,
        "frozen_science_executable": executable,
        "bindings": bindings,
        "authority": {
            "continuation_authorized": False,
            "submission_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }
    return evidence


def clean_partial_continuation_summary(
    inspection_profile: dict[str, object],
    plasma_evidence: dict[str, object],
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    *,
    frozen_e03_migration: bool,
) -> dict[str, object]:
    """Expose the authenticated available diagnostics without overstating CT proof."""

    migration_eligibility = None
    if frozen_e03_migration:
        migration_contract = plasma_evidence.get("migration_contract")
        if not isinstance(migration_contract, dict):
            raise ValueError("frozen-E03 migration lacks its exact retained contract")
        migration_eligibility = {
            "schema_version": 1,
            "policy": FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY,
            "eligible": plasma_evidence["continuation_eligible"],
            "authorizing": False,
            "authorization_effect": "none",
            "requires_exact_recommended_profile_waiver": False,
            "migration_contract": migration_contract,
            "migration_contract_sha256": compact_json_sha256(migration_contract),
        }
    return {
        "inspection": inspection_profile,
        "policy": (
            FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
            if frozen_e03_migration
            else plasma_evidence["policy"]
        ),
        "authorizing": not frozen_e03_migration,
        "authorization_effect": "none" if frozen_e03_migration else "continuation-eligible",
        "ct_divergence_claimed": not frozen_e03_migration,
        "ct_divergence_reason": (
            FROZEN_E03_CT_DIVERGENCE_REASON if frozen_e03_migration else None
        ),
        "reason": FROZEN_E03_CT_DIVERGENCE_REASON if frozen_e03_migration else None,
        "migration_eligibility": migration_eligibility,
        "available_diagnostics": {
            "mhd_history_columns": sorted(mhd),
            "user_history_columns": sorted(user),
            "plasma_continuation_evidence": plasma_evidence,
        },
    }


def validate_clean_partial_continuation_evidence(
    info: dict[str, object], inspection: dict[str, object], final_time: Decimal
) -> dict[str, object]:
    """Authenticate current schema-4 or the exact frozen-E03 migration proof."""

    label = f"{info['case_id']}/{info['segment']} clean_partial continuation evidence"
    retained_path = normalized_path(Path(info["path"]).parent / "segment_inspection.json")
    retained, evidence = read_json_file(retained_path, label)
    if evidence["mode"] != "0644" or retained != inspection:
        raise ValueError(f"{label} differs from the retained schema-4 inspection")
    frozen_e03_migration = set(inspection) == FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS
    expected_inspection_keys = (
        FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS
        if frozen_e03_migration
        else CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS
    )
    require_exact_keys(inspection, expected_inspection_keys, label)
    manifest = info["manifest"]
    if (
        inspection["schema_version"] != 4
        or inspection["execution_epoch"] != EXECUTION_EPOCH
        or inspection["manifest"] != info["path"]
        or inspection["job_id"] != manifest["job_id"]
        or inspection["case_id"] != info["case_id"]
        or inspection["segment"] != info["segment"]
        or inspection["accepted"] is not False
        or inspection["clean_for_continuation"] is not True
        or (
            not frozen_e03_migration
            and inspection["plasma_continuation_policy"] != CONTINUATION_PLASMA_POLICY
        )
        or inspection["restart_time_marker_bypass"] is not False
        or require_utc(inspection["inspected_utc"], f"{label} timestamp")
        > current_utc() + FUTURE_SKEW
    ):
        raise ValueError(f"{label} identity or disposition differs")
    checks = inspection["checks"]
    if not isinstance(checks, dict):
        raise ValueError(f"{label} checks must be an object")
    expected_checks = (
        FROZEN_E03_NO_MAX_NDIV_CHECKS
        if frozen_e03_migration
        else CONTROLLER_CLEAN_PARTIAL_CHECKS
    )
    require_exact_keys(checks, expected_checks, f"{label} checks")
    if checks != expected_checks:
        raise ValueError(f"{label} does not prove a physically clean partial endpoint")

    output_root = normalized_path(Path(info["path"]).parents[1] / "output")
    mhd_payload, mhd_profile = authenticate_controller_retained_file(
        inspection["mhd_history"], f"{label} MHD history"
    )
    user_payload, user_profile = authenticate_controller_retained_file(
        inspection["user_history"], f"{label} user history"
    )
    if (
        Path(str(mhd_profile["path"])).parent != output_root
        or Path(str(user_profile["path"])).parent != output_root
        or not Path(str(mhd_profile["path"])).name.endswith(".mhd.hst")
        or not Path(str(user_profile["path"])).name.endswith(".user.hst")
    ):
        raise ValueError(f"{label} retained history paths differ from controller output")
    mhd = parse_controller_history(mhd_payload, f"{label} MHD history")
    user = parse_controller_history(user_payload, f"{label} user history")
    expected_ranks = int(info["nodes"]) * RANKS_PER_NODE
    reproduced_plasma = (
        require_frozen_e03_no_max_ndiv_migration(
            info, inspection, mhd, user, expected_ranks
        )
        if frozen_e03_migration
        else continuation_plasma_evidence(str(info["case_id"]), mhd, user)
    )
    if (
        (
            not frozen_e03_migration
            and inspection["plasma_continuation_evidence"] != reproduced_plasma
        )
        or Decimal(str(mhd["time"][-1])) != final_time
        or inspection["final_hardwall_projection_count"] != mhd["lf_hwproj"][-1]
    ):
        raise ValueError(f"{label} plasma continuation evidence is absent, stale, or changed")

    maxima = inspection["maximum_strict_failure_counts"]
    if not isinstance(maxima, dict):
        raise ValueError(f"{label} strict-failure maxima must be an object")
    require_exact_keys(maxima, STRICT_LF_FAILURE_COLUMNS, f"{label} strict-failure maxima")
    expected_maxima = {key: max(mhd[key]) for key in STRICT_LF_FAILURE_COLUMNS}
    if maxima != expected_maxima or any(value != 0.0 for value in expected_maxima.values()):
        raise ValueError(f"{label} retains a nonzero or stale strict LF failure maximum")

    snapshots = inspection["snapshots"]
    snapshot_times = inspection["snapshot_times"]
    restarts = inspection["restarts"]
    restart_times = inspection["restart_times"]
    marker_modes = inspection["restart_time_marker_modes"]
    if not all(isinstance(value, list) and value for value in (
        snapshots, snapshot_times, restarts, restart_times, marker_modes
    )):
        raise ValueError(f"{label} retained product/time evidence is incomplete")
    if (
        len(snapshots) != len(snapshot_times)
        or len(restarts) != len(restart_times)
        or len(restarts) != len(marker_modes)
    ):
        raise ValueError(f"{label} retained product/time cardinality differs")
    for index, snapshot in enumerate(snapshots):
        authenticate_controller_product(
            snapshot, output_root / "bin", expected_ranks, f"{label} snapshot {index}"
        )
    for index, (restart, time_value, modes) in enumerate(
        zip(restarts, restart_times, marker_modes)
    ):
        authenticate_controller_product(
            restart, output_root / "rst", expected_ranks, f"{label} restart {index}"
        )
        validated = validate_restart_group(
            restart,
            expected_root=output_root / "rst",
            expected_count=expected_ranks,
            expected_time=decimal_value(time_value, f"{label} restart time {index}"),
            label=f"{label} restart {index}",
        )
        observed_modes = [
            item["loadability"]["marker_mode"] for item in validated["restart_files"]
        ]
        if modes != observed_modes:
            raise ValueError(f"{label} restart marker-mode evidence differs")
    for key, values in (("snapshot_times", snapshot_times), ("restart_times", restart_times)):
        times = [decimal_value(value, f"{label} {key}") for value in values]
        if (
            times != sorted(times)
            or any(not (info["start"] <= value <= final_time) for value in times)
            or sum(value == final_time for value in times) != 1
        ):
            raise ValueError(f"{label} {key} does not authenticate the unique endpoint")
    terminal = inspection["terminal_restart"]
    terminal_matches = [record for record, value in zip(restarts, restart_times) if decimal_value(value, label) == final_time]
    if (
        terminal_matches != [terminal]
        or decimal_value(inspection["terminal_restart_time"], f"{label} terminal restart time")
        != final_time
    ):
        raise ValueError(f"{label} terminal restart evidence differs")
    return clean_partial_continuation_summary(
        evidence,
        reproduced_plasma,
        mhd,
        user,
        frozen_e03_migration=frozen_e03_migration,
    )


def validate_r12_historical_partial_inventory(
    infos: list[dict[str, object]],
) -> list[dict[str, object]]:
    """Authenticate the exact retained R12 partial without making it a parent."""

    claimed = [
        info
        for info in infos
        if info["segment"] == R12_HISTORICAL_PARTIAL_SEGMENT
        or str(info["manifest"].get("job_id")) == R12_HISTORICAL_PARTIAL_JOB_ID
    ]
    exact = [
        info
        for info in claimed
        if info["segment"] == R12_HISTORICAL_PARTIAL_SEGMENT
        and str(info["manifest"].get("job_id")) == R12_HISTORICAL_PARTIAL_JOB_ID
    ]
    if claimed and (len(claimed) != 1 or len(exact) != 1):
        raise ValueError("R12 historical clean_partial inventory identity differs")
    if not exact:
        return []
    info = exact[0]
    manifest = info["manifest"]
    accounting = manifest.get("accounting")
    inspection = manifest.get("scientific_inspection")
    if (
        info["case_id"] != "R12"
        or info["state"] != "recorded"
        or info["index"] != 0
        or info["start"] != 0
        or info["target"] != Decimal("0.25")
        or info["nodes"] != R12_FRESH_RERUN_NODES
        or not isinstance(accounting, dict)
        or accounting.get("result") != "clean_partial"
        or not isinstance(inspection, dict)
    ):
        raise ValueError("R12 historical clean_partial inventory contract differs")
    final_time = decimal_value(
        inspection.get("final_time"), "R12 historical clean_partial final time"
    )
    if (
        decimal_value(
            inspection.get("required_time"), "R12 historical clean_partial required time"
        )
        != info["target"]
        or not (info["start"] < final_time < info["target"])
    ):
        raise ValueError("R12 historical clean_partial inventory endpoint differs")
    info["clean_partial_continuation_evidence"] = (
        validate_clean_partial_continuation_evidence(info, inspection, final_time)
    )
    terminal_restart(info)
    return exact


def validated_lineage(
    case_id: str,
    infos: list[dict[str, object]],
    historical_inventory: list[dict[str, object]] | None = None,
) -> tuple[list[dict[str, object]], dict[str, object] | None]:
    """Return one parent-authenticated lineage while preserving cancelled identities."""

    historical_inventory = historical_inventory or []
    recorded = [info for info in infos if info["state"] == "recorded"]
    active = [info for info in infos if info["state"] in {"prepared", "submitted"}]
    cancelled = [info for info in infos if info["state"] == "cancelled"]
    if len(active) > 1:
        raise ValueError(f"{case_id} has multiple active manifests")
    index_owners: dict[int, dict[str, object]] = {}
    for info in [*infos, *historical_inventory]:
        index = int(info["index"])
        if index in index_owners:
            raise ValueError(f"{case_id} reuses operational segment index {index}")
        index_owners[index] = info
    recorded.sort(key=lambda item: item["index"])
    cancelled_indexes = {int(info["index"]) for info in cancelled}
    prior: dict[str, object] | None = None
    for info in recorded:
        manifest = info["manifest"]
        accounting = manifest.get("accounting")
        inspection = manifest.get("scientific_inspection")
        if not isinstance(accounting, dict) or not isinstance(inspection, dict):
            raise ValueError(f"{case_id}/{info['segment']} recorded evidence is incomplete")
        result = accounting.get("result")
        if result not in {"accepted", "clean_partial"}:
            raise ValueError(f"{case_id} has a non-accepted recorded result")
        final_time = decimal_value(inspection.get("final_time"), "inspection final_time")
        required_time = decimal_value(inspection.get("required_time"), "inspection required_time")
        if required_time != info["target"]:
            raise ValueError(f"{case_id}/{info['segment']} inspection target differs from preparation")
        if result == "accepted" and final_time != info["target"]:
            raise ValueError(f"{case_id}/{info['segment']} accepted endpoint is not exact")
        if result == "clean_partial" and not (info["start"] < final_time < info["target"]):
            raise ValueError(f"{case_id}/{info['segment']} clean_partial endpoint is invalid")
        if result == "clean_partial":
            info["clean_partial_continuation_evidence"] = (
                validate_clean_partial_continuation_evidence(info, inspection, final_time)
            )
        if prior is None:
            expected_root_index = max(
                (int(item["index"]) for item in historical_inventory), default=-1
            ) + 1
            if (
                info["index"] != expected_root_index
                or info["start"] != 0
                or manifest["command"].get("parent_segment") is not None
            ):
                raise ValueError(
                    f"{case_id} recorded lineage does not start fresh after retained identities"
                )
        else:
            skipped = set(range(int(prior["index"]) + 1, int(info["index"])))
            if int(info["index"]) <= int(prior["index"]) or skipped != (
                skipped & cancelled_indexes
            ):
                raise ValueError(
                    f"{case_id} recorded segment indexes are not contiguous through "
                    "retained cancellations"
                )
        terminal_restart(info)
        if prior is not None:
            validate_parent_link(info, prior)
        prior = info
    active_info = active[0] if active else None
    if active_info is not None:
        prior_operational_indexes = [
            int(info["index"]) for info in recorded + cancelled + historical_inventory
        ]
        expected_index = max(prior_operational_indexes, default=-1) + 1
        if active_info["index"] != expected_index:
            raise ValueError(
                f"{case_id} active segment index does not continue retained operational "
                "identities"
            )
        if prior is None:
            if active_info["start"] != 0 or active_info["manifest"]["command"].get("parent_segment") is not None:
                raise ValueError(f"{case_id} active fresh segment is not a fresh lineage root")
        else:
            validate_parent_link(active_info, prior)
    return recorded, active_info


def recost_json_sha256(value: object) -> str:
    """Return the schema-2 recost generator's stable JSON-line digest."""

    return sha256_bytes((json.dumps(value, sort_keys=True) + "\n").encode())


def compact_json_sha256(value: object) -> str:
    """Return the qualification utility's compact canonical JSON digest."""

    return sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    )


def validate_latest_recost_publication(audit_path: Path) -> None:
    """Require one uniquely latest promoted schema-2 recost publication."""

    publications: list[tuple[datetime, Path]] = []
    accounting = expected_path("accounting")
    for candidate in sorted(accounting.glob("*_recost_evidence.json.publication_audit.json")):
        value, evidence = read_json_file(candidate, f"recost publication inventory {candidate.name}")
        if value.get("record_type") != "stage-i-recost-recommendation-publication-audit":
            continue
        require_immutable_publication_evidence(
            evidence, evidence["sha256"], f"recost publication inventory {candidate.name}"
        )
        if value.get("schema_version") != 1 or value.get("execution_epoch") != EXECUTION_EPOCH:
            raise ValueError("recost publication inventory contains an invalid schema-2 audit")
        publications.append(
            (
                require_utc(value.get("published_utc"), "recost inventory publication"),
                normalized_path(candidate),
            )
        )
    if not publications:
        raise ValueError("schema-2 recost publication inventory is empty")
    latest_time = max(timestamp for timestamp, _ in publications)
    latest = [path for timestamp, path in publications if timestamp == latest_time]
    if len(latest) != 1 or latest[0] != normalized_path(audit_path):
        raise ValueError("selected recost is not the uniquely latest promoted publication")


def validate_recost_request_review_chain(
    request_binding: object,
    recost: dict[str, object],
) -> dict[str, object]:
    """Authenticate the exact recost request review and its process assurance."""

    if not isinstance(request_binding, dict):
        raise ValueError("recost publication request binding must be an object")
    require_exact_keys(request_binding, ("path", "sha256"), "recost publication request")
    request_path = normalized_path(
        Path(require_nonempty_string(request_binding["path"], "recost request path"))
    )
    request_sha256 = require_sha256(request_binding["sha256"], "recost request SHA-256")
    provenance = recost.get("provenance")
    if (
        request_path.parent != expected_path("accounting")
        or not isinstance(provenance, dict)
        or request_sha256 != provenance.get("request_sha256")
    ):
        raise ValueError("recost request publication binding differs")
    request, request_evidence = read_json_file(request_path, "recost request")
    if (
        request_evidence["mode"] != "0644"
        or request_evidence["sha256"] != request_sha256
        or request.get("schema_version") != 2
        or request.get("record_type") != "stage-i-recost-recommendation-request"
        or request.get("execution_epoch") != EXECUTION_EPOCH
        or request.get("checkpoint") != recost["checkpoint"]
        or request.get("generated_utc") != recost["generated_utc"]
        or request.get("expires_utc") != recost["expires_utc"]
        or request.get("requested_by") != recost["requested_by"]
        or request.get("scope") != recost["scope"]
    ):
        raise ValueError("recost request bytes or promoted identity differ")

    review_path = request_path.with_name(f"{request_path.name}.independent_review.json")
    review, review_evidence = read_json_file(
        review_path, "recost request independent review"
    )
    review_sha256 = require_sha256(
        provenance.get("request_independent_review_sha256"),
        "recost request independent review SHA-256",
    )
    if review_evidence["mode"] != "0444" or review_evidence["sha256"] != review_sha256:
        raise ValueError("recost request independent review bytes or mode differ")
    require_exact_keys(
        review,
        (
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate", "scope",
        ),
        "recost request independent review",
    )
    reviewer = review["reviewer"]
    scope = review["scope"]
    if not isinstance(reviewer, dict) or not isinstance(scope, dict):
        raise ValueError("recost request independent review fields differ")
    require_exact_keys(
        reviewer,
        (
            "agent_id", "role", "declared_process_independence",
            "identity_assurance", "identity_assurance_limitation",
        ),
        "recost request reviewer",
    )
    reviewer_agent = require_nonempty_string(
        reviewer["agent_id"], "recost request reviewer agent ID"
    )
    base_scope = {"non_authorizing": True}
    reviewed = require_utc(
        review["reviewed_utc"], "recost request independent review timestamp"
    )
    generated = require_utc(recost["generated_utc"], "schema-2 recost generation timestamp")
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-recost-request-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-evidence-generation"
        or reviewer["role"] != "independent recost request reviewer"
        or reviewer["declared_process_independence"] is not True
        or reviewer["identity_assurance"] != REQUEST_REVIEW_IDENTITY_ASSURANCE
        or reviewer["identity_assurance_limitation"] != REQUEST_REVIEW_IDENTITY_LIMITATION
        or review["candidate"] != {"path": str(request_path), "sha256": request_sha256}
        or scope != base_scope
        or reviewed < generated
        or reviewed > current_utc() + FUTURE_SKEW
    ):
        raise ValueError(
            "recost request review or declared non-cryptographic process independence differs"
        )
    assurance = declared_process_independence_assurance(
        {
            "recost-request-author": require_nonempty_string(
                recost["requested_by"], "schema-2 recost requester"
            ),
            "recost-request-independent-reviewer": reviewer_agent,
        },
        "recost request independent review",
    )
    return {
        "request": request,
        "request_evidence": request_evidence,
        "independent_review": review,
        "independent_review_evidence": review_evidence,
        "independent_review_assurance": assurance,
    }


def validate_recost_publication_chain(
    recost_path: Path,
    review_path: Path,
    audit_path: Path,
) -> dict[str, object]:
    """Authenticate the sole current profile/projection/reconciliation authority."""

    recost_path = normalized_path(recost_path)
    review_path = normalized_path(review_path)
    audit_path = normalized_path(audit_path)
    if (
        recost_path.parent != expected_path("accounting")
        or RECOST_ARTIFACT_RE.fullmatch(recost_path.name) is None
        or review_path != recost_independent_review_path(recost_path)
        or audit_path != recost_publication_audit_path(recost_path)
    ):
        raise ValueError("recost publication chain paths differ from the canonical schema")
    recost, recost_evidence = read_json_file(recost_path, "promoted schema-2 recost")
    review, review_evidence = read_json_file(review_path, "recost independent review")
    audit, audit_evidence = read_json_file(audit_path, "recost publication audit")
    for evidence, label in (
        (recost_evidence, "promoted schema-2 recost"),
        (review_evidence, "recost independent review"),
        (audit_evidence, "recost publication audit"),
    ):
        require_immutable_publication_evidence(evidence, evidence["sha256"], label)

    require_exact_keys(
        recost,
        (
            "schema_version", "record_type", "checkpoint", "artifact_name",
            "execution_epoch", "generated_utc", "expires_utc", "requested_by",
            "scope", "predecessor_recost", "authority", "publication_requirements",
            "recommendations", "barrier", "budget", "storage", "ledger",
            "reservations", "manifests", "r17_readiness", "promoted_f113",
            "reconcile", "provenance",
        ),
        "schema-2 recost",
    )
    match = RECOST_ARTIFACT_RE.fullmatch(recost_path.name)
    assert match is not None
    generated = require_fresh(
        recost_evidence, recost["generated_utc"], "schema-2 recost", R17_READINESS_MAX_AGE
    )
    expires = require_utc(recost["expires_utc"], "schema-2 recost expiry")
    if (
        recost["schema_version"] != 2
        or recost["record_type"] != "stage-i-recost-recommendation-evidence"
        or recost["artifact_name"] != recost_path.name
        or recost["checkpoint"] != f"F-{match.group('number')}"
        or recost["execution_epoch"] != EXECUTION_EPOCH
        or expires <= generated
        or expires - generated > R17_READINESS_MAX_AGE
        or current_utc() > expires
        or recost["authority"]
        != {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or recost["publication_requirements"]
        != {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        }
    ):
        raise ValueError("schema-2 recost identity, lifetime, or non-authority boundary differs")
    require_nonempty_string(recost["requested_by"], "schema-2 recost requester")
    require_nonempty_string(recost["scope"], "schema-2 recost scope")

    require_exact_keys(
        review,
        (
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate", "scope",
        ),
        "recost independent review",
    )
    reviewer = review["reviewer"]
    if not isinstance(reviewer, dict):
        raise ValueError("recost reviewer must be an object")
    require_exact_keys(reviewer, ("agent_id", "independent_from_generator"), "recost reviewer")
    require_nonempty_string(reviewer["agent_id"], "recost reviewer agent ID")
    review_scope = review["scope"]
    if (
        not isinstance(review_scope, dict)
        or review_scope != {"non_authorizing": True}
    ):
        raise ValueError("recost independent review scope differs")
    reviewed = require_utc(review["reviewed_utc"], "recost review timestamp")
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-recost-recommendation-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or reviewer["independent_from_generator"] is not True
        or review["candidate"]
        != {"path": str(recost_path), "sha256": recost_evidence["sha256"]}
    ):
        raise ValueError("recost independent review differs")
    review_assurance = declared_process_independence_assurance(
        {
            "recost-requester": require_nonempty_string(
                recost["requested_by"], "schema-2 recost requester"
            ),
            "recost-independent-reviewer": str(reviewer["agent_id"]),
        },
        "recost independent review",
    )

    require_exact_keys(
        audit,
        (
            "schema_version", "record_type", "execution_epoch", "published_utc",
            "transaction_id", "artifact", "recost_recommendations",
            "independent_review", "authority", "counts", "generator",
            "scheduler_evidence", "source_bundle", "stage_i_helper", "utility",
            "forensic_copy", "publication", "generalized_publication_context",
        ),
        "recost publication audit",
    )
    published = require_utc(audit["published_utc"], "recost publication timestamp")
    if (
        audit["schema_version"] != 1
        or audit["record_type"] != "stage-i-recost-recommendation-publication-audit"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"]
        != {
            "action_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or reviewed < generated
        or published < reviewed
        or published > current_utc() + FUTURE_SKEW
        or not require_nonempty_string(audit["transaction_id"], "recost transaction ID")
        or audit["recost_recommendations"] != recost["recommendations"]
        or audit["counts"] != recost["reconcile"].get("counts")
        or audit["publication"]
        != "same-directory-link-fsync-copy-exchange-forensic-retirement-fsync"
    ):
        raise ValueError("recost publication audit or chronology differs")
    validate_declared_publication_binding(
        audit["artifact"],
        expected=recost_path,
        digest=recost_evidence["sha256"],
        mode="0444",
        label="recost publication artifact",
    )
    validate_declared_publication_binding(
        audit["independent_review"],
        expected=review_path,
        digest=review_evidence["sha256"],
        mode="0444",
        label="recost publication review",
    )
    provenance = recost["provenance"]
    assert isinstance(provenance, dict)
    expected_generator = {
        "sha256": provenance.get("generator_sha256"),
        "revision": provenance.get("generator_revision"),
        "mode": "0755",
    }
    generator = audit["generator"]
    stage_i_helper = audit["stage_i_helper"]
    source_bundle = audit["source_bundle"]
    utility = audit["utility"]
    if (
        not isinstance(generator, dict)
        or not isinstance(stage_i_helper, dict)
        or not isinstance(source_bundle, dict)
        or not isinstance(utility, dict)
    ):
        raise ValueError("recost checkpoint publication provenance must be objects")
    require_exact_keys(generator, ("path", "sha256", "mode", "revision"), "recost audit generator")
    require_exact_keys(
        stage_i_helper,
        ("path", "sha256", "committed", "reconcile_execution", "revision"),
        "recost audit Stage I helper",
    )
    require_exact_keys(
        source_bundle,
        ("path", "sha256", "mode", "verified_revisions"),
        "recost audit source bundle",
    )
    require_exact_keys(
        utility,
        ("path", "sha256", "execution", "committed"),
        "recost audit checkpoint utility",
    )
    if (
        {key: generator[key] for key in expected_generator} != expected_generator
        or stage_i_helper
        != {
            "path": str(CONTROLLER_HELPER),
            "sha256": provenance.get("stage_i_helper_sha256"),
            "committed": True,
            "reconcile_execution": "descriptor",
            "revision": provenance.get("stage_i_helper_revision"),
        }
        or source_bundle
        != {
            "path": str(promoted_source_bundle_path(str(provenance.get("stage_i_helper_revision")))),
            "sha256": provenance.get("source_bundle_sha256"),
            "mode": "0644",
            "verified_revisions": provenance.get("source_bundle_verified_revisions"),
        }
        or utility["execution"] != "authenticated-descriptor"
        or utility["committed"] is not True
        or SHA256_RE.fullmatch(str(utility["sha256"])) is None
    ):
        raise ValueError("recost checkpoint publication provenance differs")

    forensic = audit["forensic_copy"]
    if not isinstance(forensic, dict):
        raise ValueError("recost checkpoint forensic binding must be an object")
    require_exact_keys(
        forensic,
        ("path", "sha256", "mode", "links", "generalized_vectors"),
        "recost checkpoint forensic binding",
    )
    forensic_path = normalized_path(Path(str(forensic["path"])))
    if (
        forensic_path.parent
        != expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_forensics")
        or forensic["sha256"] != recost_evidence["sha256"]
        or forensic["mode"] != "0444"
        or forensic["links"] != 1
        or forensic["generalized_vectors"] != "exact-artifact-payload"
    ):
        raise ValueError("recost checkpoint forensic binding differs")
    _, forensic_evidence = read_stable_regular_file(forensic_path, "recost checkpoint forensic copy")
    require_immutable_publication_evidence(
        forensic_evidence, recost_evidence["sha256"], "recost checkpoint forensic copy"
    )

    context = audit["generalized_publication_context"]
    if not isinstance(context, dict):
        raise ValueError("recost generalized publication context must be an object")
    require_exact_keys(
        context,
        (
            "schema_version", "artifact_name", "checkpoint", "generated_utc",
            "expires_utc", "request", "artifact", "independent_review",
            "authority", "publication_requirements", "recommendations", "barrier",
            "provenance", "predecessor_recost", "reconciliation", "budget",
            "storage", "r17_readiness", "controller_enforcement",
        ),
        "recost generalized publication context",
    )
    if (
        context["schema_version"] != 2
        or context["artifact_name"] != recost["artifact_name"]
        or context["checkpoint"] != recost["checkpoint"]
        or context["generated_utc"] != recost["generated_utc"]
        or context["expires_utc"] != recost["expires_utc"]
        or context["artifact"]
        != {"basename": recost_path.name, "sha256": recost_evidence["sha256"]}
        or context["independent_review"] != audit["independent_review"]
        or context["authority"] != recost["authority"]
        or context["publication_requirements"] != recost["publication_requirements"]
        or context["recommendations"] != recost["recommendations"]
        or context["barrier"] != recost["barrier"]
        or context["provenance"] != recost["provenance"]
        or context["predecessor_recost"] != recost["predecessor_recost"]
        or context["reconciliation"] != recost["reconcile"]
        or context["budget"] != recost["budget"]
        or context["storage"] != recost["storage"]
        or context["r17_readiness"] != recost["r17_readiness"]
        or context["controller_enforcement"]
        != {
            "state": recost["recommendations"]["controller_consumption_state"],
            "launch_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("recost generalized publication context differs")
    request_review_chain = validate_recost_request_review_chain(
        context["request"], recost
    )
    validate_latest_recost_publication(audit_path)
    return {
        "artifact": recost,
        "artifact_evidence": recost_evidence,
        "independent_review": review,
        "independent_review_evidence": review_evidence,
        "independent_review_assurance": review_assurance,
        "request_review_chain": request_review_chain,
        "publication_audit": audit,
        "publication_audit_evidence": audit_evidence,
        "generated_utc": generated,
        "expires_utc": expires,
        "published_utc": published,
    }


def validate_recost_reconciliation(
    recost: dict[str, object],
    counts: dict[str, int],
    qualification_evidence: dict[str, object],
) -> None:
    """Require embedded recost reconciliation to match independently read stores."""

    report = recost.get("reconcile")
    if not isinstance(report, dict):
        raise ValueError("recost reconciliation must be an object")
    if (
        report.get("execution_epoch") != EXECUTION_EPOCH
        or report.get("root") != str(CANONICAL_ROOT)
        or report.get("consistent") is not True
        or report.get("issues") != []
        or report.get("counts") != counts
    ):
        raise ValueError("recost reconciliation differs from canonical state")
    qualification = report.get("qualification")
    if not isinstance(qualification, dict) or (
        qualification.get("state") != "approved"
        or qualification.get("path") != str(qualification_path())
        or qualification.get("sha256") != qualification_evidence["sha256"]
        or qualification.get("approved_executable_revision") != SOURCE_REVISION
        or qualification.get("approved_executable_sha256") != EXECUTABLE_SHA256
    ):
        raise ValueError("recost reconciliation qualification differs")


def validate_recorded_accounting(
    info: dict[str, object], ledger_row: dict[str, str], reservation: dict[str, object]
) -> Decimal:
    """Require one recorded manifest/reservation/ledger scheduler and budget record."""

    manifest = info["manifest"]
    accounting = manifest.get("accounting")
    if not isinstance(accounting, dict) or accounting != ledger_row:
        raise ValueError(f"{info['case_id']}/{info['segment']} accounting differs from the ledger")
    if (
        ledger_row["state"] != "COMPLETED"
        or ledger_row["exit_code"] != "0:0"
        or ledger_row["execution_epoch"] != EXECUTION_EPOCH
        or ledger_row["case_id"] != info["case_id"]
        or ledger_row["case_name"] != EXPECTED_CASES[info["case_id"]][0]
        or ledger_row["segment"] != info["segment"]
        or ledger_row["executable_revision"] != SOURCE_REVISION
        or ledger_row["executable_sha256"] != EXECUTABLE_SHA256
        or ledger_row["input_revision"] != SOURCE_REVISION
        or ledger_row["result"] not in {"accepted", "clean_partial"}
    ):
        raise ValueError(f"{info['case_id']}/{info['segment']} scheduler/provenance fields are invalid")
    if JOB_RE.fullmatch(ledger_row["job_id"]) is None or ledger_row["job_id"] != str(manifest["job_id"]):
        raise ValueError(f"{info['case_id']}/{info['segment']} recorded job ID is invalid")
    submitted = require_scheduler_timestamp(ledger_row["submitted_utc"], "ledger submitted_utc")
    completed = require_scheduler_timestamp(ledger_row["completed_utc"], "ledger completed_utc")
    if completed < submitted or completed > current_utc() + FUTURE_SKEW:
        raise ValueError(f"{info['case_id']}/{info['segment']} scheduler timestamps are invalid")
    if (
        DIGITS_RE.fullmatch(ledger_row["nodes"]) is None
        or DIGITS_RE.fullmatch(ledger_row["elapsed_seconds"]) is None
    ):
        raise ValueError("ledger nodes/elapsed fields differ from the controller schema")
    nodes = int(ledger_row["nodes"])
    if nodes != info["nodes"]:
        raise ValueError(f"{info['case_id']}/{info['segment']} ledger nodes differ")
    requested_seconds = walltime_seconds(ledger_row["requested_walltime"], "ledger walltime")
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError("ledger requested walltime exceeds controller policy")
    if (
        ledger_row["requested_walltime"] != manifest["allocation"]["requested_walltime"]
        or abs(
            decimal_value(ledger_row["reserved_node_hours"], "ledger reserved node-hours")
            - info["reserved"]
        )
        > Decimal("0.000001")
    ):
        raise ValueError("ledger requested allocation differs from the manifest")
    elapsed = int(ledger_row["elapsed_seconds"])
    if elapsed <= 0:
        raise ValueError("ledger elapsed seconds must be positive")
    actual = decimal_value(ledger_row["actual_node_hours"], "ledger actual node-hours")
    if abs(actual - Decimal(nodes * elapsed) / Decimal(3600)) > Decimal("0.000001"):
        raise ValueError("ledger actual node-hours differ from scheduler accounting")
    if (
        reservation.get("state") != "recorded"
        or reservation.get("job_id") != ledger_row["job_id"]
        or reservation.get("result") != ledger_row["result"]
        or abs(decimal_value(reservation.get("actual_node_hours"), "reservation actual") - actual)
        > Decimal("0.000001")
    ):
        raise ValueError("recorded reservation differs from ledger accounting")
    return actual


def validate_operational_job_id_bindings(
    infos: list[dict[str, object]],
    reservations_by_manifest: dict[str, dict[str, object]],
) -> None:
    """Require each non-null job ID to bind one exact manifest/reservation identity."""

    identities_by_job: dict[str, tuple[str, str, str, str]] = {}
    for info in infos:
        manifest_job = info["manifest"].get("job_id")
        reservation_job = reservations_by_manifest[info["path"]].get("job_id")
        if manifest_job != reservation_job:
            raise ValueError(
                f"{info['case_id']}/{info['segment']} operational manifest/reservation "
                "job ID differs"
            )
        if manifest_job is None:
            continue
        if not isinstance(manifest_job, str) or JOB_RE.fullmatch(manifest_job) is None:
            raise ValueError(
                f"{info['case_id']}/{info['segment']} operational job ID is invalid"
            )
        identity = (
            info["path"], info["case_id"], info["segment"], info["state"],
        )
        if manifest_job in identities_by_job:
            raise ValueError(
                "non-null job ID binds multiple operational manifest/reservation identities"
            )
        identities_by_job[manifest_job] = identity


def validate_canonical_state(
    cases: dict[str, dict[str, object]],
    authority: dict[str, object],
    promoted: dict[str, object],
) -> dict[str, object]:
    """Authenticate exact drained-barrier stores bound by the promoted recost."""

    transaction_state = discover_transaction_state()
    recost = authority["artifact"]
    assert isinstance(recost, dict)
    provenance = recost["provenance"]
    assert isinstance(provenance, dict)

    ledger, ledger_evidence = read_ledger(ledger_path())
    ledger_binding = recost.get("ledger")
    if not isinstance(ledger_binding, dict) or (
        ledger_binding.get("rows") != len(ledger)
        or ledger_binding.get("sha256") != ledger_evidence["sha256"]
        or provenance.get("ledger_sha256") != ledger_evidence["sha256"]
    ):
        raise ValueError("recost ledger binding differs from the canonical ledger")
    reservations_value, reservations_evidence = read_json_value(reservations_path(), "reservations")
    if not isinstance(reservations_value, list):
        raise ValueError("reservations store must contain a list")
    reservations = reservations_value
    reservations_binding = recost.get("reservations")
    if not isinstance(reservations_binding, dict) or (
        reservations_binding.get("rows") != len(reservations)
        or reservations_binding.get("active") != 0
        or reservations_binding.get("sha256") != reservations_evidence["sha256"]
        or provenance.get("reservations_sha256") != reservations_evidence["sha256"]
    ):
        raise ValueError("recost reservations binding differs from the canonical store")
    qualification, qualification_evidence = read_json_file(qualification_path(), "qualification")
    if qualification_evidence != promoted["qualification"]:
        raise ValueError("canonical qualification changed after recost review")

    manifest_block = recost.get("manifests")
    if not isinstance(manifest_block, dict):
        raise ValueError("recost manifest authority must be an object")
    manifest_bindings = manifest_block.get("bindings")
    if not isinstance(manifest_bindings, list) or manifest_block.get("rows") != len(manifest_bindings):
        raise ValueError("recost manifest bindings must be an exact list")
    discovered = sorted(
        normalized_path(path)
        for path in expected_path(f"runs/mks24-stage-i/{EXECUTION_EPOCH}").glob(
            "*/*/manifest/prepared_run.json"
        )
    )
    bound_paths = sorted(
        expected_path(str(item.get("path", "")))
        for item in manifest_bindings
        if isinstance(item, dict)
    )
    if discovered != bound_paths:
        raise ValueError("recost manifest bindings differ from the canonical manifest store")
    infos: list[dict[str, object]] = []
    manifest_evidence = []
    for binding in manifest_bindings:
        if not isinstance(binding, dict) or set(binding) != {"path", "sha256"}:
            raise ValueError("recost manifest binding must contain exact path/SHA-256")
        path = expected_path(str(binding["path"]))
        manifest, retained = read_json_file(path, "canonical manifest")
        if retained["sha256"] != binding["sha256"]:
            raise ValueError("canonical manifest changed after recost review")
        if manifest.get("state") not in {"recorded", "cancelled"}:
            raise ValueError("recost authority is not a drained-barrier snapshot")
        infos.append(validate_manifest(manifest, retained, cases, promoted))
        manifest_evidence.append(retained)

    by_manifest = {info["path"]: info for info in infos}
    reservations_by_manifest: dict[str, dict[str, object]] = {}
    active_infos = []
    prepared_count = 0
    active_reserved = Decimal("0")
    for reservation in reservations:
        if not isinstance(reservation, dict):
            raise ValueError("reservation record must be an object")
        columns = frozenset(reservation)
        if not RESERVATION_REQUIRED_COLUMNS.issubset(columns):
            raise ValueError("reservation record lacks required controller fields")
        if not columns.issubset(RESERVATION_REQUIRED_COLUMNS | RESERVATION_OPTIONAL_COLUMNS):
            raise ValueError("reservation record contains unsupported controller fields")
        manifest_path_value = str(normalized_path(Path(str(reservation.get("manifest", "")))))
        if manifest_path_value in reservations_by_manifest:
            raise ValueError("reservations contain a duplicate manifest")
        if manifest_path_value not in by_manifest:
            raise ValueError("reservation lacks a bound canonical manifest")
        info = by_manifest[manifest_path_value]
        reservations_by_manifest[manifest_path_value] = reservation
        info["reservation"] = reservation
        prepared_utc = require_utc(reservation.get("prepared_utc"), "reservation prepared_utc")
        if (
            prepared_utc > current_utc() + FUTURE_SKEW
            or reservation.get("prepared_utc") != info["manifest"].get("prepared_utc")
        ):
            raise ValueError(f"{info['case_id']}/{info['segment']} reservation prepared_utc is invalid")
        execution_intent = reservation.get("execution_intent_sha256")
        if execution_intent is not None and SHA256_RE.fullmatch(str(execution_intent)) is None:
            raise ValueError("reservation execution-intent digest is invalid")
        if (
            reservation.get("execution_epoch") != EXECUTION_EPOCH
            or reservation.get("case_id") != info["case_id"]
            or reservation.get("case_name") != EXPECTED_CASES[info["case_id"]][0]
            or reservation.get("segment") != info["segment"]
            or reservation.get("state") != info["state"]
            or int(reservation.get("nodes", -1)) != info["nodes"]
            or reservation.get("requested_walltime")
            != info["manifest"]["allocation"]["requested_walltime"]
            or abs(decimal_value(reservation.get("reserved_node_hours"), "reservation reserved") - info["reserved"])
            > Decimal("0.000001")
        ):
            raise ValueError(f"{info['case_id']}/{info['segment']} reservation differs from manifest")
        if info["state"] in {"prepared", "submitted"}:
            active_infos.append(info)
            active_reserved += info["reserved"]
            if info["state"] == "prepared":
                if reservation.get("job_id") is not None:
                    raise ValueError("prepared reservation unexpectedly has a job ID")
                prepared_count += 1
            elif reservation.get("job_id") != info["manifest"].get("job_id"):
                raise ValueError("submitted reservation job ID differs from manifest")
        elif info["state"] == "cancelled":
            if (
                not isinstance(reservation.get("notes"), str)
                or JOB_RE.fullmatch(str(reservation.get("job_id", ""))) is None
                or reservation.get("job_id") != info["manifest"].get("job_id")
            ):
                raise ValueError("cancelled reservation job/notes differ from manifest")
    if set(reservations_by_manifest) != set(by_manifest):
        raise ValueError("canonical manifests and reservations are not one-to-one")
    validate_operational_job_id_bindings(infos, reservations_by_manifest)
    if prepared_count > 1:
        raise ValueError("controller state contains more than one prepared packet")
    if len(active_infos) > MAX_LANES:
        raise ValueError("controller state exceeds the four-lane ceiling")
    active_cases = [info["case_id"] for info in active_infos]
    if len(active_cases) != len(set(active_cases)):
        raise ValueError("controller state contains duplicate active case lanes")
    active_nodes = sum(int(info["nodes"]) for info in active_infos)
    if active_nodes > MAX_NODES:
        raise ValueError("controller state exceeds the ten-node ceiling")
    if R17 in active_cases and len(active_infos) != 1:
        raise ValueError("R17 is not exclusive in controller state")
    if len(active_infos) > 1 and any(case_id not in CONCURRENT_CASES for case_id in active_cases):
        raise ValueError("overlapping controller lanes are not restricted to R03-R16")
    submitted_authentication = {
        str(info["manifest"]["job_id"]): authenticate_submitted_job(info)
        for info in active_infos
        if info["state"] == "submitted"
    }

    ledger_by_job: dict[str, dict[str, str]] = {}
    actual = Decimal("0")
    cumulative = Decimal("0")
    for row_index, row in enumerate(ledger, start=1):
        job = row["job_id"]
        if job in ledger_by_job:
            raise ValueError("ledger contains duplicate jobs")
        ledger_by_job[job] = row
        matches = [info for info in infos if info["state"] == "recorded" and str(info["manifest"].get("job_id")) == job]
        if len(matches) != 1:
            raise ValueError("ledger job does not bind one recorded canonical manifest")
        value = validate_recorded_accounting(matches[0], row, reservations_by_manifest[matches[0]["path"]])
        actual += value
        cumulative += value
        cumulative_tolerance = Decimal(row_index) * Decimal("0.000001")
        if (
            abs(
                decimal_value(row["cumulative_stage_i_node_hours"], "ledger cumulative")
                - cumulative
            )
            > cumulative_tolerance
        ):
            raise ValueError("ledger cumulative node-hours are not coherent")
    recorded_infos = [info for info in infos if info["state"] == "recorded"]
    if len(recorded_infos) != len(ledger):
        raise ValueError("recorded manifests and ledger rows are not one-to-one")
    if actual + active_reserved > STAGE_I_BUDGET_NODE_HOURS or actual + active_reserved > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError("canonical state exceeds a controller budget ceiling")

    if active_infos or active_reserved != 0:
        raise ValueError("schema-2 recost authority may be consumed only at a drained barrier")
    if ledger_binding.get("cumulative_stage_i_node_hours") != ledger[-1][
        "cumulative_stage_i_node_hours"
    ]:
        raise ValueError("recost ledger cumulative use differs from canonical accounting")

    counts = {
        "transactions": sum(
            len(snapshot["entries"]) for snapshot in transaction_state.values()
        ),
        "reservations": len(reservations),
        "active_reservations": len(active_infos),
        "ledger_rows": len(ledger),
        "manifests": len(infos),
    }
    actual_evidence = {
        "ledger": ledger_evidence,
        "reservations": reservations_evidence,
        "qualification": qualification_evidence,
        "manifests": manifest_evidence,
    }
    validate_recost_reconciliation(recost, counts, qualification_evidence)

    by_case: dict[str, list[dict[str, object]]] = {case_id: [] for case_id in ALL_CASES}
    for info in infos:
        by_case[info["case_id"]].append(info)
    r12_historical_inventory = validate_r12_historical_partial_inventory(by_case["R12"])
    r12_historical_paths = {str(info["path"]) for info in r12_historical_inventory}
    lineages = {}
    for case_id in ALL_CASES:
        case_inventory = r12_historical_inventory if case_id == "R12" else []
        lineage_infos = [
            info
            for info in by_case[case_id]
            if str(info["path"]) not in r12_historical_paths
        ]
        lineages[case_id] = validated_lineage(case_id, lineage_infos, case_inventory)
    lineage_digest = authenticated_lineages_sha256(lineages)
    if (
        manifest_block.get("authenticated_lineages_sha256") != lineage_digest
        or provenance.get("authenticated_lineages_sha256") != lineage_digest
    ):
        raise ValueError("recost authenticated-lineage digest differs from canonical state")
    barrier = recost.get("barrier")
    if not isinstance(barrier, dict):
        raise ValueError("recost barrier must be an object")
    recorded_segments = barrier.get("recorded_segments")
    if not isinstance(recorded_segments, list) or not recorded_segments:
        raise ValueError("recost barrier recorded segments must be nonempty")
    expected_barrier = [
        {
            "job_id": row["job_id"],
            "case_id": row["case_id"],
            "segment": row["segment"],
            "result": row["result"],
        }
        for row in ledger[-len(recorded_segments):]
    ]
    if (
        recorded_segments != expected_barrier
        or barrier.get("job_ids") != sorted(item["job_id"] for item in expected_barrier)
        or not isinstance(barrier.get("scheduler_evidence"), list)
    ):
        raise ValueError("recost barrier differs from the canonical ledger suffix")
    next_indexes = {
        case_id: max((int(info["index"]) for info in by_case[case_id]), default=-1) + 1
        for case_id in ALL_CASES
    }
    return {
        "actual_node_hours": actual,
        "active_reserved_node_hours": active_reserved,
        "active_nodes": active_nodes,
        "active_infos": active_infos,
        "counts": counts,
        "lineages": lineages,
        "case_infos": by_case,
        "historical_inventory": {"R12": r12_historical_inventory},
        "next_indexes": next_indexes,
        "authority_observed": authority["generated_utc"],
        "evidence": {
            **actual_evidence,
            "transaction_stores": transaction_state,
            "submitted_authentication": submitted_authentication,
            "historical_inventory": {
                "R12": [
                    {
                        "manifest": info["evidence"],
                        "job_id": str(info["manifest"]["job_id"]),
                        "segment": info["segment"],
                        "result": info["manifest"]["accounting"]["result"],
                        "classification": "inventory-only-non-authorizing",
                        "continuation_authorized": False,
                    }
                    for info in r12_historical_inventory
                ]
            },
        },
        "transaction_state": transaction_state,
        "submitted_scheduler": {
            job_id: {
                "scheduler": value["scheduler"],
                "batch_script_sha256": value["batch_script_sha256"],
            }
            for job_id, value in submitted_authentication.items()
        },
    }


def lineage_endpoint(lineage: list[dict[str, object]]) -> Decimal:
    """Return the authenticated endpoint of one recorded lineage."""

    if not lineage:
        return Decimal("0")
    return decimal_value(
        lineage[-1]["manifest"]["scientific_inspection"]["final_time"], "lineage endpoint"
    )


def lineage_complete(lineage: list[dict[str, object]]) -> bool:
    """Return whether one lineage terminates in accepted exact t=10 evidence."""

    if not lineage or lineage_endpoint(lineage) != Decimal("10"):
        return False
    accounting = lineage[-1]["manifest"].get("accounting")
    return isinstance(accounting, dict) and accounting.get("result") == "accepted"


def authenticated_lineages_sha256(
    lineages: dict[str, tuple[list[dict[str, object]], dict[str, object] | None]]
) -> str:
    """Reproduce the schema-2 recost authenticated-lineage digest."""

    value = {
        case_id: [
            {
                "manifest": item["path"],
                "sha256": item["evidence"]["sha256"],
                "job_id": item["manifest"]["job_id"],
                "segment": item["segment"],
                "result": item["manifest"]["accounting"]["result"],
                "final_time": item["manifest"]["scientific_inspection"]["final_time"],
            }
            for item in lineage
        ]
        for case_id, (lineage, _) in sorted(lineages.items())
        if lineage
    }
    return recost_json_sha256(value)


def build_manifest_inventory_sha256() -> str:
    """Return the recost-compatible digest of the qualified build manifest."""

    inventory = []
    for path in sorted(build_manifest_path().iterdir(), key=lambda item: item.name):
        payload, evidence = read_stable_regular_file(path, f"build manifest inventory {path.name}")
        if not payload:
            raise ValueError("qualified build manifest contains an empty file")
        inventory.append(
            {
                "name": path.name,
                "mode": evidence["mode"],
                "sha256": evidence["sha256"],
            }
        )
    return recost_json_sha256(inventory)


def validate_recost_profiles(
    recost: dict[str, object],
    state: dict[str, object],
    promoted: dict[str, object],
) -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    """Validate exact next profiles supplied by the promoted recost authority."""

    recommendations = recost.get("recommendations")
    if not isinstance(recommendations, dict):
        raise ValueError("recost recommendations must be an object")
    required = {
        "mode",
        "authorizing",
        "recommended_next_profiles",
        "bounded_concurrency",
        "controller_consumption_state",
        "non_authorizing_reason",
    }
    allowed = required | {"sole_next_segment_recommendation"}
    if not required.issubset(recommendations) or not set(recommendations).issubset(allowed):
        raise ValueError("recost recommendation schema differs")
    mode = recommendations["mode"]
    profiles_value = recommendations["recommended_next_profiles"]
    bounded = recommendations["bounded_concurrency"]
    if (
        mode not in {"sole-next-profile", "bounded-wave"}
        or recommendations["authorizing"] is not False
        or not isinstance(profiles_value, list)
        or not profiles_value
        or not isinstance(bounded, dict)
        or bounded
        != {
            "max_active_segments": MAX_LANES,
            "max_wave_nodes": sum(
                int(profile.get("nodes", 0)) for profile in profiles_value
                if isinstance(profile, dict)
            ),
            "r17_exclusive_and_last": True,
        }
        or not require_nonempty_string(
            recommendations["controller_consumption_state"],
            "recost controller-consumption state",
        )
        or not require_nonempty_string(
            recommendations["non_authorizing_reason"], "recost non-authorizing reason"
        )
    ):
        raise ValueError("recost bounded-concurrency recommendation differs")
    if mode == "sole-next-profile" and len(profiles_value) != 1:
        raise ValueError("recost sole-next-profile mode must contain one profile")
    if mode == "bounded-wave" and len(profiles_value) < 2:
        raise ValueError("recost bounded-wave mode must contain multiple profiles")
    if len(profiles_value) > MAX_LANES or int(bounded["max_wave_nodes"]) > MAX_NODES:
        raise ValueError("recost recommendation exceeds promoted concurrency ceilings")

    exact_keys = {
        "acceptance_criterion",
        "acceptance_policy",
        "athena_walltime",
        "build_manifest",
        "build_manifest_sha256",
        "case_id",
        "controller_walltime_max_seconds",
        "cpus_per_task",
        "estimated_storage_bytes",
        "executable",
        "executable_revision",
        "executable_sha256",
        "input_file",
        "input_revision",
        "input_sha256",
        "nodes",
        "output_layout",
        "segment",
        "parent_job_id",
        "parent_result",
        "parent_segment",
        "restart_file",
        "restart_file_sha256",
        "restart_time",
        "ranks_per_node",
        "recommendation_basis",
        "time_tlim_target",
        "walltime",
        "source_bundle",
        "source_bundle_sha256",
    }
    build_digest = build_manifest_inventory_sha256()
    profiles: list[dict[str, object]] = []
    by_case: dict[str, dict[str, object]] = {}
    for index, item in enumerate(profiles_value):
        if not isinstance(item, dict):
            raise ValueError(f"recost next profile {index} must be an object")
        require_exact_keys(item, exact_keys, f"recost next profile {index}")
        case_id = require_nonempty_string(item["case_id"], f"recost next profile {index} case")
        if case_id not in ALL_CASES or case_id in by_case or case_id == "R02":
            raise ValueError(f"recost next profile case {case_id} is unsupported or duplicated")
        segment_id = require_nonempty_string(item["segment"], f"recost next profile {case_id} segment")
        segment_index, start, target = parse_segment(segment_id, f"recost next profile {case_id}")
        target_alignment = (
            R12_CONTINUATION_ALIGNMENT
            if case_id == "R12"
            else STANDARD_CONTINUATION_ALIGNMENT
        )
        if (
            segment_index != state["next_indexes"][case_id]
            or decimal_value(item["time_tlim_target"], f"{case_id} recost target") != target
        ):
            raise ValueError(f"recost next profile {case_id} segment lineage differs")
        nodes, slurm_seconds, _, _ = validate_resource_values(
            case_id,
            item["nodes"],
            item["walltime"],
            item["athena_walltime"],
            Decimal(int(item["nodes"]) * walltime_seconds(item["walltime"], "recost walltime"))
            / Decimal(3600),
            f"recost next profile {case_id}",
        )
        if (
            item["ranks_per_node"] != RANKS_PER_NODE
            or item["cpus_per_task"] != CPUS_PER_TASK
            or item["controller_walltime_max_seconds"] != MAX_SEGMENT_SECONDS
            or item["output_layout"] != "rank-local"
            or item["executable"] != str(executable_path())
            or item["executable_revision"] != SOURCE_REVISION
            or item["executable_sha256"] != EXECUTABLE_SHA256
            or item["build_manifest"] != str(build_manifest_path())
            or item["build_manifest_sha256"] != build_digest
            or item["input_file"] != f"inputs/cgl_lf_paper/{EXPECTED_CASES[case_id][1]}"
            or item["input_revision"] != SOURCE_REVISION
            or item["input_sha256"] != EXPECTED_CASES[case_id][3]
            or not isinstance(item["estimated_storage_bytes"], int)
            or item["estimated_storage_bytes"] <= 0
            or not isinstance(item["recommendation_basis"], dict)
            or not item["recommendation_basis"]
            or not require_nonempty_string(item["acceptance_policy"], f"{case_id} acceptance policy")
            or not require_nonempty_string(
                item["acceptance_criterion"], f"{case_id} acceptance criterion"
            )
        ):
            raise ValueError(f"recost next profile {case_id} scientific/build profile differs")
        expected_bundle = (
            {
                "path": str(f115_source_bundle_path()),
                "sha256": R03_F115_SOURCE_BUNDLE_SHA256,
            }
            if case_id == R03 and segment_id == R03_F115_SEGMENT
            else {
                "path": promoted["source_bundle"]["path"],
                "sha256": promoted["source_bundle"]["sha256"],
            }
        )
        if (
            item["source_bundle"] != expected_bundle["path"]
            or item["source_bundle_sha256"] != expected_bundle["sha256"]
        ):
            raise ValueError(f"recost next profile {case_id} source bundle differs")
        lineage, active = state["lineages"][case_id]
        if active is not None:
            raise ValueError("drained recost recommendation collides with an active lane")
        parent_fields = (
            item["parent_job_id"],
            item["parent_result"],
            item["parent_segment"],
            item["restart_file"],
            item["restart_file_sha256"],
            item["restart_time"],
        )
        if not lineage:
            if any(value is not None for value in parent_fields) or start != 0:
                raise ValueError(f"fresh recost next profile {case_id} has parent ambiguity")
            if target != INITIAL_TARGETS[case_id]:
                raise ValueError(f"fresh recost next profile {case_id} target differs")
        else:
            parent = lineage[-1]
            restart = terminal_restart(parent)
            endpoint = lineage_endpoint(lineage)
            historical_r03_f115 = case_id == R03 and segment_id == R03_F115_SEGMENT
            if (
                abs(start - endpoint) > Decimal("0.000001")
                or target <= endpoint
                or (
                    not historical_r03_f115
                    and (target - start) % target_alignment != 0
                )
                or target - endpoint
                > (
                    Decimal("1.0")
                    if case_id == R03
                    else MAX_REVIEWED_INCREMENTS[case_id]
                )
                or item["parent_job_id"] != str(parent["manifest"]["job_id"])
                or item["parent_result"] != parent["manifest"]["accounting"]["result"]
                or item["parent_segment"] != parent["segment"]
                or item["restart_file"] != restart["restart_file"]
                or item["restart_file_sha256"] != restart["restart_sha256"]
                or decimal_value(item["restart_time"], f"{case_id} recost restart time")
                != endpoint
                or nodes != parent["nodes"]
            ):
                raise ValueError(f"continuation recost next profile {case_id} lineage differs")
        if case_id == "R12" and state["historical_inventory"]["R12"] and not lineage:
            if (
                segment_id != R12_FRESH_RERUN_SEGMENT
                or nodes != R12_FRESH_RERUN_NODES
                or nodes * int(item["ranks_per_node"]) != R12_FRESH_RERUN_RANKS
                or item["walltime"] != R12_FRESH_RERUN_WALLTIME
                or item["athena_walltime"] != R12_FRESH_RERUN_ATHENA_WALLTIME
                or any(value is not None for value in parent_fields)
                or start != 0
                or target != R12_FRESH_RERUN_TARGET
            ):
                raise ValueError("R12 fresh-rerun recost profile differs")
        profiles.append(item)
        by_case[case_id] = item

    if R17 in by_case:
        if len(profiles) != 1 or int(profiles[0]["nodes"]) != 8:
            raise ValueError("R17 recost recommendation must be exclusive on eight nodes")
        incomplete = [
            case_id for case_id in tuple(f"R{index:02d}" for index in range(2, 17))
            if not lineage_complete(state["lineages"][case_id][0])
        ]
        if incomplete:
            raise ValueError("R17 recost recommendation precedes accepted exact t=10 predecessors")
    elif recost.get("r17_readiness") is not None:
        raise ValueError("non-R17 recost retains an R17 readiness chain")
    return profiles, by_case


def validate_r12_fresh_rerun_transition(
    state: dict[str, object],
    profiles_by_case: dict[str, dict[str, object]],
) -> dict[str, object] | None:
    """Bind retained R12 s00 to inventory and require the exact fresh s01 rerun."""

    inventory = state.get("historical_inventory")
    if not isinstance(inventory, dict):
        raise ValueError("canonical state lacks historical-inventory classification")
    retained = inventory.get("R12")
    if not isinstance(retained, list):
        raise ValueError("canonical R12 historical inventory differs")
    profile = profiles_by_case.get("R12")
    if not retained:
        if profile is not None and profile.get("segment") == R12_FRESH_RERUN_SEGMENT:
            raise ValueError("R12 fresh rerun lacks its exact historical s00 inventory")
        return None
    if len(retained) != 1:
        raise ValueError("R12 historical inventory is not unique")
    historical = retained[0]
    manifest = historical.get("manifest")
    accounting = manifest.get("accounting") if isinstance(manifest, dict) else None
    if (
        not isinstance(manifest, dict)
        or historical.get("segment") != R12_HISTORICAL_PARTIAL_SEGMENT
        or str(manifest.get("job_id")) != R12_HISTORICAL_PARTIAL_JOB_ID
        or not isinstance(accounting, dict)
        or accounting.get("result") != "clean_partial"
    ):
        raise ValueError("R12 historical inventory identity differs")
    lineage, active = state["lineages"]["R12"]
    if active is not None:
        raise ValueError("drained R12 fresh-rerun transition has an active lane")
    transition_state = "fresh-rerun-required"
    fresh_lineage_root = None
    if lineage:
        root = lineage[0]
        root_manifest = root.get("manifest")
        command = root_manifest.get("command") if isinstance(root_manifest, dict) else None
        allocation = (
            root_manifest.get("allocation") if isinstance(root_manifest, dict) else None
        )
        if (
            root.get("segment") != R12_FRESH_RERUN_SEGMENT
            or root.get("index") != 1
            or root.get("start") != 0
            or root.get("target") != R12_FRESH_RERUN_TARGET
            or root.get("nodes") != R12_FRESH_RERUN_NODES
            or not isinstance(command, dict)
            or not isinstance(allocation, dict)
            or allocation.get("ranks_per_node") != RANKS_PER_NODE
            or root.get("nodes") * allocation["ranks_per_node"] != R12_FRESH_RERUN_RANKS
            or allocation.get("requested_walltime") != R12_FRESH_RERUN_WALLTIME
            or command.get("athena_walltime") != R12_FRESH_RERUN_ATHENA_WALLTIME
            or command.get("parent_segment") is not None
            or command.get("source_restart_file") is not None
            or command.get("restart_sha256") is not None
        ):
            raise ValueError("R12 scientific lineage does not begin at the fresh s01 rerun")
        transition_state = "fresh-rerun-established"
        fresh_lineage_root = root["evidence"]
    elif profile is not None:
        parent_fields = (
            profile.get("parent_job_id"),
            profile.get("parent_result"),
            profile.get("parent_segment"),
            profile.get("restart_file"),
            profile.get("restart_file_sha256"),
            profile.get("restart_time"),
        )
        _, start, target = parse_segment(
            str(profile.get("segment")), "R12 fresh-rerun profile"
        )
        if (
            profile.get("segment") != R12_FRESH_RERUN_SEGMENT
            or profile.get("nodes") != R12_FRESH_RERUN_NODES
            or profile.get("ranks_per_node") != RANKS_PER_NODE
            or profile.get("nodes") * profile["ranks_per_node"] != R12_FRESH_RERUN_RANKS
            or profile.get("walltime") != R12_FRESH_RERUN_WALLTIME
            or profile.get("athena_walltime") != R12_FRESH_RERUN_ATHENA_WALLTIME
            or any(value is not None for value in parent_fields)
            or start != 0
            or target != R12_FRESH_RERUN_TARGET
        ):
            raise ValueError("R12 fresh-rerun transition profile differs")
    return {
        "policy": R12_FRESH_RERUN_POLICY,
        "authorizing": False,
        "authorization_effect": "none",
        "continuation_authorized": False,
        "waiver_authorized": False,
        "transition_state": transition_state,
        "historical_inventory": {
            "manifest": historical["evidence"],
            "job_id": R12_HISTORICAL_PARTIAL_JOB_ID,
            "segment": R12_HISTORICAL_PARTIAL_SEGMENT,
            "result": "clean_partial",
            "classification": "inventory-only-non-authorizing",
        },
        "fresh_rerun_profile": profile if not lineage else None,
        "fresh_lineage_root": fresh_lineage_root,
    }


def scoped_measurement(
    info: dict[str, object],
) -> tuple[Decimal, dict[str, object]] | None:
    """Reproduce one scoped-v2 normalized measurement from canonical evidence."""

    case_id = str(info["case_id"])
    segment = str(info["segment"])
    manifest = info["manifest"]
    accounting = manifest.get("accounting")
    inspection = manifest.get("scientific_inspection")
    if not isinstance(accounting, dict) or not isinstance(inspection, dict):
        raise ValueError("scoped recost measurement lacks recorded accounting or inspection")
    job_id = str(manifest.get("job_id"))
    actual = recost_decimal_value(
        accounting.get("actual_node_hours"), f"scoped recost job {job_id} actual"
    )
    if actual <= 0:
        return None
    _, start, _ = parse_segment(segment, f"scoped recost job {job_id} segment")
    interval = decimal_value(
        inspection.get("final_time"), f"scoped recost job {job_id} final time"
    ) - Decimal(str(float(start)))
    if interval <= 0:
        raise ValueError(f"scoped recost job {job_id} measured interval is invalid")
    cells = resolution_cell_count(
        EXPECTED_CASES[case_id][2], f"scoped recost {case_id} resolution"
    )
    normalized_rate = actual / Decimal(cells) / interval
    return (
        normalized_rate,
        {
            "job_id": job_id,
            "case_id": case_id,
            "segment": segment,
            "actual_node_hours": format(actual, "f"),
            "observed_cells": cells,
            "observed_simulation_interval": format(interval, "f"),
            "normalized_node_hours_per_cell_per_simulation_time": format(
                normalized_rate, "f"
            ),
        },
    )


def expected_scoped_measurement_bases(
    state: dict[str, object],
) -> tuple[dict[str, object], dict[str, object]]:
    """Select the exact global and R12-local scoped-v2 measurement bases."""

    global_measurements: list[tuple[Decimal, dict[str, object]]] = []
    r12_measurements: list[tuple[Decimal, dict[str, object]]] = []
    for case_id in ALL_CASES:
        lineage, _ = state["lineages"][case_id]
        infos = lineage
        if case_id == "R12" and not infos:
            infos = state["historical_inventory"]["R12"]
        for info in infos:
            measurement = scoped_measurement(info)
            if measurement is None:
                continue
            if case_id == "R12":
                r12_measurements.append(measurement)
            else:
                global_measurements.append(measurement)
    if not global_measurements:
        raise ValueError("scoped recost lacks a positive non-R12 measurement basis")
    _, global_basis = max(
        global_measurements, key=lambda item: (item[0], str(item[1]["job_id"]))
    )
    _, r12_basis = (
        max(r12_measurements, key=lambda item: (item[0], str(item[1]["job_id"])))
        if r12_measurements
        else (Decimal("0"), global_basis)
    )
    return global_basis, r12_basis


def is_mandatory_fresh_r12_calibration_wave(
    state: dict[str, object], profiles: list[dict[str, object]]
) -> bool:
    """Return whether this drained wave contains the exact mandatory fresh R12 run."""

    inventory = state["historical_inventory"]["R12"]
    lineage, active = state["lineages"]["R12"]
    if len(inventory) != 1 or lineage or active is not None:
        return False
    historical = inventory[0]
    manifest = historical.get("manifest")
    accounting = manifest.get("accounting") if isinstance(manifest, dict) else None
    if (
        not isinstance(manifest, dict)
        or not isinstance(accounting, dict)
        or str(manifest.get("job_id")) != R12_HISTORICAL_PARTIAL_JOB_ID
        or historical.get("segment") != R12_HISTORICAL_PARTIAL_SEGMENT
        or accounting.get("result") != "clean_partial"
    ):
        return False
    r12_profiles = [profile for profile in profiles if profile.get("case_id") == "R12"]
    if len(r12_profiles) != 1:
        return False
    profile = r12_profiles[0]
    parent_fields = (
        profile.get("parent_job_id"),
        profile.get("parent_result"),
        profile.get("parent_segment"),
        profile.get("restart_file"),
        profile.get("restart_file_sha256"),
        profile.get("restart_time"),
    )
    return (
        profile.get("segment") == R12_FRESH_RERUN_SEGMENT
        and profile.get("nodes") == R12_FRESH_RERUN_NODES
        and profile.get("ranks_per_node") == RANKS_PER_NODE
        and profile.get("nodes") * profile.get("ranks_per_node") == R12_FRESH_RERUN_RANKS
        and profile.get("walltime") == R12_FRESH_RERUN_WALLTIME
        and profile.get("athena_walltime") == R12_FRESH_RERUN_ATHENA_WALLTIME
        and decimal_value(
            profile.get("time_tlim_target"), "fresh R12 calibration target"
        )
        == R12_FRESH_RERUN_TARGET
        and all(value is None for value in parent_fields)
    )


def validate_recost_budget(
    recost: dict[str, object],
    state: dict[str, object],
    profiles: list[dict[str, object]],
    cases: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Validate the recost publication as the sole current projection authority."""

    budget = recost.get("budget")
    provenance = recost.get("provenance")
    storage = recost.get("storage")
    if not isinstance(budget, dict) or not isinstance(provenance, dict):
        raise ValueError("recost budget/provenance must be objects")
    if recost_json_sha256(budget) != provenance.get("computed_projection_sha256"):
        raise ValueError("recost projection digest differs")
    require_exact_keys(budget, SCOPED_BUDGET_KEYS, "scoped-v2 recost budget")
    if budget.get("method") != NODE_HOUR_PROJECTION_METHOD:
        raise ValueError("recost node-hour projection method is not scoped-v2")

    reserved_by_case: dict[str, Decimal] = {}
    for profile in profiles:
        case_id = str(profile["case_id"])
        value = (
            Decimal(
                int(profile["nodes"])
                * walltime_seconds(profile["walltime"], "recost profile")
            )
            / Decimal(3600)
        )
        reserved_by_case[case_id] = reserved_by_case.get(case_id, Decimal("0")) + value
    reserved = sum(reserved_by_case.values(), Decimal("0"))
    actual = state["actual_node_hours"]
    committed = actual + reserved
    envelope = STAGE_I_BUDGET_NODE_HOURS
    project = PROJECT_BUDGET_NODE_HOURS
    global_basis, r12_basis = expected_scoped_measurement_bases(state)
    retained_global_basis = budget.get("measurement_basis")
    if not isinstance(retained_global_basis, dict):
        raise ValueError("scoped-v2 global measurement basis must be an object")
    require_exact_keys(
        retained_global_basis, SCOPED_MEASUREMENT_BASIS_KEYS, "scoped-v2 global measurement basis"
    )
    if retained_global_basis != global_basis:
        raise ValueError("scoped-v2 global measurement basis differs from canonical evidence")
    breakdown = budget.get("case_breakdown")
    if not isinstance(breakdown, dict) or set(breakdown) != set(ALL_CASES):
        raise ValueError("scoped-v2 recost case breakdown differs from the frozen matrix")
    remaining = Decimal("0")
    for case_id in ALL_CASES:
        item = breakdown[case_id]
        if not isinstance(item, dict):
            raise ValueError(f"scoped-v2 recost {case_id} breakdown must be an object")
        require_exact_keys(item, SCOPED_CASE_BREAKDOWN_KEYS, f"scoped-v2 recost {case_id}")
        expected_basis = r12_basis if case_id == "R12" else global_basis
        retained_basis = item.get("projection_measurement_basis")
        if not isinstance(retained_basis, dict):
            raise ValueError(f"scoped-v2 recost {case_id} measurement basis must be an object")
        require_exact_keys(
            retained_basis,
            SCOPED_MEASUREMENT_BASIS_KEYS,
            f"scoped-v2 recost {case_id} measurement basis",
        )
        if retained_basis != expected_basis:
            raise ValueError(f"scoped-v2 recost {case_id} measurement basis differs")
        cells = resolution_cell_count(
            EXPECTED_CASES[case_id][2], f"scoped-v2 recost {case_id} resolution"
        )
        lineage, _ = state["lineages"][case_id]
        progress = min(
            Decimal("1"), max(Decimal("0"), lineage_endpoint(lineage) / Decimal("10"))
        )
        remaining_time = Decimal("10") * (Decimal("1") - progress)
        observed_rate = recost_decimal_value(
            expected_basis["normalized_node_hours_per_cell_per_simulation_time"],
            f"scoped-v2 recost {case_id} normalized measurement rate",
        )
        observed_projection = observed_rate * Decimal(cells) * remaining_time
        authorized_reserved = reserved_by_case.get(case_id, Decimal("0"))
        projected_remaining = max(observed_projection, authorized_reserved)
        if (
            recost_decimal_value(
                item.get("matrix_full_case_node_hours_reference_only"),
                f"scoped-v2 recost {case_id} matrix reference",
            )
            != cases[case_id]["_estimated_node_hours"]
            or recost_decimal_value(
                item.get("authenticated_progress_fraction"),
                f"scoped-v2 recost {case_id} progress",
            )
            != progress
            or recost_decimal_value(
                item.get("remaining_simulation_time"),
                f"scoped-v2 recost {case_id} remaining time",
            )
            != remaining_time
            or item.get("projected_cells") != str(cells)
            or recost_decimal_value(
                item.get("observed_rate_projected_remaining_node_hours"),
                f"scoped-v2 recost {case_id} observed projection",
            )
            != observed_projection
            or recost_decimal_value(
                item.get("authorized_profile_reserved_node_hours"),
                f"scoped-v2 recost {case_id} authorized reservation",
            )
            != authorized_reserved
            or recost_decimal_value(
                item.get("projected_remaining_node_hours"),
                f"scoped-v2 recost {case_id} projected remainder",
            )
            != projected_remaining
        ):
            raise ValueError(f"scoped-v2 recost {case_id} projection arithmetic differs")
        remaining += projected_remaining
    total = actual + remaining
    margin = envelope - total
    if (
        recost_decimal_value(budget.get("actual_stage_i_node_hours"), "recost actual use")
        != actual
        or recost_decimal_value(
            budget.get("authorized_wave_reserved_node_hours"), "recost authorized-wave use"
        )
        != reserved
        or recost_decimal_value(
            budget.get("actual_plus_authorized_wave_node_hours"), "recost committed use"
        )
        != committed
        or recost_decimal_value(
            budget.get("computed_remaining_stage_i_node_hours"), "recost computed remaining use"
        )
        != remaining
        or recost_decimal_value(
            budget.get("computed_stage_i_total_node_hours"), "recost computed total use"
        )
        != total
        or recost_decimal_value(
            budget.get("promoted_stage_i_envelope_node_hours"), "recost Stage I envelope"
        )
        != envelope
        or recost_decimal_value(
            budget.get("project_ceiling_node_hours"), "recost project ceiling"
        )
        != project
        or recost_decimal_value(
            budget.get("computed_stage_i_margin_node_hours"),
            "recost computed Stage I margin",
            allow_negative=True,
        )
        != margin
    ):
        raise ValueError("recost projection differs from canonical accounting or ceilings")
    if committed > envelope or committed > project:
        raise ValueError("recost authorized wave exceeds a controller budget ceiling")
    if total > project:
        raise ValueError(
            "promoted recost credible remaining-campaign projection exceeds project ceiling"
        )
    if total > envelope and not is_mandatory_fresh_r12_calibration_wave(state, profiles):
        raise ValueError(
            "above-envelope scoped-v2 recost is allowed only for the exact mandatory "
            "fresh R12 calibration wave"
        )
    if not isinstance(storage, dict):
        raise ValueError("recost storage projection must be an object")
    for key in (
        "available_bytes",
        "retained_stage_i_bytes",
        "required_safety_bytes",
        "projected_authorized_wave_growth_bytes",
        "headroom_after_authorized_wave_and_safety_bytes",
    ):
        if isinstance(storage.get(key), bool) or not isinstance(storage.get(key), int):
            raise ValueError("recost storage projection contains a non-integer value")
    if (
        storage["headroom_after_authorized_wave_and_safety_bytes"] < 0
        or available_storage_bytes(CANONICAL_ROOT)
        < storage["required_safety_bytes"] + storage["projected_authorized_wave_growth_bytes"]
    ):
        raise ValueError("recost storage projection no longer has required live headroom")
    return budget


def campaign_projection(
    cases: dict[str, dict[str, object]],
    state: dict[str, object],
    additional_reserved: dict[str, Decimal] | None = None,
) -> dict[str, object]:
    """Independently project the full remaining campaign from matrix estimates."""

    reserved_by_case: dict[str, Decimal] = {}
    for info in state["active_infos"]:
        case_id = str(info["case_id"])
        reserved_by_case[case_id] = reserved_by_case.get(case_id, Decimal("0")) + info["reserved"]
    for case_id, value in (additional_reserved or {}).items():
        reserved_by_case[case_id] = reserved_by_case.get(case_id, Decimal("0")) + value
    remaining = Decimal("0")
    breakdown = {}
    for case_id in ALL_CASES:
        estimate = cases[case_id]["_estimated_node_hours"]
        if not isinstance(estimate, Decimal):
            raise ValueError(f"matrix {case_id} estimate is not authenticated")
        endpoint = lineage_endpoint(state["lineages"][case_id][0])
        progress = min(Decimal("1"), max(Decimal("0"), endpoint / Decimal("10")))
        scaled = estimate * (Decimal("1") - progress)
        reserved = reserved_by_case.get(case_id, Decimal("0"))
        projected = max(scaled, reserved)
        remaining += projected
        breakdown[case_id] = {
            "matrix_full_case_node_hours": decimal_text(estimate),
            "authenticated_progress_fraction": decimal_text(progress),
            "scaled_remaining_node_hours": decimal_text(scaled),
            "committed_profile_reserved_node_hours": decimal_text(reserved),
            "projected_remaining_node_hours": decimal_text(projected),
        }
    total = state["actual_node_hours"] + remaining
    return {
        "method": (
            "authenticated ledger actuals plus frozen-matrix full-case estimates scaled "
            "by byte-authenticated terminal-lineage progress and floored by committed "
            "active/planned reservations"
        ),
        "actual_stage_i_node_hours": decimal_text(state["actual_node_hours"]),
        "computed_remaining_stage_i_node_hours": decimal_text(remaining),
        "computed_stage_i_total_node_hours": decimal_text(total),
        "promoted_stage_i_envelope_node_hours": decimal_text(STAGE_I_BUDGET_NODE_HOURS),
        "project_ceiling_node_hours": decimal_text(PROJECT_BUDGET_NODE_HOURS),
        "computed_stage_i_margin_node_hours": decimal_text(STAGE_I_BUDGET_NODE_HOURS - total),
        "case_breakdown": breakdown,
    }


def projection_total(projection: dict[str, object]) -> Decimal:
    """Return the exact independently computed total from one projection."""

    return decimal_value(
        projection.get("computed_stage_i_total_node_hours"), "campaign projected total"
    )


def require_projection_within_budget(projection: dict[str, object], label: str) -> None:
    """Require a credible full-campaign projection to fit both ceilings."""

    total = projection_total(projection)
    if total > STAGE_I_BUDGET_NODE_HOURS or total > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError(f"{label} credible remaining-campaign projection exceeds budget")


def validate_authority_projection(
    authority: dict[str, object], projection: dict[str, object]
) -> None:
    """Bind the independent allocation authority to reproduced campaign arithmetic."""

    if authority["artifact"].get("projection") != projection:
        raise ValueError("allocation authority projection differs from authenticated campaign state")
    require_projection_within_budget(projection, "allocation authority")


def validate_active_reviewed_profile(
    case_id: str,
    lineage: list[dict[str, object]],
    active: dict[str, object] | None,
    profile: dict[str, object],
) -> None:
    """Require an active lane to be the exact next segment in its reviewed profile."""

    if active is None:
        return
    endpoint = lineage_endpoint(lineage)
    expected_target = min(Decimal("10"), endpoint + profile["next_increment"])
    manifest = active["manifest"]
    if (
        active["index"] != len(lineage)
        or abs(active["start"] - endpoint) > Decimal("0.000001")
        or active["target"] != expected_target
        or active["nodes"] != int(profile["nodes"])
        or manifest["allocation"]["requested_walltime"] != profile["walltime"]
        or manifest["command"]["athena_walltime"] != profile["athena_walltime"]
    ):
        raise ValueError(f"active {case_id} lane differs from its fresh reviewed profile")


def require_immutable_publication_evidence(
    evidence: dict[str, object], digest: str, label: str
) -> None:
    """Require one exact digest-bound immutable single-link publication."""

    if (
        evidence.get("sha256") != digest
        or evidence.get("mode") != "0444"
        or evidence.get("nlink") != 1
    ):
        raise ValueError(f"{label} is not the exact immutable published artifact")


def validate_declared_publication_binding(
    binding: object,
    *,
    expected: Path,
    digest: str,
    mode: str,
    label: str,
) -> None:
    """Require one audit-declared path/hash/mode/link binding."""

    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding must be an object")
    require_exact_keys(binding, ("path", "sha256", "mode", "links"), f"{label} binding")
    declared = Path(str(binding["path"]))
    declared = (
        expected_path(declared.as_posix())
        if not declared.is_absolute()
        else normalized_path(declared)
    )
    if binding["sha256"] != digest or binding["mode"] != mode or binding["links"] != 1:
        raise ValueError(f"{label} binding differs")
    if declared != normalized_path(expected):
        raise ValueError(f"{label} path differs")


def validate_f118_current_source_authority(binding: object) -> dict[str, object]:
    """Authenticate production F118 and immutable historical F116/F115."""

    if not isinstance(binding, dict):
        raise ValueError("recost F118 current-source authority binding must be an object")
    require_exact_keys(
        binding,
        (
            "checkpoint", "evidence", "provenance_review", "plasma_review",
            "publication_audit", "final_source_bundle",
        ),
        "recost F118 current-source authority binding",
    )
    if binding["checkpoint"] != "F-118":
        raise ValueError("recost current-source authority is not F-118")

    def recost_binding(value: object, expected: Path, label: str) -> str:
        if not isinstance(value, dict):
            raise ValueError(f"{label} binding must be an object")
        require_exact_keys(value, ("path", "sha256"), f"{label} binding")
        relative = Path(require_nonempty_string(value["path"], f"{label} path"))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"{label} path must be canonical-root relative")
        if expected_path(relative.as_posix()) != normalized_path(expected):
            raise ValueError(f"{label} path differs")
        return require_sha256(value["sha256"], f"{label} SHA-256")

    artifact_path = f118_current_source_authority_path()
    provenance_review_path = f118_provenance_security_review_path()
    plasma_review_path = f118_plasma_scientific_review_path()
    audit_path = f118_publication_audit_path()
    artifact_sha256 = recost_binding(binding["evidence"], artifact_path, "F118 evidence")
    provenance_review_sha256 = recost_binding(
        binding["provenance_review"], provenance_review_path, "F118 provenance review"
    )
    plasma_review_sha256 = recost_binding(
        binding["plasma_review"], plasma_review_path, "F118 plasma review"
    )
    audit_sha256 = recost_binding(binding["publication_audit"], audit_path, "F118 audit")

    evidence, evidence_file = read_json_file(artifact_path, "F118 current-source authority")
    provenance_review, provenance_review_file = read_json_file(
        provenance_review_path, "F118 provenance/security review"
    )
    plasma_review, plasma_review_file = read_json_file(
        plasma_review_path, "F118 plasma/scientific review"
    )
    audit, audit_file = read_json_file(audit_path, "F118 publication audit")
    for retained, digest, label in (
        (evidence_file, artifact_sha256, "F118 current-source authority"),
        (provenance_review_file, provenance_review_sha256, "F118 provenance/security review"),
        (plasma_review_file, plasma_review_sha256, "F118 plasma/scientific review"),
        (audit_file, audit_sha256, "F118 publication audit"),
    ):
        require_immutable_publication_evidence(retained, digest, label)

    require_exact_keys(
        evidence,
        (
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        ),
        "F118 current-source authority",
    )
    generated = require_utc(evidence["generated_utc"], "F118 generation time")
    scope = evidence["scope"]
    predecessors = evidence["predecessor_authorities"]
    implementation = evidence["implementation"]
    catalog = evidence["source_archive_catalog"]
    if not all(isinstance(value, dict) for value in (scope, predecessors, implementation, catalog)):
        raise ValueError("F118 scope, predecessors, implementation, and catalog must be objects")
    require_exact_keys(
        scope, ("relationship", "summary", "preserves", "does_not_authorize"), "F118 scope"
    )
    require_exact_keys(predecessors, ("historical_f116",), "F118 predecessor authorities")
    require_exact_keys(
        implementation,
        (
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
        ),
        "F118 implementation",
    )
    require_exact_keys(catalog, ("before", "after"), "F118 source-archive catalog")
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != "F-118"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or generated > current_utc() + FUTURE_SKEW
        or scope["relationship"] != "current-source-selection-only-supersession"
        or not require_nonempty_string(scope["summary"], "F118 scope summary")
        or scope["preserves"] != F118_SCOPE_PRESERVES
        or scope["does_not_authorize"] != F118_SCOPE_DOES_NOT_AUTHORIZE
        or evidence["authorization"] != F118_AUTHORIZATION
        or evidence["validation"] != F118_VALIDATION_CLAIMS
        or evidence["publication_requirements"] != F118_PUBLICATION_REQUIREMENTS
    ):
        raise ValueError("F118 identity, source-selection-only scope, or authority differs")

    historical_f116_bindings = predecessors["historical_f116"]
    if not isinstance(historical_f116_bindings, dict):
        raise ValueError("F118 historical F116 predecessor binding must be an object")
    require_exact_keys(
        historical_f116_bindings,
        ("evidence", "publication_audit", "provenance_review", "plasma_review"),
        "F118 historical F116 predecessor binding",
    )
    historical_f116_expected = {
        "evidence": f116_current_source_authority_path(),
        "publication_audit": f116_publication_audit_path(),
        "provenance_review": f116_provenance_security_review_path(),
        "plasma_review": f116_plasma_scientific_review_path(),
    }
    historical_f116_loaded: dict[str, tuple[dict[str, object], dict[str, object]]] = {}
    historical_f116_digests: dict[str, str] = {}
    for key, path in historical_f116_expected.items():
        value = historical_f116_bindings[key]
        if not isinstance(value, dict):
            raise ValueError(f"F118 historical F116 {key} binding must be an object")
        require_exact_keys(value, ("path", "sha256"), f"F118 historical F116 {key} binding")
        digest = require_sha256(value["sha256"], f"F118 historical F116 {key} SHA-256")
        if (
            value["path"] != path.relative_to(CANONICAL_ROOT).as_posix()
            or (
                CANONICAL_ROOT
                == Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
                and digest != F116_CANONICAL_SHA256[key]
            )
        ):
            raise ValueError(f"F118 historical F116 {key} binding differs")
        retained_value, retained_file = read_json_file(path, f"historical F116 {key}")
        require_immutable_publication_evidence(retained_file, digest, f"historical F116 {key}")
        historical_f116_loaded[key] = retained_value, retained_file
        historical_f116_digests[f"{key}_sha256"] = digest
    historical_f116_evidence = historical_f116_loaded["evidence"][0]
    require_exact_keys(
        historical_f116_evidence,
        (
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        ),
        "historical F116 evidence",
    )
    historical_f116_predecessors = historical_f116_evidence["predecessor_authorities"]
    historical_f116_implementation = historical_f116_evidence["implementation"]
    if (
        not isinstance(historical_f116_predecessors, dict)
        or not isinstance(historical_f116_implementation, dict)
    ):
        raise ValueError("historical F116 predecessors or implementation differ")
    require_exact_keys(
        historical_f116_predecessors, ("historical_f115",),
        "historical F116 predecessor authorities",
    )
    require_exact_keys(
        historical_f116_implementation,
        ("publisher", "committed_tools", "intermediate_36140_bundle", "current_source_bundle"),
        "historical F116 implementation",
    )
    if (
        historical_f116_evidence["schema_version"] != 1
        or historical_f116_evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or historical_f116_evidence["checkpoint"] != "F-116"
        or historical_f116_evidence["execution_epoch"] != EXECUTION_EPOCH
        or historical_f116_evidence["authorization"] != F116_AUTHORIZATION
    ):
        raise ValueError("historical F116 evidence identity or authority differs")
    historical = historical_f116_predecessors["historical_f115"]
    if not isinstance(historical, dict):
        raise ValueError("F116 historical F115 predecessor binding must be an object")
    require_exact_keys(
        historical,
        ("evidence", "publication_audit", "provenance_review", "plasma_review"),
        "F116 historical F115 predecessor binding",
    )
    historical_expected = {
        "evidence": (r03_f115_path(), R03_F115_SHA256),
        "publication_audit": (
            r03_f115_publication_audit_path(),
            R03_F115_PUBLICATION_AUDIT_SHA256,
        ),
        "provenance_review": (
            r03_f115_provenance_security_review_path(),
            R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256,
        ),
        "plasma_review": (
            r03_f115_plasma_scientific_review_path(),
            R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
        ),
    }
    historical_digests = {}
    for key, (path, digest) in historical_expected.items():
        value = historical[key]
        if not isinstance(value, dict):
            raise ValueError(f"F116 historical F115 {key} binding must be an object")
        require_exact_keys(value, ("path", "sha256"), f"F116 historical F115 {key} binding")
        relative = path.relative_to(CANONICAL_ROOT).as_posix()
        if value != {"path": relative, "sha256": digest}:
            raise ValueError(f"F116 historical F115 {key} binding differs")
        _, retained = read_stable_regular_file(path, f"F116 historical F115 {key}")
        require_immutable_publication_evidence(retained, digest, f"F116 historical F115 {key}")
        historical_digests[f"{key}_sha256"] = digest

    historical_f116_audit = historical_f116_loaded["publication_audit"][0]
    require_exact_keys(
        historical_f116_audit,
        (
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f115_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        ),
        "historical F116 publication audit",
    )
    historical_f116_reviews = historical_f116_audit["independent_reviews"]
    if not isinstance(historical_f116_reviews, dict):
        raise ValueError("historical F116 publication reviews must be an object")
    require_exact_keys(
        historical_f116_reviews,
        (
            "reviews_bind_exact_published_f116_sha256",
            "provenance_security", "plasma_scientific_continuation",
        ),
        "historical F116 publication reviews",
    )
    if (
        historical_f116_audit["checkpoint"] != "F-116"
        or historical_f116_audit["execution_epoch"] != EXECUTION_EPOCH
        or historical_f116_audit["historical_f115_authority"] != historical_digests
        or historical_f116_audit["authority_and_enforcement"] != F116_AUTHORIZATION
        or historical_f116_reviews["reviews_bind_exact_published_f116_sha256"]
        != historical_f116_digests["evidence_sha256"]
    ):
        raise ValueError("historical F116 publication audit differs")
    for key, audit_key in (
        ("provenance_review", "provenance_security"),
        ("plasma_review", "plasma_scientific_continuation"),
    ):
        validate_declared_publication_binding(
            historical_f116_reviews[audit_key],
            expected=historical_f116_expected[key],
            digest=historical_f116_digests[f"{key}_sha256"],
            mode="0444",
            label=f"historical F116 {key}",
        )

    def bundle_declaration(value: object, *, current: bool, label: str) -> dict[str, object]:
        if not isinstance(value, dict):
            raise ValueError(f"{label} must be an object")
        keys = {
            "path", "sha256", "complete_history", "head", "advertised_tip",
            "verified_revisions", "selected_as_current",
        }
        keys |= {"candidate_path", "subject"} if current else {"role"}
        require_exact_keys(value, keys, label)
        relative = Path(require_nonempty_string(value["path"], f"{label} path"))
        head = require_git_revision(value["head"], f"{label} head")
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative.parent.as_posix() != "source-archives"
            or relative.name != f"athenak-feature-cgl-through-{head[:9]}.bundle"
        ):
            raise ValueError(f"{label} path does not bind its head")
        revisions_value = value["verified_revisions"]
        if not isinstance(revisions_value, list) or not revisions_value:
            raise ValueError(f"{label} verified revisions must be a nonempty list")
        revisions = [
            require_git_revision(revision, f"{label} verified revision")
            for revision in revisions_value
        ]
        tip = value["advertised_tip"]
        if not isinstance(tip, dict):
            raise ValueError(f"{label} advertised tip must be an object")
        require_exact_keys(tip, ("revision", "name"), f"{label} advertised tip")
        tip_name = require_nonempty_string(tip["name"], f"{label} advertised tip name")
        branch_tip = (
            tip_name.startswith("refs/heads/")
            and len(tip_name) > len("refs/heads/")
            and not any(character.isspace() for character in tip_name)
        )
        if (
            len(set(revisions)) != len(revisions)
            or value["complete_history"] is not True
            or tip["revision"] != head
            or not (tip_name == "HEAD" or branch_tip)
        ):
            raise ValueError(f"{label} history or advertised tip differs")
        if current:
            candidate = Path(
                require_nonempty_string(value["candidate_path"], f"{label} candidate path")
            )
            require_nonempty_string(value["subject"], f"{label} subject")
            if not candidate.is_absolute() or value["selected_as_current"] is not True:
                raise ValueError(f"{label} current-selection declaration differs")
        elif (
            value["selected_as_current"] is not False
            or value["role"] != "retained-non-current-bridge"
        ):
            raise ValueError(f"{label} retained-bridge declaration differs")
        return {
            "path": relative.as_posix(),
            "path_object": expected_path(relative.as_posix()),
            "sha256": require_sha256(value["sha256"], f"{label} SHA-256"),
            "head": head,
            "tip_name": tip_name,
            "verified_revisions": revisions,
        }

    historical_f116_current = bundle_declaration(
        historical_f116_implementation["current_source_bundle"],
        current=True,
        label="historical F116 current source bundle",
    )
    bridge = bundle_declaration(
        implementation["intermediate_36140_bundle"],
        current=False,
        label="F118 retained bridge bundle",
    )
    final = bundle_declaration(
        implementation["current_source_bundle"],
        current=True,
        label="F118 current source bundle",
    )
    predecessor = implementation["predecessor_current_source_bundle"]
    if not isinstance(predecessor, dict):
        raise ValueError("F118 predecessor current source bundle must be an object")
    expected_predecessor = dict(
        historical_f116_implementation["current_source_bundle"]
    )
    expected_predecessor.pop("candidate_path")
    expected_predecessor["selected_as_current"] = False
    expected_predecessor["role"] = "retained-non-current-predecessor"
    require_exact_keys(
        predecessor,
        tuple(expected_predecessor),
        "F118 predecessor current source bundle",
    )
    if predecessor != expected_predecessor:
        raise ValueError("F118 predecessor current bundle differs from immutable F116")
    if (
        bridge["path"] != F116_BRIDGE_SOURCE_BUNDLE
        or bridge["head"] != F116_BRIDGE_REVISION
        or bridge["sha256"] != F116_BRIDGE_SHA256
        or bridge["tip_name"] != "refs/heads/feature/cgl-landau-fluid"
        or final["path_object"] != promoted_source_bundle_path(str(final["head"]))
    ):
        raise ValueError("F118 bridge or final source-bundle identity differs")

    promoted_revision = str(final["head"])
    recost_bundle_binding = binding["final_source_bundle"]
    if not isinstance(recost_bundle_binding, dict):
        raise ValueError("recost F118 final source bundle binding must be an object")
    require_exact_keys(
        recost_bundle_binding,
        ("path", "sha256", "verified_revisions"),
        "recost F118 final source bundle binding",
    )
    if recost_bundle_binding != {
        "path": final["path"],
        "sha256": final["sha256"],
        "verified_revisions": final["verified_revisions"],
    }:
        raise ValueError("recost current-source authority does not bind the F118 final bundle")

    repository = controller_repository_path()
    head = git_read_only(repository, ["rev-parse", "--verify", "HEAD^{commit}"])
    try:
        live_head = head.stdout.decode("ascii").strip()
    except UnicodeDecodeError as error:
        raise ValueError("live promoted source-authority HEAD is not ASCII") from error
    if head.returncode or live_head != promoted_revision:
        raise ValueError("F118 promoted revision is not live repository HEAD")

    tools = implementation["committed_tools"]
    if not isinstance(tools, list) or len(tools) != len(F118_REQUIRED_TOOLS):
        raise ValueError("F118 committed-tool vector differs")
    if [item.get("path") for item in tools if isinstance(item, dict)] != sorted(
        F118_REQUIRED_TOOLS
    ):
        raise ValueError("F118 committed-tool vector is not in production deterministic order")
    retained_tools: dict[str, dict[str, str]] = {}
    for item in tools:
        if not isinstance(item, dict):
            raise ValueError("F118 committed-tool record must be an object")
        require_exact_keys(item, ("path", "revision", "sha256", "mode"), "F118 committed tool")
        path = require_nonempty_string(item["path"], "F118 committed-tool path")
        if path in retained_tools or path not in F118_REQUIRED_TOOLS:
            raise ValueError("F118 committed-tool vector contains unsupported or duplicate paths")
        identity = {
            "path": path,
            "revision": require_git_revision(item["revision"], f"F118 committed tool {path} revision"),
            "sha256": require_sha256(item["sha256"], f"F118 committed tool {path} SHA-256"),
            "mode": require_nonempty_string(item["mode"], f"F118 committed tool {path} mode"),
        }
        if identity["revision"] != promoted_revision or identity["mode"] != F118_REQUIRED_TOOLS[path]:
            raise ValueError(f"F118 committed tool identity differs: {path}")
        committed = git_read_only(repository, ["show", f"{promoted_revision}:{path}"])
        _, live = read_stable_regular_file(repository / path, f"F118 live committed tool {path}")
        tree = git_read_only(repository, ["ls-tree", promoted_revision, "--", path])
        expected_git_mode = "100755" if identity["mode"] == "0755" else "100644"
        try:
            tree_line = tree.stdout.decode("ascii").strip()
        except UnicodeDecodeError as error:
            raise ValueError(f"F118 committed tool mode is not ASCII: {path}") from error
        if (
            committed.returncode
            or sha256_bytes(committed.stdout) != identity["sha256"]
            or live["sha256"] != identity["sha256"]
            or live["mode"] != identity["mode"]
            or tree.returncode
            or not tree_line.startswith(f"{expected_git_mode} blob ")
            or not tree_line.endswith(f"\t{path}")
        ):
            raise ValueError(f"F118 committed/live tool bytes or mode differ: {path}")
        retained_tools[path] = identity
    if set(retained_tools) != set(F118_REQUIRED_TOOLS):
        raise ValueError("F118 committed-tool vector differs")
    publisher = implementation["publisher"]
    if publisher != retained_tools[F116_PUBLISHER_RELATIVE]:
        raise ValueError("F118 publisher binding differs from the committed-tool vector")

    helper = retained_tools["scripts/frontier/cgl_lf_stage_i.py"]
    generator = retained_tools["scripts/frontier/cgl_lf_stage_i_recost.py"]
    matrix_relative = "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
    committed_matrix = git_read_only(repository, ["show", f"{promoted_revision}:{matrix_relative}"])
    if committed_matrix.returncode or sha256_bytes(committed_matrix.stdout) != MATRIX_SHA256:
        raise ValueError("F118 final source does not preserve the frozen Stage I matrix")
    matrix = {"path": matrix_relative, "revision": promoted_revision, "sha256": MATRIX_SHA256}

    required_revisions = set(F116_PRODUCTION_REQUIRED_REVISIONS) | {promoted_revision}
    if not required_revisions.issubset(final["verified_revisions"]):
        raise ValueError("F118 final source bundle omits required retained history")
    _, bridge_evidence = read_stable_regular_file(
        bridge["path_object"], "F116 retained bridge source bundle"
    )
    _, bundle_evidence = read_stable_regular_file(
        final["path_object"], "F116 final source bundle"
    )
    bridge_helper = git_read_only(
        repository, ["show", f"{F116_BRIDGE_REVISION}:scripts/frontier/cgl_lf_stage_i.py"]
    )
    if bridge_helper.returncode:
        raise ValueError("F116 bridge revision lacks the Stage I helper")
    bridge_validation = validate_source_bundle_coverage(
        bridge["path_object"],
        bundle_sha256=str(bridge["sha256"]),
        controller_revision=F116_BRIDGE_REVISION,
        controller_sha256=sha256_bytes(bridge_helper.stdout),
        required_revisions=bridge["verified_revisions"],
        expected_evidence=bridge_evidence,
        allow_branch_ref=True,
        label="F116 retained bridge source bundle",
    )
    bundle_validation = validate_source_bundle_coverage(
        final["path_object"],
        bundle_sha256=str(final["sha256"]),
        controller_revision=promoted_revision,
        controller_sha256=helper["sha256"],
        required_revisions=final["verified_revisions"],
        expected_evidence=bundle_evidence,
        allow_branch_ref=True,
        label="F116 final source bundle",
    )
    if (
        bridge_validation["advertised_heads"]
        != [{"revision": bridge["head"], "name": bridge["tip_name"]}]
        or bundle_validation["advertised_heads"]
        != [{"revision": final["head"], "name": final["tip_name"]}]
    ):
        raise ValueError("F116 source-bundle advertised-tip declaration differs")

    before_catalog = catalog["before"]
    after_catalog = catalog["after"]
    if not isinstance(before_catalog, dict) or not isinstance(after_catalog, dict):
        raise ValueError("F116 before/after source-archive catalogs must be objects")
    require_exact_keys(
        before_catalog,
        (
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "predecessor_current_source_bundle_listed_exactly_once",
            "final_bundle_listed", "corrupt_c7_listed", "historical_f115_preserved",
        ),
        "F118 predecessor source-archive catalog",
    )
    require_exact_keys(
        after_catalog,
        (
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "predecessor_current_source_bundle_listed_exactly_once",
            "final_bundle_listed_exactly_once", "corrupt_c7_listed",
            "historical_f115_preserved", "historical_f116_preserved",
            "all_prior_checksum_entries_preserved", "sole_current_source_bundle",
        ),
        "F118 current source-archive catalog",
    )
    for value, label in (
        (before_catalog["readme_sha256"], "F116 predecessor README SHA-256"),
        (before_catalog["sha256sums_sha256"], "F116 predecessor SHA256SUMS SHA-256"),
        (after_catalog["readme_sha256"], "F116 current README SHA-256"),
        (after_catalog["sha256sums_sha256"], "F116 current SHA256SUMS SHA-256"),
    ):
        require_sha256(value, label)
    readme_path = expected_path("source-archives/README.md")
    sums_path = expected_path("source-archives/SHA256SUMS")
    _, readme_evidence = read_stable_regular_file(readme_path, "F116 source-archive README")
    sums_payload, sums_evidence = read_stable_regular_file(sums_path, "F116 source-archive SHA256SUMS")
    if (
        before_catalog["bridge_listed_exactly_once"] is not True
        or before_catalog["predecessor_current_source_bundle_listed_exactly_once"] is not True
        or before_catalog["final_bundle_listed"] is not False
        or before_catalog["corrupt_c7_listed"] is not False
        or before_catalog["historical_f115_preserved"] is not True
        or after_catalog["readme_sha256"] != readme_evidence["sha256"]
        or after_catalog["sha256sums_sha256"] != sums_evidence["sha256"]
        or after_catalog["bridge_listed_exactly_once"] is not True
        or after_catalog["predecessor_current_source_bundle_listed_exactly_once"] is not True
        or after_catalog["final_bundle_listed_exactly_once"] is not True
        or after_catalog["corrupt_c7_listed"] is not False
        or after_catalog["historical_f115_preserved"] is not True
        or after_catalog["historical_f116_preserved"] is not True
        or after_catalog["all_prior_checksum_entries_preserved"] is not True
        or after_catalog["sole_current_source_bundle"] != final["path"]
        or readme_evidence["mode"] != "0644"
        or sums_evidence["mode"] != "0644"
    ):
        raise ValueError("F118 source-archive catalog declaration differs")
    try:
        sums_rows = [
            tuple(line.split("  ", 1))
            for line in sums_payload.decode("utf-8").splitlines()
        ]
    except UnicodeDecodeError as error:
        raise ValueError("F116 source-archive SHA256SUMS is not UTF-8") from error
    if any(
        len(row) != 2
        or SHA256_RE.fullmatch(row[0]) is None
        or not row[1]
        for row in sums_rows
    ):
        raise ValueError("F116 source-archive SHA256SUMS contains malformed rows")
    sums = {name: digest for digest, name in sums_rows}
    if len(sums) != len(sums_rows) or sums.get(Path(str(bridge["path"])).name) != bridge["sha256"]:
        raise ValueError("F116 source-archive SHA256SUMS bridge binding differs")
    if (
        sums.get(Path(str(final["path"])).name) != final["sha256"]
        or sums.get(Path(str(predecessor["path"])).name) != predecessor["sha256"]
        or sums.get(Path(R03_F115_SOURCE_BUNDLE).name) != R03_F115_SOURCE_BUNDLE_SHA256
        or F116_CORRUPT_C7_NAME in sums
    ):
        raise ValueError("F118 source-archive SHA256SUMS current/historical authority differs")

    verified_expected = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": final["sha256"],
        "final_head": final["head"],
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }
    reviewers = []
    candidate_paths = []
    reviewed_times = []
    for value, kind, decision, label in (
        (
            provenance_review,
            "provenance-security",
            "approved-for-publication",
            "F118 provenance/security review",
        ),
        (
            plasma_review,
            "plasma-scientific-continuation",
            "approved",
            "F118 plasma/scientific review",
        ),
    ):
        require_exact_keys(
            value,
            (
                "schema_version", "record_type", "checkpoint", "execution_epoch",
                "review_kind", "decision", "reviewed_candidate", "published_f118",
                "reviewer", "reviewed_utc", "findings", "limitations", "verified",
            ),
            label,
        )
        reviewer = value["reviewer"]
        if not isinstance(reviewer, dict):
            raise ValueError(f"{label} reviewer must be an object")
        require_exact_keys(reviewer, ("agent_id", "identity"), f"{label} reviewer")
        reviewer_id = require_nonempty_string(reviewer["agent_id"], f"{label} reviewer ID")
        require_nonempty_string(reviewer["identity"], f"{label} reviewer identity")
        candidate = value["reviewed_candidate"]
        if not isinstance(candidate, dict):
            raise ValueError(f"{label} reviewed candidate must be an object")
        require_exact_keys(candidate, ("path", "sha256"), f"{label} reviewed candidate")
        require_nonempty_string(candidate["path"], f"{label} reviewed-candidate path")
        reviewed = require_utc(value["reviewed_utc"], f"{label} timestamp")
        if (
            value["schema_version"] != 1
            or value["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or value["checkpoint"] != "F-118"
            or value["execution_epoch"] != EXECUTION_EPOCH
            or value["review_kind"] != kind
            or value["decision"] != decision
            or candidate["sha256"] != artifact_sha256
            or value["published_f118"] != {"path": str(artifact_path), "sha256": artifact_sha256}
            or value["verified"] != verified_expected
            or not isinstance(value["findings"], list)
            or not value["findings"]
            or any(not isinstance(item, str) or not item for item in value["findings"])
            or not isinstance(value["limitations"], list)
            or not value["limitations"]
            or any(not isinstance(item, str) or not item for item in value["limitations"])
            or reviewed < generated
            or reviewed > current_utc() + FUTURE_SKEW
        ):
            raise ValueError(f"{label} identity, independence, or verification differs")
        reviewers.append(reviewer_id)
        candidate_paths.append(candidate["path"])
        reviewed_times.append(reviewed)
    if len(set(reviewers)) != 2 or len(set(candidate_paths)) != 1:
        raise ValueError("F118 independent reviews must have distinct reviewers and one candidate")
    review_assurance = declared_process_independence_assurance(
        {
            "F118-provenance-security-reviewer": reviewers[0],
            "F118-plasma-scientific-reviewer": reviewers[1],
        },
        "F118 independent reviews",
    )

    require_exact_keys(
        audit,
        (
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f116_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        ),
        "F118 publication audit",
    )
    reviews = audit["independent_reviews"]
    audit_catalog = audit["source_archive_catalog"]
    if not isinstance(reviews, dict) or not isinstance(audit_catalog, dict):
        raise ValueError("F116 audit reviews/catalog must be objects")
    require_exact_keys(
        reviews,
        (
            "reviews_bind_exact_published_f118_sha256",
            "provenance_security", "plasma_scientific_continuation",
        ),
        "F118 audit reviews",
    )
    require_exact_keys(
        audit_catalog,
        (
            "readme", "sha256sums", "bridge_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger", "sole_current_source_bundle",
        ),
        "F118 audit source-archive catalog",
    )
    validate_declared_publication_binding(
        audit["artifact"],
        expected=artifact_path,
        digest=artifact_sha256,
        mode="0444",
        label="F118 audit artifact",
    )
    validate_declared_publication_binding(
        reviews["provenance_security"],
        expected=provenance_review_path,
        digest=provenance_review_sha256,
        mode="0444",
        label="F118 audit provenance review",
    )
    validate_declared_publication_binding(
        reviews["plasma_scientific_continuation"],
        expected=plasma_review_path,
        digest=plasma_review_sha256,
        mode="0444",
        label="F118 audit plasma review",
    )
    validate_declared_publication_binding(
        audit_catalog["readme"],
        expected=readme_path,
        digest=readme_evidence["sha256"],
        mode="0644",
        label="F118 audit source-archive README",
    )
    validate_declared_publication_binding(
        audit_catalog["sha256sums"],
        expected=sums_path,
        digest=sums_evidence["sha256"],
        mode="0644",
        label="F118 audit source-archive SHA256SUMS",
    )
    expected_bridge_audit = {
        "path": str(bridge["path_object"]),
        "sha256": bridge["sha256"],
        "mode": "0644",
        "links": 1,
        "head": bridge["head"],
        "role": "retained-non-current-bridge",
        "selected_as_current": False,
    }
    expected_final_audit = {
        "path": str(final["path_object"]),
        "sha256": final["sha256"],
        "mode": "0644",
        "links": 1,
        "head": final["head"],
        "selected_as_current": True,
    }
    expected_predecessor_audit = {
        "path": str(expected_path(str(predecessor["path"]))),
        "sha256": predecessor["sha256"],
        "mode": "0644",
        "links": 1,
        "head": predecessor["head"],
        "role": "retained-non-current-predecessor",
        "selected_as_current": False,
    }
    published = require_utc(audit["published_utc"], "F118 publication time")
    if (
        audit["schema_version"] != 1
        or audit["record_type"] != "stage-i-current-source-authority-supersession-publication-audit"
        or audit["checkpoint"] != "F-118"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["publication"] != F118_PUBLICATION_METHOD
        or audit["historical_f116_authority"] != historical_f116_digests
        or audit["authority_and_enforcement"] != F118_AUTHORIZATION
        or reviews["reviews_bind_exact_published_f118_sha256"] != artifact_sha256
        or audit_catalog["bridge_bundle"] != expected_bridge_audit
        or audit_catalog["predecessor_current_source_bundle"] != expected_predecessor_audit
        or audit_catalog["current_source_bundle"] != expected_final_audit
        or audit_catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
        or audit_catalog["sole_current_source_bundle"] != str(final["path_object"])
        or published < max(reviewed_times)
        or published > current_utc() + FUTURE_SKEW
    ):
        raise ValueError("F118 publication audit identity, chronology, or authority differs")

    return {
        "controller_helper": {**helper, "evidence": read_stable_regular_file(CONTROLLER_HELPER, "F118 promoted helper")[1]},
        "generator": generator,
        "matrix": matrix,
        "source_bundle": {
            "path": str(final["path_object"]),
            "sha256": final["sha256"],
            "verified_revisions": final["verified_revisions"],
            "evidence": bundle_evidence,
            "independent_validation": bundle_validation,
        },
        "authority_evidence": {
            "artifact": evidence_file,
            "provenance_review": provenance_review_file,
            "plasma_review": plasma_review_file,
            "publication_audit": audit_file,
        },
        "independent_review_assurance": review_assurance,
    }


def validate_f115_independent_review(
    value: dict[str, object],
    evidence: dict[str, object],
    *,
    path: Path,
    digest: str,
    kind: str,
    decision: str,
    label: str,
) -> str:
    """Require one immutable independent review of the exact published F-115 bytes."""

    require_immutable_publication_evidence(evidence, digest, label)
    reviewer = value.get("reviewer")
    published = value.get("published_f115")
    candidate = value.get("reviewed_candidate")
    if (
        value.get("schema_version") != 1
        or value.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-independent-review"
        or value.get("checkpoint") != "F-115"
        or value.get("execution_epoch") != EXECUTION_EPOCH
        or value.get("review_kind") != kind
        or value.get("decision") != decision
        or not isinstance(reviewer, dict)
        or not isinstance(published, dict)
        or not isinstance(candidate, dict)
        or published
        != {"path": str(r03_f115_path()), "sha256": R03_F115_SHA256}
        or candidate.get("sha256") != R03_F115_SHA256
    ):
        raise ValueError(f"{label} identity/exact-byte binding differs")
    if kind == "provenance-security":
        reviewer_authorized = (
            reviewer == F115_PROVENANCE_REVIEWER
            and value.get("scope") == list(F115_PROVENANCE_REVIEW_SCOPE)
        )
    else:
        verified = value.get("verified")
        reviewer_authorized = (
            reviewer == F115_PLASMA_REVIEWER
            and isinstance(verified, dict)
            and verified.get("authorization_limitations")
            == F115_PLASMA_AUTHORIZATION_LIMITATIONS
        )
    if not reviewer_authorized:
        raise ValueError(f"{label} reviewer lacks explicit F-115 authority")
    reviewed = require_utc(value.get("reviewed_utc"), f"{label} reviewed_utc")
    if reviewed > current_utc() + FUTURE_SKEW:
        raise ValueError(f"{label} review time is implausibly future")
    if normalized_path(Path(str(evidence["path"]))) != normalized_path(path):
        raise ValueError(f"{label} path differs")
    return str(reviewer["agent_id"])


def validate_f115_s02_manifest(
    info: dict[str, object],
    parent: dict[str, object],
    sole: dict[str, object],
) -> None:
    """Require any retained F-115 s02 intent/result to preserve the exact profile."""

    manifest = info["manifest"]
    allocation = manifest["allocation"]
    command = manifest["command"]
    if (
        info["segment"] != R03_F115_SEGMENT
        or info["index"] != 2
        or abs(info["start"] - Decimal(str(sole["restart_time"]))) > Decimal("0.000001")
        or info["target"] != Decimal("0.5")
        or info["state"] not in {"prepared", "submitted", "recorded"}
        or info["nodes"] != 1
        or allocation.get("ranks_per_node") != RANKS_PER_NODE
        or allocation.get("cpus_per_task") != CPUS_PER_TASK
        or allocation.get("requested_walltime") != sole["walltime"]
        or allocation.get("requested_seconds") != 7200
        or command.get("athena_walltime") != sole["athena_walltime"]
        or command.get("production_utility")
        != {
            "committed": True,
            "path": str(CONTROLLER_HELPER),
            "revision": R03_F115_CONTROLLER_REVISION,
            "sha256": R03_F115_CONTROLLER_SHA256,
        }
        or command.get("source_bundle")
        != {
            "path": str(f115_source_bundle_path()),
            "sha256": R03_F115_SOURCE_BUNDLE_SHA256,
            "verified_revisions": [SOURCE_REVISION, R03_F115_CONTROLLER_REVISION],
        }
    ):
        raise ValueError("retained R03 F-115 s02 helper/source/resource profile differs")
    validate_parent_link(info, parent)
    if info["state"] == "recorded":
        accounting = manifest.get("accounting")
        inspection = manifest.get("scientific_inspection")
        if (
            not isinstance(accounting, dict)
            or accounting.get("result") != "accepted"
            or not isinstance(inspection, dict)
            or decimal_value(inspection.get("final_time"), "R03 F-115 s02 final time")
            != Decimal("0.5")
        ):
            raise ValueError("recorded R03 F-115 s02 is not an exact accepted endpoint")


def validate_r03_authorization(
    lineage: list[dict[str, object]],
    active: dict[str, object] | None,
    case_infos: list[dict[str, object]],
) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    """Validate the permanent published F-115 anchor and any exact R03 s02 intent."""

    s00 = [
        info for info in case_infos
        if info["state"] == "recorded" and info["segment"] == "s00_rankio_t0_t0p5"
    ]
    if len(s00) != 1:
        raise ValueError("R03 F-115 requires the exact retained clean-partial s00 parent")
    parent = s00[0]
    parent_endpoint = decimal_value(
        parent["manifest"]["scientific_inspection"]["final_time"],
        "R03 F-115 parent endpoint",
    )
    path = r03_f115_path()
    authorization, evidence = read_json_file(path, "R03 sole-next authorization")
    require_immutable_publication_evidence(evidence, R03_F115_SHA256, "R03 F-115 authority")
    f115_bundle_validation = validate_f115_source_bundle_coverage()
    if (
        authorization.get("schema_version") != 1
        or authorization.get("execution_epoch") != EXECUTION_EPOCH
        or authorization.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-evidence"
        or authorization.get("checkpoint") != "F-115"
    ):
        raise ValueError("R03 sole-next authorization identity is invalid")
    authorization_block = authorization.get("authorization")
    if not isinstance(authorization_block, dict):
        raise ValueError("R03 F-115 authorization block is malformed")
    sole = authorization_block.get("sole_next_segment_profile")
    if not isinstance(sole, dict):
        raise ValueError("R03 sole-next authorization lacks an exact profile")
    if not lineage:
        raise ValueError("R03 sole-next authorization cannot authorize an absent lineage")
    restart = terminal_restart(parent)
    expected = {
        "athena_walltime": "01:50:00",
        "case_id": R03,
        "cpus_per_task": CPUS_PER_TASK,
        "executable": str(executable_path()),
        "executable_revision": SOURCE_REVISION,
        "executable_sha256": EXECUTABLE_SHA256,
        "matrix": str(frozen_matrix_path()),
        "matrix_sha256": MATRIX_SHA256,
        "nodes": 1,
        "override": "time/tlim=0.5",
        "parent_job_id": str(parent["manifest"]["job_id"]),
        "parent_result": parent["manifest"]["accounting"]["result"],
        "parent_segment": parent["segment"],
        "ranks_per_node": RANKS_PER_NODE,
        "restart_file": restart["restart_file"],
        "restart_time": float(parent_endpoint),
        "segment": R03_F115_SEGMENT,
        "source_bundle": str(f115_source_bundle_path()),
        "source_bundle_sha256": R03_F115_SOURCE_BUNDLE_SHA256,
        "source_dir": str(FROZEN_SOURCE),
        "time_tlim_target": 0.5,
        "walltime": "02:00:00",
    }
    if sole != expected:
        raise ValueError(
            "R03 F-115 sole-next profile differs from exact authenticated "
            "parent/resources/provenance"
        )
    index, start, target = parse_segment(str(sole["segment"]), "R03 sole-next segment")
    if index != 2 or abs(start - parent_endpoint) > Decimal("0.000001"):
        raise ValueError("R03 sole-next authorization does not continue its lineage")
    if decimal_value(sole["time_tlim_target"], "R03 authorized target") != target:
        raise ValueError("R03 sole-next target differs from its segment")
    validate_resource_values(
        R03, 1, sole["walltime"], sole["athena_walltime"],
        Decimal(walltime_seconds(sole["walltime"], "R03 walltime")) / Decimal(3600),
        "R03 sole-next profile",
    )
    cancelled = [
        info for info in case_infos
        if info["state"] == "cancelled" and info["segment"] == R03_F114_SEGMENT
    ]
    if len(cancelled) != 1 or int(cancelled[0]["index"]) != 1:
        raise ValueError("R03 F-115 requires the retained cancelled s01 operational identity")
    cancelled_info = cancelled[0]
    cancelled_manifest = cancelled_info["manifest"]
    if (
        abs(cancelled_info["start"] - parent_endpoint) > Decimal("0.000001")
        or cancelled_info["target"] != Decimal("0.5")
        or cancelled_info["nodes"] != 1
        or cancelled_manifest["allocation"]["requested_walltime"] != "02:00:00"
        or cancelled_manifest["command"]["athena_walltime"] != "01:50:00"
        or cancelled_manifest["command"].get("source_bundle")
        != {
            "path": str(expected_path(R03_F114_SOURCE_BUNDLE)),
            "sha256": R03_F114_SOURCE_BUNDLE_SHA256,
            "verified_revisions": [SOURCE_REVISION, R03_F114_CONTROLLER_REVISION],
        }
    ):
        raise ValueError("R03 F-115 retained cancelled s01 provenance/resources differ")
    validate_parent_link(cancelled_info, parent)
    cancellation = authorization.get("cancelled_submission")
    if not isinstance(cancellation, dict):
        raise ValueError("R03 F-115 cancellation evidence is missing")
    cancellation_evidence = cancellation.get("evidence")
    cancelled_reservation = cancelled_info.get("reservation")
    cancelled_job_id = str(cancelled_info["manifest"].get("job_id", ""))
    if (
        cancellation.get("state") != "CANCELLED"
        or cancellation.get("exit_code") != "0:0"
        or cancellation.get("allocated_nodes") != 0
        or cancellation.get("elapsed_seconds") != 0
        or cancellation.get("reusable") is not False
        or JOB_RE.fullmatch(cancelled_job_id) is None
        or str(cancellation.get("job_id")) != cancelled_job_id
        or not isinstance(cancelled_reservation, dict)
        or str(cancelled_reservation.get("job_id")) != cancelled_job_id
        or not isinstance(cancellation_evidence, dict)
    ):
        raise ValueError("R03 F-115 cancellation/no-reuse boundary differs")
    validate_declared_publication_binding(
        cancellation_evidence.get("live_cancelled_manifest"),
        expected=Path(str(cancelled_info["evidence"]["path"])),
        digest=str(cancelled_info["evidence"]["sha256"]),
        mode=str(cancelled_info["evidence"]["mode"]),
        label="R03 F-115 retained cancelled s01 manifest",
    )
    expected_supersession = {
        "segment": {
            "from": R03_F114_SEGMENT,
            "reason": "s01 is an immutable cancelled no-start manifest and directory",
            "to": R03_F115_SEGMENT,
        },
        "source_bundle": {
            "from": str(expected_path(R03_F114_SOURCE_BUNDLE)),
            "to": str(f115_source_bundle_path()),
        },
        "source_bundle_sha256": {
            "from": R03_F114_SOURCE_BUNDLE_SHA256,
            "to": R03_F115_SOURCE_BUNDLE_SHA256,
        },
    }
    if (
        authorization_block.get("supersedes_f114_profile_only_where_explicitly_listed")
        != expected_supersession
    ):
        raise ValueError("R03 F-115 supersession scope differs")
    scope = authorization.get("scope")
    requirements = authorization.get("publication_requirements")
    if (
        not isinstance(scope, dict)
        or scope.get("relationship")
        != "source-bundle-binding-and-cancelled-segment-identity-supersession"
        or not isinstance(requirements, list)
        or not any("Do not reuse cancelled job" in str(item) for item in requirements)
        or not any("do not invoke direct sbatch" in str(item) for item in requirements)
    ):
        raise ValueError("R03 F-115 scope/publication requirements differ")

    implementation = authorization.get("implementation")
    implementation_bundle = (
        implementation.get("source_bundle") if isinstance(implementation, dict) else None
    )
    implementation_helper = (
        implementation.get("stage_i_helper") if isinstance(implementation, dict) else None
    )
    if implementation_bundle != {
        "complete_history": True,
        "head": R03_F115_CONTROLLER_REVISION,
        "links": 1,
        "mode": "0644",
        "path": f115_source_bundle_path().relative_to(CANONICAL_ROOT).as_posix(),
        "sha256": R03_F115_SOURCE_BUNDLE_SHA256,
        "verified_revisions": list(F115_SOURCE_BUNDLE_REQUIRED_REVISIONS),
    }:
        raise ValueError("R03 F-115 source-bundle implementation declaration differs")
    if not isinstance(implementation_helper, dict) or (
        implementation_helper.get("links") != 1
        or implementation_helper.get("mode") != "0644"
        or implementation_helper.get("sha256") != R03_F115_CONTROLLER_SHA256
    ):
        raise ValueError("R03 F-115 helper implementation declaration differs")

    s02_infos = [info for info in case_infos if info["segment"] == R03_F115_SEGMENT]
    if len(s02_infos) > 1:
        raise ValueError("R03 F-115 s02 operational identity is duplicated")
    if s02_infos:
        validate_f115_s02_manifest(s02_infos[0], parent, sole)
    if active is not None and active["segment"] != sole["segment"]:
        raise ValueError("active R03 lane differs from the sole reviewed authorization")
    if active is not None:
        active_manifest = active["manifest"]
        if (
            active["nodes"] != 1
            or active_manifest["allocation"]["requested_walltime"] != sole["walltime"]
            or active_manifest["command"]["athena_walltime"] != sole["athena_walltime"]
            or active_manifest["command"]["source_restart_file"] != sole["restart_file"]
        ):
            raise ValueError("active R03 resources/restart differ from the sole authorization")
    audit_path = r03_f115_publication_audit_path()
    audit, audit_evidence = read_json_file(audit_path, "R03 F-115 publication audit")
    require_immutable_publication_evidence(
        audit_evidence, R03_F115_PUBLICATION_AUDIT_SHA256, "R03 F-115 publication audit"
    )
    if (
        audit.get("schema_version") != 1
        or audit.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-publication-audit"
        or audit.get("checkpoint") != "F-115"
        or audit.get("execution_epoch") != EXECUTION_EPOCH
        or audit.get("publication")
        != "atomic-write-fsync-rename-fsync-under-canonical-stage-i-lock"
    ):
        raise ValueError("R03 F-115 publication audit identity differs")
    validate_declared_publication_binding(
        audit.get("artifact"),
        expected=path,
        digest=R03_F115_SHA256,
        mode="0444",
        label="R03 F-115 publication audit artifact",
    )
    authority = audit.get("authority_and_enforcement")
    if not isinstance(authority, dict):
        raise ValueError("R03 F-115 audit authority/enforcement block is missing")
    require_exact_keys(
        authority,
        (
            "authorization_kind", "direct_sbatch_authorized", "enforcement_chain",
            "f115_authority", "reuse_cancelled_job_or_s01_authorized",
            "shared_root_acknowledgement_authorized_by_f115",
            "shared_root_acknowledgement_requires_separate_exact_isolation_review_after_prepare",
            "sole_next_segment_profile",
        ),
        "R03 F-115 audit authority/enforcement",
    )
    if (
        authority["authorization_kind"] != "procedural-pre-prepare-sole-profile-authority"
        or authority["direct_sbatch_authorized"] is not False
        or authority["reuse_cancelled_job_or_s01_authorized"] is not False
        or authority["shared_root_acknowledgement_authorized_by_f115"] is not False
        or authority[
            "shared_root_acknowledgement_requires_separate_exact_isolation_review_after_prepare"
        ]
        is not True
        or authority["sole_next_segment_profile"] != sole
        or not isinstance(authority["enforcement_chain"], list)
    ):
        raise ValueError("R03 F-115 audit over-authorizes or changes the sole profile")
    validate_declared_publication_binding(
        authority["f115_authority"],
        expected=path,
        digest=R03_F115_SHA256,
        mode="0444",
        label="R03 F-115 audit authority",
    )

    review_specs = {
        "provenance_security": (
            r03_f115_provenance_security_review_path(),
            R03_F115_PROVENANCE_SECURITY_REVIEW_SHA256,
            "provenance-security",
            "approved-for-publication",
        ),
        "plasma_scientific_continuation": (
            r03_f115_plasma_scientific_review_path(),
            R03_F115_PLASMA_SCIENTIFIC_REVIEW_SHA256,
            "plasma-scientific-continuation",
            "approved",
        ),
    }
    review_bindings = audit.get("independent_reviews")
    if not isinstance(review_bindings, dict):
        raise ValueError("R03 F-115 audit independent-review bindings are missing")
    require_exact_keys(
        review_bindings,
        (
            "plasma_scientific_continuation", "provenance_security",
            "reviews_bind_exact_published_f115_sha256",
        ),
        "R03 F-115 audit independent reviews",
    )
    if review_bindings["reviews_bind_exact_published_f115_sha256"] != R03_F115_SHA256:
        raise ValueError("R03 F-115 independent reviews do not bind the published bytes")
    retained_reviews: dict[str, dict[str, object]] = {}
    reviewers = []
    for key, (review_path, digest, kind, decision) in review_specs.items():
        review, review_evidence = read_json_file(review_path, f"R03 F-115 {key} review")
        validate_declared_publication_binding(
            review_bindings[key],
            expected=review_path,
            digest=digest,
            mode="0444",
            label=f"R03 F-115 {key} review",
        )
        reviewers.append(
            validate_f115_independent_review(
                review,
                review_evidence,
                path=review_path,
                digest=digest,
                kind=kind,
                decision=decision,
                label=f"R03 F-115 {key} review",
            )
        )
        retained_reviews[key] = review_evidence
    if len(set(reviewers)) != len(reviewers):
        raise ValueError("R03 F-115 independent reviews do not have distinct reviewers")
    review_assurance = declared_process_independence_assurance(
        {
            "F115-provenance-security-reviewer": reviewers[0],
            "F115-plasma-scientific-reviewer": reviewers[1],
        },
        "R03 F-115 independent reviews",
    )

    reproducible = audit.get("reproducible_implementation_authority")
    if not isinstance(reproducible, dict):
        raise ValueError("R03 F-115 audit reproducible implementation authority is missing")
    bundle = reproducible.get("authoritative_source_bundle")
    helper = reproducible.get("committed_stage_i_helper")
    if not isinstance(bundle, dict) or not isinstance(helper, dict):
        raise ValueError("R03 F-115 audit reproducible implementation bindings are malformed")
    validate_declared_publication_binding(
        {
            key: bundle.get(key) for key in ("path", "sha256", "mode", "links")
        },
        expected=f115_source_bundle_path(),
        digest=R03_F115_SOURCE_BUNDLE_SHA256,
        mode="0644",
        label="R03 F-115 authoritative source bundle",
    )
    if (
        bundle.get("complete_history") is not True
        or bundle.get("head") != R03_F115_CONTROLLER_REVISION
        or helper
        != {
            "path": "scripts/frontier/cgl_lf_stage_i.py",
            "sha256": R03_F115_CONTROLLER_SHA256,
        }
    ):
        raise ValueError("R03 F-115 reproducible controller provenance differs")
    catalog = audit.get("source_archive_catalog")
    if not isinstance(catalog, dict) or (
        catalog.get("corrupt_c7_absent_from_active_checksum_ledger") is not True
        or catalog.get("corrupt_c7_retained_as_incident_evidence") is not True
        or catalog.get("full_active_checksum_ledger") != "passed"
        or catalog.get("new_bundle_present_exactly_once") is not True
    ):
        raise ValueError("R03 F-115 source-archive catalog claims differ")
    if (
        require_utc(audit.get("published_utc"), "R03 F-115 publication time")
        > current_utc() + FUTURE_SKEW
        or require_utc(audit.get("audit_generated_utc"), "R03 F-115 audit time")
        > current_utc() + FUTURE_SKEW
    ):
        raise ValueError("R03 F-115 publication audit is implausibly future")
    current_sole = None if s02_infos else sole
    return current_sole, {
        "artifact": evidence,
        "publication_audit": audit_evidence,
        "independent_reviews": retained_reviews,
        "independent_review_assurance": review_assurance,
        "source_bundle_validation": f115_bundle_validation,
        "planner_scope": {
            "planning_only": True,
            "direct_sbatch_authorized": False,
            "shared_root_acknowledgement_authorized": False,
            "requires_separate_exact_isolation_review_after_prepare": True,
        },
        "continuation": restart,
        "s02_manifest": s02_infos[0]["evidence"] if s02_infos else None,
    }


def available_storage_bytes(path: Path) -> int:
    """Return currently available bytes at one filesystem path."""

    profile = os.statvfs(path)
    return int(profile.f_bavail) * int(profile.f_frsize)


def parse_scheduler_evidence(payload: bytes, job_id: str) -> dict[str, str]:
    """Parse one exact retained eight-column top-level Slurm allocation row."""

    try:
        rows = [
            row for row in csv.reader(payload.decode("utf-8").splitlines(), delimiter="|")
            if row and row[0] == job_id
        ]
    except UnicodeDecodeError as error:
        raise ValueError("retained scheduler evidence is not UTF-8") from error
    for row in rows:
        while row and not row[-1]:
            row.pop()
    if len(rows) != 1 or len(rows[0]) != 8:
        raise ValueError(f"expected one eight-column scheduler row for {job_id}")
    keys = (
        "job_id", "job_name", "state", "exit_code", "nodes", "elapsed_seconds",
        "submitted_utc", "completed_utc",
    )
    result = dict(zip(keys, rows[0]))
    result["state"] = result["state"].split()[0].split("+")[0]
    return result


def parse_account_scheduler_timestamp(value: object, label: str) -> datetime:
    """Require one offset-qualified Slurm accounting timestamp."""

    text = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(text)
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include an explicit UTC offset")
    return parsed


def parse_optional_account_scheduler_timestamp(
    value: object, label: str
) -> datetime | None:
    """Parse one optional account-wide Slurm timestamp as UTC."""

    if value in ACCOUNT_MISSING_TIMESTAMPS:
        return None
    return parse_account_scheduler_timestamp(value, label).astimezone(timezone.utc)


def parse_account_scheduler_evidence(payload: bytes) -> list[dict[str, object]]:
    """Parse one complete all-users top-level Slurm account query."""

    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("R17 account scheduler evidence is not UTF-8") from error
    if not lines or lines[0] != ACCOUNT_SCHEDULER_HEADER:
        raise ValueError("R17 account scheduler evidence has an invalid header")
    records = []
    seen = set()
    for line in lines[1:]:
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != len(ACCOUNT_SACCT_FIELDS):
            raise ValueError("R17 account scheduler evidence has invalid columns")
        (
            job_id, job_name, state_value, exit_code, nodes_text, elapsed_text,
            submit, start, end, partition, account, owner,
        ) = fields
        state_parts = state_value.split()
        state = state_parts[0].split("+")[0] if state_parts else ""
        if (
            ACCOUNT_JOB_RE.fullmatch(job_id) is None
            or job_id in seen
            or not job_name
            or not state
            or not exit_code
            or not partition
            or account.casefold() != ACCOUNT.casefold()
            or not owner
        ):
            raise ValueError("R17 account scheduler evidence has invalid job identity")
        try:
            nodes = int(nodes_text)
            elapsed = int(elapsed_text)
        except ValueError as error:
            raise ValueError("R17 account scheduler evidence has invalid numeric values") from error
        submit_time = parse_optional_account_scheduler_timestamp(
            submit, "R17 account scheduler submit time"
        )
        start_time = parse_optional_account_scheduler_timestamp(
            start, "R17 account scheduler start time"
        )
        end_time = parse_optional_account_scheduler_timestamp(
            end, "R17 account scheduler end time"
        )
        if (
            nodes < 0
            or elapsed < 0
            or submit_time is None
            or (
                start_time is None
                and (
                    elapsed != 0
                    or state in {"RUNNING", "COMPLETING", "SUSPENDED", "COMPLETED"}
                    or (end_time is not None and submit_time > end_time)
                )
            )
            or (start_time is not None and nodes <= 0)
            or (start_time is not None and submit_time > start_time)
            or (
                start_time is not None
                and end_time is not None
                and (
                    end_time < start_time
                    or abs((end_time - start_time).total_seconds() - elapsed) > 1.0
                )
            )
        ):
            raise ValueError("R17 account scheduler evidence has invalid time ordering")
        seen.add(job_id)
        records.append(
            {
                "job_id": job_id,
                "job_name": job_name,
                "state": state,
                "exit_code": exit_code,
                "nodes": nodes,
                "elapsed_seconds": elapsed,
                "submit_utc": submit,
                "start_utc": start,
                "end_utc": end,
                "partition": partition,
                "account": account,
                "owner": owner,
            }
        )
    if not records:
        raise ValueError("R17 account scheduler evidence is empty")
    return sorted(records, key=lambda record: str(record["job_id"]))


def r17_account_scheduler_query_contract(start: datetime, end: datetime) -> dict[str, object]:
    """Reproduce the producer's exact all-users execution-interval query."""

    if start.tzinfo is None or end.tzinfo is None or not start < end:
        raise ValueError("R17 account scheduler query interval is invalid")
    query_start = start - timedelta(seconds=1)
    query_end = end + timedelta(seconds=1)
    return {
        "account": ACCOUNT,
        "all_users": True,
        "allocations_only": True,
        "expanded_arrays": True,
        "start_utc": query_start.astimezone(timezone.utc).isoformat(timespec="seconds"),
        "end_utc": query_end.astimezone(timezone.utc).isoformat(timespec="seconds"),
        "start_argument": query_start.strftime("%Y-%m-%dT%H:%M:%S"),
        "end_argument": query_end.strftime("%Y-%m-%dT%H:%M:%S"),
        "start_scheduler_offset": query_start.strftime("%z"),
        "end_scheduler_offset": query_end.strftime("%z"),
        "fields": list(ACCOUNT_SACCT_FIELDS),
    }


def read_r17_relative_artifact(
    binding: object, label: str, mode: str, expected_relative: Path
) -> tuple[bytes, dict[str, object]]:
    """Read one exact root-relative R17 artifact through its retained binding."""

    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding must be an object")
    require_exact_keys(binding, ("path", "sha256"), f"{label} binding")
    relative = Path(require_nonempty_string(binding["path"], f"{label} path"))
    if relative.is_absolute() or ".." in relative.parts or relative != expected_relative:
        raise ValueError(f"{label} path differs")
    digest = require_sha256(binding["sha256"], f"{label} SHA-256")
    payload, evidence = read_stable_regular_file(expected_path(relative.as_posix()), label)
    if evidence["sha256"] != digest or evidence["mode"] != mode:
        raise ValueError(f"{label} digest or mode differs")
    return payload, evidence


def validate_r17_account_exclusivity_artifacts(
    qualification: dict[str, object], qualification_reviewed: datetime
) -> tuple[dict[str, object], dict[str, dict[str, object]]]:
    """Reproduce the retained all-users R17 execution-exclusivity proof."""

    job_id = require_nonempty_string(qualification.get("job_id"), "R17 qualification job ID")
    if JOB_RE.fullmatch(job_id) is None:
        raise ValueError("R17 qualification job ID is invalid")
    scheduler_payload, scheduler_profile = read_r17_relative_artifact(
        qualification.get("scheduler_evidence"),
        "R17 qualification scheduler evidence",
        "0444",
        Path(f"accounting/{job_id}.r17_qualification.sacct.txt"),
    )
    scheduler = parse_scheduler_evidence(scheduler_payload, job_id)
    account_payload, account_profile = read_r17_relative_artifact(
        qualification.get("account_scheduler_evidence"),
        "R17 account scheduler evidence",
        "0444",
        Path(f"accounting/{job_id}.r17_qualification.account.sacct.txt"),
    )
    exclusivity_payload, exclusivity_profile = read_r17_relative_artifact(
        qualification.get("account_exclusivity_evidence"),
        "R17 account exclusivity evidence",
        "0444",
        Path(
            "accounting/"
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
        ),
    )
    try:
        exclusivity = json.loads(
            exclusivity_payload.decode("utf-8"),
            object_pairs_hook=duplicate_rejecting_object,
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError("R17 account exclusivity evidence is not unambiguous UTF-8 JSON") from error
    if not isinstance(exclusivity, dict):
        raise ValueError("R17 account exclusivity evidence must be an object")
    require_exact_keys(
        exclusivity,
        (
            "schema_version",
            "record_type",
            "execution_epoch",
            "measured_utc",
            "query_contract",
            "visibility_contract",
            "raw_account_scheduler_sha256",
            "qualification_job",
            "qualification_job_sha256",
            "account_jobs",
            "account_jobs_sha256",
            "overlapping_job_ids",
            "exclusive_entire_execution_interval",
        ),
        "R17 account exclusivity evidence",
    )
    records = parse_account_scheduler_evidence(account_payload)
    target_jobs = [record for record in records if record["job_id"] == job_id]
    if len(target_jobs) != 1:
        raise ValueError("R17 account scheduler evidence lacks one exact qualification job")
    target = target_jobs[0]
    start = parse_account_scheduler_timestamp(target["start_utc"], "R17 account target start")
    end = parse_account_scheduler_timestamp(target["end_utc"], "R17 account target end")
    measured = require_utc(exclusivity["measured_utc"], "R17 account exclusivity measurement")
    visibility = exclusivity["visibility_contract"]
    if not isinstance(visibility, dict):
        raise ValueError("R17 account exclusivity visibility contract must be an object")
    require_exact_keys(
        visibility,
        ("private_data", "all_users_job_visibility"),
        "R17 account exclusivity visibility contract",
    )
    private_data = require_nonempty_string(
        visibility["private_data"], "R17 account exclusivity PrivateData"
    )
    private_settings = {
        item.strip().casefold() for item in private_data.split(",") if item.strip()
    }
    target_scheduler = {
        key: target[key]
        for key in (
            "job_id",
            "job_name",
            "state",
            "exit_code",
            "nodes",
            "elapsed_seconds",
            "submit_utc",
            "start_utc",
            "end_utc",
            "partition",
            "account",
        )
    }
    overlapping = []
    for record in records:
        record_start = parse_optional_account_scheduler_timestamp(
            record["start_utc"], "R17 account scheduler overlap start"
        )
        record_end = parse_optional_account_scheduler_timestamp(
            record["end_utc"], "R17 account scheduler overlap end"
        )
        if record_start is not None and record_start < end.astimezone(timezone.utc) and (
            record_end is None or record_end > start.astimezone(timezone.utc)
        ):
            overlapping.append(str(record["job_id"]))
    expected = {
        "schema_version": 1,
        "record_type": "stage-i-r17-account-exclusivity-evidence",
        "execution_epoch": EXECUTION_EPOCH,
        "measured_utc": exclusivity["measured_utc"],
        "query_contract": r17_account_scheduler_query_contract(start, end),
        "visibility_contract": {
            "private_data": private_data,
            "all_users_job_visibility": True,
        },
        "raw_account_scheduler_sha256": sha256_bytes(account_payload),
        "qualification_job": target_scheduler,
        "qualification_job_sha256": compact_json_sha256(target_scheduler),
        "account_jobs": records,
        "account_jobs_sha256": compact_json_sha256(records),
        "overlapping_job_ids": [job_id],
        "exclusive_entire_execution_interval": True,
    }
    if (
        target["job_name"] != scheduler["job_name"]
        or target["state"] != scheduler["state"]
        or target["exit_code"] != scheduler["exit_code"]
        or target["nodes"] != int(scheduler["nodes"])
        or target["elapsed_seconds"] != int(scheduler["elapsed_seconds"])
        or target["submit_utc"] != scheduler["submitted_utc"]
        or target["end_utc"] != scheduler["completed_utc"]
        or target["end_utc"] != qualification.get("completed_utc")
        or target["state"] != qualification.get("state")
        or target["exit_code"] != qualification.get("exit_code")
        or target["nodes"] != qualification.get("nodes")
        or target["partition"] != PARTITION
        or visibility["all_users_job_visibility"] is not True
        or not private_settings
        or any(item in {"all", "jobs"} for item in private_settings)
        or measured < end.astimezone(timezone.utc)
        or measured > qualification_reviewed
        or measured > current_utc() + FUTURE_SKEW
        or overlapping != [job_id]
        or exclusivity != expected
    ):
        raise ValueError("R17 account-wide execution exclusivity evidence differs")
    return exclusivity, {
        "scheduler_evidence": scheduler_profile,
        "account_scheduler_evidence": account_profile,
        "account_exclusivity_evidence": exclusivity_profile,
    }


def validate_r17_rank_inventory(value: object, label: str) -> str:
    """Authenticate exactly one retained 0644 file for each of the 64 R17 ranks."""

    if not isinstance(value, list) or len(value) != 64:
        raise ValueError(f"{label} must contain exactly one file for each of 64 ranks")
    ranks = set()
    paths = set()
    for index, binding in enumerate(value):
        if not isinstance(binding, dict):
            raise ValueError(f"{label} binding {index} must be an object")
        require_exact_keys(binding, ("path", "sha256"), f"{label} binding {index}")
        relative = Path(require_nonempty_string(binding["path"], f"{label} binding {index} path"))
        rank_parts = [
            part for part in relative.parts if re.fullmatch(r"rank_[0-9]{8}", part)
        ]
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or len(rank_parts) != 1
            or relative.as_posix() in paths
        ):
            raise ValueError(f"{label} binding {index} is not one unique rank-local file")
        require_sha256(binding["sha256"], f"{label} binding {index} SHA-256")
        payload, evidence = read_stable_regular_file(
            expected_path(relative.as_posix()), f"{label} file {index}"
        )
        if (
            not payload
            or evidence["sha256"] != binding["sha256"]
            or evidence["mode"] != "0644"
        ):
            raise ValueError(f"{label} binding {index} bytes or mode differ")
        paths.add(relative.as_posix())
        ranks.add(rank_parts[0])
    if ranks != {f"rank_{rank:08d}" for rank in range(64)}:
        raise ValueError(f"{label} must contain exactly one file for each of 64 ranks")
    if value != sorted(value, key=lambda item: str(item["path"])):
        raise ValueError(f"{label} is not in deterministic rank-local path order")
    return recost_json_sha256(value)


def validate_r17_meshblock_decomposition_proof(
    value: object, terminal_output_inventory_sha256: str
) -> None:
    """Require a complete balanced 64-rank proof of the R17 logical mesh."""

    label = "R17 meshblock decomposition proof"
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    require_exact_keys(
        value,
        (
            "schema_version",
            "record_type",
            "resolution",
            "mesh_shape",
            "meshblock_shape",
            "logical_meshblock_grid",
            "logical_meshblocks",
            "ranks",
            "meshblocks_per_rank",
            "complete_block_rank_inventory",
            "complete_block_rank_inventory_sha256",
            "terminal_rank_local_output_inventory_sha256",
            "checks",
        ),
        label,
    )
    rank_records = value["complete_block_rank_inventory"]
    if (
        value["schema_version"] != 1
        or value["record_type"] != "stage-i-r17-decomposition-evidence"
        or value["resolution"] != "384x384x768"
        or value["mesh_shape"] != [384, 384, 768]
        or value["meshblock_shape"] != [32, 32, 64]
        or value["logical_meshblock_grid"] != [12, 12, 12]
        or value["logical_meshblocks"] != 1728
        or value["ranks"] != 64
        or value["meshblocks_per_rank"] != 27
        or not isinstance(rank_records, list)
        or len(rank_records) != 64
    ):
        raise ValueError(f"{label} dimensions or counts differ")

    observed_locations = []
    for expected_rank, record in enumerate(rank_records):
        if not isinstance(record, dict):
            raise ValueError(f"{label} rank record {expected_rank} must be an object")
        require_exact_keys(
            record,
            ("rank", "rank_name", "logical_meshblocks"),
            f"{label} rank record {expected_rank}",
        )
        locations = record["logical_meshblocks"]
        if (
            record["rank"] != expected_rank
            or record["rank_name"] != f"rank_{expected_rank:08d}"
            or not isinstance(locations, list)
            or len(locations) != 27
            or locations != sorted(locations)
        ):
            raise ValueError(f"{label} does not retain exactly 27 meshblocks per rank")
        rank_locations = set()
        for location in locations:
            if (
                not isinstance(location, list)
                or len(location) != 4
                or any(type(coordinate) is not int for coordinate in location)
                or location[3] != 0
            ):
                raise ValueError(f"{label} contains an invalid logical location")
            logical = tuple(location)
            if logical in rank_locations:
                raise ValueError(f"{label} duplicates a logical location within one rank")
            rank_locations.add(logical)
            observed_locations.append(logical)

    expected_locations = {
        (lx1, lx2, lx3, 0)
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    }
    if (
        value["complete_block_rank_inventory_sha256"] != compact_json_sha256(rank_records)
        or value["terminal_rank_local_output_inventory_sha256"]
        != terminal_output_inventory_sha256
        or value["checks"]
        != {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        }
        or len(observed_locations) != 1728
        or len(set(observed_locations)) != 1728
        or set(observed_locations) != expected_locations
    ):
        raise ValueError(f"{label} does not authenticate the exact 1728 logical meshblocks")


def r17_parameter_contract() -> dict[str, str]:
    """Return the producer's exact frozen R17 qualification parameter contract."""

    return {
        "job/basename": "qualification_R17_n08",
        "mesh/ix1_bc": "periodic",
        "mesh/ix2_bc": "periodic",
        "mesh/ix3_bc": "periodic",
        "mesh/nghost": "2",
        "mesh/nx1": "384",
        "mesh/nx2": "384",
        "mesh/nx3": "768",
        "mesh/ox1_bc": "periodic",
        "mesh/ox2_bc": "periodic",
        "mesh/ox3_bc": "periodic",
        "mesh/x1max": "1.0",
        "mesh/x1min": "0.0",
        "mesh/x2max": "1.0",
        "mesh/x2min": "0.0",
        "mesh/x3max": "2.0",
        "mesh/x3min": "0.0",
        "meshblock/nx1": "32",
        "meshblock/nx2": "32",
        "meshblock/nx3": "64",
        "mhd/backup_limiters": "false",
        "mhd/cgl_firehose_threshold": "parallel",
        "mhd/cgl_heat_flux": "landau_fluid",
        "mhd/cgl_heat_flux_integrator": "sts",
        "mhd/cgl_lf_record_pressure_work": "true",
        "mhd/cgl_lf_strict_admissibility": "true",
        "mhd/eos": "cgl",
        "mhd/firehose_limiter": "true",
        "mhd/lf_coefficient_mode": "local",
        "mhd/lf_k_parallel": "6.283185307179586",
        "mhd/limiter_hardwall": "true",
        "mhd/limiter_nu_coll": "1.0e10",
        "mhd/mirror_limiter": "true",
        "mhd/passive": "false",
        "mhd/reconstruct": "plm",
        "mhd/rsolver": "hlle",
        "output1/dt": "0.02",
        "output1/file_type": "hst",
        "output2/dt": "0.25",
        "output2/file_type": "bin",
        "output2/single_file_per_rank": "true",
        "output2/variable": "mhd_w_bcc",
        "output3/dt": "1.0",
        "output3/file_type": "rst",
        "output3/single_file_per_rank": "true",
        "problem/beta0": "10.0",
        "problem/paper_mode": "turbulence",
        "problem/passive_delta": "false",
        "problem/pgen_name": "cgl_lf_paper",
        "problem/user_hist": "true",
        "time/evolution": "dynamic",
        "time/integrator": "rk2",
        "time/sts_integrator": "rkl2",
        "time/tlim": "0.25",
        "turb_driving/dedt": "0.32",
        "turb_driving/driving_type": "1",
        "turb_driving/projection_policy": "mks24_alfvenic_perpendicular",
        "turb_driving/record_injected_work": "true",
        "turb_driving/tcorr": "2.0",
    }


def read_r17_bound_artifact(
    binding: object, label: str, *, mode: str | None = None
) -> tuple[bytes, dict[str, object], Path]:
    """Authenticate one canonical-root-relative producer evidence binding."""

    if not isinstance(binding, dict):
        raise ValueError(f"{label} binding must be an object")
    require_exact_keys(binding, ("path", "sha256"), f"{label} binding")
    relative_text = require_nonempty_string(binding["path"], f"{label} path")
    relative = Path(relative_text)
    if (
        relative.is_absolute()
        or ".." in relative.parts
        or relative.as_posix() != relative_text
    ):
        raise ValueError(f"{label} path must be canonical-root relative")
    payload, evidence = read_stable_regular_file(expected_path(relative_text), label)
    if evidence["sha256"] != binding["sha256"] or (
        mode is not None and evidence["mode"] != mode
    ):
        raise ValueError(f"{label} bytes or mode differ")
    return payload, evidence, relative


def parse_bound_json(payload: bytes, label: str) -> dict[str, object]:
    """Parse one unambiguous bound JSON object."""

    try:
        value = json.loads(
            payload.decode("utf-8"), object_pairs_hook=duplicate_rejecting_object
        )
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise ValueError(f"{label} is not unambiguous UTF-8 JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must contain an object")
    return value


def validate_r17_physics_contract(scientific: dict[str, object]) -> dict[str, object]:
    """Independently enforce the qualification producer's R17 physics gates."""

    measurements = scientific.get("physics_measurements")
    checks = scientific.get("checks")
    if not isinstance(measurements, dict) or not isinstance(checks, dict):
        raise ValueError("R17 scientific evidence lacks complete physics measurements/checks")
    require_exact_keys(
        measurements,
        (
            "finite_rank_outputs", "mass_relative_drift_max",
            "mhd_user_mass_mismatch_max", "lf_bad_counts_total",
            "normalized_ct_divb_max", "normalized_ct_divb_threshold",
            "normalized_ct_divb_below_threshold",
        ),
        "R17 physics measurements",
    )
    require_exact_keys(checks, R17_SCIENTIFIC_CHECKS, "R17 scientific checks")

    def producer_float_text(key: str) -> float:
        value = measurements[key]
        if not isinstance(value, str):
            raise ValueError(f"R17 physics measurement {key} is not producer text")
        try:
            parsed = float(value)
        except ValueError as error:
            raise ValueError(f"R17 physics measurement {key} is invalid") from error
        if not math.isfinite(parsed) or format(parsed, ".17g") != value:
            raise ValueError(f"R17 physics measurement {key} is not canonical finite text")
        return parsed

    mass_drift = producer_float_text("mass_relative_drift_max")
    mass_mismatch = producer_float_text("mhd_user_mass_mismatch_max")
    normalized_divb = producer_float_text("normalized_ct_divb_max")
    if (
        measurements["finite_rank_outputs"] != 64
        or isinstance(measurements["lf_bad_counts_total"], bool)
        or measurements["lf_bad_counts_total"] != 0
        or measurements["normalized_ct_divb_threshold"]
        != R17_MAX_NORMALIZED_CT_DIVB_TEXT
        or measurements["normalized_ct_divb_below_threshold"] is not True
        or mass_drift < 0.0
        or mass_drift > CONTINUATION_MASS_TOLERANCE
        or mass_mismatch < 0.0
        or mass_mismatch > CONTINUATION_MASS_TOLERANCE
        or normalized_divb < 0.0
        or normalized_divb >= float(R17_MAX_NORMALIZED_CT_DIVB_TEXT)
        or checks != R17_SCIENTIFIC_CHECKS
    ):
        raise ValueError("R17 physics measurements or scientific checks fail policy")
    return measurements


def validate_r17_frozen_science_build_contract(
    value: object, profile: dict[str, object]
) -> dict[str, object]:
    """Authenticate the producer's exact frozen R17 science/build contract."""

    if not isinstance(value, dict):
        raise ValueError("R17 frozen science/build contract must be an object")
    require_exact_keys(
        value,
        (
            "case_id", "case_name", "profile_class", "resolution", "mesh_shape",
            "meshblock_shape", "target_time", "scientific_policy", "run_basename",
            "source_revision", "source_bundle_sha256", "matrix_sha256", "input_sha256",
            "provenance_sha256", "execution_intent_sha256",
            "execution_contract_sha256", "parameter_contract",
            "parameter_contract_sha256", "executable_revision", "executable_sha256",
            "build_manifest_inventory_sha256",
        ),
        "R17 frozen science/build contract",
    )
    parameter_contract = r17_parameter_contract()
    exact = {
        "case_id": R17,
        "case_name": EXPECTED_CASES[R17][0],
        "profile_class": "scale_separation_384x384x768",
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "target_time": 0.25,
        "scientific_policy": "active_hardwall",
        "run_basename": "qualification_R17_n08",
        "source_revision": SOURCE_REVISION,
        "matrix_sha256": MATRIX_SHA256,
        "input_sha256": EXPECTED_CASES[R17][3],
        "parameter_contract": parameter_contract,
        "parameter_contract_sha256": compact_json_sha256(parameter_contract),
        "executable_revision": SOURCE_REVISION,
        "executable_sha256": profile["executable_sha256"],
        "build_manifest_inventory_sha256": profile["build_manifest_sha256"],
    }
    if any(value.get(key) != expected for key, expected in exact.items()):
        raise ValueError("R17 frozen science/build contract differs")
    for key in (
        "source_bundle_sha256", "provenance_sha256", "execution_intent_sha256",
        "execution_contract_sha256",
    ):
        require_sha256(value[key], f"R17 frozen contract {key}")
    return value


def validate_r17_operational_qualification_contract(
    qualification: dict[str, object],
    profile: dict[str, object],
    qualification_path: Path,
    qualification_profile: dict[str, object],
    qualification_review: dict[str, object],
    qualification_reviewed: datetime,
) -> tuple[dict[str, object], dict[str, dict[str, object]], set[str]]:
    """Authenticate the full exact schema-2 contract emitted by the R17 producer."""

    raw_keys = {
        "schema_version", "record_type", "execution_epoch", "completed_utc",
        "measured_utc", "measured_by", "job_id", "state", "exit_code", "nodes",
        "ranks", "prepared_wave", "qualification_evidence", "scientific_evidence",
        "scheduler_evidence", "account_scheduler_evidence",
        "account_exclusivity_evidence", "executable_sha256",
        "build_manifest_inventory", "build_manifest_inventory_sha256",
        "rank_local_outputs", "rank_local_output_inventory_sha256",
        "rank_local_restarts", "rank_local_restart_inventory_sha256",
        "decomposition_evidence", "restart_load_evidence",
        "physics_validation_evidence", "frozen_science_build_contract",
        "independent_review_contract", "authority",
    }
    require_exact_keys(qualification, raw_keys, "R17 operational qualification")
    expected_qualification = Path(
        "accounting/"
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_operational_qualification.json"
    )
    completed = parse_account_scheduler_timestamp(
        qualification["completed_utc"], "R17 qualification completion"
    ).astimezone(timezone.utc)
    measured = require_utc(
        qualification["measured_utc"], "R17 qualification measurement"
    )
    measured_by = require_nonempty_string(
        qualification["measured_by"], "R17 qualification measurement author"
    )
    if (
        qualification_path != expected_path(expected_qualification.as_posix())
        or qualification_profile["mode"] != "0444"
        or qualification["schema_version"] != 2
        or qualification["record_type"] != "stage-i-r17-operational-qualification"
        or qualification["execution_epoch"] != EXECUTION_EPOCH
        or qualification["state"] != "COMPLETED"
        or qualification["exit_code"] != "0:0"
        or qualification["nodes"] != 8
        or qualification["ranks"] != 64
        or qualification["executable_sha256"] != profile["executable_sha256"]
        or measured < completed
        or measured > qualification_reviewed
        or qualification["authority"]
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("R17 operational qualification identity or authority differs")
    frozen = validate_r17_frozen_science_build_contract(
        qualification["frozen_science_build_contract"], profile
    )

    expected_review_contract = {
        "required": True,
        "path": (
            expected_qualification.with_name(
                f"{expected_qualification.name}.independent_review.json"
            ).as_posix()
        ),
        "mode": "0444",
        "schema_version": 1,
        "record_type": "stage-i-r17-operational-qualification-independent-review",
        "execution_epoch": EXECUTION_EPOCH,
        "decision": "approved",
        "candidate_path": str(qualification_path),
        "candidate_sha256_required": True,
        "reviewer_must_differ_from": [measured_by],
        "reviewed_after_utc": qualification["measured_utc"],
    }
    if (
        qualification["independent_review_contract"] != expected_review_contract
        or qualification_reviewed < measured
        or qualification_review.get("reviewer") == measured_by
    ):
        raise ValueError("R17 independent-review contract differs")

    profiles: dict[str, dict[str, object]] = {}
    wave_payload, profiles["prepared_wave"], _ = read_r17_bound_artifact(
        qualification["prepared_wave"], "R17 prepared wave"
    )
    evidence_payload, profiles["qualification_evidence"], _ = read_r17_bound_artifact(
        qualification["qualification_evidence"], "R17 qualification evidence"
    )
    scientific_payload, profiles["scientific_evidence"], _ = read_r17_bound_artifact(
        qualification["scientific_evidence"], "R17 scientific evidence"
    )
    wave = parse_bound_json(wave_payload, "R17 prepared wave")
    evidence = parse_bound_json(evidence_payload, "R17 qualification evidence")
    scientific = parse_bound_json(scientific_payload, "R17 scientific evidence")
    provenance = wave.get("provenance")
    packets = [
        packet
        for wave_record in wave.get("waves", [])
        if isinstance(wave_record, dict)
        for packet in wave_record.get("packets", [])
        if isinstance(packet, dict)
        and isinstance(packet.get("execution_intent"), dict)
        and packet["execution_intent"].get("case_id") == R17
    ] if isinstance(wave.get("waves"), list) else []
    if not isinstance(provenance, dict) or len(packets) != 1:
        raise ValueError("R17 prepared wave does not retain one exact R17 packet")
    packet = packets[0]
    intent = packet["execution_intent"]
    source = provenance.get("source")
    source_bundle = provenance.get("source_bundle")
    matrix = provenance.get("matrix")
    executable = provenance.get("executable")
    build = provenance.get("build_manifest")
    if (
        wave.get("provenance_sha256") != frozen["provenance_sha256"]
        or compact_json_sha256(provenance) != frozen["provenance_sha256"]
        or not all(isinstance(item, dict) for item in (
            source, source_bundle, matrix, executable, build
        ))
        or source.get("revision") != frozen["source_revision"]
        or source_bundle.get("sha256") != frozen["source_bundle_sha256"]
        or matrix.get("sha256") != frozen["matrix_sha256"]
        or executable.get("revision") != frozen["executable_revision"]
        or executable.get("sha256") != frozen["executable_sha256"]
        or build.get("inventory_sha256") != frozen["build_manifest_inventory_sha256"]
        or packet.get("execution_intent_sha256") != frozen["execution_intent_sha256"]
        or intent.get("execution_intent_sha256") != frozen["execution_intent_sha256"]
        or intent.get("execution_contract_sha256") != frozen["execution_contract_sha256"]
        or any(intent.get(key) != frozen[key] for key in (
            "case_id", "case_name", "profile_class", "target_time",
            "scientific_policy", "run_basename",
        ))
        or not isinstance(intent.get("input"), dict)
        or intent["input"].get("sha256") != frozen["input_sha256"]
    ):
        raise ValueError("R17 prepared wave differs from frozen science/build contract")

    require_exact_keys(
        evidence,
        (
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        ),
        "R17 qualification evidence",
    )
    prepared_evidence = evidence["prepared_wave"]
    if (
        not isinstance(prepared_evidence, dict)
        or set(prepared_evidence) != {"path", "sha256", "size_bytes"}
        or evidence["schema_version"] != 2
        or evidence["record_type"] != "cgl_lf_stage_i_qualification_evidence"
        or evidence["project_root"] != str(CANONICAL_ROOT)
        or evidence["case_id"] != R17
        or evidence["target_time"] != 0.25
        or prepared_evidence["path"] != str(expected_path(qualification["prepared_wave"]["path"]))
        or prepared_evidence["sha256"] != qualification["prepared_wave"]["sha256"]
    ):
        raise ValueError("R17 qualification evidence binding differs")

    output_sha256 = validate_r17_rank_inventory(
        qualification["rank_local_outputs"], "R17 qualification terminal outputs"
    )
    restart_sha256 = validate_r17_rank_inventory(
        qualification["rank_local_restarts"], "R17 qualification terminal restarts"
    )
    if (
        qualification["rank_local_output_inventory_sha256"] != output_sha256
        or qualification["rank_local_restart_inventory_sha256"] != restart_sha256
    ):
        raise ValueError("R17 terminal rank-local inventory digest differs")
    validate_r17_meshblock_decomposition_proof(
        qualification["decomposition_evidence"], output_sha256
    )
    build_inventory = qualification["build_manifest_inventory"]
    if (
        not isinstance(build_inventory, list)
        or build_inventory != sorted(build_inventory, key=lambda item: str(item.get("name", "")))
        or qualification["build_manifest_inventory_sha256"]
        != recost_json_sha256(build_inventory)
        or qualification["build_manifest_inventory_sha256"] != profile["build_manifest_sha256"]
    ):
        raise ValueError("R17 build-manifest inventory digest differs")
    for item in build_inventory:
        if not isinstance(item, dict):
            raise ValueError("R17 build-manifest inventory record must be an object")
        require_exact_keys(item, ("name", "mode", "sha256"), "R17 build-manifest inventory record")
        name = require_nonempty_string(item["name"], "R17 build-manifest file name")
        payload, retained = read_stable_regular_file(
            build_manifest_path() / name, f"R17 build-manifest file {name}"
        )
        if (
            not payload
            or item["mode"] != "0644"
            or retained["mode"] != "0644"
            or retained["sha256"] != item["sha256"]
        ):
            raise ValueError("R17 build-manifest inventory bytes or mode differ")
    if (
        scientific.get("schema_version") != 6
        or scientific.get("record_type")
        != "cgl_lf_stage_i_qualification_scientific_evidence"
        or scientific.get("case_id") != R17
        or scientific.get("nodes") != 8
        or scientific.get("execution_intent_sha256") != frozen["execution_intent_sha256"]
        or scientific.get("execution_contract_sha256") != frozen["execution_contract_sha256"]
        or scientific.get("terminal_rank_local_outputs") != qualification["rank_local_outputs"]
        or scientific.get("terminal_rank_local_output_inventory_sha256") != output_sha256
        or scientific.get("terminal_rank_local_restarts") != qualification["rank_local_restarts"]
        or scientific.get("terminal_rank_local_restart_inventory_sha256") != restart_sha256
        or scientific.get("r17_decomposition") != qualification["decomposition_evidence"]
        or scientific.get("accepted_for_operational_qualification") is not True
        or scientific.get("accepted_for_profile_selection") is not False
    ):
        raise ValueError("R17 scientific evidence differs from qualification contract")
    physics_measurements = validate_r17_physics_contract(scientific)

    validation_records = {}
    expected_measurements = {
        "restart_load_evidence": {
            "rank_local_restart_inventory_sha256": restart_sha256,
            "loaded_rank_count": 64,
            "load_state": "COMPLETED",
            "load_exit_code": "0:0",
        },
        "physics_validation_evidence": {
            "rank_local_output_inventory_sha256": output_sha256,
            **physics_measurements,
        },
    }
    record_types = {
        "restart_load_evidence": "stage-i-r17-restart-load-evidence",
        "physics_validation_evidence": "stage-i-r17-physics-validation-evidence",
    }
    measurement_authors = set()
    for key in ("restart_load_evidence", "physics_validation_evidence"):
        payload, profiles[key], relative = read_r17_bound_artifact(
            qualification[key], f"R17 {key.replace('_', ' ')}", mode="0444"
        )
        expected_relative = Path(
            "accounting/"
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_{key}.json"
        )
        record = parse_bound_json(payload, f"R17 {key.replace('_', ' ')}")
        require_exact_keys(
            record,
            (
                "schema_version", "record_type", "execution_epoch", "measured_utc",
                "measured_by", "job_id", "executable_sha256",
                "build_manifest_inventory_sha256", "passed", "measurements",
            ),
            f"R17 {key.replace('_', ' ')}",
        )
        author = require_nonempty_string(record["measured_by"], f"R17 {key} author")
        measurement_authors.add(author)
        if (
            relative != expected_relative
            or record["schema_version"] != 1
            or record["record_type"] != record_types[key]
            or record["execution_epoch"] != EXECUTION_EPOCH
            or record["measured_utc"] != qualification["measured_utc"]
            or author != measured_by
            or record["job_id"] != qualification["job_id"]
            or record["executable_sha256"] != qualification["executable_sha256"]
            or record["build_manifest_inventory_sha256"]
            != qualification["build_manifest_inventory_sha256"]
            or record["passed"] is not True
            or record["measurements"] != expected_measurements[key]
        ):
            raise ValueError("R17 retained validation evidence differs")
        validation_records[key] = {
            **qualification[key],
            "measured_utc": measured.isoformat(),
            "measured_by": author,
            "measurements": record["measurements"],
        }
    return validation_records, profiles, measurement_authors


def validate_recost_r17_readiness(
    recost: dict[str, object],
    authority: dict[str, object],
    state: dict[str, object],
    profile: dict[str, object],
    budget: dict[str, object],
) -> dict[str, object]:
    """Require the sole strong recost-compatible R17 readiness publication chain."""

    readiness = recost.get("r17_readiness")
    provenance = recost.get("provenance")
    if not isinstance(readiness, dict) or not isinstance(provenance, dict):
        raise ValueError("R17 recost recommendation lacks strong readiness publication")
    expanded = dict(readiness)
    qualification_expanded = expanded.pop("operational_qualification_evidence", None)
    publication_chain = expanded.pop("publication_chain", None)
    path = r17_readiness_path()
    value, evidence = read_json_file(path, "R17 readiness evidence")
    readiness_sha256 = require_sha256(
        provenance.get("r17_readiness_evidence_sha256"), "R17 readiness evidence SHA-256"
    )
    require_immutable_publication_evidence(evidence, readiness_sha256, "R17 readiness evidence")
    if value != expanded:
        raise ValueError("embedded R17 readiness differs from its promoted artifact")
    require_exact_keys(
        value,
        (
            "schema_version", "record_type", "execution_epoch", "root", "generated_utc",
            "expires_utc", "reviewed_by", "predecessor_lineages_sha256",
            "storage_evidence_sha256", "computed_projection_sha256", "executable_sha256",
            "build_manifest_sha256", "required_retained_bytes", "nodes", "ranks",
            "storage_ready", "node_hour_ready", "rank_64_ready",
            "operational_qualification", "operational_qualification_review",
            "current_source_authority",
        ),
        "R17 readiness evidence",
    )
    generated = require_fresh(evidence, value["generated_utc"], "R17 readiness", R17_READINESS_MAX_AGE)
    expires = require_utc(value["expires_utc"], "R17 readiness expiry")
    lineage_sha = authenticated_lineages_sha256(state["lineages"])
    if (
        value["schema_version"] != 1
        or value["record_type"] != "stage-i-r17-readiness"
        or value["execution_epoch"] != EXECUTION_EPOCH
        or value["root"] != str(CANONICAL_ROOT)
        or expires <= generated
        or expires - generated > R17_READINESS_MAX_AGE
        or current_utc() > expires
        or expires < authority["expires_utc"]
        or value["predecessor_lineages_sha256"] != lineage_sha
        or value["storage_evidence_sha256"] != provenance.get("storage_evidence_sha256")
        or value["computed_projection_sha256"] != provenance.get("computed_projection_sha256")
        or value["current_source_authority"] != provenance.get("source_authority")
        or value["executable_sha256"] != profile["executable_sha256"]
        or value["build_manifest_sha256"] != profile["build_manifest_sha256"]
        or value["required_retained_bytes"] < R17_REQUIRED_RETENTION_BYTES
        or value["nodes"] != 8
        or value["ranks"] != 64
        or any(
            value[key] is not True
            for key in ("storage_ready", "node_hour_ready", "rank_64_ready")
        )
        or available_storage_bytes(CANONICAL_ROOT) < value["required_retained_bytes"]
        or recost_json_sha256(budget) != value["computed_projection_sha256"]
    ):
        raise ValueError("R17 readiness reviewed bindings or 24-hour launch window differ")
    readiness_reviewer = require_nonempty_string(value["reviewed_by"], "R17 readiness reviewer")

    review_path = recost_independent_review_path(path)
    audit_path = recost_publication_audit_path(path)
    review, review_evidence = read_json_file(review_path, "R17 readiness independent review")
    audit, audit_evidence = read_json_file(audit_path, "R17 readiness publication audit")
    require_immutable_publication_evidence(
        review_evidence, review_evidence["sha256"], "R17 readiness independent review"
    )
    require_immutable_publication_evidence(
        audit_evidence, audit_evidence["sha256"], "R17 readiness publication audit"
    )
    require_exact_keys(
        review,
        (
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate",
        ),
        "R17 readiness independent review",
    )
    reviewed = require_utc(review["reviewed_utc"], "R17 readiness review timestamp")
    require_exact_keys(
        audit,
        (
            "schema_version", "record_type", "execution_epoch", "published_utc",
            "artifact", "independent_review", "authority",
        ),
        "R17 readiness publication audit",
    )
    published = require_utc(audit["published_utc"], "R17 readiness publication timestamp")
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-r17-readiness-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or review["reviewer"] != readiness_reviewer
        or review["candidate"] != {"path": str(path), "sha256": readiness_sha256}
        or audit["schema_version"] != 1
        or audit["record_type"] != "stage-i-r17-readiness-publication-audit"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"]
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or reviewed < generated
        or published < reviewed
        or authority["generated_utc"] < published
        or current_utc() - reviewed > R17_READINESS_MAX_AGE
        or current_utc() - published > R17_READINESS_MAX_AGE
    ):
        raise ValueError("R17 readiness review/publication chain or freshness differs")
    validate_declared_publication_binding(
        audit["artifact"],
        expected=path,
        digest=readiness_sha256,
        mode="0444",
        label="R17 readiness publication artifact",
    )
    validate_declared_publication_binding(
        audit["independent_review"],
        expected=review_path,
        digest=review_evidence["sha256"],
        mode="0444",
        label="R17 readiness publication review",
    )
    expected_chain = {
        "independent_review_sha256": review_evidence["sha256"],
        "publication_audit_sha256": audit_evidence["sha256"],
        "reviewed_by": readiness_reviewer,
        "published_utc": published.isoformat(),
    }
    if publication_chain != expected_chain:
        raise ValueError("embedded R17 readiness publication chain differs")

    def read_relative_binding(binding: object, label: str, mode: str) -> tuple[dict[str, object], dict[str, object]]:
        if not isinstance(binding, dict) or set(binding) != {"path", "sha256"}:
            raise ValueError(f"{label} binding differs")
        relative = Path(require_nonempty_string(binding["path"], f"{label} path"))
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"{label} path must be canonical-root relative")
        retained_path = expected_path(relative.as_posix())
        retained, retained_evidence = read_json_file(retained_path, label)
        if retained_evidence["sha256"] != binding["sha256"] or retained_evidence["mode"] != mode:
            raise ValueError(f"{label} bytes/mode differ")
        return retained, retained_evidence

    qualification, qualification_evidence = read_relative_binding(
        value["operational_qualification"], "R17 operational qualification", "0444"
    )
    qualification_review, qualification_review_evidence = read_relative_binding(
        value["operational_qualification_review"],
        "R17 operational qualification independent review",
        "0444",
    )
    if not isinstance(qualification_expanded, dict):
        raise ValueError("recost omits expanded R17 operational qualification")
    raw_qualification_keys = {
        "schema_version", "record_type", "execution_epoch", "completed_utc",
        "measured_utc", "measured_by", "job_id", "state", "exit_code", "nodes",
        "ranks", "prepared_wave", "qualification_evidence", "scientific_evidence",
        "scheduler_evidence", "account_scheduler_evidence",
        "account_exclusivity_evidence", "executable_sha256",
        "build_manifest_inventory", "build_manifest_inventory_sha256",
        "rank_local_outputs", "rank_local_output_inventory_sha256",
        "rank_local_restarts", "rank_local_restart_inventory_sha256",
        "decomposition_evidence", "restart_load_evidence",
        "physics_validation_evidence", "frozen_science_build_contract",
        "independent_review_contract", "authority",
    }
    expanded_qualification_keys = raw_qualification_keys | {
        "authenticated_rank_local_outputs",
        "authenticated_rank_local_restarts",
        "authenticated_decomposition_evidence",
        "authenticated_account_exclusivity_evidence",
        "authenticated_validation_evidence",
        "reviewed_by",
    }
    require_exact_keys(
        qualification, raw_qualification_keys, "R17 operational qualification"
    )
    require_exact_keys(
        qualification_expanded,
        expanded_qualification_keys,
        "expanded R17 operational qualification",
    )
    require_exact_keys(
        qualification_review,
        (
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate",
        ),
        "R17 operational qualification independent review",
    )
    expanded_base = {
        key: item
        for key, item in qualification_expanded.items()
        if key
        not in {
            "authenticated_rank_local_outputs",
            "authenticated_rank_local_restarts",
            "authenticated_decomposition_evidence",
            "authenticated_account_exclusivity_evidence",
            "authenticated_validation_evidence",
            "reviewed_by",
        }
    }
    qualification_reviewer = require_nonempty_string(
        qualification_review["reviewer"], "R17 operational qualification reviewer"
    )
    qualification_reviewed = require_utc(
        qualification_review["reviewed_utc"], "R17 operational qualification review timestamp"
    )
    (
        authenticated_validation_evidence,
        qualification_profiles,
        measurement_authors,
    ) = validate_r17_operational_qualification_contract(
        qualification,
        profile,
        normalized_path(Path(str(qualification_evidence["path"]))),
        qualification_evidence,
        qualification_review,
        qualification_reviewed,
    )
    (
        authenticated_account_exclusivity,
        account_exclusivity_profiles,
    ) = validate_r17_account_exclusivity_artifacts(qualification, qualification_reviewed)
    if (
        qualification_expanded["authenticated_decomposition_evidence"]
        != qualification["decomposition_evidence"]
    ):
        raise ValueError("expanded R17 decomposition evidence differs")
    if (
        qualification_expanded["authenticated_account_exclusivity_evidence"]
        != authenticated_account_exclusivity
    ):
        raise ValueError("expanded R17 account exclusivity evidence differs")
    if (
        expanded_base != qualification
        or qualification_expanded["authenticated_rank_local_outputs"]
        != qualification["rank_local_outputs"]
        or qualification_expanded["authenticated_rank_local_restarts"]
        != qualification["rank_local_restarts"]
        or qualification_expanded["authenticated_validation_evidence"]
        != authenticated_validation_evidence
        or qualification_review.get("schema_version") != 1
        or qualification_review.get("record_type")
        != "stage-i-r17-operational-qualification-independent-review"
        or qualification_review.get("execution_epoch") != EXECUTION_EPOCH
        or qualification_review.get("decision") != "approved"
        or qualification_review.get("candidate")
        != {
            "path": qualification_evidence["path"],
            "sha256": qualification_evidence["sha256"],
        }
        or qualification_expanded.get("reviewed_by") != qualification_reviewer
        or qualification_reviewer == readiness_reviewer
        or qualification_reviewer in measurement_authors
        or qualification_reviewed > current_utc() + FUTURE_SKEW
    ):
        raise ValueError("R17 operational qualification/review chain differs")
    review_assurance = declared_process_independence_assurance(
        {
            "R17-readiness-reviewer": readiness_reviewer,
            "R17-operational-qualification-reviewer": qualification_reviewer,
            **{
                f"R17-measurement-author-{index}": agent
                for index, agent in enumerate(sorted(measurement_authors))
            },
        },
        "R17 qualification/readiness review chain",
    )
    return {
        "artifact": evidence,
        "independent_review": review_evidence,
        "publication_audit": audit_evidence,
        "operational_qualification": qualification_evidence,
        "operational_qualification_review": qualification_review_evidence,
        **qualification_profiles,
        **account_exclusivity_profiles,
        "independent_review_assurance": review_assurance,
        "live_available_bytes": available_storage_bytes(CANONICAL_ROOT),
    }


def acceptance_prose(case_id: str, target: Decimal, codes: tuple[str, ...]) -> str:
    """Return expanded exact-boundary acceptance prose."""

    clauses = " ".join(f"{code}: {POLICY_TEXT[code]}" for code in codes)
    return (
        f"Accept {case_id} only at exact t={decimal_text(target)} under policies "
        f"{'+'.join(codes)}. {clauses} Any endpoint, provenance, product-inventory, "
        "scheduler-accounting, or policy failure blocks acceptance and successor "
        "packet planning."
    )


def build_packet(
    case_id: str,
    case: dict[str, object],
    profile: dict[str, object],
    next_index: int,
    start: Decimal,
    static: dict[str, object],
    continuation: dict[str, object] | None,
    authorization: dict[str, object] | None = None,
) -> dict[str, object]:
    """Build one exact command-free, explicitly non-launching packet specification."""

    segment = require_nonempty_string(profile["segment"], f"{case_id} recost profile segment")
    segment_index, segment_start, target = parse_segment(segment, f"{case_id} recost profile segment")
    if (
        segment_index != next_index
        or abs(segment_start - start) > Decimal("0.000001")
        or decimal_value(profile["time_tlim_target"], f"{case_id} recost profile target") != target
    ):
        raise ValueError(f"{case_id} recost profile differs from packet lineage")
    codes = policy_codes(case_id)
    if authorization is not None:
        sole = authorization.get("sole_next_segment_profile")
        if not isinstance(sole, dict):
            raise ValueError(f"{case_id} authorization lacks a sole-next profile")
        authorized_segment = require_nonempty_string(
            sole.get("segment"), f"{case_id} authorized segment"
        )
        authorized_index, authorized_start, authorized_target = parse_segment(
            authorized_segment, f"{case_id} authorized segment"
        )
        if (
            authorized_index != next_index
            or abs(authorized_start - start) > Decimal("0.000001")
            or authorized_target != target
        ):
            raise ValueError(f"{case_id} authorized segment differs from the packet lineage")
        segment = authorized_segment
    if continuation is None and start != 0:
        raise ValueError(f"{case_id} continuation lacks authenticated restart provenance")
    production_provenance = {
        key: value for key, value in static.items() if key != "inputs"
    }
    production_provenance["frozen_input"] = static["inputs"][case_id]
    if authorization is not None:
        production_provenance["controller_helper"] = {
            "path": str(CONTROLLER_HELPER),
            "revision": R03_F115_CONTROLLER_REVISION,
            "sha256": R03_F115_CONTROLLER_SHA256,
        }
        production_provenance["source_bundle"] = {
            "path": str(f115_source_bundle_path()),
            "sha256": R03_F115_SOURCE_BUNDLE_SHA256,
            "verified_revisions": list(F115_SOURCE_BUNDLE_REQUIRED_REVISIONS),
            "independent_validation": authorization["source_bundle_validation"],
        }
    return {
        "case_id": case_id,
        "case_name": case["name"],
        "input": case["input"],
        "resolution": case["resolution"],
        "model_role": model_role(case_id),
        "physics_policy_codes": list(codes),
        "acceptance_criterion": acceptance_prose(case_id, target, codes),
        "authorization_state": (
            "planning_only_f115_profile_reference_non_authorizing"
            if authorization is not None
            else "planning_only_reviewed_recost_profile_non_authorizing"
        ),
        "authorization_evidence": authorization,
        "lineage": {
            "kind": "fresh" if start == 0 else "continuation",
            "segment": segment,
            "start_time": decimal_text(start),
            "target_time": decimal_text(target),
            "restart_required": start != 0,
            "continuation_provenance": continuation,
        },
        "allocation": {
            "nodes": profile["nodes"],
            "ranks_per_node": RANKS_PER_NODE,
            "total_ranks": profile["nodes"] * RANKS_PER_NODE,
            "cpus_per_task": CPUS_PER_TASK,
            "walltime": profile["walltime"],
            "athena_walltime": profile["athena_walltime"],
            "estimated_storage_bytes": profile["estimated_storage_bytes"],
            "recommendation_basis": profile["recommendation_basis"],
        },
        "runtime_parameters": {"time_tlim": decimal_text(target)},
        "production_provenance": production_provenance,
    }


def packet_reserved_node_hours(packet: dict[str, object]) -> Decimal:
    """Return one packet's exact requested node-hours."""

    allocation = packet["allocation"]
    return Decimal(int(allocation["nodes"]) * walltime_seconds(
        allocation["walltime"], f"{packet['case_id']} packet walltime"
    )) / Decimal(3600)


def build_wave_plan(
    matrix: dict[str, object],
    matrix_evidence: dict[str, object],
    authority: dict[str, object],
) -> dict[str, object]:
    """Build one evidence-bound deterministic read-only advisory diagnostic."""

    cases, input_evidence = validate_matrix(matrix, matrix_evidence)
    recost = authority["artifact"]
    assert isinstance(recost, dict)
    static = validate_static_provenance(recost, matrix_evidence, input_evidence)
    state = validate_canonical_state(cases, authority, static)
    lineages = state["lineages"]

    r03_lineage, r03_active = lineages[R03]
    r03_pending_f115_sole, r03_authority = validate_r03_authorization(
        r03_lineage, r03_active, state["case_infos"][R03]
    )
    profiles, profiles_by_case = validate_recost_profiles(recost, state, static)
    r12_fresh_rerun = validate_r12_fresh_rerun_transition(state, profiles_by_case)
    budget = validate_recost_budget(recost, state, profiles, cases)
    predecessor_state = {
        case_id: lineage_complete(lineages[case_id][0])
        for case_id in tuple(f"R{index:02d}" for index in range(2, 17))
    }
    r17_lineage, r17_active = lineages[R17]
    if r17_lineage or r17_active:
        if not all(predecessor_state.values()):
            raise ValueError("R17 started before every authenticated R02-R16 lineage reached exact t=10")
    if R17 not in profiles_by_case and (r17_lineage or r17_active):
        raise ValueError("R17 has started; a lower-resolution recost recommendation is forbidden")

    r17_readiness = (
        validate_recost_r17_readiness(
            recost, authority, state, profiles_by_case[R17], budget
        )
        if R17 in profiles_by_case
        else None
    )
    packets = []
    for profile_value in profiles:
        case_id = str(profile_value["case_id"])
        lineage, active = lineages[case_id]
        if active is not None or lineage_complete(lineage):
            raise ValueError(f"recost recommends occupied or complete case {case_id}")
        endpoint = lineage_endpoint(lineage)
        continuation = terminal_restart(lineage[-1]) if lineage else None
        authorization = None
        if case_id == R03 and r03_pending_f115_sole is not None:
            if r03_authority is None:
                raise ValueError("recost R03 recommendation lacks exact F115 sole-next authority")
            compatible = {
                "case_id": profile_value["case_id"],
                "segment": profile_value["segment"],
                "nodes": profile_value["nodes"],
                "walltime": profile_value["walltime"],
                "athena_walltime": profile_value["athena_walltime"],
                "cpus_per_task": profile_value["cpus_per_task"],
                "ranks_per_node": profile_value["ranks_per_node"],
                "executable": profile_value["executable"],
                "executable_revision": profile_value["executable_revision"],
                "executable_sha256": profile_value["executable_sha256"],
                "input_sha256": profile_value["input_sha256"],
                "parent_job_id": profile_value["parent_job_id"],
                "parent_result": profile_value["parent_result"],
                "parent_segment": profile_value["parent_segment"],
                "restart_file": profile_value["restart_file"],
                "restart_time": profile_value["restart_time"],
                "source_bundle": profile_value["source_bundle"],
                "source_bundle_sha256": profile_value["source_bundle_sha256"],
                "time_tlim_target": profile_value["time_tlim_target"],
            }
            if any(
                compatible.get(key) != value
                for key, value in r03_pending_f115_sole.items()
                if key in compatible
            ):
                raise ValueError("recost R03 profile differs from immutable F115 sole-next profile")
            continuation = r03_authority["continuation"]
            authorization = {
                **{
                    key: value for key, value in r03_authority.items()
                    if key != "continuation"
                },
                "sole_next_segment_profile": r03_pending_f115_sole,
            }
        packets.append(
            build_packet(
                case_id,
                cases[case_id],
                profile_value,
                int(state["next_indexes"][case_id]),
                endpoint,
                static,
                continuation,
                authorization,
            )
        )

    wave_status = "r17_exclusive" if R17 in profiles_by_case else "drained_barrier_wave"
    total_nodes = sum(int(packet["allocation"]["nodes"]) for packet in packets)
    total_lanes = len(packets)
    planned_reserved_node_hours = sum(
        (packet_reserved_node_hours(packet) for packet in packets), Decimal("0")
    )
    total_committed_node_hours = state["actual_node_hours"] + planned_reserved_node_hours
    if (
        total_committed_node_hours > STAGE_I_BUDGET_NODE_HOURS
        or total_committed_node_hours > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError("planned wave exceeds a controller budget ceiling")
    if wave_status in {"r17_active", "r17_exclusive"}:
        if total_lanes != 1 or total_nodes != 8:
            raise ValueError("R17 plan is not exclusive on eight nodes")
    elif total_lanes > MAX_LANES or total_nodes > MAX_NODES:
        raise ValueError("planned wave exceeds reviewed concurrency ceilings")
    packet_cases = [packet["case_id"] for packet in packets]
    if len(packet_cases) != len(set(packet_cases)):
        raise ValueError("planned wave contains duplicate case lanes")

    core = {
        "schema_version": 4,
        "record_type": "cgl_lf_stage_i_read_only_wave_diagnostic",
        "execution_epoch": EXECUTION_EPOCH,
        "authority": (
            "planning-only advisory packet specifications; no prepare, shared-root "
            "acknowledgement, check-submit, submit, or launch authority"
        ),
        "read_only": True,
        "disclosures": [
            "This diagnostic never prepares, acknowledges shared-root isolation, "
            "check-submits, submits, cancels, or records a job.",
            "Published F-115 is procedural pre-prepare authority for only its exact R03 "
            "s02 profile; this planner merely authenticates and references that authority "
            "and does not extend or exercise it.",
            "F-115 does not authorize shared-root acknowledgement; a separate exact "
            "isolation review is required after prepare before such acknowledgement may "
            "be considered.",
            "Live Slurm authentication is an instantaneous observation and is rechecked "
            "before emission; the production controller must authenticate state again.",
            "Restart bytes, independent profile publication, transactions, and the full "
            "remaining-campaign budget are authenticated here, but this output grants no "
            "execution permission.",
            "The replacement source bundle authenticates the committed controller "
            "implementation at its bound head; the qualified executable and frozen "
            "scientific inputs remain separately bound to their unchanged source revision "
            "and bytes.",
            "The planner states future scientific review criteria but does not authenticate "
            "a CT/divergence metric, set a numerical CT threshold, or accept a run.",
            "Every emitted packet requires separate controller/checkpoint consumption of "
            "its underlying authority before launch.",
            "The promoted schema-2 recost artifact and its independent review/publication "
            "audit are the sole current profile, projection, and reconciliation authority.",
            "The exact F-118 evidence/reviews/publication-audit chain is the sole current "
            "helper/source-bundle authority; F-115 remains historical R03 s02 authority only.",
            "Independent-review separation is authenticated only as exact retained process, "
            "role, and agent declarations; artifact digests do not cryptographically "
            "authenticate reviewer identity.",
            "Retained R12 s00/job 4766856 is authenticated only as historical inventory. "
            "It cannot authorize continuation or a waiver; R12 resumes only as fresh "
            "s01_rankio_t0_t0p12 with no parent or restart.",
        ],
        "evidence": {
            "matrix": matrix_evidence,
            "recost_authority": {
                key: value.isoformat() if isinstance(value, datetime) else value
                for key, value in authority.items()
            },
            "canonical_state": state["evidence"],
            "r03_f115_authority": r03_authority,
            "r12_fresh_rerun_transition": r12_fresh_rerun,
            "r17_readiness": r17_readiness,
        },
        "policy": {
            "max_lanes": MAX_LANES,
            "max_nodes": MAX_NODES,
            "max_segment_seconds": MAX_SEGMENT_SECONDS,
            "shutdown_margin_seconds": SHUTDOWN_MARGIN_SECONDS,
            "ranks_per_node": RANKS_PER_NODE,
            "cpus_per_task": CPUS_PER_TASK,
            "recost_bounded_concurrency": recost["recommendations"]["bounded_concurrency"],
            "r17_exclusive_last": True,
        },
        "physics_policies": POLICY_TEXT,
        "observed_state": {
            "updated_utc": authority["generated_utc"].isoformat(),
            "reconciliation_counts": state["counts"],
            "actual_node_hours": decimal_text(state["actual_node_hours"]),
            "active_reserved_node_hours": "0",
            "planned_reserved_node_hours": decimal_text(planned_reserved_node_hours),
            "total_committed_node_hours": decimal_text(total_committed_node_hours),
            "credible_remaining_campaign_projection": budget,
            "active_lanes": [],
            "predecessors_accepted_exact_t10": predecessor_state,
        },
        "wave": {
            "status": wave_status,
            "active_lane_count": 0,
            "active_nodes": 0,
            "planned_lane_count": len(packets),
            "planned_nodes": sum(int(packet["allocation"]["nodes"]) for packet in packets),
            "total_lane_count": total_lanes,
            "total_nodes": total_nodes,
            "packets": packets,
            "deferred": [],
        },
    }
    core["evidence"]["global_toctou_boundary"] = revalidate_global_boundary(
        [
            matrix_evidence,
            authority,
            static,
            state,
            r03_authority,
            r12_fresh_rerun,
            r17_readiness,
            packets,
        ],
        state["transaction_state"],
        state["submitted_scheduler"],
    )
    return {**core, "plan_sha256": sha256_bytes(canonical_json(core))}


def plan_from_paths(
    recost_path_value: Path,
    review_path_value: Path,
    audit_path_value: Path,
) -> dict[str, object]:
    """Read all exact evidence and build one non-authorizing advisory diagnostic."""

    matrix, matrix_evidence = read_json_file(frozen_matrix_path(), "matrix")
    authority = validate_recost_publication_chain(
        recost_path_value, review_path_value, audit_path_value
    )
    return build_wave_plan(matrix, matrix_evidence, authority)


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--recost-artifact", required=True, type=Path)
    command.add_argument("--recost-independent-review", required=True, type=Path)
    command.add_argument("--recost-publication-audit", required=True, type=Path)
    return command


def main(argv: list[str] | None = None) -> int:
    """Read exact canonical evidence and print one command-free JSON plan."""

    args = parser().parse_args(argv)
    try:
        plan = plan_from_paths(
            args.recost_artifact,
            args.recost_independent_review,
            args.recost_publication_audit,
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Stage I wave planner failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(plan, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
