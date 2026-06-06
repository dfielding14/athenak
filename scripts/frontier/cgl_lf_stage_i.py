#!/usr/bin/env python3
"""Prepare and account Frontier MKS24 Stage I production segments.

This utility is separate from ``cgl_lf_frontier.py`` because that launcher is
restricted to non-production ``debug`` QOS work.  Stage I jobs use Frontier's
``batch`` partition with its default production QOS, are prepared and
submitted one segment at a time, and must belong to the tracked mapped-case
manifest.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import ctypes
import csv
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation
import errno
import fcntl
from functools import wraps
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import pwd
import re
import shlex
import shutil
import stat
import struct
import subprocess
import sys
import tempfile
import uuid


ROOT_DIR = Path(__file__).resolve().parents[2]
PRODUCTION_UTILITY_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i.py")
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
TRUSTED_GROUP_WRITABLE_PROJECT_BOUNDARY = Path(
    "/lustre/orion/ast207/proj-shared"
)
TRUSTED_GROUP_WRITABLE_PROJECT_MODE = 0o2770
TRUSTED_GROUP_WRITABLE_PROJECT_UID = 0
TRUSTED_GROUP_WRITABLE_PROJECT_GID = 31114
DURABLE_DIRECTORY_MODE = 0o755
SACCT = Path("/usr/bin/sacct")
SBATCH = Path("/usr/bin/sbatch")
SCONTROL = Path("/usr/bin/scontrol")
SQUEUE = Path("/usr/bin/squeue")
GIT = Path("/usr/lib/git/git")
SYSTEM_PYTHON = Path("/usr/bin/python3.11")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
DEFAULT_MATRIX = ROOT_DIR / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
ACCOUNT = "AST207"
PARTITION = "batch"
CGL_JOB_NAME_PREFIX = "cgl_"
PRODUCTION_QOS = "normal (Frontier default; no -q directive)"
PROJECT_BUDGET_NODE_HOURS = 4000.0
HISTORICAL_DEBUG_NODE_HOURS = 0.851670
HISTORICAL_E01_STAGE_I_NODE_HOURS = 9.962778
HISTORICAL_E02_PIPELINE_NODE_HOURS = 15.628610
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
AUTHORIZED_CASE_IDS = frozenset(f"R{number:02d}" for number in range(2, 18))
R17_CASE_ID = "R17"
R17_PREDECESSOR_CASE_IDS = tuple(
    f"R{number:02d}" for number in range(2, 17)
)
R17_READINESS_RELATIVE = Path(
    "accounting/mks24_stage_i_E03_forcing_policy_R17_readiness_evidence.json"
)
F116_CURRENT_SOURCE_AUTHORITY_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F116_current_source_authority_supersession_evidence.json"
)
F116_PROVENANCE_REVIEW_RELATIVE = Path(
    f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json"
)
F116_PLASMA_REVIEW_RELATIVE = Path(
    f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json"
)
F116_PUBLICATION_AUDIT_RELATIVE = Path(
    f"{F116_CURRENT_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json"
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
F118_PROVENANCE_REVIEW_RELATIVE = Path(
    f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json"
)
F118_PLASMA_REVIEW_RELATIVE = Path(
    f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json"
)
F118_PUBLICATION_AUDIT_RELATIVE = Path(
    f"{F118_CURRENT_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json"
)
SHARED_ROOT_STALE_CAMPAIGN_ID = "beta25-accel05-gamma10001-purecgl-256"
SHARED_ROOT_STALE_JOB_ID = "4743106"
SHARED_ROOT_STALE_JOB_NAME = "b25_a05_g10001_pcgl_s00"
SHARED_ROOT_CLEARANCE_REFRESH_AUDIT_PATTERN = re.compile(
    rf"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F(?P<checkpoint>[0-9]+)_"
    r"shared_root_isolation_clearance_refresh\.json\.publication_audit\.json"
)
SHARED_ROOT_CLEARANCE_REFRESH_AUDIT_REVIEW_PATTERN = re.compile(
    rf"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F(?P<checkpoint>[0-9]+)_"
    r"shared_root_isolation_clearance_refresh\.json\.publication_audit\.json"
    r"\.independent_review\.json"
)
SHARED_ROOT_CLEARANCE_SUPERSESSION_RELATIVE = Path(
    "accounting/"
    f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_"
    "standing_shared_root_isolation_clearance_supersession.json"
)
SHARED_ROOT_CLEARANCE_LEGACY_RELATIVES = tuple(
    Path("accounting") / name
    for name in (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance.json.independent_review.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance.json.publication_audit.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance.json.publication_audit.json.independent_review.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance_supersession.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance_supersession.json.independent_review.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance_supersession.json.publication_audit.json",
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_standing_shared_root_isolation_clearance_supersession.json.publication_audit.json.independent_review.json",
    )
)
SHARED_ROOT_CLEARANCE_SOURCE_AUTHORITY_RELATIVES = (
    F118_CURRENT_SOURCE_AUTHORITY_RELATIVE,
    F118_PROVENANCE_REVIEW_RELATIVE,
    F118_PLASMA_REVIEW_RELATIVE,
    F118_PUBLICATION_AUDIT_RELATIVE,
)
F116_SOURCE_AUTHORITY_TRANSACTION_ID_PATTERN = re.compile(
    r"\d{4}-\d{2}-\d{2}T\d{6}\+0000-[0-9a-f]{32}"
)
F116_SOURCE_AUTHORITY_STAGING_SUFFIX = ".staging"
F116_SOURCE_AUTHORITY_RETIRED_SUFFIX = ".retired"
F116_SOURCE_AUTHORITY_FORENSIC_ENTRY_PATTERN = re.compile(
    r"\.cgl-source-authority-retired-(?:directory-)?[0-9a-f]{32}\.forensic"
)
F116_SOURCE_AUTHORITY_TRANSACTION_PAYLOADS = {
    "bundle": ("final.bundle", "0644"),
    "evidence": ("evidence.json", "0444"),
    "provenance_review": ("provenance_review.json", "0444"),
    "plasma_review": ("plasma_review.json", "0444"),
    "audit": ("audit.json", "0444"),
    "readme_before": ("README.before", "0644"),
    "sha256sums_before": ("SHA256SUMS.before", "0644"),
    "readme_after": ("README.after", "0644"),
    "sha256sums_after": ("SHA256SUMS.after", "0644"),
}
F116_SOURCE_AUTHORITY_JOURNAL_RECOVERY_NAMES = frozenset({
    ".journal.json.recovery.tmp",
    ".journal.json.recovery.alternate.tmp",
})
R03_F115_SOURCE_AUTHORITY_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F115_source_bundle_recovery_supersession_evidence.json"
)
R03_F115_SOURCE_AUTHORITY_SHA256 = (
    "cb50beb064678a9446ac33801a0023547d06c432d8c59b6bd3fbe34b11cf0391"
)
R03_F115_PUBLICATION_AUDIT_SHA256 = (
    "5923e3872b1d4a84a147bcd1d81fddc20ee79b72bfc410683d083039781f5b1f"
)
R03_F115_PROVENANCE_REVIEW_SHA256 = (
    "6fcd19f9267f36332742f1f103968098216bd8e6b42fa9821964ab8df704bacd"
)
R03_F115_PLASMA_REVIEW_SHA256 = (
    "a78357ed90e593809b1a82a641d2b40d651b16569ec943fc1781940e569b8440"
)
R03_F115_MANIFEST_RELATIVE = Path(
    "runs/mks24-stage-i/E03-forcing-policy/R03/"
    "s02_rankio_t0p312823_t0p5/manifest/prepared_run.json"
)
R03_F115_MANIFEST_SHA256 = (
    "305414d48f3cf4dffe50c889c78a99f84c73791138e353968fdf07758147a7c1"
)
R03_F115_SEGMENT = "s02_rankio_t0p312823_t0p5"
R03_F115_JOB_ID = "4766828"
R03_F115_CONTROLLER_REVISION = "5834a91e448a69ec0df5d011b7be3fe786666806"
R03_F115_CONTROLLER_SHA256 = (
    "0a6b30a60ba70faf32dae722c4472e96acb38031d8c1a4b52b816fc93e58e53a"
)
R03_F115_SOURCE_BUNDLE_RELATIVE = Path(
    "source-archives/athenak-feature-cgl-through-5834a91e4.bundle"
)
R03_F115_SOURCE_BUNDLE_SHA256 = (
    "a6aa40f8f3350d65022be6d898d5575a60185b05bba7b72f184c3f2ebd401ec8"
)
R03_F114_CONTROLLER_REVISION = "c7e4fa30ea7162e4d5ce070a45dfec5b56ca2052"
QUALIFIED_SOURCE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
F116_BRIDGE_REVISION = "36140ea825cb853b298714c27720440fdab60b9e"
F116_BRIDGE_SHA256 = (
    "2c2f57a166877387244dd5bb6bdf87beb12492ea075a7431939b78e5df7307a0"
)
F116_BRIDGE_NAME = "athenak-feature-cgl-through-36140ea82.bundle"
F116_CORRUPT_C7_NAME = "athenak-feature-cgl-through-c7e4fa30e.bundle"
F116_PRODUCTION_REQUIRED_REVISIONS = frozenset({
    QUALIFIED_SOURCE_REVISION,
    R03_F114_CONTROLLER_REVISION,
    "b0d3e8d526d8f3b4e5977333000db24c09d4eab1",
    "38aedd2a65c3c11855721858f5dbcce20bae11e4",
    R03_F115_CONTROLLER_REVISION,
    "469d38a841ef25d5c044713071adccd55ff90fef",
    "e1f4f4a0b62c3649b4d80a25a188b991f856ebbe",
    F116_BRIDGE_REVISION,
})
F116_REQUIRED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
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
F116_PRESERVES = [
    "The immutable F-115 evidence, reviews, publication audit, and historical R03 s02 authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F116_DOES_NOT_AUTHORIZE = [
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
F118_PRESERVES = [
    "The immutable F-116 evidence, reviews, publication audit, selected source bundle, and nested F-115 historical authority.",
    "The qualified executable, frozen source revision, inputs, matrix, restart lineages, targets, resources, qualification, and Stage I budget policy.",
    "Every prior active source-archive checksum-ledger entry and the corrupt-C7 incident-evidence exclusion.",
]
F118_DOES_NOT_AUTHORIZE = F116_DOES_NOT_AUTHORIZE
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
F118_REQUIRED_TOOLS = F116_REQUIRED_TOOLS
LEGACY_F114_RECOST_NAME = (
    "mks24_stage_i_E03_forcing_policy_R03_s00_clean_partial_recost_evidence.json"
)
LEGACY_F114_RECOST_SHA256 = (
    "27dd154b8594693afe0308f4de63a7faaed6daf81a396ab31166bb70a8c43b86"
)
LEGACY_F114_PUBLICATION_AUDIT_SHA256 = (
    "3bf50ef1359f7f1798da945185d42b61488cb98779d5a6942839c03e3eb6bc73"
)
LEGACY_F114_GENERATED_UTC = "2026-06-05T02:44:51.687220+00:00"
LEGACY_F114_PUBLISHED_UTC = "2026-06-05T03:03:11+00:00"
F117_RECOST_ARTIFACT_NAME = (
    "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json"
)
F119_RECOST_ARTIFACT_NAME = (
    "mks24_stage_i_E03_forcing_policy_F119_recost_evidence.json"
)
F119_LEGACY_F114_BOOTSTRAP = (
    "exact-retained-legacy-F114-after-authenticated-F117-failed-attempt"
)
F117_FAILED_ATTEMPT_RELATIVES = {
    "packet": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_recost_draft_packet.json"
    ),
    "request": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_recost_request.json"
    ),
    "reconciliation": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_reconciliation_evidence.json"
    ),
    "storage": Path(
        "accounting/mks24_stage_i_E03_forcing_policy_F117_storage_evidence.json"
    ),
}
F117_FAILED_ATTEMPT_SHA256 = {
    "packet": "8abceb6b2e3031a21b87f95b19053f2bf2b35f30d642c7a86bc4d9c6614071ad",
    "request": "5dea85b34b4f09cb3a926927ad73b4443e2f8ce4fdd2b2a1b1f84b70294bc632",
    "reconciliation": "677547889992d700e05559e322d95e7df6078a736bb2fbdff27ada38410c83e6",
    "storage": "637e43a3cd12f3051e1eae1d3102a2106d0445e26ac500b8946f8fea2064f218",
}
F117_FORBIDDEN_PROMOTION_RELATIVES = (
    Path(f"accounting/{F117_RECOST_ARTIFACT_NAME}"),
    Path(f"accounting/{F117_RECOST_ARTIFACT_NAME}.staged"),
    Path(f"accounting/{F117_RECOST_ARTIFACT_NAME}.independent_review.json"),
    Path(f"accounting/{F117_RECOST_ARTIFACT_NAME}.publication_audit.json"),
    Path(
        "accounting/"
        "mks24_stage_i_E03_forcing_policy_F117_recost_request.json.independent_review.json"
    ),
)
R17_RECOST_PUBLICATION_AUDIT_PATTERN = re.compile(
    r"mks24_stage_i_E03_forcing_policy_F([0-9]+)_recost_evidence"
    r"\.json\.publication_audit\.json"
)
R17_RECOST_PUBLICATION_METHOD = (
    "same-directory-link-fsync-copy-exchange-forensic-retirement-fsync"
)
R17_AUTHORIZATION_MAX_LIFETIME = timedelta(hours=24)
R17_AUTHORIZATION_FUTURE_SKEW = timedelta(minutes=5)
R17_MINIMUM_RETAINED_BYTES = 958_271_710_272
R17_MAX_NORMALIZED_CT_DIVB = Decimal("1e-12")
R17_FROZEN_MATRIX_SHA256 = (
    "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
)
R17_FROZEN_INPUT_SHA256 = (
    "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04"
)
R17_FROZEN_PARAMETER_CONTRACT = {
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
CONTINUATION_MASS_TOLERANCE = 1.0e-12
CONTINUATION_MAX_NORMALIZED_CT_DIVB = 1.0e-12
CONTINUATION_ACTIVITY_ABSOLUTE_GT = 1.0e-6
CONTINUATION_ACTIVITY_NORMALIZED_GT = 1.0e-8
CONTINUATION_PLASMA_POLICY = "stage-i-clean-partial-continuation-v2"
FROZEN_E03_CONTINUATION_MIGRATION_POLICY = (
    "stage-i-frozen-e03-exact-retained-continuation-migration-v1"
)
FROZEN_E03_CT_DIVERGENCE_REASON = (
    "the qualified historical E03 executable did not retain normalized CT "
    "divB in its exact user-history schema"
)
R03_ACCEPTED_CONTINUATION_JOB_ID = "4766828"
R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID = "4766856"
R12_FRESH_RERUN_SEGMENT = "s01_rankio_t0_t0p12"
R12_FRESH_RERUN_NODES = 4
R12_FRESH_RERUN_RANKS_PER_NODE = 8
R12_FRESH_RERUN_TOTAL_RANKS = 32
R12_FRESH_RERUN_WALLTIME = "02:00:00"
R12_FRESH_RERUN_ATHENA_WALLTIME = "01:50:00"
R12_FRESH_RERUN_TARGET = 0.12
FROZEN_E03_EXECUTABLE = {
    "revision": "9e07542281e4e6d125582f253df3ad2e3b8b154d",
    "sha256": "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c",
}
FROZEN_E03_CONTINUATION_MIGRATIONS = {
    ("R03", "4762472"): {
        "segment": "s00_rankio_t0_t0p5",
        "manifest": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/manifest/prepared_run.json"
            ),
            "mode": 0o644,
            "size_bytes": 34030,
            "sha256": (
                "ad3428b51945d65bed2b9ddeb337f4972bdc5b6ec07ac0cd6105ae6368055444"
            ),
        },
        "inspection": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R03/"
                "s00_rankio_t0_t0p5/manifest/segment_inspection.json"
            ),
            "mode": 0o644,
            "size_bytes": 23551,
            "sha256": (
                "24de23eaac2805eb7a0aa6e1c86ee6f1f133b573c9815548052fc618c5d3a643"
            ),
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
            "sha256": (
                "7cb77ecf6d756b5a64c5bb02d917a5e75482df94b67f678cc34a68e05299f5af"
            ),
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
            "sha256": (
                "766a79819ad14653042f578af6aa7f0ab8feba2abf5ca640fc4009a528001551"
            ),
        },
        "independent_validation": {
            "path": "accounting/4762472.stage_i.independent_validation.json",
            "mode": 0o644,
            "size_bytes": 24565,
            "sha256": (
                "da15fc8a74013fd0f13276b3421ccc827ce0319c350cf265b4548297bfa1185e"
            ),
            "schema_version": 1,
            "record_type": "stage-i-binary-aware-clean-partial-validation",
        },
        "independent_review": None,
    },
    ("R12", "4766856"): {
        "segment": "s00_rankio_t0_t0p25",
        "manifest": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/manifest/prepared_run.json"
            ),
            "mode": 0o644,
            "size_bytes": 78976,
            "sha256": (
                "ae1bc8256b3713dcf70ec99dde05a60d9d018650766e80ba29d49214c5e9bd19"
            ),
        },
        "inspection": {
            "path": (
                "runs/mks24-stage-i/E03-forcing-policy/R12/"
                "s00_rankio_t0_t0p25/manifest/segment_inspection.json"
            ),
            "mode": 0o644,
            "size_bytes": 66333,
            "sha256": (
                "f901e51978a2e07d90728a10e8ca0f2c5abc393020a0c3ce8a4064ee5c914ab0"
            ),
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
            "sha256": (
                "19f2fd19c003f644c9a53f0ba7bed5a754de46dd905d374eaf6cc65661e6d9b5"
            ),
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
            "sha256": (
                "06507bd49e5b10f9db2c75f8858d70187b5ea7b0ce9a6e139db35aecd45d88aa"
            ),
        },
        "independent_validation": {
            "path": "accounting/4766856.stage_i.independent_validation.json",
            "mode": 0o444,
            "size_bytes": 278070,
            "sha256": (
                "1d69d1b8f6cb35928f2334336399f3f7eee07aa511bf9ef36490670fe4da27e5"
            ),
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
            "sha256": (
                "45c20708e5619bcd46b3ec16f82f517c8acee415c1d2002ec0d43386c89d94f9"
            ),
        },
    },
}
CALLER_CHILD_ENVIRONMENT_EXACT_DENY = frozenset({
    "BASH_ENV", "BASHOPTS", "CDPATH", "ENV", "GLOBIGNORE", "IFS",
    "LIBPATH", "PERL5LIB", "PERLLIB", "PROMPT_COMMAND", "PS4",
    "PYTHONHOME", "PYTHONPATH", "RUBYLIB", "RUBYOPT", "SHELLOPTS",
    "SHLIB_PATH", "ZDOTDIR",
})
CALLER_CHILD_ENVIRONMENT_PREFIX_DENY = (
    "BASH_FUNC_", "DYLD_", "GIT_", "LD_", "PYTHON", "SBATCH_", "SLURM_",
)
SELF_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_CONTROLLER_DESCRIPTOR"
PYTHON_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_CONTROLLER_PYTHON_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_STAGE_I_CONTROLLER_SOURCE"
REPOSITORY_ROOT_ENV = "_CGL_LF_STAGE_I_CONTROLLER_REPOSITORY_ROOT"
RENAME_NOREPLACE = 1
RENAME_EXCHANGE = 2
RENAMEAT2_UNSUPPORTED_ERRNOS = frozenset({
    errno.EINVAL,
    errno.ENOSYS,
    getattr(errno, "ENOTSUP", errno.EINVAL),
    getattr(errno, "EOPNOTSUPP", errno.EINVAL),
})
R17_PROFILE_KEYS = frozenset({
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
})
R17_SOLE_PROFILE_KEYS = (
    "athena_walltime",
    "case_id",
    "nodes",
    "parent_job_id",
    "parent_result",
    "parent_segment",
    "restart_file",
    "restart_time",
    "segment",
    "source_bundle",
    "source_bundle_sha256",
    "time_tlim_target",
    "walltime",
)
CONCURRENT_CASE_IDS = frozenset(f"R{number:02d}" for number in range(3, 17))
MAX_ACTIVE_STAGE_I_SEGMENTS = 4
MAX_ACTIVE_STAGE_I_NODES = 10
CASE_NODE_PROFILES = {
    "R02": frozenset({1}),
    "R03": frozenset({1}),
    **{
        f"R{number:02d}": frozenset({1, 2, 4})
        for number in range(4, 16)
    },
    "R16": frozenset({1, 2}),
    R17_CASE_ID: frozenset({8}),
}
REQUIRED_CASE_FINAL_TIME = 10.0
COMPLETED_R16_NODE_HOURS = 6.145556
COMPLETED_R02_STANDARD_LAYOUT_PILOT_NODE_HOURS = 0.473333
COMPLETED_R17_HIGH_RESOLUTION_PILOT_NODE_HOURS = 4.235556
MEASURED_STAGE_I_RESERVED_NODE_HOURS = 900.0
PROMOTED_STAGE_I_RESERVED_NODE_HOURS = 1400.0
CURRENT_STAGE_I_RESERVED_NODE_HOURS = PROMOTED_STAGE_I_RESERVED_NODE_HOURS
MAX_SEGMENT_SECONDS = 2 * 60 * 60
# ParameterInput::LoadFromFile accepts a terminator found in its eleventh 4 KiB
# chunk and then seeks one byte past the marker for the trailing newline.
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
QUALIFIED_RESTART_BINARY_ABIS = {
    (
        "9e07542281e4e6d125582f253df3ad2e3b8b154d",
        "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c",
    ): {
        "mesh_time_offset_after_parameter_dump": 232,
        "mesh_time_format": "<d",
        "allowed_marker_modes": frozenset({
            "full_precision", "legacy_default_precision",
        }),
    },
}
HISTORICAL_SUBMITTED_R03_UTILITY_TRANSITION = {
    "project_root": str(DEFAULT_ROOT),
    "state": "submitted",
    "job_id": "4762472",
    "case_id": "R03",
    "segment": "s00_rankio_t0_t0p5",
    "production_utility_revision": "dbe7e50045bfe0c3a99ce0ba0b21e8735bbf1103",
    "production_utility_sha256": "54ec671bb45aa27735a174d40b4b2e6009070716346ea09699bbe62421bbfada",
    "source_bundle_sha256": "b8437f066f8391a696efaaaf0de531430a9dac27c95dfc38aa7328dd49fb19fe",
    "executable_revision": "9e07542281e4e6d125582f253df3ad2e3b8b154d",
    "executable_sha256": "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c",
}
CANONICAL_SOURCE_BUNDLE_RECOVERY = {
    "job_id": "4766485",
    "manifest_relative": (
        "runs/mks24-stage-i/E03-forcing-policy/R03/"
        "s01_rankio_t0p312823_t0p5/manifest/prepared_run.json"
    ),
    "source_bundle_relative": (
        "source-archives/athenak-feature-cgl-through-c7e4fa30e.bundle"
    ),
    "source_bundle_expected_sha256": (
        "d94d559108470157f07981c9da8fec343a992128c0c4333df87181b81f0c505e"
    ),
    "source_bundle_observed_sha256": (
        "10fbc6a736efc36054c594b2ef1f0f33700586366c7de483d00d37c1ca28a0c8"
    ),
    "incident_relative": (
        "accounting/mks24_stage_i_E03_forcing_policy_source_bundle_incidents/"
        "20260605T035710Z-c7e4fa30e-corruption/incident.json"
    ),
    "incident_sha256": (
        "7e017a4eed5d3a056f7a8dc0b24c131f5a2a1ab89093878d3d83a2473a4992fd"
    ),
    "authorization_sha256": (
        "b9218128173b11d245eeced1e7115a88a5a865508d9a6bc1e024177b9dbc4fa1"
    ),
    "publication_audit_sha256": (
        "7f3f9c402a90cd1381061345a0b781c6e8e7a4363d67bb4c6623ac77e5818400"
    ),
}
# Allow scheduler timestamp formatting and host-clock skew around the persisted
# pre-sbatch ambiguity barrier, but never an unrelated later submission.
SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS = 5 * 60
EXPECTED_RANKS_PER_NODE = 8
EXPECTED_CPUS_PER_TASK = 7
SEGMENT_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,28}")
JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
GIT_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
R17_RECOST_SEGMENT_PATTERN = re.compile(
    r"s(?P<index>[0-9]+)_rankio_t(?P<start>[0-9]+(?:p[0-9]+)?)"
    r"_t(?P<target>[0-9]+(?:p[0-9]+)?)"
)
R17_SCOPED_NODE_HOUR_METHOD = "observed-stage-i-scoped-node-hour-rate-v2"
BATCH_SCRIPT_DIGEST_PLACEHOLDER = "0" * 64
BATCH_SCRIPT_DIGEST_PATTERN = re.compile(
    r"(?m)^BATCH_SCRIPT_SHA256=([0-9a-f]{64})$"
)
LEDGER_NODE_HOUR_TOLERANCE = 5.0e-7 + 1.0e-12
LEDGER_CUMULATIVE_NODE_HOUR_TOLERANCE = 5.0e-7 + 1.0e-12
LEDGER_COLUMNS = (
    "execution_epoch",
    "job_id",
    "submitted_utc",
    "completed_utc",
    "case_id",
    "case_name",
    "segment",
    "state",
    "exit_code",
    "nodes",
    "requested_walltime",
    "elapsed_seconds",
    "reserved_node_hours",
    "actual_node_hours",
    "cumulative_stage_i_node_hours",
    "executable_revision",
    "executable_sha256",
    "input_revision",
    "input_file",
    "output_dir",
    "result",
    "notes",
)
RESERVATION_REQUIRED_COLUMNS = frozenset({
    "execution_epoch",
    "manifest",
    "case_id",
    "case_name",
    "segment",
    "nodes",
    "requested_walltime",
    "reserved_node_hours",
    "state",
    "prepared_utc",
})
RESERVATION_OPTIONAL_COLUMNS = frozenset({
    "execution_intent_sha256",
    "job_id",
    "actual_node_hours",
    "result",
    "notes",
})
TRANSACTION_KINDS = frozenset({
    "prepared", "submit_pending", "submitted", "submit_cleared",
    "recorded", "cancelled", "cancelled_submitted",
})
TRANSACTION_COMMON_COLUMNS = frozenset({
    "schema_version",
    "execution_epoch",
    "transaction_id",
    "kind",
    "created_utc",
    "manifest_path",
    "prior_reservations",
    "prior_reservations_sha256",
})
TRANSACTION_PAYLOAD_COLUMNS = frozenset({
    "manifest",
    "reservations",
    "ledger_row",
})
NONTERMINAL_STATES = {
    "PENDING",
    "RUNNING",
    "CONFIGURING",
    "COMPLETING",
    "SUSPENDED",
}
SHARED_ROOT_ACTIVE_STATES = {
    "prepared",
    "submitted",
    "pending",
    "running",
}
STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
_ACTIVE_ROOT_LOCKS: dict[
    Path, tuple[object, int, int, os.stat_result, tuple[tuple[str, object], ...]]
] = {}


def utc_now() -> str:
    """Return a deterministic UTC timestamp."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def fsync_directory(path: Path) -> None:
    """Persist directory-entry updates after atomic replacement or removal."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_file(path: Path) -> None:
    """Persist one retained file after a copied or appended payload."""

    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def directory_security_identity(profile: os.stat_result) -> tuple[int, int, int]:
    """Return stable directory mode and ownership metadata."""

    return stat.S_IMODE(profile.st_mode), profile.st_uid, profile.st_gid


def directory_content_identity(profile: os.stat_result) -> tuple[object, ...]:
    """Return directory metadata that changes with child namespace content."""

    return tuple(
        getattr(profile, field)
        for field in (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
            "st_size", "st_mtime_ns", "st_ctime_ns",
        )
    )


def require_directory_descriptor_binding(path: Path, descriptor: int,
                                         label: str,
                                         expected: os.stat_result | None = None,
                                         ) -> os.stat_result:
    """Require an open directory to remain bound to its exact pathname."""

    require_owner_symlink_free_path(path, label)
    opened = os.fstat(descriptor)
    try:
        named = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} pathname changed: {path}") from error
    if (
        not stat.S_ISDIR(opened.st_mode)
        or not stat.S_ISDIR(named.st_mode)
        or (named.st_dev, named.st_ino) != (opened.st_dev, opened.st_ino)
    ):
        raise ValueError(f"{label} pathname changed: {path}")
    if directory_security_identity(named) != directory_security_identity(opened):
        raise ValueError(f"{label} pathname security metadata changed: {path}")
    if expected is not None and (
        (opened.st_dev, opened.st_ino) != (expected.st_dev, expected.st_ino)
        or directory_security_identity(opened)
        != directory_security_identity(expected)
    ):
        raise ValueError(f"{label} security metadata changed: {path}")
    return opened


def require_bound_directory_entry(parent_fd: int, name: str, descriptor: int,
                                  expected: os.stat_result, label: str) -> None:
    """Require one direct child to retain an exact directory inode and profile."""

    try:
        named = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        opened = os.fstat(descriptor)
    except OSError as error:
        raise ValueError(f"{label} directory entry changed") from error
    if (
        not stat.S_ISDIR(named.st_mode)
        or not stat.S_ISDIR(opened.st_mode)
        or (named.st_dev, named.st_ino)
        != (expected.st_dev, expected.st_ino)
        or (opened.st_dev, opened.st_ino)
        != (expected.st_dev, expected.st_ino)
        or directory_security_identity(named)
        != directory_security_identity(expected)
        or directory_security_identity(opened)
        != directory_security_identity(expected)
    ):
        raise ValueError(f"{label} directory entry changed")


def require_durable_directory_profile(profile: os.stat_result,
                                      parent_profile: os.stat_result,
                                      path: Path) -> None:
    """Require the exact intended public profile for a durable directory."""

    parent_mode = stat.S_IMODE(parent_profile.st_mode)
    expected_mode = DURABLE_DIRECTORY_MODE
    expected_gid = os.getegid()
    if parent_mode & stat.S_ISGID:
        expected_mode |= stat.S_ISGID
        expected_gid = parent_profile.st_gid
    if (
        not stat.S_ISDIR(profile.st_mode)
        or profile.st_uid != os.geteuid()
        or profile.st_gid != expected_gid
        or stat.S_IMODE(profile.st_mode) != expected_mode
    ):
        raise ValueError(
            f"directory creation outcome has unexpected owner/group/mode; "
            f"expected {os.geteuid()}:{expected_gid} {expected_mode:04o}: {path}"
        )


def stable_bound_directory_entries(path: Path, descriptor: int,
                                   expected: os.stat_result,
                                   label: str) -> tuple[list[str], os.stat_result]:
    """Return one stable descriptor-bound directory listing or reject churn."""

    for _attempt in range(3):
        require_directory_descriptor_binding(path, descriptor, label, expected)
        before = os.fstat(descriptor)
        entries = sorted(os.listdir(descriptor))
        after = os.fstat(descriptor)
        require_directory_descriptor_binding(path, descriptor, label, expected)
        if directory_content_identity(before) == directory_content_identity(after):
            return entries, after
    raise ValueError(f"{label} contents changed repeatedly while listing")


def regular_file_stable_identity(profile: os.stat_result) -> tuple[object, ...]:
    """Return exact inode metadata that must remain stable across rename."""

    return tuple(
        getattr(profile, field)
        for field in (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
            "st_size", "st_mtime_ns", "st_ctime_ns",
        )
    )


def regular_file_rename_identity(profile: os.stat_result) -> tuple[object, ...]:
    """Return inode metadata that remains stable across an authenticated rename."""

    return tuple(
        getattr(profile, field)
        for field in (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
            "st_size", "st_mtime_ns",
        )
    )


def regular_file_hardlink_identity(profile: os.stat_result) -> tuple[object, ...]:
    """Return inode metadata stable while a retained hard-link alias exists."""

    return tuple(
        getattr(profile, field)
        for field in (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_size",
            "st_mtime_ns",
        )
    )


def descriptor_bytes(descriptor: int) -> bytes:
    """Read exact bytes without changing one retained descriptor offset."""

    blocks = []
    offset = 0
    while True:
        block = os.pread(descriptor, 1024 * 1024, offset)
        if not block:
            return b"".join(blocks)
        blocks.append(block)
        offset += len(block)


def descriptor_sha256(descriptor: int) -> str:
    """Hash exact descriptor bytes without changing its retained offset."""

    digest = hashlib.sha256()
    offset = 0
    while True:
        block = os.pread(descriptor, 1024 * 1024, offset)
        if not block:
            return digest.hexdigest()
        digest.update(block)
        offset += len(block)


class RegularFileBinding:
    """Hold one exact inode and, when readable, its exact retained bytes."""

    def __init__(self, directory_fd: int, name: str, label: str, *,
                 allow_unreadable: bool = False,
                 retain_payload: bool = True,
                 expected_payload: bytes | None = None,
                 expected_sha256: str | None = None):
        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as error:
            raise ValueError(f"{label} target is unavailable") from error
        if not stat.S_ISREG(named.st_mode):
            raise ValueError(f"{label} target is not a regular file")
        readable = True
        try:
            descriptor = os.open(
                name,
                os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
        except OSError as error:
            if not allow_unreadable or error.errno not in {errno.EACCES, errno.EPERM}:
                raise ValueError(f"{label} target changed") from error
            readable = False
            try:
                descriptor = os.open(
                    name,
                    getattr(os, "O_PATH", os.O_RDONLY)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=directory_fd,
                )
            except OSError as path_error:
                raise ValueError(f"{label} target changed") from path_error
        try:
            opened = os.fstat(descriptor)
            if (
                not stat.S_ISREG(opened.st_mode)
                or regular_file_stable_identity(opened)
                != regular_file_stable_identity(named)
            ):
                raise ValueError(f"{label} target changed")
            payload = (
                descriptor_bytes(descriptor)
                if readable and retain_payload
                else None
            )
            digest = (
                hashlib.sha256(payload).hexdigest()
                if payload is not None
                else descriptor_sha256(descriptor) if readable else None
            )
            if expected_payload is not None and payload != expected_payload:
                raise ValueError(f"{label} exact bytes differ")
            if expected_sha256 is not None and digest != expected_sha256:
                raise ValueError(f"{label} exact bytes differ")
        except BaseException:
            os.close(descriptor)
            raise
        self.descriptor = descriptor
        self.profile = opened
        self.payload = payload
        self.digest = digest
        self.label = label

    def inode_is_bound(self, directory_fd: int, name: str) -> bool:
        """Return whether a direct child still names this exact retained inode."""

        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            opened = os.fstat(self.descriptor)
        except OSError:
            return False
        return (
            stat.S_ISREG(named.st_mode)
            and stat.S_ISREG(opened.st_mode)
            and (named.st_dev, named.st_ino)
            == (self.profile.st_dev, self.profile.st_ino)
            == (opened.st_dev, opened.st_ino)
        )

    def assert_inode(self, label: str | None = None) -> None:
        """Require the retained descriptor profile and exact bytes to be unchanged."""

        retained_label = label or self.label
        try:
            opened = os.fstat(self.descriptor)
        except OSError as error:
            raise ValueError(f"{retained_label} descriptor changed") from error
        if (
            regular_file_stable_identity(opened)
            != regular_file_stable_identity(self.profile)
            or (
                self.payload is not None
                and descriptor_bytes(self.descriptor) != self.payload
            )
            or (
                self.digest is not None
                and descriptor_sha256(self.descriptor) != self.digest
            )
        ):
            raise ValueError(f"{retained_label} exact bytes or inode metadata changed")

    def assert_bound(self, directory_fd: int, name: str,
                     label: str | None = None) -> None:
        """Require one direct child and retained descriptor to remain exact."""

        retained_label = label or self.label
        self.assert_inode(retained_label)
        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        except OSError as error:
            raise ValueError(f"{retained_label} target changed") from error
        if regular_file_stable_identity(named) != regular_file_stable_identity(
            self.profile
        ):
            raise ValueError(f"{retained_label} target changed")

    def accept_bound_rename(self, directory_fd: int, name: str,
                            label: str | None = None) -> None:
        """Accept only rename-induced ctime drift while retaining exact bytes."""

        retained_label = label or self.label
        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            opened = os.fstat(self.descriptor)
        except OSError as error:
            raise ValueError(f"{retained_label} target changed") from error
        if (
            not stat.S_ISREG(named.st_mode)
            or not stat.S_ISREG(opened.st_mode)
            or regular_file_rename_identity(named)
            != regular_file_rename_identity(self.profile)
            or regular_file_rename_identity(opened)
            != regular_file_rename_identity(self.profile)
            or (
                self.payload is not None
                and descriptor_bytes(self.descriptor) != self.payload
            )
            or (
                self.digest is not None
                and descriptor_sha256(self.descriptor) != self.digest
            )
        ):
            raise ValueError(f"{retained_label} exact bytes or inode metadata changed")
        self.profile = opened
        self.assert_bound(directory_fd, name, retained_label)

    def assert_bound_hardlink(self, directory_fd: int, name: str,
                              expected_links: int,
                              label: str | None = None) -> None:
        """Require one exact retained inode while a hard-link move is staged."""

        retained_label = label or self.label
        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            opened = os.fstat(self.descriptor)
        except OSError as error:
            raise ValueError(f"{retained_label} target changed") from error
        if (
            not stat.S_ISREG(named.st_mode)
            or not stat.S_ISREG(opened.st_mode)
            or named.st_nlink != expected_links
            or opened.st_nlink != expected_links
            or regular_file_hardlink_identity(named)
            != regular_file_hardlink_identity(self.profile)
            or regular_file_hardlink_identity(opened)
            != regular_file_hardlink_identity(self.profile)
            or (
                self.payload is not None
                and descriptor_bytes(self.descriptor) != self.payload
            )
            or (
                self.digest is not None
                and descriptor_sha256(self.descriptor) != self.digest
            )
        ):
            raise ValueError(f"{retained_label} exact bytes or inode metadata changed")

    def accept_bound_hardlink_unlink(self, directory_fd: int, name: str,
                                     label: str | None = None) -> None:
        """Accept only ctime/link-count drift after exact alias removal."""

        retained_label = label or self.label
        try:
            named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            opened = os.fstat(self.descriptor)
        except OSError as error:
            raise ValueError(f"{retained_label} target changed") from error
        if (
            not stat.S_ISREG(named.st_mode)
            or not stat.S_ISREG(opened.st_mode)
            or named.st_nlink != 1
            or opened.st_nlink != 1
            or regular_file_hardlink_identity(named)
            != regular_file_hardlink_identity(self.profile)
            or regular_file_hardlink_identity(opened)
            != regular_file_hardlink_identity(self.profile)
            or (
                self.payload is not None
                and descriptor_bytes(self.descriptor) != self.payload
            )
            or (
                self.digest is not None
                and descriptor_sha256(self.descriptor) != self.digest
            )
        ):
            raise ValueError(f"{retained_label} exact bytes or inode metadata changed")
        self.profile = opened
        self.assert_bound(directory_fd, name, retained_label)

    def close(self) -> None:
        """Close the retained inode descriptor."""

        os.close(self.descriptor)


def open_regular_file_binding(directory_fd: int, name: str, label: str, *,
                              allow_absent: bool = False,
                              allow_unreadable: bool = False,
                              retain_payload: bool = True,
                              expected_payload: bytes | None = None,
                              expected_sha256: str | None = None,
                              ) -> RegularFileBinding | None:
    """Open one exact direct-child binding, optionally allowing absence."""

    try:
        os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    except FileNotFoundError:
        if allow_absent:
            return None
        raise ValueError(f"{label} target is unavailable") from None
    return RegularFileBinding(
        directory_fd,
        name,
        label,
        allow_unreadable=allow_unreadable,
        retain_payload=retain_payload,
        expected_payload=expected_payload,
        expected_sha256=expected_sha256,
    )


def require_requested_regular_file_profile(binding: RegularFileBinding,
                                           mode: int, label: str) -> None:
    """Require an exact owner-controlled requested publication profile."""

    binding_error = None
    try:
        binding.assert_inode(label)
    except ValueError as error:
        binding_error = error
    profile = os.fstat(binding.descriptor)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or stat.S_IMODE(profile.st_mode) != stat.S_IMODE(mode)
    ):
        raise ValueError(
            f"{label} must be an owner-controlled regular "
            f"{stat.S_IMODE(mode):04o} single-link file"
        ) from binding_error
    if binding_error is not None:
        raise binding_error


def durably_authenticate_bound_file(directory_fd: int, name: str,
                                    binding: RegularFileBinding,
                                    label: str) -> None:
    """Persist and reauthenticate one exact bound regular file."""

    binding.assert_bound(directory_fd, name, label)
    os.fsync(binding.descriptor)
    os.fsync(directory_fd)
    binding.assert_bound(directory_fd, name, label)


def copy_bound_file_payload(source: RegularFileBinding,
                            destination_descriptor: int) -> None:
    """Copy exact bytes from a retained descriptor without pathname reopening."""

    offset = 0
    while offset < source.profile.st_size:
        block = os.pread(
            source.descriptor,
            min(1024 * 1024, source.profile.st_size - offset),
            offset,
        )
        if not block:
            raise ValueError("copy source ended before its retained size")
        written = 0
        while written < len(block):
            count = os.write(destination_descriptor, block[written:])
            if count <= 0:
                raise ValueError("copy destination stopped accepting bytes")
            written += count
        offset += len(block)
    if os.pread(source.descriptor, 1, offset):
        raise ValueError("copy source grew beyond its retained size")


def regular_file_profile(directory_fd: int, name: str,
                         label: str) -> os.stat_result | None:
    """Return one descriptor-authenticated regular-file profile."""

    try:
        named = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    except FileNotFoundError:
        return None
    if not stat.S_ISREG(named.st_mode):
        raise ValueError(f"{label} target is not a regular file")
    try:
        descriptor = os.open(
            name,
            getattr(os, "O_PATH", os.O_RDONLY) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
    except OSError as error:
        raise ValueError(f"{label} target changed") from error
    try:
        opened = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (
        not stat.S_ISREG(opened.st_mode)
        or regular_file_stable_identity(opened)
        != regular_file_stable_identity(named)
    ):
        raise ValueError(f"{label} target changed")
    return opened


def profile_identity(profile: os.stat_result) -> tuple[int, int]:
    """Return the immutable filesystem identity of one retained inode."""

    return profile.st_dev, profile.st_ino


def require_bound_regular_entry(directory_fd: int, name: str,
                                expected: os.stat_result, label: str) -> None:
    """Require one direct-child name to retain an authenticated inode."""

    retained = regular_file_profile(directory_fd, name, label)
    if (
        retained is None
        or regular_file_stable_identity(retained)
        != regular_file_stable_identity(expected)
    ):
        raise ValueError(f"{label} target changed")


def require_regular_entry_absent(directory_fd: int, name: str, label: str) -> None:
    """Require one direct-child name to remain absent."""

    try:
        os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
    except FileNotFoundError:
        return
    raise ValueError(f"{label} target reappeared")


class MutationAuthorityError(ValueError):
    """Report that the held lock or an original mutation parent is no longer valid."""


def is_canonical_mutation_path(path: Path) -> bool:
    """Return whether one absolute target is beneath the canonical Stage I root."""

    canonical = DEFAULT_ROOT.expanduser().absolute()
    absolute = path.expanduser().absolute()
    try:
        absolute.relative_to(canonical)
    except ValueError:
        return False
    return True


def require_metadata_mutation_boundary(path: Path, parent: Path, parent_fd: int,
                                       parent_profile: os.stat_result,
                                       label: str) -> None:
    """Reauthenticate the lock and original parent immediately around mutation."""

    try:
        require_canonical_mutation_lock(path)
        retained = require_directory_descriptor_binding(
            parent, parent_fd, f"{label} parent", parent_profile
        )
        if is_canonical_mutation_path(path) and (
            retained.st_uid != os.geteuid()
            or stat.S_IMODE(retained.st_mode) & 0o022
        ):
            raise ValueError(
                f"{label} canonical parent is not owner-controlled or is "
                f"group/world-writable: {parent}"
            )
    except MutationAuthorityError:
        raise
    except ValueError as error:
        raise MutationAuthorityError(str(error)) from error


def fsync_and_reauthenticate_entries(directory_fd: int, names: tuple[str, ...],
                                     label: str) -> None:
    """Durably authenticate the current descriptor-relative namespace state."""

    before_error = None
    directory_before = None
    entries_before = None
    try:
        directory_before = os.fstat(directory_fd)
        entries_before = tuple(
            regular_file_profile(directory_fd, name, f"{label} entry {name}")
            for name in names
        )
    except BaseException as error:
        before_error = error
    os.fsync(directory_fd)
    try:
        directory_after = os.fstat(directory_fd)
        entries_after = tuple(
            regular_file_profile(directory_fd, name, f"{label} entry {name}")
            for name in names
        )
    except BaseException as error:
        raise ValueError(f"{label} state changed after directory fsync") from error
    if before_error is not None:
        raise ValueError(f"{label} state changed before directory fsync") from before_error
    if directory_before is None or entries_before is None:
        raise ValueError(f"{label} pre-fsync state is unavailable")
    if (
        not stat.S_ISDIR(directory_before.st_mode)
        or not stat.S_ISDIR(directory_after.st_mode)
        or (directory_after.st_dev, directory_after.st_ino)
        != (directory_before.st_dev, directory_before.st_ino)
        or directory_security_identity(directory_after)
        != directory_security_identity(directory_before)
    ):
        raise ValueError(f"{label} directory descriptor changed")
    for name, before, after in zip(names, entries_before, entries_after):
        if (before is None) != (after is None) or (
            before is not None
            and after is not None
            and regular_file_stable_identity(before)
            != regular_file_stable_identity(after)
        ):
            raise ValueError(f"{label} entry {name} changed during directory fsync")


class Renameat2Unsupported(OSError):
    """Report an explicit filesystem/kernel rejection of nonzero rename flags."""


class LustrePosixRenameUnsupported(ValueError):
    """Report that the reviewed descriptor-relative POSIX fallback is unavailable."""


def require_direct_child_name(name: str, label: str) -> None:
    """Require one normalized descriptor-relative direct-child name."""

    if not name or name in {".", ".."} or "/" in name:
        raise ValueError(f"{label} contains an invalid direct-child name")


def fsync_directory_descriptors(*descriptors: int) -> None:
    """Persist each distinct retained directory descriptor exactly once."""

    retained = set()
    for descriptor in descriptors:
        profile = os.fstat(descriptor)
        identity = profile.st_dev, profile.st_ino
        if identity not in retained:
            os.fsync(descriptor)
            retained.add(identity)


def renameat2_between(source_directory_fd: int, source: str,
                      target_directory_fd: int, target: str,
                      flags: int, label: str) -> None:
    """Perform one raw descriptor-relative Linux renameat2 operation."""

    require_direct_child_name(source, label)
    require_direct_child_name(target, label)
    if flags not in {RENAME_NOREPLACE, RENAME_EXCHANGE}:
        raise ValueError(f"{label} requested unsupported renameat2 flags: {flags}")

    try:
        operation = ctypes.CDLL(None, use_errno=True).renameat2
    except AttributeError as error:
        raise Renameat2Unsupported(
            errno.ENOSYS,
            "descriptor-relative renameat2 is unsupported",
            f"{source} -> {target}",
        ) from error
    operation.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    operation.restype = ctypes.c_int
    if operation(
        source_directory_fd,
        os.fsencode(source),
        target_directory_fd,
        os.fsencode(target),
        flags,
    ) != 0:
        retained_errno = ctypes.get_errno()
        if retained_errno in RENAMEAT2_UNSUPPORTED_ERRNOS:
            raise Renameat2Unsupported(
                retained_errno,
                os.strerror(retained_errno),
                f"{source} -> {target}",
            )
        raise OSError(
            retained_errno, os.strerror(retained_errno), f"{source} -> {target}"
        )


def renameat2(directory_fd: int, source: str, target: str,
              flags: int, label: str) -> None:
    """Perform one same-directory Linux renameat2 operation."""

    renameat2_between(directory_fd, source, directory_fd, target, flags, label)


def link_bound_noreplace_between(source_directory_fd: int, source: str,
                                 target_directory_fd: int, target: str,
                                 expected: RegularFileBinding, label: str,
                                 mutation_guard) -> None:
    """Lustre-compatible no-replace move using a durable hard-link commit."""

    mutation_guard()
    source_is_expected = expected.inode_is_bound(source_directory_fd, source)
    target_is_expected = expected.inode_is_bound(target_directory_fd, target)
    if not source_is_expected:
        raise ValueError(f"{label} Lustre hard-link publication source changed")
    if expected.profile.st_nlink not in {1, 2}:
        raise ValueError(
            f"{label} Lustre hard-link publication requires one recoverable source"
        )
    link_error = None
    if target_is_expected:
        expected.assert_bound_hardlink(
            source_directory_fd, source, 2, f"{label} retained linked source"
        )
        expected.assert_bound_hardlink(
            target_directory_fd, target, 2, f"{label} retained linked target"
        )
    else:
        if expected.profile.st_nlink != 1:
            raise ValueError(
                f"{label} Lustre hard-link publication lost its recovery alias"
            )
        require_regular_entry_absent(target_directory_fd, target, label)
        try:
            os.link(
                source,
                target,
                src_dir_fd=source_directory_fd,
                dst_dir_fd=target_directory_fd,
                follow_symlinks=False,
            )
        except BaseException as error:
            link_error = error
        source_is_expected = expected.inode_is_bound(source_directory_fd, source)
        target_is_expected = expected.inode_is_bound(target_directory_fd, target)
    if not target_is_expected:
        if source_is_expected:
            try:
                expected.assert_bound(source_directory_fd, source, label)
            except BaseException as state_error:
                raise ValueError(
                    f"{label} Lustre hard-link publication has an unsupported "
                    f"namespace state"
                ) from state_error
            try:
                require_regular_entry_absent(target_directory_fd, target, label)
            except ValueError as collision:
                raise FileExistsError(
                    errno.EEXIST, os.strerror(errno.EEXIST), target
                ) from collision
            if link_error is not None:
                raise ValueError(
                    f"{label} Lustre hard-link publication failed before mutation"
                ) from link_error
        raise ValueError(
            f"{label} Lustre hard-link publication state is ambiguous"
        ) from link_error

    if source_is_expected:
        expected.assert_bound_hardlink(
            source_directory_fd, source, 2, f"{label} linked source"
        )
        expected.assert_bound_hardlink(
            target_directory_fd, target, 2, f"{label} linked target"
        )
        fsync_directory_descriptors(source_directory_fd, target_directory_fd)
        expected.assert_bound_hardlink(
            source_directory_fd, source, 2, f"{label} linked source"
        )
        expected.assert_bound_hardlink(
            target_directory_fd, target, 2, f"{label} linked target"
        )
        unlink_error = None
        try:
            os.unlink(source, dir_fd=source_directory_fd)
        except BaseException as error:
            unlink_error = error
        fsync_directory_descriptors(source_directory_fd, target_directory_fd)
        source_is_expected = expected.inode_is_bound(source_directory_fd, source)
        target_is_expected = expected.inode_is_bound(target_directory_fd, target)
        if source_is_expected and target_is_expected:
            expected.assert_bound_hardlink(
                source_directory_fd, source, 2, f"{label} retained linked source"
            )
            expected.assert_bound_hardlink(
                target_directory_fd, target, 2, f"{label} retained linked target"
            )
            expected.profile = os.fstat(expected.descriptor)
            raise ValueError(
                f"{label} Lustre two-link recovery state was durably preserved"
            ) from unlink_error
        if source_is_expected or not target_is_expected:
            raise ValueError(
                f"{label} Lustre hard-link publication has an unsupported "
                f"post-unlink namespace state"
            ) from unlink_error

    require_regular_entry_absent(
        source_directory_fd, source, f"{label} retired linked source"
    )
    expected.accept_bound_hardlink_unlink(target_directory_fd, target, label)
    fsync_directory_descriptors(source_directory_fd, target_directory_fd)
    mutation_guard()
    require_regular_entry_absent(
        source_directory_fd, source, f"{label} retired linked source"
    )
    expected.assert_bound(target_directory_fd, target, label)


def unlink_bound_hardlink_alias(directory_fd: int, alias: str, retained: str,
                                expected: RegularFileBinding, label: str,
                                mutation_guard) -> None:
    """Remove one exact recovery alias, accepting a classified post-unlink error."""

    mutation_guard()
    expected.assert_bound_hardlink(directory_fd, alias, 2, f"{label} alias")
    expected.assert_bound_hardlink(directory_fd, retained, 2, label)
    unlink_error = None
    try:
        os.unlink(alias, dir_fd=directory_fd)
    except BaseException as error:
        unlink_error = error
    os.fsync(directory_fd)
    mutation_guard()
    alias_is_expected = expected.inode_is_bound(directory_fd, alias)
    retained_is_expected = expected.inode_is_bound(directory_fd, retained)
    if alias_is_expected and retained_is_expected:
        expected.assert_bound_hardlink(directory_fd, alias, 2, f"{label} alias")
        expected.assert_bound_hardlink(directory_fd, retained, 2, label)
        raise ValueError(
            f"{label} recovery alias remains in a durable two-link state"
        ) from unlink_error
    if alias_is_expected or not retained_is_expected:
        raise ValueError(
            f"{label} recovery alias has an unsupported post-unlink state"
        ) from unlink_error
    require_regular_entry_absent(directory_fd, alias, f"{label} removed alias")
    expected.accept_bound_hardlink_unlink(directory_fd, retained, label)


def rename_bound_noreplace(directory_fd: int, source: str, target: str,
                           expected: RegularFileBinding, label: str,
                           mutation_guard) -> None:
    """Publish exact bytes and preserve any ambiguous forward state durably."""

    mutation_guard()
    expected.assert_bound(directory_fd, source, label)
    require_regular_entry_absent(directory_fd, target, label)
    operation_error = None
    try:
        try:
            renameat2(directory_fd, source, target, RENAME_NOREPLACE, label)
        except Renameat2Unsupported:
            link_bound_noreplace_between(
                directory_fd,
                source,
                directory_fd,
                target,
                expected,
                label,
                mutation_guard,
            )
    except BaseException as error:
        try:
            fsync_and_reauthenticate_entries(
                directory_fd, (source, target), f"{label} ambiguous publication"
            )
        except BaseException as durable_error:
            raise ValueError(
                f"{label} ambiguous rename state is not durably authenticated"
            ) from durable_error
        if expected.inode_is_bound(directory_fd, source):
            expected.assert_bound(directory_fd, source, label)
            if isinstance(error, FileExistsError):
                raise ValueError(f"{label} target already exists") from error
            raise ValueError(f"{label} rename failed before publication") from error
        if (
            expected.inode_is_bound(directory_fd, target)
            and regular_file_profile(
                directory_fd, source, f"{label} ambiguous source"
            ) is None
        ):
            try:
                expected.accept_bound_rename(directory_fd, target, label)
            except BaseException as binding_error:
                raise ValueError(
                    f"{label} ambiguous rename state was durably preserved"
                ) from binding_error
            operation_error = error
        else:
            raise ValueError(
                f"{label} ambiguous rename state was durably preserved"
            ) from error
    if operation_error is None:
        try:
            expected.accept_bound_rename(directory_fd, target, label)
            require_regular_entry_absent(directory_fd, source, label)
            os.fsync(directory_fd)
            mutation_guard()
            expected.assert_bound(directory_fd, target, label)
            require_regular_entry_absent(directory_fd, source, label)
        except BaseException as error:
            operation_error = error
    if operation_error is not None:
        try:
            fsync_and_reauthenticate_entries(
                directory_fd, (source, target), f"{label} failed publication"
            )
        except BaseException as durable_error:
            raise ValueError(
                f"{label} post-rename state is not durably authenticated"
            ) from durable_error
        error = operation_error
        if isinstance(error, MutationAuthorityError):
            raise error
        try:
            mutation_guard()
        except MutationAuthorityError as authority_error:
            raise authority_error from error
        raise ValueError(
            f"{error}; {label} durable forward recovery state was preserved"
        ) from error


def exchange_bound_entries(directory_fd: int, source: str, target: str,
                           source_expected: RegularFileBinding,
                           target_expected: RegularFileBinding, label: str,
                           mutation_guard) -> None:
    """Exchange exact bytes and preserve any ambiguous forward state durably."""

    mutation_guard()
    source_expected.assert_bound(directory_fd, source, label)
    target_expected.assert_bound(directory_fd, target, f"{label} predecessor")
    operation_error = None
    try:
        try:
            renameat2(directory_fd, source, target, RENAME_EXCHANGE, label)
        except Renameat2Unsupported as error:
            raise LustrePosixRenameUnsupported(
                f"{label} occupied-target exchange is unsupported after "
                "renameat2 rejection"
            ) from error
    except BaseException as error:
        try:
            fsync_and_reauthenticate_entries(
                directory_fd, (source, target), f"{label} ambiguous exchange"
            )
        except BaseException as durable_error:
            raise ValueError(
                f"{label} ambiguous exchange state is not durably authenticated"
            ) from durable_error
        pretransition = (
            source_expected.inode_is_bound(directory_fd, source)
            and target_expected.inode_is_bound(directory_fd, target)
        )
        posttransition = (
            source_expected.inode_is_bound(directory_fd, target)
            and target_expected.inode_is_bound(directory_fd, source)
        )
        if pretransition:
            source_expected.assert_bound(directory_fd, source, label)
            target_expected.assert_bound(
                directory_fd, target, f"{label} predecessor"
            )
            if isinstance(error, LustrePosixRenameUnsupported):
                raise error
            raise ValueError(f"{label} exchange failed before mutation") from error
        if posttransition:
            try:
                source_expected.accept_bound_rename(directory_fd, target, label)
                target_expected.accept_bound_rename(
                    directory_fd, source, f"{label} predecessor"
                )
            except BaseException as binding_error:
                raise ValueError(
                    f"{label} ambiguous exchange state was durably preserved"
                ) from binding_error
            operation_error = error
        else:
            raise ValueError(
                f"{label} ambiguous exchange state was durably preserved"
            ) from error
    if operation_error is None:
        try:
            source_expected.accept_bound_rename(directory_fd, target, label)
            target_expected.accept_bound_rename(
                directory_fd, source, f"{label} predecessor"
            )
            os.fsync(directory_fd)
            mutation_guard()
            source_expected.assert_bound(directory_fd, target, label)
            target_expected.assert_bound(
                directory_fd, source, f"{label} predecessor"
            )
        except BaseException as error:
            operation_error = error
    if operation_error is not None:
        try:
            fsync_and_reauthenticate_entries(
                directory_fd, (source, target), f"{label} failed exchange"
            )
        except BaseException as durable_error:
            raise ValueError(
                f"{label} post-exchange state is not durably authenticated"
            ) from durable_error
        error = operation_error
        if isinstance(error, MutationAuthorityError):
            raise error
        try:
            mutation_guard()
        except MutationAuthorityError as authority_error:
            raise authority_error from error
        raise ValueError(
            f"{error}; {label} durable forward recovery state was preserved"
        ) from error


def quarantine_bound_predecessor(source_directory_fd: int, source: str,
                                 expected: RegularFileBinding,
                                 source_directory: Path,
                                 source_directory_profile: os.stat_result,
                                 mutation_path: Path, label: str, *,
                                 retained: tuple[
                                     tuple[str, RegularFileBinding, str], ...
                                 ] = ()) -> Path:
    """Forensically retire exact bytes and preserve ambiguous forward state."""

    forensic_parent = source_directory.parent
    forensic_identity = hashlib.sha256(
        (
            f"{source_directory.absolute()}\0{source}\0"
            f"{expected.digest or hashlib.sha256(expected.payload or b'').hexdigest()}"
        ).encode("utf-8")
    ).hexdigest()
    forensic_name = (
        f".cgl_lf_stage_i_replaced_{forensic_identity}.forensic"
    )
    forensic_fd = os.open(
        forensic_parent,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        forensic_parent_profile = require_directory_descriptor_binding(
            forensic_parent, forensic_fd, "metadata forensic parent"
        )

        def mutation_guard() -> None:
            require_metadata_mutation_boundary(
                mutation_path,
                source_directory,
                source_directory_fd,
                source_directory_profile,
                label,
            )
            require_metadata_mutation_boundary(
                forensic_parent / forensic_name,
                forensic_parent,
                forensic_fd,
                forensic_parent_profile,
                f"{label} forensic",
            )

        def persist_retirement_state(state_label: str) -> None:
            errors = []
            for durable_fd, names, durable_label in (
                (
                    source_directory_fd,
                    (source, *(name for name, _binding, _label in retained)),
                    f"{state_label} source",
                ),
                (
                    forensic_fd,
                    (forensic_name,),
                    f"{state_label} forensic",
                ),
            ):
                try:
                    fsync_and_reauthenticate_entries(
                        durable_fd, names, durable_label
                    )
                except BaseException as error:
                    errors.append(error)
            if errors:
                raise ValueError(
                    f"{state_label} is not durably authenticated"
                ) from errors[0]

        mutation_guard()
        expected.assert_bound(source_directory_fd, source, label)
        for name, binding, retained_label in retained:
            binding.assert_bound(source_directory_fd, name, retained_label)
        operation_error = None
        try:
            try:
                renameat2_between(
                    source_directory_fd,
                    source,
                    forensic_fd,
                    forensic_name,
                    RENAME_NOREPLACE,
                    label,
                )
            except Renameat2Unsupported as error:
                def linked_retirement_guard() -> None:
                    mutation_guard()
                    for name, binding, retained_label in retained:
                        binding.assert_bound(
                            source_directory_fd, name, retained_label
                        )

                link_bound_noreplace_between(
                    source_directory_fd,
                    source,
                    forensic_fd,
                    forensic_name,
                    expected,
                    label,
                    linked_retirement_guard,
                )
        except BaseException as error:
            try:
                persist_retirement_state(f"{label} ambiguous retirement")
            except BaseException as durable_error:
                raise ValueError(
                    f"{label} ambiguous retirement state is not durably authenticated"
                ) from durable_error
            if expected.inode_is_bound(source_directory_fd, source):
                if expected.inode_is_bound(forensic_fd, forensic_name):
                    expected.assert_bound_hardlink(
                        source_directory_fd,
                        source,
                        2,
                        f"{label} linked source",
                    )
                    expected.assert_bound_hardlink(
                        forensic_fd,
                        forensic_name,
                        2,
                        f"{label} linked forensic recovery",
                    )
                    raise ValueError(
                        f"{label} durable two-link recovery state "
                        "was preserved"
                    ) from error
                expected.assert_bound(source_directory_fd, source, label)
                for name, binding, retained_label in retained:
                    binding.assert_bound(source_directory_fd, name, retained_label)
                if isinstance(error, FileExistsError):
                    raise ValueError(
                        f"{label} forensic target already exists"
                    ) from error
                raise ValueError(f"{label} retirement failed before mutation") from error
            if (
                expected.inode_is_bound(forensic_fd, forensic_name)
                and regular_file_profile(
                    source_directory_fd, source, f"{label} ambiguous source"
                ) is None
            ):
                try:
                    expected.accept_bound_rename(forensic_fd, forensic_name, label)
                    for name, binding, retained_label in retained:
                        binding.assert_bound(
                            source_directory_fd, name, retained_label
                        )
                except BaseException as binding_error:
                    raise ValueError(
                        f"{label} ambiguous retirement state was durably preserved"
                    ) from binding_error
                operation_error = error
            else:
                raise ValueError(
                    f"{label} ambiguous retirement state was durably preserved"
                ) from error
        if operation_error is None:
            try:
                expected.accept_bound_rename(forensic_fd, forensic_name, label)
                require_regular_entry_absent(
                    source_directory_fd, source, f"{label} retired source"
                )
                for name, binding, retained_label in retained:
                    binding.assert_bound(source_directory_fd, name, retained_label)
                os.fsync(source_directory_fd)
                os.fsync(forensic_fd)
                mutation_guard()
                expected.assert_bound(forensic_fd, forensic_name, label)
                require_regular_entry_absent(
                    source_directory_fd, source, f"{label} retired source"
                )
                for name, binding, retained_label in retained:
                    binding.assert_bound(source_directory_fd, name, retained_label)
            except BaseException as error:
                operation_error = error
        if operation_error is not None:
            try:
                persist_retirement_state(f"{label} failed retirement")
            except BaseException as durable_error:
                raise ValueError(
                    f"{label} post-retirement state is not durably authenticated"
                ) from durable_error
            error = operation_error
            if isinstance(error, MutationAuthorityError):
                raise error
            try:
                mutation_guard()
            except MutationAuthorityError as authority_error:
                raise authority_error from error
            raise ValueError(
                f"{error}; {label} durable forensic recovery state was preserved"
            ) from error
    finally:
        os.close(forensic_fd)
    return forensic_parent / forensic_name


METADATA_TEMPORARY_SUFFIX = ".cgl-lf-stage-i.tmp"
METADATA_PREDECESSOR_RECOVERY_SUFFIX = ".predecessor-recovery"
TRANSACTION_METADATA_TEMPORARY_PATTERN = re.compile(
    rf"\.(?P<target>[^/]+\.json){re.escape(METADATA_TEMPORARY_SUFFIX)}"
)
LEGACY_TRANSACTION_METADATA_TEMPORARY_PATTERN = re.compile(
    r"\.(?P<target>[^/]+\.json)\.[0-9]+\.[0-9a-f]{32}\.tmp"
)


def metadata_temporary_name(target: str) -> str:
    """Return the sole controller-owned atomic temporary for one direct child."""

    if Path(target).name != target or target in {"", ".", ".."}:
        raise ValueError(f"metadata target name is invalid: {target!r}")
    return f".{target}{METADATA_TEMPORARY_SUFFIX}"


def metadata_predecessor_recovery_name(temporary_name: str) -> str:
    """Return the deterministic predecessor recovery for one metadata target."""

    require_direct_child_name(temporary_name, "metadata temporary")
    return f"{temporary_name}{METADATA_PREDECESSOR_RECOVERY_SUFFIX}"


def metadata_temporary_public_target(temporary_name: str) -> str | None:
    """Return the public target encoded by one deterministic controller temporary."""

    if (
        not temporary_name.startswith(".")
        or not temporary_name.endswith(METADATA_TEMPORARY_SUFFIX)
    ):
        return None
    target = temporary_name[1:-len(METADATA_TEMPORARY_SUFFIX)]
    try:
        require_direct_child_name(target, "metadata temporary target")
    except ValueError:
        return None
    return target


def require_owner_controlled_metadata_recovery(
    directory_fd: int,
    name: str,
    binding: RegularFileBinding,
    label: str,
    *,
    mode: int | None = None,
    expected_payload: bytes | None = None,
) -> None:
    """Require one exact single-link deterministic metadata recovery file."""

    binding.assert_bound(directory_fd, name, label)
    profile = os.fstat(binding.descriptor)
    if (
        profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or stat.S_IMODE(profile.st_mode) & 0o022
        or (mode is not None and stat.S_IMODE(profile.st_mode) != mode)
        or binding.payload is None
        or (
            expected_payload is not None
            and binding.payload != expected_payload
        )
    ):
        raise ValueError(f"{label} is not an exact owner-controlled recovery file")


def unlink_bound_deterministic_recovery(
    directory_fd: int,
    name: str,
    expected: RegularFileBinding,
    label: str,
    mutation_guard,
) -> None:
    """Durably remove one exact deterministic recovery file."""

    mutation_guard()
    require_owner_controlled_metadata_recovery(
        directory_fd, name, expected, label
    )
    unlink_error = None
    try:
        os.unlink(name, dir_fd=directory_fd)
    except BaseException as error:
        unlink_error = error
    os.fsync(directory_fd)
    mutation_guard()
    if expected.inode_is_bound(directory_fd, name):
        expected.assert_bound(directory_fd, name, label)
        raise ValueError(f"{label} remains after durable unlink") from unlink_error
    require_regular_entry_absent(directory_fd, name, f"{label} removed recovery")


def write_bound_metadata_recovery(
    directory_fd: int,
    name: str,
    payload: bytes,
    mode: int,
    label: str,
    mutation_guard,
) -> RegularFileBinding:
    """Create one exact deterministic recovery file without namespace replacement."""

    mutation_guard()
    require_regular_entry_absent(directory_fd, name, label)
    descriptor = None
    try:
        descriptor = os.open(
            name,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
            0o600,
            dir_fd=directory_fd,
        )
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            if written <= 0:
                raise ValueError(f"{label} write made no progress")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    except BaseException as error:
        if descriptor is not None:
            try:
                os.fsync(descriptor)
            except BaseException:
                pass
        try:
            os.fsync(directory_fd)
        except BaseException:
            pass
        raise error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    os.fsync(directory_fd)
    mutation_guard()
    retained = open_regular_file_binding(
        directory_fd,
        name,
        label,
        expected_payload=payload,
    )
    if retained is None:
        raise ValueError(f"{label} is unavailable after creation")
    try:
        require_owner_controlled_metadata_recovery(
            directory_fd,
            name,
            retained,
            label,
            mode=mode,
            expected_payload=payload,
        )
    except BaseException:
        retained.close()
        raise
    return retained


def rewrite_bound_metadata_in_place(
    directory_fd: int,
    target: str,
    current: RegularFileBinding,
    successor_name: str,
    successor: RegularFileBinding,
    predecessor_name: str,
    predecessor: RegularFileBinding,
    mode: int,
    label: str,
    mutation_guard,
) -> None:
    """Complete one journaled metadata successor on the bound public inode."""

    mutation_guard()
    require_owner_controlled_metadata_recovery(
        directory_fd,
        successor_name,
        successor,
        f"{label} successor",
        mode=mode,
    )
    require_owner_controlled_metadata_recovery(
        directory_fd,
        predecessor_name,
        predecessor,
        f"{label} predecessor recovery",
        mode=0o400,
    )
    current.assert_bound(directory_fd, target, label)
    current_profile = os.fstat(current.descriptor)
    current_mode = stat.S_IMODE(current_profile.st_mode)
    if (
        current_profile.st_uid != os.geteuid()
        or current_profile.st_nlink != 1
        or current_mode & 0o022
        or current.payload is None
        or (
            current_mode != 0o600
            and current.payload != predecessor.payload
        )
    ):
        raise ValueError(f"{label} public inode is not a recoverable transition state")

    descriptor = None
    try:
        if current_mode != 0o600:
            os.fchmod(current.descriptor, 0o600)
            os.fsync(current.descriptor)
            os.fsync(directory_fd)
            mutation_guard()
        descriptor = os.open(
            target,
            os.O_RDWR | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        opened = os.fstat(descriptor)
        if (
            (opened.st_dev, opened.st_ino)
            != (current_profile.st_dev, current_profile.st_ino)
            or opened.st_uid != os.geteuid()
            or opened.st_nlink != 1
            or stat.S_IMODE(opened.st_mode) != 0o600
        ):
            raise ValueError(f"{label} public inode changed before in-place rewrite")
        successor.assert_bound(
            directory_fd, successor_name, f"{label} successor"
        )
        predecessor.assert_bound(
            directory_fd, predecessor_name, f"{label} predecessor recovery"
        )
        mutation_guard()
        os.ftruncate(descriptor, 0)
        os.lseek(descriptor, 0, os.SEEK_SET)
        payload = successor.payload
        if payload is None:
            raise ValueError(f"{label} successor bytes are unavailable")
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            if written <= 0:
                raise ValueError(f"{label} in-place rewrite made no progress")
            offset += written
        os.fsync(descriptor)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
        os.fsync(directory_fd)
        mutation_guard()
        successor.assert_bound(
            directory_fd, successor_name, f"{label} successor"
        )
        predecessor.assert_bound(
            directory_fd, predecessor_name, f"{label} predecessor recovery"
        )
    except BaseException as error:
        if descriptor is not None:
            try:
                os.fsync(descriptor)
            except BaseException:
                pass
        try:
            os.fsync(directory_fd)
            successor.assert_bound(
                directory_fd, successor_name, f"{label} successor"
            )
            predecessor.assert_bound(
                directory_fd, predecessor_name, f"{label} predecessor recovery"
            )
        except BaseException as durable_error:
            raise ValueError(
                f"{label} in-place recovery state is not durably authenticated"
            ) from durable_error
        raise ValueError(
            f"{error}; {label} durable in-place forward recovery state was preserved"
        ) from error
    finally:
        if descriptor is not None:
            os.close(descriptor)

    completed = open_regular_file_binding(
        directory_fd,
        target,
        f"{label} completed successor",
        expected_payload=successor.payload,
    )
    if completed is None:
        raise ValueError(f"{label} completed successor is unavailable")
    try:
        require_requested_regular_file_profile(
            completed, mode, f"{label} completed successor"
        )
        durably_authenticate_bound_file(
            directory_fd, target, completed, f"{label} completed successor"
        )
    finally:
        completed.close()


def transaction_metadata_temporary_target(name: str) -> str | None:
    """Return the journal target encoded by one recognized temporary name."""

    for pattern in (
        TRANSACTION_METADATA_TEMPORARY_PATTERN,
        LEGACY_TRANSACTION_METADATA_TEMPORARY_PATTERN,
    ):
        match = pattern.fullmatch(name)
        if match is not None:
            return str(match.group("target"))
    return None


def retire_metadata_temporary(parent: Path, parent_fd: int, temporary_name: str,
                              parent_profile: os.stat_result, label: str,
                              ) -> Path | None:
    """Forensically retire one authenticated crash remnant without pathname unlink."""

    require_directory_descriptor_binding(
        parent, parent_fd, f"{label} parent", parent_profile
    )
    predecessor_name = metadata_predecessor_recovery_name(temporary_name)
    predecessor = open_regular_file_binding(
        parent_fd,
        predecessor_name,
        f"{label} predecessor recovery",
        allow_absent=True,
        allow_unreadable=True,
    )
    temporary = open_regular_file_binding(
        parent_fd,
        temporary_name,
        label,
        allow_absent=True,
        allow_unreadable=True,
    )
    if temporary is None:
        if predecessor is not None:
            predecessor.close()
            raise ValueError(f"{label} predecessor recovery lacks its successor")
        return None
    try:
        public_target = metadata_temporary_public_target(temporary_name)
        if predecessor is not None:
            target = None
            try:
                if public_target is None:
                    raise ValueError(
                        f"{label} in-place recovery lacks its public target identity"
                    )
                require_owner_controlled_metadata_recovery(
                    parent_fd,
                    temporary_name,
                    temporary,
                    f"{label} successor",
                )
                require_owner_controlled_metadata_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"{label} predecessor recovery",
                    mode=0o400,
                )
                target = open_regular_file_binding(
                    parent_fd,
                    public_target,
                    f"{label} public target",
                    allow_unreadable=True,
                )
                if target is None:
                    raise ValueError(
                        f"{label} in-place recovery public target is unavailable"
                    )

                def in_place_guard() -> None:
                    for name, retained_label in (
                        (temporary_name, label),
                        (predecessor_name, f"{label} predecessor recovery"),
                        (public_target, f"{label} public target"),
                    ):
                        require_metadata_mutation_boundary(
                            parent / name,
                            parent,
                            parent_fd,
                            parent_profile,
                            retained_label,
                        )

                final_mode = stat.S_IMODE(temporary.profile.st_mode)
                rewrite_bound_metadata_in_place(
                    parent_fd,
                    public_target,
                    target,
                    temporary_name,
                    temporary,
                    predecessor_name,
                    predecessor,
                    final_mode,
                    label,
                    in_place_guard,
                )
                unlink_bound_deterministic_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"{label} predecessor recovery",
                    in_place_guard,
                )
                predecessor.close()
                predecessor = None
                unlink_bound_deterministic_recovery(
                    parent_fd,
                    temporary_name,
                    temporary,
                    label,
                    in_place_guard,
                )
                return None
            finally:
                if target is not None:
                    target.close()

        if (
            temporary.profile.st_uid != os.geteuid()
            or temporary.profile.st_nlink not in {1, 2}
            or stat.S_IMODE(temporary.profile.st_mode) & 0o022
        ):
            raise ValueError(f"{label} is not an authenticated controller temporary")
        if temporary.profile.st_nlink == 2:
            public_target = metadata_temporary_public_target(temporary_name)
            if (
                public_target is None
                or not temporary.inode_is_bound(parent_fd, public_target)
            ):
                raise ValueError(
                    f"{label} two-link publication lacks its deterministic "
                    "authenticated public target"
                )

            def recovery_guard() -> None:
                require_metadata_mutation_boundary(
                    parent / temporary_name,
                    parent,
                    parent_fd,
                    parent_profile,
                    label,
                )
                require_metadata_mutation_boundary(
                    parent / public_target,
                    parent,
                    parent_fd,
                    parent_profile,
                    f"{label} public target",
                )

            unlink_bound_hardlink_alias(
                parent_fd,
                temporary_name,
                public_target,
                temporary,
                f"{label} completed hard-link publication",
                recovery_guard,
            )
            return None
        if public_target is not None:
            completed = open_regular_file_binding(
                parent_fd,
                public_target,
                f"{label} completed public target",
                allow_absent=True,
                allow_unreadable=True,
            )
            if completed is not None:
                try:
                    if (
                        temporary.payload is not None
                        and completed.payload == temporary.payload
                    ):
                        final_mode = stat.S_IMODE(temporary.profile.st_mode)
                        require_owner_controlled_metadata_recovery(
                            parent_fd,
                            temporary_name,
                            temporary,
                            label,
                            mode=final_mode,
                        )
                        require_requested_regular_file_profile(
                            completed,
                            final_mode,
                            f"{label} completed public target",
                        )

                        def completed_guard() -> None:
                            for name, retained_label in (
                                (temporary_name, label),
                                (public_target, f"{label} completed public target"),
                            ):
                                require_metadata_mutation_boundary(
                                    parent / name,
                                    parent,
                                    parent_fd,
                                    parent_profile,
                                    retained_label,
                                )

                        unlink_bound_deterministic_recovery(
                            parent_fd,
                            temporary_name,
                            temporary,
                            label,
                            completed_guard,
                        )
                        durably_authenticate_bound_file(
                            parent_fd,
                            public_target,
                            completed,
                            f"{label} completed public target",
                        )
                        return None
                finally:
                    completed.close()
        return quarantine_bound_predecessor(
            parent_fd,
            temporary_name,
            temporary,
            parent,
            parent_profile,
            parent / temporary_name,
            label,
        )
    finally:
        if predecessor is not None:
            predecessor.close()
        temporary.close()


def validate_canonical_root_lock_binding(
    root: Path,
    stream: object,
    root_fd: int | None = None,
    root_profile: os.stat_result | None = None,
    authority_profile: tuple[tuple[str, object], ...] | None = None,
) -> None:
    """Require the held lock descriptor and pathname to remain identical."""

    require_canonical_root_authority(root, authority_profile)
    if root_fd is not None:
        require_directory_descriptor_binding(
            root, root_fd, "canonical Stage I root", root_profile
        )
    lock_path = canonical_root_lock_path(root)
    descriptor = stream.fileno()
    opened = os.fstat(descriptor)
    require_canonical_root_lock_profile(opened, lock_path)
    require_owner_symlink_free_path(lock_path, "Stage I lock")
    try:
        named = (
            os.stat(lock_path.name, dir_fd=root_fd, follow_symlinks=False)
            if root_fd is not None
            else lock_path.lstat()
        )
    except OSError as error:
        raise ValueError(
            f"Stage I lock path changed during canonical mutation: {lock_path}"
        ) from error
    if (
        not stat.S_ISREG(named.st_mode)
        or (named.st_dev, named.st_ino) != (opened.st_dev, opened.st_ino)
        or named.st_uid != opened.st_uid
        or named.st_nlink != opened.st_nlink
        or named.st_size != opened.st_size
        or stat.S_IMODE(named.st_mode) != stat.S_IMODE(opened.st_mode)
    ):
        raise ValueError(
            f"Stage I lock path changed during canonical mutation: {lock_path}"
        )


def require_canonical_mutation_lock(path: Path) -> None:
    """Reauthenticate the canonical mutation lock for one target path."""

    canonical = DEFAULT_ROOT.expanduser().absolute()
    absolute = path.expanduser().absolute()
    try:
        absolute.relative_to(canonical)
    except ValueError:
        return
    require_canonical_root_authority(canonical)
    active = _ACTIVE_ROOT_LOCKS.get(canonical)
    if active is None:
        raise ValueError(f"canonical Stage I mutation lacks its root lock: {absolute}")
    validate_canonical_root_lock_binding(
        canonical, active[0], active[2], active[3], active[4]
    )


def copy_file(source: Path, destination: Path) -> None:
    """Publish one descriptor-bound exact copy without clobbering a raced target.

    Source and destination parents remain bound through publication.  A failed
    or ambiguous copy leaves either the prior public namespace or an exact
    fsynced temporary/public recovery state that a retry can classify.
    """

    source = source.expanduser().absolute()
    destination = destination.expanduser().absolute()
    require_canonical_mutation_lock(destination)
    source_parent_fd = os.open(
        source.parent,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    destination_parent_fd = os.open(
        destination.parent,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    source_binding = None
    destination_binding = None
    replacement = None
    temporary_name = metadata_temporary_name(destination.name)
    try:
        source_parent_profile = require_directory_descriptor_binding(
            source.parent, source_parent_fd, "copy source parent"
        )
        destination_parent_profile = require_directory_descriptor_binding(
            destination.parent, destination_parent_fd, "copy destination parent"
        )
        source_binding = open_regular_file_binding(
            source_parent_fd,
            source.name,
            f"copy source {source}",
            retain_payload=False,
        )
        if source_binding is None or source_binding.digest is None:
            raise ValueError(f"copy source is unavailable: {source}")
        source_mode = stat.S_IMODE(source_binding.profile.st_mode)
        if source_mode & 0o022:
            raise ValueError(
                f"copy source mode would publish a group/world-writable "
                f"destination: {source}"
            )

        def source_guard() -> None:
            require_directory_descriptor_binding(
                source.parent,
                source_parent_fd,
                "copy source parent",
                source_parent_profile,
            )
            source_binding.assert_bound(
                source_parent_fd, source.name, f"copy source {source}"
            )

        def destination_guard() -> None:
            require_metadata_mutation_boundary(
                destination,
                destination.parent,
                destination_parent_fd,
                destination_parent_profile,
                "copy",
            )

        source_guard()
        destination_guard()
        retire_metadata_temporary(
            destination.parent,
            destination_parent_fd,
            temporary_name,
            destination_parent_profile,
            f"copy temporary {destination}",
        )
        source_guard()
        destination_guard()
        destination_binding = open_regular_file_binding(
            destination_parent_fd,
            destination.name,
            f"copy destination {destination}",
            allow_absent=True,
            retain_payload=False,
        )
        if destination_binding is not None:
            if destination_binding.digest != source_binding.digest:
                raise ValueError(f"copy destination already exists with different bytes: {destination}")
            require_requested_regular_file_profile(
                destination_binding, source_mode, f"copy destination {destination}"
            )
            durably_authenticate_bound_file(
                destination_parent_fd,
                destination.name,
                destination_binding,
                f"copy destination {destination}",
            )
            source_guard()
            destination_guard()
            return
        try:
            descriptor = os.open(
                temporary_name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o600,
                dir_fd=destination_parent_fd,
            )
            try:
                copy_bound_file_payload(source_binding, descriptor)
                os.fchmod(descriptor, source_mode)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        except BaseException as error:
            partial = open_regular_file_binding(
                destination_parent_fd,
                temporary_name,
                f"copy partial temporary {destination}",
                allow_absent=True,
                retain_payload=False,
            )
            try:
                if partial is None:
                    source_guard()
                    destination_guard()
                    raise ValueError(
                        f"copy failed before temporary creation: {destination}"
                    ) from error
                if (
                    partial.profile.st_uid != os.geteuid()
                    or partial.profile.st_nlink != 1
                    or stat.S_IMODE(partial.profile.st_mode) & 0o022
                ):
                    raise ValueError(
                        f"copy temporary is not owner-controlled: {destination}"
                    ) from error
                durably_authenticate_bound_file(
                    destination_parent_fd,
                    temporary_name,
                    partial,
                    f"copy partial temporary {destination}",
                )
                source_guard()
                destination_guard()
            finally:
                if partial is not None:
                    partial.close()
            raise ValueError(
                f"{error}; copy durable temporary recovery state was preserved: "
                f"{destination}"
            ) from error
        source_guard()
        replacement = open_regular_file_binding(
            destination_parent_fd,
            temporary_name,
            f"copy replacement {destination}",
            retain_payload=False,
            expected_sha256=source_binding.digest,
        )
        if replacement is None:
            raise ValueError(f"copy replacement is unavailable: {destination}")

        def mutation_guard() -> None:
            source_guard()
            destination_guard()
            require_requested_regular_file_profile(
                replacement, source_mode, f"copy replacement {destination}"
            )

        rename_bound_noreplace(
            destination_parent_fd,
            temporary_name,
            destination.name,
            replacement,
            f"copy publication {destination}",
            mutation_guard,
        )
        mutation_guard()
        replacement.assert_bound(
            destination_parent_fd,
            destination.name,
            f"copy destination {destination}",
        )
    finally:
        if replacement is not None:
            replacement.close()
        if destination_binding is not None:
            destination_binding.close()
        if source_binding is not None:
            source_binding.close()
        os.close(destination_parent_fd)
        os.close(source_parent_fd)


def mkdir_durable(path: Path) -> None:
    """Create and durably classify each descriptor-relative directory component.

    POSIX ``mkdirat`` does not return the created inode descriptor.  Each child
    is therefore opened immediately from its bound parent, then the child and
    parent are fsynced and the public pathname is reauthenticated.  Authority
    loss preserves that durable forward state; this helper never rolls it back.
    Creation requests 0755 and requires exact inherited 2755 beneath a setgid
    parent, so ambient umask can never broaden the resulting profile.
    """

    path = path.expanduser().absolute()
    require_canonical_mutation_lock(path)
    missing = []
    current = path
    while True:
        try:
            profile = os.lstat(current)
        except FileNotFoundError:
            missing.append(current)
            current = current.parent
            continue
        except OSError as error:
            raise ValueError(f"directory creation path is unavailable: {current}") from error
        if not stat.S_ISDIR(profile.st_mode) or stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"directory creation path is not a directory: {current}")
        break
    targets = list(reversed(missing)) or [path]
    for target in targets:
        parent = target.parent
        parent_fd = os.open(
            parent,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
        child_fd = None
        try:
            parent_profile = require_directory_descriptor_binding(
                parent, parent_fd, "directory creation parent"
            )

            def mutation_guard() -> None:
                require_metadata_mutation_boundary(
                    target,
                    parent,
                    parent_fd,
                    parent_profile,
                    "directory creation",
                )

            mutation_guard()
            operation_error = None
            try:
                os.mkdir(target.name, DURABLE_DIRECTORY_MODE, dir_fd=parent_fd)
            except FileExistsError:
                pass
            except BaseException as error:
                operation_error = error
            try:
                child_fd = os.open(
                    target.name,
                    os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=parent_fd,
                )
            except BaseException as error:
                try:
                    child_fd = os.open(
                        target.name,
                        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                        dir_fd=parent_fd,
                    )
                except BaseException as retry_error:
                    os.fsync(parent_fd)
                    mutation_guard()
                    if operation_error is not None:
                        raise ValueError(
                            f"directory creation failed before a durable outcome: {target}"
                        ) from operation_error
                    raise ValueError(
                        f"directory creation outcome is unavailable: {target}"
                    ) from retry_error
            child_profile = os.fstat(child_fd)
            try:
                require_bound_directory_entry(
                    parent_fd,
                    target.name,
                    child_fd,
                    child_profile,
                    "directory creation outcome",
                )
                os.fsync(child_fd)
                os.fsync(parent_fd)
                require_durable_directory_profile(
                    child_profile, parent_profile, target
                )
                mutation_guard()
                require_directory_descriptor_binding(
                    target, child_fd, "directory creation outcome", child_profile
                )
            except BaseException as error:
                try:
                    os.fsync(child_fd)
                    os.fsync(parent_fd)
                    require_bound_directory_entry(
                        parent_fd,
                        target.name,
                        child_fd,
                        child_profile,
                        "directory creation recovery outcome",
                    )
                except BaseException as durable_error:
                    raise ValueError(
                        f"directory creation outcome is not durably authenticated: "
                        f"{target}"
                    ) from durable_error
                if isinstance(error, MutationAuthorityError):
                    raise error
                try:
                    mutation_guard()
                except MutationAuthorityError as authority_error:
                    raise authority_error from error
                raise ValueError(
                    f"{error}; directory creation durable forward recovery state "
                    f"was preserved: {target}"
                ) from error
        finally:
            if child_fd is not None:
                os.close(child_fd)
            os.close(parent_fd)
    require_canonical_mutation_lock(path)


def write_text(path: Path, value: str, mode: int | None = None, *,
               expected_predecessor: bytes | None = None) -> None:
    """Durably publish text with recoverable existing-target transitions."""

    if mode is not None and mode != stat.S_IMODE(mode):
        raise ValueError(f"metadata publication mode is invalid: {mode!r}")
    payload = value.encode("utf-8")
    path = path.expanduser().absolute()
    require_canonical_mutation_lock(path)
    parent = path.parent
    try:
        parent_fd = os.open(
            parent, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
        )
    except OSError as error:
        raise ValueError(f"metadata parent directory is unavailable: {parent}") from error
    temporary_name = metadata_temporary_name(path.name)
    predecessor_name = metadata_predecessor_recovery_name(temporary_name)
    initial = None
    replacement = None
    predecessor = None
    try:
        parent_profile = require_directory_descriptor_binding(
            parent, parent_fd, "metadata parent"
        )

        def mutation_guard() -> None:
            require_metadata_mutation_boundary(
                path, parent, parent_fd, parent_profile, "metadata"
            )

        def authenticate_completed_target(final_mode: int | None) -> None:
            completed = open_regular_file_binding(
                parent_fd,
                path.name,
                f"metadata completed successor {path}",
                expected_payload=payload,
            )
            if completed is None:
                raise ValueError(f"metadata completed successor is unavailable: {path}")
            try:
                if final_mode is None:
                    require_owner_controlled_metadata_recovery(
                        parent_fd,
                        path.name,
                        completed,
                        f"metadata completed successor {path}",
                    )
                else:
                    require_requested_regular_file_profile(
                        completed,
                        final_mode,
                        f"metadata completed successor {path}",
                    )
                durably_authenticate_bound_file(
                    parent_fd,
                    path.name,
                    completed,
                    f"metadata completed successor {path}",
                )
            finally:
                completed.close()

        mutation_guard()
        initial = open_regular_file_binding(
            parent_fd,
            path.name,
            "metadata predecessor",
            allow_absent=True,
            allow_unreadable=True,
        )
        replacement = open_regular_file_binding(
            parent_fd,
            temporary_name,
            f"metadata replacement {path}",
            allow_absent=True,
            allow_unreadable=True,
        )
        predecessor = open_regular_file_binding(
            parent_fd,
            predecessor_name,
            f"metadata predecessor recovery {path}",
            allow_absent=True,
            allow_unreadable=True,
        )

        if (
            replacement is not None
            and replacement.profile.st_nlink == 2
            and initial is not None
            and replacement.inode_is_bound(parent_fd, path.name)
        ):
            if predecessor is not None:
                raise ValueError(
                    f"metadata hard-link publication has an unexpected "
                    f"predecessor recovery: {path}"
                )
            replacement.close()
            replacement = None
            initial.close()
            initial = None
            retire_metadata_temporary(
                parent,
                parent_fd,
                temporary_name,
                parent_profile,
                f"metadata temporary {path}",
            )
            initial = open_regular_file_binding(
                parent_fd,
                path.name,
                "metadata predecessor",
                allow_absent=True,
                allow_unreadable=True,
            )

        if initial is not None and initial.payload == payload and (
            mode is None or stat.S_IMODE(initial.profile.st_mode) == mode
        ):
            authenticate_completed_target(mode)
            if predecessor is not None:
                require_owner_controlled_metadata_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"metadata predecessor recovery {path}",
                    mode=0o400,
                    expected_payload=expected_predecessor,
                )
                unlink_bound_deterministic_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"metadata predecessor recovery {path}",
                    mutation_guard,
                )
                predecessor.close()
                predecessor = None
                authenticate_completed_target(mode)
            if replacement is not None:
                if replacement.payload == payload:
                    require_owner_controlled_metadata_recovery(
                        parent_fd,
                        temporary_name,
                        replacement,
                        f"metadata replacement {path}",
                        mode=mode or stat.S_IMODE(initial.profile.st_mode),
                        expected_payload=payload,
                    )
                    unlink_bound_deterministic_recovery(
                        parent_fd,
                        temporary_name,
                        replacement,
                        f"metadata replacement {path}",
                        mutation_guard,
                    )
                else:
                    quarantine_bound_predecessor(
                        parent_fd,
                        temporary_name,
                        replacement,
                        parent,
                        parent_profile,
                        path,
                        f"metadata predecessor {path}",
                        retained=((
                            path.name,
                            initial,
                            f"metadata replacement {path}",
                        ),),
                    )
                replacement.close()
                replacement = None
                authenticate_completed_target(mode)
            return

        if predecessor is not None:
            predecessor_valid = True
            try:
                require_owner_controlled_metadata_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"metadata predecessor recovery {path}",
                    mode=0o400,
                    expected_payload=expected_predecessor,
                )
            except ValueError:
                predecessor_valid = False
            target_is_complete_predecessor = (
                initial is not None
                and initial.payload is not None
                and (
                    expected_predecessor is None
                    or initial.payload == expected_predecessor
                )
            )
            if not predecessor_valid and target_is_complete_predecessor:
                unlink_bound_deterministic_recovery(
                    parent_fd,
                    predecessor_name,
                    predecessor,
                    f"metadata incomplete predecessor recovery {path}",
                    mutation_guard,
                )
                predecessor.close()
                predecessor = None
            elif not predecessor_valid:
                raise ValueError(
                    f"metadata predecessor recovery is not exact: {path}"
                )

        if predecessor is not None:
            if initial is None or replacement is None:
                raise ValueError(
                    f"metadata in-place recovery state is incomplete: {path}"
                )
            final_mode = mode or stat.S_IMODE(replacement.profile.st_mode)
            require_owner_controlled_metadata_recovery(
                parent_fd,
                temporary_name,
                replacement,
                f"metadata replacement {path}",
                mode=final_mode,
                expected_payload=payload,
            )
            rewrite_bound_metadata_in_place(
                parent_fd,
                path.name,
                initial,
                temporary_name,
                replacement,
                predecessor_name,
                predecessor,
                final_mode,
                f"metadata replacement {path}",
                mutation_guard,
            )
            authenticate_completed_target(final_mode)
            unlink_bound_deterministic_recovery(
                parent_fd,
                predecessor_name,
                predecessor,
                f"metadata predecessor recovery {path}",
                mutation_guard,
            )
            predecessor.close()
            predecessor = None
            authenticate_completed_target(final_mode)
            unlink_bound_deterministic_recovery(
                parent_fd,
                temporary_name,
                replacement,
                f"metadata replacement {path}",
                mutation_guard,
            )
            replacement.close()
            replacement = None
            authenticate_completed_target(final_mode)
            return

        if initial is not None:
            require_owner_controlled_metadata_recovery(
                parent_fd,
                path.name,
                initial,
                f"metadata predecessor {path}",
                expected_payload=expected_predecessor,
            )
            final_mode = mode or stat.S_IMODE(initial.profile.st_mode)
        else:
            if expected_predecessor is not None:
                raise ValueError(f"metadata predecessor is unavailable: {path}")
            final_mode = mode or 0o600

        if replacement is not None:
            replacement_is_exact = True
            try:
                require_owner_controlled_metadata_recovery(
                    parent_fd,
                    temporary_name,
                    replacement,
                    f"metadata replacement {path}",
                    mode=final_mode,
                    expected_payload=payload,
                )
            except ValueError:
                replacement_is_exact = False
            if not replacement_is_exact:
                unlink_bound_deterministic_recovery(
                    parent_fd,
                    temporary_name,
                    replacement,
                    f"metadata incomplete replacement {path}",
                    mutation_guard,
                )
                replacement.close()
                replacement = None

        if replacement is None:
            replacement = write_bound_metadata_recovery(
                parent_fd,
                temporary_name,
                payload,
                final_mode,
                f"metadata replacement {path}",
                mutation_guard,
            )

        mutation_guard()
        if initial is None:
            rename_bound_noreplace(
                parent_fd,
                temporary_name,
                path.name,
                replacement,
                f"metadata publication {path}",
                mutation_guard,
            )
            authenticate_completed_target(final_mode)
            return

        try:
            exchange_bound_entries(
                parent_fd,
                temporary_name,
                path.name,
                replacement,
                initial,
                f"metadata replacement {path}",
                mutation_guard,
            )
        except LustrePosixRenameUnsupported:
            if initial.payload is None:
                raise ValueError(
                    f"metadata predecessor bytes are unavailable: {path}"
                ) from None
            predecessor = write_bound_metadata_recovery(
                parent_fd,
                predecessor_name,
                initial.payload,
                0o400,
                f"metadata predecessor recovery {path}",
                mutation_guard,
            )
            rewrite_bound_metadata_in_place(
                parent_fd,
                path.name,
                initial,
                temporary_name,
                replacement,
                predecessor_name,
                predecessor,
                final_mode,
                f"metadata replacement {path}",
                mutation_guard,
            )
            authenticate_completed_target(final_mode)
            unlink_bound_deterministic_recovery(
                parent_fd,
                predecessor_name,
                predecessor,
                f"metadata predecessor recovery {path}",
                mutation_guard,
            )
            predecessor.close()
            predecessor = None
            authenticate_completed_target(final_mode)
            unlink_bound_deterministic_recovery(
                parent_fd,
                temporary_name,
                replacement,
                f"metadata replacement {path}",
                mutation_guard,
            )
            replacement.close()
            replacement = None
            authenticate_completed_target(final_mode)
            return

        try:
            quarantine_bound_predecessor(
                parent_fd,
                temporary_name,
                initial,
                parent,
                parent_profile,
                path,
                f"metadata predecessor {path}",
                retained=((
                    path.name,
                    replacement,
                    f"metadata replacement {path}",
                ),),
            )
        except BaseException as error:
            try:
                fsync_and_reauthenticate_entries(
                    parent_fd,
                    (temporary_name, path.name),
                    f"metadata replacement {path} failed predecessor retirement",
                )
            except BaseException as durable_error:
                raise ValueError(
                    f"metadata replacement recovery state is not "
                    f"durably authenticated: {path}"
                ) from durable_error
            if isinstance(error, MutationAuthorityError):
                raise
            try:
                mutation_guard()
            except MutationAuthorityError as authority_error:
                raise authority_error from error
            raise ValueError(
                f"{error}; metadata replacement durable forward recovery "
                f"state was preserved: {path}"
            ) from error
        mutation_guard()
        replacement.assert_bound(
            parent_fd, path.name, f"metadata replacement {path}"
        )
    finally:
        if predecessor is not None:
            predecessor.close()
        if replacement is not None:
            replacement.close()
        if initial is not None:
            initial.close()
        os.close(parent_fd)


def write_json(path: Path, value: object, mode: int | None = None) -> None:
    """Atomically and durably write stable JSON metadata."""

    write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n", mode=mode)


def stable_json_sha256(value: object) -> str:
    """Return the checksum produced by ``write_json`` for one value."""

    return hashlib.sha256(
        (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    ).hexdigest()


def unlink_durable(path: Path) -> None:
    """Retire one retained entry to a no-clobber forensic inode."""

    path = path.expanduser().absolute()
    require_canonical_mutation_lock(path)
    parent_fd = os.open(
        path.parent, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    )
    initial = None
    try:
        parent_profile = require_directory_descriptor_binding(
            path.parent, parent_fd, "unlink parent"
        )
        initial = open_regular_file_binding(
            parent_fd,
            path.name,
            "unlink",
            allow_absent=True,
            allow_unreadable=True,
        )
        if initial is None:
            raise FileNotFoundError(path)
        quarantine_bound_predecessor(
            parent_fd,
            path.name,
            initial,
            path.parent,
            parent_profile,
            path,
            f"retired entry {path}",
        )
    finally:
        if initial is not None:
            initial.close()
        os.close(parent_fd)


def append_ledger_row(path: Path, row: dict[str, object],
                      prior_rows: list[dict[str, str]]) -> None:
    """Idempotently publish a complete ledger while its journal remains retained.

    The complete prospective CSV is atomically exchanged with the retained
    baseline.  Pre-publication failure leaves the old ledger and any surviving
    temporary is recoverable; ambiguous post-publication retry recognizes only
    the exact prospective ledger and never duplicates the row.
    """

    if frozenset(row) != frozenset(LEDGER_COLUMNS):
        raise ValueError("ledger append row has invalid columns")
    path = path.expanduser().absolute()
    retained_row = {column: str(row[column]) for column in LEDGER_COLUMNS}
    prior_actual = sum(float(item["actual_node_hours"]) for item in prior_rows)
    validate_ledger_cumulative_fields(
        retained_row, prior_actual, "ledger append row"
    )
    current_rows = read_ledger({"ledger": path})
    prospective_rows = [*prior_rows, retained_row]
    if current_rows == prospective_rows:
        expected = ledger_csv_text(prospective_rows).encode("utf-8")
        parent_fd = os.open(
            path.parent,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
        binding = None
        try:
            parent_profile = require_directory_descriptor_binding(
                path.parent, parent_fd, "ledger parent"
            )
            require_metadata_mutation_boundary(
                path, path.parent, parent_fd, parent_profile, "ledger"
            )
            binding = open_regular_file_binding(
                parent_fd,
                path.name,
                "ledger retry result",
                expected_payload=expected,
            )
            if binding is None:
                raise ValueError("ledger retry result is unavailable")
            require_requested_regular_file_profile(
                binding, 0o644, "ledger retry result"
            )
            durably_authenticate_bound_file(
                parent_fd, path.name, binding, "ledger retry result"
            )
            require_metadata_mutation_boundary(
                path, path.parent, parent_fd, parent_profile, "ledger"
            )
        finally:
            if binding is not None:
                binding.close()
            os.close(parent_fd)
        return
    if current_rows != prior_rows:
        raise ValueError("ledger append baseline differs from retained ledger")
    write_text(
        path,
        ledger_csv_text(prospective_rows),
        mode=0o644,
        expected_predecessor=ledger_csv_text(prior_rows).encode("utf-8"),
    )


def ledger_csv_text(rows: list[dict[str, object]]) -> str:
    """Serialize one complete canonical Stage I ledger."""

    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=LEDGER_COLUMNS)
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def canonical_root_lock_path(root: Path) -> Path:
    """Return the cooperative Stage I mutation lock location."""

    return root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock"


def require_canonical_root_lock_profile(profile: os.stat_result,
                                        lock_path: Path) -> None:
    """Require one cooperative canonical-root lock profile."""

    mode = stat.S_IMODE(profile.st_mode)
    if not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"Stage I lock is not a regular file: {lock_path}")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"Stage I lock owner differs: {lock_path}")
    if profile.st_nlink != 1:
        raise ValueError(f"Stage I lock link count differs: {lock_path}")
    if profile.st_size != 0:
        raise ValueError(f"Stage I lock must be empty: {lock_path}")
    if mode & ~0o644:
        raise ValueError(f"Stage I lock mode is too permissive: {lock_path}")


def preserve_created_canonical_lock_evidence(root_fd: int, lock_name: str,
                                             descriptor: int) -> None:
    """Durably retain a bound created lock after pre-yield authority failure.

    Once canonical authority is invalid, no cleanup namespace mutation is
    permitted.  The retained lock is explicit fail-closed recovery evidence.
    """

    try:
        opened = os.fstat(descriptor)
        named = os.stat(lock_name, dir_fd=root_fd, follow_symlinks=False)
        require_canonical_root_lock_profile(opened, Path(lock_name))
        if regular_file_stable_identity(named) != regular_file_stable_identity(opened):
            raise ValueError("created Stage I lock pathname changed before preservation")
        os.fsync(descriptor)
        os.fsync(root_fd)
        retained = os.stat(lock_name, dir_fd=root_fd, follow_symlinks=False)
        if regular_file_stable_identity(retained) != regular_file_stable_identity(
            os.fstat(descriptor)
        ):
            raise ValueError("created Stage I lock changed during preservation")
    except BaseException as error:
        try:
            os.fsync(root_fd)
        except BaseException:
            pass
        raise ValueError(
            "created Stage I lock recovery evidence is not durably authenticated"
        ) from error


@contextmanager
def canonical_root_lock(root: Path):
    """Take a nonblocking reentrant lock for canonical-root mutations."""

    resolved = root.expanduser().absolute()
    if resolved != DEFAULT_ROOT.expanduser().absolute():
        yield
        return
    authority_profile = require_canonical_root_authority(resolved)
    active = _ACTIVE_ROOT_LOCKS.get(resolved)
    if active is not None:
        stream, depth, root_fd, root_profile, authority_profile = active
        validate_canonical_root_lock_binding(
            resolved, stream, root_fd, root_profile, authority_profile
        )
        _ACTIVE_ROOT_LOCKS[resolved] = (
            stream, depth + 1, root_fd, root_profile, authority_profile
        )
        try:
            yield
        finally:
            try:
                validate_canonical_root_lock_binding(
                    resolved, stream, root_fd, root_profile, authority_profile
                )
            finally:
                _ACTIVE_ROOT_LOCKS[resolved] = (
                    stream, depth, root_fd, root_profile, authority_profile
                )
        return
    if not resolved.is_dir():
        raise ValueError(f"canonical Stage I root is unavailable: {resolved}")
    lock_path = canonical_root_lock_path(resolved)
    lock_name = lock_path.name
    try:
        root_fd = os.open(
            resolved, os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
        )
    except OSError as error:
        raise ValueError(f"canonical Stage I root is unavailable: {resolved}") from error
    try:
        root_profile = require_directory_descriptor_binding(
            resolved, root_fd, "canonical Stage I root"
        )
        require_canonical_root_authority(resolved, authority_profile)
    except BaseException:
        os.close(root_fd)
        raise
    created = False
    setup_complete = False
    descriptor = None
    stream = None
    try:
        previous = os.umask(0)
        try:
            try:
                descriptor = os.open(
                    lock_name,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                    0o644,
                    dir_fd=root_fd,
                )
                created = True
            except FileExistsError:
                try:
                    descriptor = os.open(
                        lock_name,
                        os.O_RDWR | os.O_NOFOLLOW,
                        dir_fd=root_fd,
                    )
                except OSError as error:
                    if error.errno == errno.ELOOP:
                        raise ValueError(
                            f"Stage I lock must not be a symlink: {lock_path}"
                        ) from error
                    raise
        finally:
            os.umask(previous)
        require_directory_descriptor_binding(
            resolved, root_fd, "canonical Stage I root", root_profile
        )
        require_canonical_root_authority(resolved, authority_profile)
        profile = os.fstat(descriptor)
        require_canonical_root_lock_profile(profile, lock_path)
        if created:
            os.fsync(root_fd)
        stream = os.fdopen(descriptor, "a+", encoding="utf-8")
        descriptor = None
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            raise ValueError(
                f"another Stage I mutation holds {lock_path}"
            ) from error
        validate_canonical_root_lock_binding(
            resolved, stream, root_fd, root_profile, authority_profile
        )
        if stat.S_IMODE(os.fstat(stream.fileno()).st_mode) != 0o644:
            os.fchmod(stream.fileno(), 0o644)
            os.fsync(stream.fileno())
            os.fsync(root_fd)
            validate_canonical_root_lock_binding(
                resolved, stream, root_fd, root_profile, authority_profile
            )
        _ACTIVE_ROOT_LOCKS[resolved] = (
            stream, 1, root_fd, root_profile, authority_profile
        )
        try:
            validate_canonical_root_lock_binding(
                resolved, stream, root_fd, root_profile, authority_profile
            )
            setup_complete = True
            yield
        finally:
            try:
                validate_canonical_root_lock_binding(
                    resolved, stream, root_fd, root_profile, authority_profile
                )
            finally:
                del _ACTIVE_ROOT_LOCKS[resolved]
                fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    except BaseException:
        if created and not setup_complete:
            cleanup_descriptor = (
                stream.fileno() if stream is not None else descriptor
            )
            if cleanup_descriptor is not None:
                preserve_created_canonical_lock_evidence(
                    root_fd, lock_name, cleanup_descriptor
                )
        raise
    finally:
        if stream is not None:
            stream.close()
        elif descriptor is not None:
            os.close(descriptor)
        os.close(root_fd)


def locked_root_action(function):
    """Lock a mutating action whose root is supplied directly."""

    @wraps(function)
    def wrapped(args: argparse.Namespace):
        root = require_root(Path(args.root), args.allow_local_root)
        with canonical_root_lock(root):
            return function(args)

    return wrapped


def locked_manifest_action(function):
    """Lock a mutating action whose root is retained in its manifest."""

    @wraps(function)
    def wrapped(args: argparse.Namespace):
        manifest = read_manifest(Path(args.manifest).expanduser().resolve())
        root = require_root(
            Path(str(manifest["project_root"])),
            getattr(args, "allow_local_root", False),
        )
        with canonical_root_lock(root):
            retained = read_manifest(Path(args.manifest).expanduser().resolve())
            if Path(str(retained.get("project_root", ""))).resolve() != root:
                raise ValueError("manifest project root changed while acquiring its lock")
            return function(args)

    return wrapped


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_walltime(value: str) -> int:
    """Parse an HH:MM:SS allocation walltime."""

    match = re.fullmatch(r"(\d{2}):([0-5]\d):([0-5]\d)", value)
    if match is None:
        raise ValueError(f"walltime must use HH:MM:SS: {value}")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def node_hours(nodes: int, seconds: int) -> float:
    """Compute allocated node hours."""

    return nodes * seconds / 3600.0


def validate_ledger_numeric_fields(row: dict[str, object],
                                   label: str) -> tuple[float, float, float]:
    """Require finite, non-negative retained Stage I accounting values."""

    try:
        nodes_text = str(row["nodes"])
        elapsed_text = str(row["elapsed_seconds"])
        if (
            re.fullmatch(r"[0-9]+", nodes_text) is None
            or re.fullmatch(r"[0-9]+", elapsed_text) is None
        ):
            raise ValueError
        nodes = int(nodes_text)
        elapsed_seconds = int(elapsed_text)
        requested_seconds = parse_walltime(str(row["requested_walltime"]))
        decimal_fields = {
            key: str(row[key])
            for key in (
                "reserved_node_hours",
                "actual_node_hours",
                "cumulative_stage_i_node_hours",
            )
        }
        if any(
            re.fullmatch(r"(?:0|[1-9][0-9]*)\.[0-9]{6}", value) is None
            for value in decimal_fields.values()
        ):
            raise ValueError
        reserved_node_hours = float(decimal_fields["reserved_node_hours"])
        actual_node_hours = float(decimal_fields["actual_node_hours"])
        cumulative_node_hours = float(decimal_fields["cumulative_stage_i_node_hours"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label} has invalid numeric fields") from error
    if (
        nodes < 1
        or elapsed_seconds < 0
        or not math.isfinite(reserved_node_hours)
        or reserved_node_hours <= 0.0
        or not math.isfinite(actual_node_hours)
        or actual_node_hours < 0.0
        or not math.isfinite(cumulative_node_hours)
        or cumulative_node_hours < 0.0
        or abs(
            reserved_node_hours - node_hours(nodes, requested_seconds)
        ) > LEDGER_NODE_HOUR_TOLERANCE
        or abs(actual_node_hours - node_hours(nodes, elapsed_seconds))
        > LEDGER_NODE_HOUR_TOLERANCE
    ):
        raise ValueError(f"{label} has invalid numeric fields")
    return (
        actual_node_hours,
        cumulative_node_hours,
        node_hours(nodes, elapsed_seconds),
    )


def validate_ledger_cumulative_fields(row: dict[str, object],
                                      prior_actual_node_hours: float,
                                      label: str) -> float:
    """Require one retained cumulative total to follow its preceding rows."""

    actual, cumulative, derived_actual = validate_ledger_numeric_fields(row, label)
    expected = prior_actual_node_hours + derived_actual
    if (
        not math.isfinite(expected)
        or abs(cumulative - expected) > LEDGER_CUMULATIVE_NODE_HOUR_TOLERANCE
    ):
        raise ValueError(f"{label} has inconsistent cumulative node-hours")
    if (
        expected > CURRENT_STAGE_I_RESERVED_NODE_HOURS
        or expected > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError(f"{label} exceeds an accounting ceiling")
    return prior_actual_node_hours + actual


def require_safe_segment(value: str) -> str:
    """Require a path- and Slurm-safe retained segment identifier."""

    if value == "analysis" or SEGMENT_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "--segment must contain 1-29 ASCII letters, digits, underscores, "
            "or hyphens, must begin with a letter or digit, and must not use "
            "the reserved analysis namespace"
        )
    return value


def require_numeric_job_id(value: str) -> str:
    """Require a top-level numeric Slurm allocation identifier."""

    if JOB_ID_PATTERN.fullmatch(value) is None:
        raise ValueError("--job-id must be a positive numeric Slurm job ID")
    return value


def require_root(root: Path, allow_local_root: bool) -> Path:
    """Require the declared project run root outside offline validation."""

    resolved = root.expanduser().resolve()
    if resolved != DEFAULT_ROOT.expanduser().resolve() and not allow_local_root:
        raise ValueError(
            f"Stage I root must be {DEFAULT_ROOT}; "
            "use --allow-local-root only for offline validation"
        )
    return resolved


def is_offline_local_root(root: Path, allow_local_root: bool) -> bool:
    """Return whether relaxed validation is permitted for a nonproject fixture."""

    return (
        allow_local_root
        and root.expanduser().resolve() != DEFAULT_ROOT.expanduser().resolve()
    )


def require_beneath_root(path: Path, root: Path, label: str,
                         allow_local_root: bool) -> None:
    """Keep submitted products in the project filesystem."""

    if is_offline_local_root(root, allow_local_root):
        return
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath {root}: {path}") from error


def require_current_epoch(manifest: dict[str, object], label: str) -> None:
    """Reject archival or untagged manifests in the current production path."""

    if manifest.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(
            f"{label} execution epoch is not {EXECUTION_EPOCH}: "
            f"{manifest.get('execution_epoch')!r}"
        )


def require_authorized_case(case_id: str) -> None:
    """Limit execution to the frozen mapped Stage I matrix."""

    if case_id not in AUTHORIZED_CASE_IDS:
        raise ValueError(
            f"{EXECUTION_EPOCH} Stage I is authorized only for mapped matrix cases "
            "R02-R17"
        )


def require_case_node_count(case_id: str, nodes: int) -> None:
    """Bind each mapped case to its reviewed Frontier allocation shape."""

    require_authorized_case(case_id)
    allowed = CASE_NODE_PROFILES[case_id]
    if nodes not in allowed:
        choices = "/".join(str(value) for value in sorted(allowed))
        raise ValueError(
            f"{case_id} canonical Stage I preparation requires --nodes={choices}"
        )


def retained_case_has_started(paths: dict[str, Path], case_id: str) -> bool:
    """Return whether one mapped case has retained any segment manifest."""

    require_authorized_case(case_id)
    return any(
        (paths["runs"] / case_id).glob("*/manifest/prepared_run.json")
    )


def require_r17_last(paths: dict[str, Path], case_id: str) -> None:
    """Require every lower-cost mapped case to complete before R17 starts."""

    if case_id != R17_CASE_ID:
        if retained_case_has_started(paths, R17_CASE_ID):
            raise ValueError(
                f"{R17_CASE_ID} has started; later {case_id} preparation is forbidden"
            )
        return
    incomplete = []
    for predecessor in R17_PREDECESSOR_CASE_IDS:
        lineage = accepted_case_lineage(paths, predecessor)
        if not lineage:
            incomplete.append(predecessor)
            continue
        final_time = float(lineage[-1]["scientific_inspection"]["final_time"])
        if (
            not math.isfinite(final_time)
            or final_time < REQUIRED_CASE_FINAL_TIME - 1.0e-10
        ):
            incomplete.append(predecessor)
    if incomplete:
        raise ValueError(
            "R17 must remain last; accepted t=10 lineages are incomplete for: "
            + ", ".join(incomplete)
        )


def require_prepare_case_policy(paths: dict[str, Path], case_id: str,
                                nodes: int, offline_local_root: bool) -> None:
    """Apply mapped allocation and ordering gates only to canonical production."""

    if offline_local_root:
        return
    require_case_node_count(case_id, nodes)
    require_r17_last(paths, case_id)


def require_retained_r17_policy(paths: dict[str, Path],
                                reservations: list[dict[str, object]]) -> None:
    """Keep the final high-resolution lane last across every retained transition."""

    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    r17_retained = retained_case_has_started(paths, R17_CASE_ID) or any(
        reservation.get("case_id") == R17_CASE_ID
        for reservation in reservations
        if isinstance(reservation, dict)
    )
    if not r17_retained:
        return
    require_r17_last(paths, R17_CASE_ID)
    lower_active = [
        str(reservation.get("case_id"))
        for reservation in active_reservations(reservations)
        if reservation.get("case_id") != R17_CASE_ID
    ]
    if lower_active:
        raise ValueError(
            f"{R17_CASE_ID} has started; active lower-resolution lanes are forbidden: "
            + ", ".join(sorted(lower_active))
        )


def parse_utc_timestamp(value: object, label: str) -> datetime:
    """Parse one retained timestamp and normalize it to UTC."""

    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} is not an ISO-8601 timestamp: {value!r}") from error
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=datetime.now().astimezone().tzinfo)
    return parsed.astimezone(timezone.utc)


def require_r17_exact_keys(value: object, keys: set[str] | frozenset[str],
                           label: str) -> dict[str, object]:
    """Require one exact R17 evidence object."""

    if not isinstance(value, dict) or set(value) != set(keys):
        found = sorted(value) if isinstance(value, dict) else type(value).__name__
        raise ValueError(
            f"{label} schema differs; expected {sorted(keys)}, found {found}"
        )
    return value


def require_r17_nonempty_string(value: object, label: str) -> str:
    """Require one trimmed nonempty R17 evidence string."""

    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a nonempty trimmed string")
    return value


def require_r17_sha256(value: object, label: str) -> str:
    """Require one lowercase SHA-256 binding."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256")
    return value


def require_r17_integer(value: object, label: str, minimum: int = 0) -> int:
    """Require one non-boolean bounded integer."""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def require_r17_decimal(value: object, label: str) -> Decimal:
    """Require one finite decimal string."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be a finite decimal string")
    try:
        retained = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"{label} must be a finite decimal string") from error
    if not retained.is_finite():
        raise ValueError(f"{label} must be a finite decimal string")
    return retained


def validate_r17_scoped_budget_projection(
    budget: dict[str, object],
    current_rows: list[dict[str, str]],
) -> None:
    """Require a final, internally coherent v2 scoped node-hour projection."""

    if budget.get("method") != R17_SCOPED_NODE_HOUR_METHOD:
        raise ValueError("latest promoted recost budget method is stale")

    basis_keys = {
        "job_id", "case_id", "segment", "actual_node_hours", "observed_cells",
        "observed_simulation_interval",
        "normalized_node_hours_per_cell_per_simulation_time",
    }

    def validate_basis(
        value: object, label: str,
    ) -> tuple[dict[str, object], Decimal, int]:
        basis = require_r17_exact_keys(value, basis_keys, label)
        job_id = require_r17_nonempty_string(basis["job_id"], f"{label} job ID")
        case_id = require_r17_nonempty_string(basis["case_id"], f"{label} case ID")
        segment = require_r17_nonempty_string(basis["segment"], f"{label} segment")
        segment_match = R17_RECOST_SEGMENT_PATTERN.fullmatch(segment)
        if (
            JOB_ID_PATTERN.fullmatch(job_id) is None
            or case_id not in AUTHORIZED_CASE_IDS
            or segment_match is None
        ):
            raise ValueError(f"{label} identity is invalid")
        try:
            segment_start = Decimal(segment_match.group("start").replace("p", "."))
            segment_target = Decimal(segment_match.group("target").replace("p", "."))
        except InvalidOperation as error:
            raise ValueError(f"{label} segment interval is invalid") from error
        actual = require_r17_decimal(
            basis["actual_node_hours"], f"{label} actual node-hours"
        )
        cells = require_r17_integer(
            basis["observed_cells"], f"{label} observed cells", 1
        )
        interval = require_r17_decimal(
            basis["observed_simulation_interval"], f"{label} observed interval"
        )
        rate = require_r17_decimal(
            basis["normalized_node_hours_per_cell_per_simulation_time"],
            f"{label} normalized rate",
        )
        matching_rows = [
            row for row in current_rows
            if isinstance(row, dict) and row.get("job_id") == job_id
        ]
        if (
            segment_start < 0
            or segment_target <= segment_start
            or actual <= 0
            or interval <= 0
            or rate <= 0
            or len(matching_rows) != 1
            or matching_rows[0].get("case_id") != case_id
            or matching_rows[0].get("segment") != segment
            or require_r17_decimal(
                matching_rows[0].get("actual_node_hours"),
                f"{label} ledger actual node-hours",
            )
            != actual
            or rate != actual / Decimal(cells) / interval
        ):
            raise ValueError(f"{label} is not a credible current-ledger measurement")
        return basis, rate, int(segment_match.group("index"))

    global_basis, global_rate, _ = validate_basis(
        budget.get("measurement_basis"),
        "latest promoted recost global measurement basis",
    )
    if global_basis["case_id"] == "R12":
        raise ValueError(
            "latest promoted recost global measurement basis must be non-R12"
        )

    breakdown = budget.get("case_breakdown")
    if not isinstance(breakdown, dict) or set(breakdown) != set(AUTHORIZED_CASE_IDS):
        raise ValueError("latest promoted recost scoped case breakdown differs")
    case_keys = {
        "matrix_full_case_node_hours_reference_only",
        "authenticated_progress_fraction", "remaining_simulation_time",
        "projected_cells", "projection_measurement_basis",
        "observed_rate_projected_remaining_node_hours",
        "authorized_profile_reserved_node_hours", "projected_remaining_node_hours",
    }
    total_authorized = Decimal("0")
    total_remaining = Decimal("0")
    required_time = Decimal(str(REQUIRED_CASE_FINAL_TIME))
    for case_id in sorted(AUTHORIZED_CASE_IDS):
        case = require_r17_exact_keys(
            breakdown[case_id], case_keys, f"latest promoted recost {case_id} budget"
        )
        if case_id == "R12":
            basis, rate, segment_index = validate_basis(
                case["projection_measurement_basis"],
                "latest promoted recost R12 measurement basis",
            )
            if (
                basis["case_id"] != "R12"
                or basis["job_id"] == R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID
                or segment_index < 1
            ):
                raise ValueError(
                    "latest promoted recost R12 basis is not a fresh final measurement"
                )
        else:
            if case["projection_measurement_basis"] != global_basis:
                raise ValueError(
                    f"latest promoted recost {case_id} basis does not reuse "
                    "the global basis"
                )
            rate = global_rate

        reference = require_r17_decimal(
            case["matrix_full_case_node_hours_reference_only"],
            f"latest promoted recost {case_id} matrix reference",
        )
        progress = require_r17_decimal(
            case["authenticated_progress_fraction"],
            f"latest promoted recost {case_id} progress",
        )
        remaining_time = require_r17_decimal(
            case["remaining_simulation_time"],
            f"latest promoted recost {case_id} remaining time",
        )
        projected_cells_text = case["projected_cells"]
        if (
            not isinstance(projected_cells_text, str)
            or re.fullmatch(r"[1-9][0-9]*", projected_cells_text) is None
        ):
            raise ValueError(
                f"latest promoted recost {case_id} projected cells are invalid"
            )
        projected_cells = Decimal(projected_cells_text)
        observed = require_r17_decimal(
            case["observed_rate_projected_remaining_node_hours"],
            f"latest promoted recost {case_id} observed-rate projection",
        )
        authorized = require_r17_decimal(
            case["authorized_profile_reserved_node_hours"],
            f"latest promoted recost {case_id} authorized reserve",
        )
        projected = require_r17_decimal(
            case["projected_remaining_node_hours"],
            f"latest promoted recost {case_id} projected remaining node-hours",
        )
        if (
            reference < 0
            or progress < 0
            or progress > 1
            or remaining_time != required_time * (Decimal("1") - progress)
            or observed != rate * projected_cells * remaining_time
            or authorized < 0
            or projected != max(observed, authorized)
            or (
                case_id in R17_PREDECESSOR_CASE_IDS
                and (
                    abs(progress - Decimal("1")) > Decimal("1e-10")
                    or abs(remaining_time) > Decimal("1e-9")
                )
            )
        ):
            raise ValueError(
                f"latest promoted recost {case_id} scoped projection is invalid"
            )
        total_authorized += authorized
        total_remaining += projected

    if (
        total_authorized
        != require_r17_decimal(
            budget.get("authorized_wave_reserved_node_hours"),
            "latest promoted recost scoped authorized reserve",
        )
        or total_remaining
        != require_r17_decimal(
            budget.get("computed_remaining_stage_i_node_hours"),
            "latest promoted recost scoped remaining node-hours",
        )
    ):
        raise ValueError("latest promoted recost scoped projection totals differ")


def require_r17_utc(value: object, label: str) -> datetime:
    """Require one timezone-aware UTC timestamp."""

    if not isinstance(value, str) or not value:
        raise ValueError(f"{label} must be a UTC ISO-8601 timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be a UTC ISO-8601 timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must include a UTC timezone")
    return parsed.astimezone(timezone.utc)


def require_r17_scheduler_utc(value: object, label: str) -> datetime:
    """Parse retained sacct time, treating its timezone-free form as UTC."""

    retained = require_r17_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(retained.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    if parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")
    return parsed.astimezone(timezone.utc)


def require_symlink_free_path(path: Path, label: str) -> None:
    """Require every existing component of one absolute path to be non-symlink."""

    absolute = path.absolute()
    if (
        not absolute.is_absolute()
        or absolute != Path(os.path.normpath(str(absolute)))
        or ".." in absolute.parts
    ):
        raise ValueError(f"{label} path must be absolute and normalized: {path}")
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current /= part
        try:
            profile = os.lstat(current)
        except OSError as error:
            raise ValueError(f"{label} path component is unavailable: {current}") from error
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symbolic link: {current}")


def require_owner_symlink_free_path(path: Path, label: str) -> None:
    """Require every path component to be root- or controller-owned and non-symlink."""

    absolute = path.absolute()
    allowed_owners = {0, os.geteuid()}
    current = Path(absolute.anchor)
    for part in absolute.parts[1:]:
        current /= part
        try:
            profile = os.lstat(current)
        except OSError as error:
            raise ValueError(f"{label} path component is unavailable: {current}") from error
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symbolic link: {current}")
        if profile.st_uid not in allowed_owners:
            raise ValueError(
                f"{label} path is not owner-controlled; component owner differs: {current}"
            )


def canonical_authority_identity(profile: os.stat_result) -> tuple[object, ...]:
    """Return the exact public-namespace identity retained by the root lock."""

    return tuple(
        getattr(profile, field)
        for field in (
            "st_dev", "st_ino", "st_mode", "st_uid", "st_gid",
        )
    )


def require_canonical_root_authority(
    root: Path,
    expected: tuple[tuple[str, object], ...] | None = None,
) -> tuple[tuple[str, object], ...]:
    """Bind the canonical public namespace through its trusted project boundary.

    A root-owned sticky system temporary directory is accepted for offline
    fixtures.  Production may additionally traverse only Frontier's exact
    ``/lustre/orion/ast207/proj-shared`` 2770 root:31114 boundary.  The
    canonical root and every controller-owned component remain non-writable
    by group or world.
    """

    absolute = root.expanduser().absolute()
    allowed_owners = {0, os.geteuid()}
    current = Path(absolute.anchor)
    components = [current]
    authority = []
    for part in absolute.parts[1:]:
        current /= part
        components.append(current)
    for current in components:
        try:
            profile = os.lstat(current)
        except OSError as error:
            raise ValueError(
                f"canonical Stage I root component is unavailable: {current}"
            ) from error
        mode = stat.S_IMODE(profile.st_mode)
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(
                f"canonical Stage I root contains a symbolic link: {current}"
            )
        if not stat.S_ISDIR(profile.st_mode) or profile.st_uid not in allowed_owners:
            raise ValueError(
                f"canonical Stage I root is not owner-controlled: {current}"
            )
        root_owned_sticky_ancestor = (
            current != absolute
            and profile.st_uid == 0
            and bool(mode & stat.S_ISVTX)
        )
        trusted_project_boundary = (
            absolute == DEFAULT_ROOT.expanduser().absolute()
            and current == TRUSTED_GROUP_WRITABLE_PROJECT_BOUNDARY
            and mode == TRUSTED_GROUP_WRITABLE_PROJECT_MODE
            and profile.st_uid == TRUSTED_GROUP_WRITABLE_PROJECT_UID
            and profile.st_gid == TRUSTED_GROUP_WRITABLE_PROJECT_GID
        )
        if mode & 0o022 and not (
            root_owned_sticky_ancestor or trusted_project_boundary
        ):
            raise ValueError(
                f"canonical Stage I root or ancestor is group/world-writable: {current}"
            )
        authority.append((str(current), canonical_authority_identity(profile)))
    retained = tuple(authority)
    if expected is not None and retained != expected:
        raise ValueError("canonical Stage I public namespace authority changed")
    return retained


def require_trusted_system_executable_path(path: Path, label: str) -> None:
    """Require a root-owned executable beneath immutable root-owned directories."""

    absolute = path.absolute()
    if (
        not absolute.is_absolute()
        or absolute != Path(os.path.normpath(str(absolute)))
        or ".." in absolute.parts
    ):
        raise ValueError(f"{label} path must be absolute and normalized: {path}")
    current = Path(absolute.anchor)
    for index, part in enumerate(absolute.parts[1:], start=1):
        current /= part
        try:
            profile = os.lstat(current)
        except OSError as error:
            raise ValueError(f"{label} path component is unavailable: {current}") from error
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symbolic link: {current}")
        if profile.st_uid != 0 or stat.S_IMODE(profile.st_mode) & 0o022:
            raise ValueError(
                f"{label} path is not immutable and root-owned: {current}"
            )
        leaf = index == len(absolute.parts) - 1
        if (not leaf and not stat.S_ISDIR(profile.st_mode)) or (
            leaf
            and (
                not stat.S_ISREG(profile.st_mode)
                or not stat.S_IMODE(profile.st_mode) & 0o111
            )
        ):
            raise ValueError(f"{label} path profile differs: {current}")


def trusted_directory_entries(path: Path, label: str, *,
                              allow_absent: bool = False) -> list[str]:
    """List one owner-controlled, symlink-free directory through its descriptor."""

    absolute = path.absolute()
    require_owner_symlink_free_path(absolute.parent, label)
    try:
        absolute.lstat()
    except FileNotFoundError:
        if allow_absent:
            return []
        raise ValueError(f"{label} is unavailable: {absolute}") from None
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {absolute}") from error
    require_owner_symlink_free_path(absolute, label)
    try:
        descriptor = os.open(
            absolute,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {absolute}") from error
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(before.st_mode)
            or before.st_uid != os.geteuid()
            or stat.S_IMODE(before.st_mode) & 0o022
        ):
            raise ValueError(
                f"{label} must be an owner-controlled non-writable directory: "
                f"{absolute}"
            )
        entries = sorted(os.listdir(descriptor))
        after = os.fstat(descriptor)
    except OSError as error:
        raise ValueError(f"{label} cannot be listed: {absolute}") from error
    finally:
        os.close(descriptor)
    require_owner_symlink_free_path(absolute, label)
    try:
        named = absolute.lstat()
    except OSError as error:
        raise ValueError(f"{label} pathname changed while listing: {absolute}") from error
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_gid", "st_nlink",
    )
    if (
        any(getattr(before, field) != getattr(after, field) for field in stable_fields)
        or (named.st_dev, named.st_ino) != (before.st_dev, before.st_ino)
        or not stat.S_ISDIR(named.st_mode)
        or directory_security_identity(named)
        != directory_security_identity(before)
        or named.st_uid != os.geteuid()
        or stat.S_IMODE(named.st_mode) & 0o022
    ):
        raise ValueError(f"{label} pathname changed while listing: {absolute}")
    return entries


def require_empty_transaction_directory(path: Path, label: str) -> None:
    """Require an owner-controlled transaction directory to be absent or empty."""

    if trusted_directory_entries(path, label, allow_absent=True):
        raise ValueError(f"{label} requires recovery before controller use")


def read_r17_evidence_bytes(path: Path, label: str, *,
                            mode: int, links: int = 1,
                            owner_controlled: bool = False,
                            symlink_free: bool = False) -> bytes:
    """Read one exact immutable R17 evidence file."""

    if symlink_free:
        if owner_controlled:
            require_owner_symlink_free_path(path, label)
        else:
            require_symlink_free_path(path, label)
    try:
        descriptor = os.open(
            path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        )
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {path}") from error
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat.S_IMODE(before.st_mode) != mode
            or before.st_nlink != links
            or (owner_controlled and before.st_uid != os.geteuid())
        ):
            ownership = "an owner-controlled" if owner_controlled else "a"
            raise ValueError(
                f"{label} must be {ownership} regular {mode:04o} "
                f"single-link file: {path}"
            )
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        after = os.fstat(descriptor)
    except OSError as error:
        raise ValueError(f"{label} cannot be read: {path}") from error
    finally:
        os.close(descriptor)
    if symlink_free:
        if owner_controlled:
            require_owner_symlink_free_path(path, label)
        else:
            require_symlink_free_path(path, label)
    try:
        named = path.lstat()
    except OSError as error:
        raise ValueError(f"{label} pathname changed while reading: {path}") from error
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if (
        any(getattr(before, field) != getattr(after, field) for field in stable_fields)
        or (named.st_dev, named.st_ino) != (before.st_dev, before.st_ino)
        or stat.S_IMODE(named.st_mode) != mode
        or named.st_nlink != links
        or (owner_controlled and named.st_uid != os.geteuid())
    ):
        raise ValueError(f"{label} pathname changed while reading: {path}")
    return b"".join(blocks)


def read_f116_authority_bytes(path: Path, label: str, *,
                              mode: int, links: int = 1) -> bytes:
    """Read one owner-controlled F116 authority or source-catalog file."""

    return read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        links=links,
        owner_controlled=True,
        symlink_free=True,
    )


def read_f116_authority_json(path: Path, label: str, *,
                             mode: int, links: int = 1
                             ) -> tuple[dict[str, object], str]:
    """Read one owner-controlled F116 JSON authority publication."""

    retained = read_f116_authority_bytes(path, label, mode=mode, links=links)
    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value, hashlib.sha256(retained).hexdigest()


def read_r17_evidence_json(path: Path, label: str, *,
                           mode: int, links: int = 1
                           ) -> tuple[dict[str, object], str]:
    """Read one exact immutable R17 JSON publication."""

    retained = read_r17_evidence_bytes(path, label, mode=mode, links=links)
    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value, hashlib.sha256(retained).hexdigest()


def read_controller_publication_json(path: Path, label: str, *,
                                     mode: int, links: int = 1
                                     ) -> tuple[dict[str, object], str]:
    """Read one owner-controlled, symlink-free controller publication."""

    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        links=links,
        owner_controlled=True,
        symlink_free=True,
    )
    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object: {path}")
    return value, hashlib.sha256(retained).hexdigest()


def require_r17_root_binding(root: Path, value: object, label: str, *,
                             mode: int) -> tuple[Path, str, bytes]:
    """Authenticate one exact root-relative R17 input binding."""

    binding = require_r17_exact_keys(value, {"path", "sha256"}, label)
    relative_text = require_r17_nonempty_string(binding["path"], f"{label} path")
    relative = Path(relative_text)
    if relative.is_absolute() or ".." in relative.parts or relative.as_posix() != relative_text:
        raise ValueError(f"{label} path must be normalized and root-relative")
    path = root / relative
    try:
        if path.resolve().relative_to(root.resolve()) != relative:
            raise ValueError(f"{label} path escapes or aliases the Stage I root")
    except ValueError as error:
        raise ValueError(f"{label} path escapes or aliases the Stage I root") from error
    expected = require_r17_sha256(binding["sha256"], f"{label} SHA-256")
    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        owner_controlled=True,
        symlink_free=True,
    )
    if hashlib.sha256(retained).hexdigest() != expected:
        raise ValueError(f"{label} checksum differs")
    return path, expected, retained


def require_r17_declared_publication(value: object, path: Path, digest: str,
                                     label: str) -> None:
    """Require one exact immutable publication binding."""

    binding = require_r17_exact_keys(
        value, {"path", "sha256", "mode", "links"}, label
    )
    if binding != {
        "path": str(path),
        "sha256": digest,
        "mode": "0444",
        "links": 1,
    }:
        raise ValueError(f"{label} differs from exact retained publication")


def require_source_authority_revision(value: object, label: str) -> str:
    """Require one full lowercase Git revision in source-authority evidence."""

    revision = require_r17_nonempty_string(value, label)
    if GIT_REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError(f"{label} must be a lowercase full Git revision")
    return revision


def caller_independent_child_environment(*, strip_scheduler: bool) -> dict[str, str]:
    """Strip caller-controlled interpreter, loader, Git, and scheduler routing."""

    prefixes = CALLER_CHILD_ENVIRONMENT_PREFIX_DENY
    if not strip_scheduler:
        prefixes = tuple(
            prefix for prefix in prefixes if prefix not in {"SBATCH_", "SLURM_"}
        )
    return {
        key: value
        for key, value in os.environ.items()
        if key not in CALLER_CHILD_ENVIRONMENT_EXACT_DENY
        and not key.startswith(prefixes)
    }


def hardened_git_environment() -> dict[str, str]:
    """Return a config-isolated Git environment independent of the caller."""

    environment = {
        "GIT_CONFIG_COUNT": "6",
        "GIT_CONFIG_KEY_0": "core.hooksPath",
        "GIT_CONFIG_VALUE_0": "/dev/null",
        "GIT_CONFIG_KEY_1": "core.fsmonitor",
        "GIT_CONFIG_VALUE_1": "false",
        "GIT_CONFIG_KEY_2": "core.attributesFile",
        "GIT_CONFIG_VALUE_2": "/dev/null",
        "GIT_CONFIG_KEY_3": "diff.external",
        "GIT_CONFIG_VALUE_3": "",
        "GIT_CONFIG_KEY_4": "core.pager",
        "GIT_CONFIG_VALUE_4": "cat",
        "GIT_CONFIG_KEY_5": "protocol.file.allow",
        "GIT_CONFIG_VALUE_5": "always",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_SYSTEM": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "GIT_OPTIONAL_LOCKS": "0",
        "GIT_TERMINAL_PROMPT": "0",
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "XDG_CONFIG_HOME": "/nonexistent",
    }
    return environment


def hardened_python_environment() -> dict[str, str]:
    """Return an isolated environment for any controller Python child."""

    return {
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONNOUSERSITE": "1",
        "PYTHONSAFEPATH": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def open_authenticated_system_python_descriptor() -> int:
    """Open and authenticate the fixed reviewed controller interpreter."""

    require_trusted_system_executable_path(
        SYSTEM_PYTHON, "controller Python interpreter"
    )
    try:
        descriptor = os.open(
            SYSTEM_PYTHON, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        )
    except OSError as error:
        raise ValueError("controller Python interpreter is unavailable") from error
    opened = os.fstat(descriptor)
    named = SYSTEM_PYTHON.lstat()
    if (
        not stat.S_ISREG(opened.st_mode)
        or opened.st_uid != 0
        or stat.S_IMODE(opened.st_mode) & 0o022
        or not stat.S_IMODE(opened.st_mode) & 0o111
        or (named.st_dev, named.st_ino) != (opened.st_dev, opened.st_ino)
        or named.st_uid != 0
        or stat.S_IMODE(named.st_mode) & 0o022
    ):
        os.close(descriptor)
        raise ValueError(
            "controller Python interpreter is not an authenticated "
            f"root-owned executable: {SYSTEM_PYTHON}"
        )
    return descriptor


def authenticated_python_binary() -> str:
    """Return the fixed authenticated controller interpreter path."""

    descriptor = open_authenticated_system_python_descriptor()
    os.close(descriptor)
    return str(SYSTEM_PYTHON)


def require_controller_source_profile(profile: os.stat_result, label: str) -> None:
    """Require one exact owner-controlled controller source inode."""

    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid not in {0, os.geteuid()}
        or stat.S_IMODE(profile.st_mode) & 0o022
        or profile.st_nlink != 1
    ):
        raise ValueError(f"{label} is not an owner-controlled single-link file")


def sha256_descriptor(descriptor: int) -> str:
    """Return a digest of exact retained descriptor bytes."""

    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    while block := os.read(descriptor, 1024 * 1024):
        digest.update(block)
    os.lseek(descriptor, 0, os.SEEK_SET)
    return digest.hexdigest()


def initial_controller_source_path() -> Path:
    """Return the named controller source for initial descriptor reexecution."""

    source = Path(__file__).absolute()
    if source.is_symlink():
        raise ValueError("controller source must not be a symbolic link")
    return source.resolve(strict=True)


def retained_controller_source_path() -> Path:
    """Return the named controller source retained across authenticated reexec."""

    if SELF_SOURCE_ENV in os.environ:
        return inherited_reexec_path(SELF_SOURCE_ENV, "source path")
    return initial_controller_source_path()


def require_reexec_source_relationship(source: Path, repository: Path) -> None:
    """Bind private reexec metadata to the exact controller repository path."""

    if source != repository / PRODUCTION_UTILITY_RELATIVE:
        raise ValueError("authenticated controller source/repository relationship differs")


def inherited_reexec_path(name: str, label: str) -> Path:
    """Read one private reexec path after authenticating the self descriptor."""

    retained = os.environ.get(name)
    if retained is None:
        raise ValueError(f"authenticated controller reexecution lacks {label}")
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(retained))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(
            f"authenticated controller {label} is not normalized and absolute"
        )
    require_owner_symlink_free_path(path, f"authenticated controller {label}")
    return path


def controller_reexec_environment(
    descriptor: int, python_descriptor: int, source: Path, repository: Path
) -> dict[str, str]:
    """Return the complete caller-independent controller reexec environment."""

    environment = hardened_python_environment()
    environment.update({
        SELF_DESCRIPTOR_ENV: str(descriptor),
        PYTHON_DESCRIPTOR_ENV: str(python_descriptor),
        SELF_SOURCE_ENV: str(source),
        REPOSITORY_ROOT_ENV: str(repository),
    })
    return environment


def require_authenticated_reexec_runtime() -> None:
    """Require the inherited controller process to use isolated fixed Python."""

    if os.environ.get("HOME") != "/nonexistent":
        raise ValueError("authenticated controller reexecution HOME is not sanitized")
    flags = sys.flags
    if (
        flags.isolated != 1
        or flags.ignore_environment != 1
        or flags.no_user_site != 1
        or flags.no_site != 1
        or flags.dont_write_bytecode != 1
    ):
        raise ValueError("authenticated controller reexecution is not isolated")
    retained = os.environ.get(PYTHON_DESCRIPTOR_ENV)
    if retained is None or re.fullmatch(r"[0-9]+", retained) is None:
        raise ValueError("authenticated controller Python descriptor marker is invalid")
    expected = os.fstat(int(retained))
    require_controller_source_profile(expected, "authenticated controller Python")
    if expected.st_uid != 0 or not stat.S_IMODE(expected.st_mode) & 0o111:
        raise ValueError("authenticated controller Python descriptor is not root-owned")
    running = os.stat("/proc/self/exe")
    if profile_identity(expected) != profile_identity(running):
        raise ValueError("authenticated controller Python descriptor is not this process")
    fixed = open_authenticated_system_python_descriptor()
    try:
        if profile_identity(os.fstat(fixed)) != profile_identity(expected):
            raise ValueError("authenticated controller interpreter differs")
    finally:
        os.close(fixed)


def require_initial_controller_runtime() -> None:
    """Reject every initial startup before fixed isolated Python is authenticated."""

    flags = sys.flags
    if (
        flags.isolated != 1
        or flags.ignore_environment != 1
        or flags.no_user_site != 1
        or flags.no_site != 1
        or flags.dont_write_bytecode != 1
    ):
        raise ValueError(
            "controller invocation requires "
            f"{SYSTEM_PYTHON} -I -S -B {PRODUCTION_UTILITY_RELATIVE}"
        )
    expected = open_authenticated_system_python_descriptor()
    try:
        running = os.stat("/proc/self/exe")
        if profile_identity(os.fstat(expected)) != profile_identity(running):
            raise ValueError(
                "controller invocation is not using the fixed "
                "authenticated controller interpreter"
            )
    finally:
        os.close(expected)


def reexec_authenticated_controller(
    descriptor: int, source: Path, repository: Path, argv: list[str]
) -> None:
    """Reexecute retained controller bytes through fixed isolated Python."""

    python_descriptor = open_authenticated_system_python_descriptor()
    os.set_inheritable(descriptor, True)
    os.set_inheritable(python_descriptor, True)
    os.execve(
        f"/proc/self/fd/{python_descriptor}",
        [
            str(SYSTEM_PYTHON), "-I", "-S", "-B",
            f"/proc/self/fd/{descriptor}", *argv,
        ],
        controller_reexec_environment(
            descriptor, python_descriptor, source, repository
        ),
    )
    raise AssertionError("authenticated controller reexecution unexpectedly returned")


def authenticate_controller_runtime(argv: list[str]) -> None:
    """Reexecute and authenticate exact controller bytes before any action."""

    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        if re.fullmatch(r"[0-9]+", inherited) is None:
            raise ValueError("controller descriptor marker is invalid")
        descriptor = int(inherited)
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError(
                "controller descriptor marker is not attached to this execution"
            )
        opened = os.fstat(descriptor)
        require_controller_source_profile(opened, "authenticated controller descriptor")
        require_authenticated_reexec_runtime()
        source = inherited_reexec_path(SELF_SOURCE_ENV, "source path")
        repository = inherited_reexec_path(REPOSITORY_ROOT_ENV, "repository root")
        require_reexec_source_relationship(source, repository)
        named = os.open(source, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        try:
            require_controller_source_profile(
                os.fstat(named), "retained controller source"
            )
            if (
                profile_identity(os.fstat(named)) != profile_identity(opened)
                or sha256_descriptor(named) != sha256_descriptor(descriptor)
            ):
                raise ValueError(
                    "retained controller source differs after reexecution"
                )
        finally:
            os.close(named)
        return

    private = (
        PYTHON_DESCRIPTOR_ENV, SELF_SOURCE_ENV, REPOSITORY_ROOT_ENV
    )
    if any(name in os.environ for name in private):
        raise ValueError(
            "private controller reexecution path is forbidden without an "
            "authenticated descriptor"
        )
    require_initial_controller_runtime()
    source = initial_controller_source_path()
    repository = source.parents[2].resolve(strict=True)
    require_reexec_source_relationship(source, repository)
    require_owner_symlink_free_path(source, "controller source")
    descriptor = os.open(source, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        require_controller_source_profile(
            os.fstat(descriptor), "retained controller source"
        )
        reexec_authenticated_controller(descriptor, source, repository, argv)
    finally:
        os.close(descriptor)


def authenticated_git_binary() -> str:
    """Return the exact root-owned absolute Git executable."""

    if not GIT.is_absolute() or ".." in GIT.parts:
        raise ValueError(f"Git binary path is not absolute and normalized: {GIT}")
    require_trusted_system_executable_path(GIT, "Git binary")
    try:
        descriptor = os.open(GIT, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ValueError(f"Git binary is unavailable: {GIT}") from error
    try:
        opened = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    try:
        named = GIT.lstat()
    except OSError as error:
        raise ValueError(f"Git binary pathname changed: {GIT}") from error
    if (
        not stat.S_ISREG(opened.st_mode)
        or opened.st_uid != 0
        or stat.S_IMODE(opened.st_mode) & 0o022
        or not stat.S_IMODE(opened.st_mode) & 0o111
        or (named.st_dev, named.st_ino) != (opened.st_dev, opened.st_ino)
        or named.st_uid != 0
        or stat.S_IMODE(named.st_mode) & 0o022
    ):
        raise ValueError(f"Git binary is not an authenticated root-owned executable: {GIT}")
    return str(GIT)


def controller_git_worktree_prefix(repository: Path) -> list[str]:
    """Bind Git worktree queries to explicit absolute repository paths."""

    worktree = repository.expanduser().absolute()
    git_dir = worktree / ".git"
    if (
        not worktree.is_absolute()
        or worktree != Path(os.path.normpath(str(worktree)))
        or not git_dir.is_dir()
    ):
        raise ValueError(f"Git worktree is unavailable or noncanonical: {repository}")
    return ["--git-dir", str(git_dir), "--work-tree", str(worktree)]


def source_authority_git(arguments: list[str], *,
                         text: bool = False) -> subprocess.CompletedProcess:
    """Run one F116 repository query without inherited Git state or replacements."""

    return subprocess.run(
        [
            authenticated_git_binary(), "--no-replace-objects",
            *controller_git_worktree_prefix(ROOT_DIR),
            *arguments,
        ],
        check=False,
        capture_output=True,
        text=text,
        env=hardened_git_environment(),
    )


def source_authority_committed_sha256(revision: str, relative: str,
                                      label: str) -> str:
    """Return the digest of one repository file at an authenticated revision."""

    retained = source_authority_git(["show", f"{revision}:{relative}"])
    if retained.returncode != 0:
        raise ValueError(f"{label} committed bytes are unavailable")
    return hashlib.sha256(retained.stdout).hexdigest()


def source_authority_committed_mode(revision: str, relative: str,
                                    label: str) -> str:
    """Return the published octal mode of one repository file."""

    retained = source_authority_git(
        ["ls-tree", revision, "--", relative], text=True
    )
    if retained.returncode != 0:
        raise ValueError(f"{label} committed mode is unavailable")
    match = re.fullmatch(r"(100644|100755) blob [0-9a-f]{40}\t(.+)\n?", retained.stdout)
    if match is None or match.group(2) != relative:
        raise ValueError(f"{label} committed mode differs")
    return "0755" if match.group(1) == "100755" else "0644"


def source_authority_commit_subject(revision: str) -> str:
    """Return the exact UTF-8 subject of one authenticated source revision."""

    retained = source_authority_git(["show", "-s", "--format=%s", revision])
    if retained.returncode != 0:
        raise ValueError("F116 final HEAD subject is unavailable")
    try:
        subject = retained.stdout.decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as error:
        raise ValueError("F116 final HEAD subject is not UTF-8") from error
    if not subject:
        raise ValueError("F116 final HEAD subject is empty")
    return subject


def source_authority_revision_is_ancestor(revision: str,
                                          promoted_revision: str) -> bool:
    """Return whether one declared source revision is in promoted HEAD history."""

    retained = source_authority_git(
        ["merge-base", "--is-ancestor", revision, promoted_revision]
    )
    if retained.returncode not in (0, 1):
        raise ValueError("F116 final source bundle ancestry cannot be authenticated")
    return retained.returncode == 0


def require_historical_f115_source_authority(
    paths: dict[str, Path],
) -> dict[str, object]:
    """Require the exact published F115 chain and immutable recorded R03/s02 manifest."""

    root = paths["root"]
    expected = {
        "evidence": (
            R03_F115_SOURCE_AUTHORITY_RELATIVE,
            R03_F115_SOURCE_AUTHORITY_SHA256,
        ),
        "publication_audit": (
            Path(f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json"),
            R03_F115_PUBLICATION_AUDIT_SHA256,
        ),
        "provenance_review": (
            Path(f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json"),
            R03_F115_PROVENANCE_REVIEW_SHA256,
        ),
        "plasma_review": (
            Path(f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json"),
            R03_F115_PLASMA_REVIEW_SHA256,
        ),
    }
    loaded: dict[str, dict[str, object]] = {}
    for key, (relative, digest) in expected.items():
        payload = read_f116_authority_bytes(
            root / relative, f"historical F115 {key}", mode=0o444
        )
        if hashlib.sha256(payload).hexdigest() != digest:
            raise ValueError(f"historical F115 {key} checksum differs")
        try:
            value = json.loads(payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError(f"historical F115 {key} is not valid JSON") from error
        if not isinstance(value, dict):
            raise ValueError(f"historical F115 {key} must be an object")
        loaded[key] = value
    if (
        loaded["evidence"].get("schema_version") != 1
        or loaded["evidence"].get("record_type")
        != "stage-i-source-bundle-recovery-supersession-evidence"
        or loaded["evidence"].get("checkpoint") != "F-115"
        or loaded["evidence"].get("execution_epoch") != EXECUTION_EPOCH
        or loaded["publication_audit"].get("schema_version") != 1
        or loaded["publication_audit"].get("record_type")
        != "stage-i-source-bundle-recovery-supersession-publication-audit"
        or loaded["publication_audit"].get("checkpoint") != "F-115"
        or loaded["publication_audit"].get("execution_epoch") != EXECUTION_EPOCH
        or any(
            loaded[key].get("schema_version") != 1
            or loaded[key].get("record_type")
            != "stage-i-source-bundle-recovery-supersession-independent-review"
            or loaded[key].get("checkpoint") != "F-115"
            or loaded[key].get("execution_epoch") != EXECUTION_EPOCH
            for key in ("provenance_review", "plasma_review")
        )
    ):
        raise ValueError("historical F115 publication chain identity differs")

    manifest_path = root / R03_F115_MANIFEST_RELATIVE
    manifest_payload = read_f116_authority_bytes(
        manifest_path, "historical F115 R03/s02 manifest", mode=0o644
    )
    if hashlib.sha256(manifest_payload).hexdigest() != R03_F115_MANIFEST_SHA256:
        raise ValueError("historical F115 R03/s02 manifest checksum differs")
    try:
        manifest = json.loads(manifest_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("historical F115 R03/s02 manifest is not valid JSON") from error
    run = manifest.get("run") if isinstance(manifest, dict) else None
    command = manifest.get("command") if isinstance(manifest, dict) else None
    accounting = manifest.get("accounting") if isinstance(manifest, dict) else None
    inspection = (
        manifest.get("scientific_inspection") if isinstance(manifest, dict) else None
    )
    utility = command.get("production_utility") if isinstance(command, dict) else None
    bundle = command.get("source_bundle") if isinstance(command, dict) else None
    try:
        final_time = float(inspection["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("historical F115 R03/s02 manifest endpoint is invalid") from error
    if (
        manifest.get("schema_version") != 3
        or manifest.get("execution_epoch") != EXECUTION_EPOCH
        or manifest.get("state") != "recorded"
        or manifest.get("job_id") != R03_F115_JOB_ID
        or manifest.get("project_root") != str(root)
        or not isinstance(run, dict)
        or run.get("case_id") != "R03"
        or run.get("segment") != R03_F115_SEGMENT
        or not isinstance(utility, dict)
        or utility.get("committed") is not True
        or utility.get("revision") != R03_F115_CONTROLLER_REVISION
        or utility.get("sha256") != R03_F115_CONTROLLER_SHA256
        or not isinstance(bundle, dict)
        or bundle.get("path") != str(root / R03_F115_SOURCE_BUNDLE_RELATIVE)
        or bundle.get("sha256") != R03_F115_SOURCE_BUNDLE_SHA256
        or bundle.get("verified_revisions")
        != [QUALIFIED_SOURCE_REVISION, R03_F115_CONTROLLER_REVISION]
        or not isinstance(accounting, dict)
        or accounting.get("result") != "accepted"
        or accounting.get("state") != "COMPLETED"
        or not isinstance(inspection, dict)
        or not math.isfinite(final_time)
        or abs(final_time - 0.5) > 1.0e-12
        or "current_source_authority" in command
    ):
        raise ValueError("historical F115 R03/s02 manifest identity or result differs")
    return {
        "scope": "historical-R03-s02-only",
        "evidence_sha256": R03_F115_SOURCE_AUTHORITY_SHA256,
        "publication_audit_sha256": R03_F115_PUBLICATION_AUDIT_SHA256,
        "provenance_review_sha256": R03_F115_PROVENANCE_REVIEW_SHA256,
        "plasma_review_sha256": R03_F115_PLASMA_REVIEW_SHA256,
    }


def source_authority_bundle_advertised_tip(
    path: Path, expected_sha256: str | None = None,
) -> tuple[str, str]:
    """Authenticate one stable complete-history bundle and return its sole tip."""

    require_symlink_free_path(path, "F116 source bundle")
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ValueError(f"F116 final source bundle is unavailable: {path}") from error
    try:
        before = os.fstat(descriptor)
        digest = hashlib.sha256()
        header = bytearray()
        header_too_large = False
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            digest.update(block)
            if b"\n\n" not in header and not header_too_large:
                header.extend(block)
                header_too_large = len(header) > 1024 * 1024
        after = os.fstat(descriptor)
    except OSError as error:
        raise ValueError(f"F116 final source bundle cannot be read: {path}") from error
    finally:
        os.close(descriptor)
    try:
        named = path.lstat()
    except OSError as error:
        raise ValueError(f"F116 final source bundle pathname changed: {path}") from error
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if (
        not stat.S_ISREG(before.st_mode)
        or stat.S_IMODE(before.st_mode) != 0o644
        or before.st_uid != os.geteuid()
        or before.st_nlink != 1
        or any(getattr(before, field) != getattr(after, field) for field in stable_fields)
        or (named.st_dev, named.st_ino) != (before.st_dev, before.st_ino)
        or stat.S_IMODE(named.st_mode) != 0o644
        or named.st_uid != os.geteuid()
        or named.st_nlink != 1
    ):
        raise ValueError(
            "F116 source bundle changed or is not an owner-controlled "
            "0644 single-link file"
        )
    require_symlink_free_path(path, "F116 source bundle")
    if expected_sha256 is not None and digest.hexdigest() != require_r17_sha256(
        expected_sha256, "F116 final source bundle expected SHA-256"
    ):
        raise ValueError("F116 final source bundle checksum differs")
    if header_too_large or b"\n\n" not in header:
        raise ValueError("F116 final source bundle header is missing or oversized")
    try:
        lines = bytes(header).split(b"\n\n", 1)[0].decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("F116 final source bundle header is not UTF-8") from error
    if not lines or lines[0] not in {"# v2 git bundle", "# v3 git bundle"}:
        raise ValueError("F116 final source bundle header is invalid")
    advertised = []
    for line in lines[1:]:
        if line.startswith("@"):
            continue
        if line.startswith("-"):
            raise ValueError("F116 final source bundle is not complete-history")
        match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
        if match is None:
            raise ValueError("F116 final source bundle advertised reference is malformed")
        advertised.append((match.group(1), match.group(2)))
    if len(advertised) != 1:
        raise ValueError("F116 final source bundle must advertise exactly one tip")
    return advertised[0]


def validate_f116_source_archive_checksum_ledger(
    root: Path, ledger_lines: list[str],
) -> dict[str, str]:
    """Authenticate every file named by the active F116 checksum ledger."""

    if not ledger_lines:
        raise ValueError("F116 source-archive checksum ledger is empty")
    ledger: dict[str, str] = {}
    for line in ledger_lines:
        match = re.fullmatch(r"([0-9a-f]{64})  ([A-Za-z0-9_.-]+)", line)
        if match is None or match.group(2) in ledger:
            raise ValueError("F116 source-archive checksum ledger is malformed")
        expected, name = match.groups()
        path = root / "source-archives" / name
        payload = read_f116_authority_bytes(
            path, f"F116 active source archive {name}", mode=0o644
        )
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f"F116 active source archive checksum differs: {name}")
        ledger[name] = expected
    return ledger


def require_authenticated_f118_source_authority_staging(
    source_transactions: Path,
) -> None:
    """Classify post-commit source-authority containers as non-authoritative debris.

    The immutable public F118 authority is authenticated before this function
    runs.  Retained transaction contents therefore have no authority and are
    never opened or interpreted here; only the transaction-root and its
    direct-child container profiles remain security boundaries.
    """

    source_transactions = source_transactions.absolute()
    label = "F118 source-authority recovery-debris root"
    require_owner_symlink_free_path(source_transactions.parent, label)
    try:
        source_transactions.lstat()
    except FileNotFoundError:
        return
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {source_transactions}") from error
    require_owner_symlink_free_path(source_transactions, label)
    try:
        descriptor = os.open(
            source_transactions,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {source_transactions}") from error
    retained_transactions: list[tuple[str, int, os.stat_result]] = []
    retained_forensics: list[tuple[str, os.stat_result]] = []
    regular_forensics: list[tuple[str, os.stat_result]] = []
    try:
        root_profile = require_directory_descriptor_binding(
            source_transactions, descriptor, label
        )
        if (
            root_profile.st_uid != os.geteuid()
            or stat.S_IMODE(root_profile.st_mode) & 0o022
        ):
            raise ValueError(f"{label} is not an owner-controlled directory")
        transaction_names, _ = stable_bound_directory_entries(
            source_transactions, descriptor, root_profile, label
        )
        for transaction_name in transaction_names:
            if F116_SOURCE_AUTHORITY_FORENSIC_ENTRY_PATTERN.fullmatch(
                transaction_name
            ) is not None:
                forensic_label = f"{label} forensic entry {transaction_name}"
                try:
                    forensic_profile = os.stat(
                        transaction_name,
                        dir_fd=descriptor,
                        follow_symlinks=False,
                    )
                except OSError as error:
                    raise ValueError(f"{forensic_label} is unavailable") from error
                forensic_mode = stat.S_IMODE(forensic_profile.st_mode)
                if stat.S_ISREG(forensic_profile.st_mode):
                    if (
                        forensic_profile.st_uid != os.geteuid()
                        or forensic_profile.st_nlink not in {1, 2}
                        or forensic_mode & 0o022
                    ):
                        raise ValueError(f"{forensic_label} profile differs")
                    regular_forensics.append(
                        (transaction_name, forensic_profile)
                    )
                elif stat.S_ISDIR(forensic_profile.st_mode):
                    if (
                        forensic_profile.st_uid != os.geteuid()
                        or forensic_mode & 0o022
                    ):
                        raise ValueError(f"{forensic_label} profile differs")
                else:
                    raise ValueError(f"{forensic_label} has an invalid type")
                retained_forensics.append((transaction_name, forensic_profile))
                continue
            transaction_id = transaction_name
            if transaction_name.endswith(F116_SOURCE_AUTHORITY_STAGING_SUFFIX):
                transaction_id = transaction_name.removesuffix(
                    F116_SOURCE_AUTHORITY_STAGING_SUFFIX
                )
            elif transaction_name.endswith(F116_SOURCE_AUTHORITY_RETIRED_SUFFIX):
                transaction_id = transaction_name.removesuffix(
                    F116_SOURCE_AUTHORITY_RETIRED_SUFFIX
                )
            if (
                F116_SOURCE_AUTHORITY_TRANSACTION_ID_PATTERN.fullmatch(transaction_id)
                is None
            ):
                raise ValueError(f"{label} contains a malformed transaction")
            transaction_label = "non-authoritative F116 source-authority transaction debris"
            try:
                transaction_descriptor = os.open(
                    transaction_name,
                    os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except OSError as error:
                raise ValueError(f"{transaction_label} is unavailable") from error
            try:
                transaction_profile = os.fstat(transaction_descriptor)
                retained_mode = stat.S_IMODE(transaction_profile.st_mode)
                if (
                    not stat.S_ISDIR(transaction_profile.st_mode)
                    or transaction_profile.st_uid != os.geteuid()
                    or retained_mode not in {0o700, 0o700 | stat.S_ISGID}
                ):
                    raise ValueError(f"{transaction_label} profile differs")
                require_bound_directory_entry(
                    descriptor,
                    transaction_name,
                    transaction_descriptor,
                    transaction_profile,
                    transaction_label,
                )
                retained_transactions.append(
                    (transaction_name, transaction_descriptor, transaction_profile)
                )
                transaction_descriptor = None
            finally:
                if transaction_descriptor is not None:
                    os.close(transaction_descriptor)
        for forensic_name, forensic_profile in regular_forensics:
            if forensic_profile.st_nlink == 2 and sum(
                profile_identity(other_profile) == profile_identity(forensic_profile)
                for _other_name, other_profile in regular_forensics
            ) != 2:
                raise ValueError(
                    f"{label} forensic entry {forensic_name} has an external "
                    "retained hardlink"
                )
        retained_names, _ = stable_bound_directory_entries(
            source_transactions, descriptor, root_profile, label
        )
        if retained_names != transaction_names:
            raise ValueError(f"{label} changed during authentication")
        for transaction_name, transaction_descriptor, transaction_profile in (
            retained_transactions
        ):
            require_bound_directory_entry(
                descriptor,
                transaction_name,
                transaction_descriptor,
                transaction_profile,
                "non-authoritative F118 source-authority transaction debris",
            )
        for forensic_name, forensic_profile in retained_forensics:
            try:
                retained_profile = os.stat(
                    forensic_name,
                    dir_fd=descriptor,
                    follow_symlinks=False,
                )
            except OSError as error:
                raise ValueError(
                    f"{label} forensic entry {forensic_name} changed"
                ) from error
            if stat.S_ISREG(forensic_profile.st_mode):
                if (
                    not stat.S_ISREG(retained_profile.st_mode)
                    or regular_file_stable_identity(retained_profile)
                    != regular_file_stable_identity(forensic_profile)
                ):
                    raise ValueError(
                        f"{label} forensic entry {forensic_name} changed"
                    )
            elif (
                not stat.S_ISDIR(retained_profile.st_mode)
                or profile_identity(retained_profile) != profile_identity(forensic_profile)
                or directory_security_identity(retained_profile)
                != directory_security_identity(forensic_profile)
            ):
                raise ValueError(f"{label} forensic entry {forensic_name} changed")
    finally:
        for _name, transaction_descriptor, _profile in retained_transactions:
            try:
                os.close(transaction_descriptor)
            except OSError:
                pass
        os.close(descriptor)


def require_current_source_authority_for_prepare(
    paths: dict[str, Path],
    *,
    bundle_provenance: dict[str, object] | None,
    utility_provenance: dict[str, object],
    matrix_path: Path,
    input_revision: str,
    offline_local_root: bool,
    now: datetime | None = None,
) -> dict[str, object] | None:
    """Authenticate fixed F118 authority and immutable F116/F115 history."""

    if offline_local_root:
        return None
    root = paths["root"]
    now = datetime.now(timezone.utc) if now is None else now.astimezone(timezone.utc)
    historical_f115 = require_historical_f115_source_authority(paths)
    f115_digests = {
        key: historical_f115[key]
        for key in (
            "evidence_sha256", "publication_audit_sha256",
            "provenance_review_sha256", "plasma_review_sha256",
        )
    }
    f115_bindings = {
        "evidence": {
            "path": R03_F115_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "sha256": R03_F115_SOURCE_AUTHORITY_SHA256,
        },
        "publication_audit": {
            "path": f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.publication_audit.json",
            "sha256": R03_F115_PUBLICATION_AUDIT_SHA256,
        },
        "provenance_review": {
            "path": f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.provenance_security_review.json",
            "sha256": R03_F115_PROVENANCE_REVIEW_SHA256,
        },
        "plasma_review": {
            "path": f"{R03_F115_SOURCE_AUTHORITY_RELATIVE}.plasma_scientific_review.json",
            "sha256": R03_F115_PLASMA_REVIEW_SHA256,
        },
    }
    fixed_paths = {
        "evidence": root / F118_CURRENT_SOURCE_AUTHORITY_RELATIVE,
        "provenance_review": root / F118_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / F118_PLASMA_REVIEW_RELATIVE,
        "publication_audit": root / F118_PUBLICATION_AUDIT_RELATIVE,
    }
    loaded: dict[str, tuple[dict[str, object], str]] = {
        key: read_f116_authority_json(path, f"F118 {key}", mode=0o444)
        for key, path in fixed_paths.items()
    }
    evidence, evidence_sha = loaded["evidence"]
    evidence = require_r17_exact_keys(
        evidence,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        },
        "F118 current-source authority",
    )
    generated = require_r17_utc(evidence["generated_utc"], "F118 generation time")
    scope = require_r17_exact_keys(
        evidence["scope"],
        {"relationship", "summary", "preserves", "does_not_authorize"},
        "F118 scope",
    )
    predecessors = require_r17_exact_keys(
        evidence["predecessor_authorities"], {"historical_f116"},
        "F118 predecessor authorities",
    )
    implementation = require_r17_exact_keys(
        evidence["implementation"],
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "predecessor_current_source_bundle",
            "current_source_bundle",
        },
        "F118 implementation",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or evidence["checkpoint"] != "F-118"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or scope["relationship"] != "current-source-selection-only-supersession"
        or not require_r17_nonempty_string(scope["summary"], "F118 scope summary")
        or scope["preserves"] != F118_PRESERVES
        or scope["does_not_authorize"] != F118_DOES_NOT_AUTHORIZE
        or evidence["authorization"] != F118_AUTHORIZATION
        or evidence["validation"] != F118_VALIDATION_CLAIMS
        or evidence["publication_requirements"] != F118_PUBLICATION_REQUIREMENTS
        or generated > now + R17_AUTHORIZATION_FUTURE_SKEW
    ):
        raise ValueError("F118 identity, historical F116 scope, or authority differs")

    historical_bindings = require_r17_exact_keys(
        predecessors["historical_f116"],
        {"evidence", "publication_audit", "provenance_review", "plasma_review"},
        "F118 historical F116 bindings",
    )
    historical_paths = {
        "evidence": root / F116_CURRENT_SOURCE_AUTHORITY_RELATIVE,
        "publication_audit": root / F116_PUBLICATION_AUDIT_RELATIVE,
        "provenance_review": root / F116_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / F116_PLASMA_REVIEW_RELATIVE,
    }
    historical_loaded: dict[str, tuple[dict[str, object], str]] = {}
    for key, path in historical_paths.items():
        binding = require_r17_exact_keys(
            historical_bindings[key], {"path", "sha256"}, f"historical F116 {key}"
        )
        digest = require_r17_sha256(binding["sha256"], f"historical F116 {key}")
        if (
            binding["path"] != path.relative_to(root).as_posix()
            or (root == DEFAULT_ROOT and digest != F116_CANONICAL_SHA256[key])
        ):
            raise ValueError(f"historical F116 {key} binding differs")
        historical_loaded[key] = read_f116_authority_json(
            path, f"historical F116 {key}", mode=0o444
        )
        if historical_loaded[key][1] != digest:
            raise ValueError(f"historical F116 {key} digest differs")
    historical_evidence, historical_evidence_sha = historical_loaded["evidence"]
    historical_evidence = require_r17_exact_keys(
        historical_evidence,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "scope", "predecessor_authorities", "implementation",
            "source_archive_catalog", "authorization", "validation",
            "publication_requirements",
        },
        "historical F116 evidence",
    )
    historical_predecessors = require_r17_exact_keys(
        historical_evidence["predecessor_authorities"],
        {"historical_f115"},
        "historical F116 predecessor authorities",
    )
    historical_implementation = require_r17_exact_keys(
        historical_evidence["implementation"],
        {
            "publisher", "committed_tools", "intermediate_36140_bundle",
            "current_source_bundle",
        },
        "historical F116 implementation",
    )
    historical_current = require_r17_exact_keys(
        historical_implementation["current_source_bundle"],
        {
            "candidate_path", "path", "sha256", "complete_history", "head",
            "advertised_tip", "verified_revisions", "selected_as_current", "subject",
        },
        "historical F116 current source bundle",
    )
    if (
        historical_evidence["schema_version"] != 1
        or historical_evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or historical_evidence["checkpoint"] != "F-116"
        or historical_evidence["execution_epoch"] != EXECUTION_EPOCH
        or historical_predecessors["historical_f115"] != f115_bindings
        or historical_evidence["authorization"] != F116_AUTHORIZATION
        or historical_current["complete_history"] is not True
        or historical_current["selected_as_current"] is not True
    ):
        raise ValueError("historical F116 evidence or nested F115 binding differs")
    historical_audit = require_r17_exact_keys(
        historical_loaded["publication_audit"][0],
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f115_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        },
        "historical F116 publication audit",
    )
    if (
        historical_audit["checkpoint"] != "F-116"
        or historical_audit["historical_f115_authority"] != f115_digests
        or historical_audit["authority_and_enforcement"] != F116_AUTHORIZATION
    ):
        raise ValueError("historical F116 publication audit differs")
    historical_artifact = require_r17_exact_keys(
        historical_audit["artifact"], {"path", "sha256", "mode", "links"},
        "historical F116 audit artifact",
    )
    if historical_artifact != {
        "path": str(historical_paths["evidence"]),
        "sha256": historical_evidence_sha,
        "mode": "0444",
        "links": 1,
    }:
        raise ValueError("historical F116 audit artifact binding differs")

    def bundle_declaration(value: object, label: str, *,
                           current: bool) -> dict[str, object]:
        keys = {
            "path", "sha256", "complete_history", "head", "advertised_tip",
            "verified_revisions", "selected_as_current",
        }
        keys |= {"candidate_path", "subject"} if current else {"role"}
        retained = require_r17_exact_keys(
            value, keys, label
        )
        relative_text = require_r17_nonempty_string(retained["path"], f"{label} path")
        relative = Path(relative_text)
        head = require_source_authority_revision(retained["head"], f"{label} head")
        revisions_value = retained["verified_revisions"]
        if not isinstance(revisions_value, list) or not revisions_value:
            raise ValueError(f"{label} verified revisions must be a nonempty list")
        revisions = [
            require_source_authority_revision(item, f"{label} verified revision")
            for item in revisions_value
        ]
        tip = require_r17_exact_keys(
            retained["advertised_tip"], {"revision", "name"},
            f"{label} advertised tip",
        )
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative.as_posix() != relative_text
            or relative.parent != Path("source-archives")
            or relative.name
            != f"athenak-feature-cgl-through-{head[:9]}.bundle"
            or len(set(revisions)) != len(revisions)
            or tip["revision"] != head
            or retained["complete_history"] is not True
        ):
            raise ValueError(f"{label} identity or complete-history declaration differs")
        if current:
            candidate_text = require_r17_nonempty_string(
                retained["candidate_path"], f"{label} candidate path"
            )
            candidate = Path(candidate_text)
            if (
                not candidate.is_absolute()
                or ".." in candidate.parts
                or str(candidate) != candidate_text
                or tip["name"] != "HEAD"
                or retained["selected_as_current"] is not True
                or not require_r17_nonempty_string(retained["subject"], f"{label} subject")
            ):
                raise ValueError(f"{label} is not the exact selected final-HEAD bundle")
        elif (
            tip["name"] != "refs/heads/feature/cgl-landau-fluid"
            or retained["selected_as_current"] is not False
            or retained["role"] != "retained-non-current-bridge"
        ):
            raise ValueError(f"{label} role or branch-ref declaration differs")
        return {
            **retained,
            "relative": relative,
            "sha256": require_r17_sha256(retained["sha256"], f"{label} SHA-256"),
            "head": head,
            "revisions": revisions,
        }

    bridge = bundle_declaration(
        implementation["intermediate_36140_bundle"],
        "F118 retained bridge bundle",
        current=False,
    )
    final = bundle_declaration(
        implementation["current_source_bundle"],
        "F118 current source bundle",
        current=True,
    )
    predecessor_bundle = require_r17_exact_keys(
        implementation["predecessor_current_source_bundle"],
        {
            "path", "sha256", "complete_history", "head", "advertised_tip",
            "verified_revisions", "selected_as_current", "role", "subject",
        },
        "F118 predecessor current source bundle",
    )
    expected_predecessor = dict(historical_current)
    expected_predecessor.pop("candidate_path")
    expected_predecessor["selected_as_current"] = False
    expected_predecessor["role"] = "retained-non-current-predecessor"
    if predecessor_bundle != expected_predecessor:
        raise ValueError("F118 predecessor current bundle differs from immutable F116")
    promoted_revision = str(final["head"])
    final_revisions = list(final["revisions"])
    selected_revisions = (
        bundle_provenance.get("verified_revisions")
        if isinstance(bundle_provenance, dict) else None
    )
    required_revisions = {
        *F116_PRODUCTION_REQUIRED_REVISIONS,
        promoted_revision,
        require_source_authority_revision(input_revision, "current input revision"),
    }
    if isinstance(selected_revisions, list) and all(
        isinstance(revision, str)
        and GIT_REVISION_PATTERN.fullmatch(revision) is not None
        for revision in selected_revisions
    ):
        required_revisions.update(selected_revisions)
    if (
        bridge["head"] != F116_BRIDGE_REVISION
        or bridge["sha256"] != F116_BRIDGE_SHA256
        or Path(str(bridge["relative"])).name != F116_BRIDGE_NAME
        or final["subject"] != source_authority_commit_subject(promoted_revision)
        or not required_revisions.issubset(set(final_revisions))
        or any(
            not source_authority_revision_is_ancestor(revision, promoted_revision)
            for revision in final_revisions
        )
    ):
        raise ValueError("F118 bridge or final source history differs")

    tools_value = implementation["committed_tools"]
    if not isinstance(tools_value, list) or len(tools_value) != len(F118_REQUIRED_TOOLS):
        raise ValueError("F118 committed-tool vector differs")
    tools: dict[str, dict[str, object]] = {}
    for retained, (relative, expected_mode) in zip(
        tools_value, sorted(F118_REQUIRED_TOOLS.items())
    ):
        tool = require_r17_exact_keys(
            retained, {"path", "revision", "sha256", "mode"},
            f"F118 committed tool {relative}",
        )
        digest = require_r17_sha256(tool["sha256"], f"F118 committed tool {relative}")
        if (
            tool["path"] != relative
            or tool["revision"] != promoted_revision
            or tool["mode"] != expected_mode
            or source_authority_committed_sha256(
                promoted_revision, relative, f"F118 committed tool {relative}"
            ) != digest
            or source_authority_committed_mode(
                promoted_revision, relative, f"F118 committed tool {relative}"
            ) != expected_mode
        ):
            raise ValueError(f"F118 committed tool differs: {relative}")
        tools[relative] = tool
    publisher = require_r17_exact_keys(
        implementation["publisher"], {"path", "revision", "sha256", "mode"},
        "F118 publisher",
    )
    if publisher != tools["scripts/frontier/cgl_lf_stage_i_source_authority.py"]:
        raise ValueError("F118 publisher binding differs from committed-tool vector")
    helper = tools[PRODUCTION_UTILITY_RELATIVE.as_posix()]
    matrix_relative = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
    if (
        utility_provenance.get("path") != str(Path(__file__).resolve())
        or utility_provenance.get("revision") != promoted_revision
        or utility_provenance.get("sha256") != helper["sha256"]
        or utility_provenance.get("committed") is not True
        or input_revision != QUALIFIED_SOURCE_REVISION
        or source_authority_committed_sha256(
            input_revision,
            matrix_relative.as_posix(),
            "F118 frozen scientific Stage I matrix",
        ) != sha256(matrix_path)
    ):
        raise ValueError("F118 tooling authority or frozen scientific source differs")

    bundle_relative = Path(str(final["relative"]))
    bundle_path = root / bundle_relative
    bundle_sha = str(final["sha256"])
    if (
        not isinstance(bundle_provenance, dict)
        or bundle_provenance.get("path") != str(bundle_path)
        or bundle_provenance.get("sha256") != bundle_sha
        or not isinstance(selected_revisions, list)
        or len(set(selected_revisions)) != len(selected_revisions)
        or any(
            not isinstance(revision, str)
            or GIT_REVISION_PATTERN.fullmatch(revision) is None
            for revision in selected_revisions
        )
        or any(revision not in final_revisions for revision in selected_revisions)
    ):
        raise ValueError("F118 selected final source bundle differs")
    verified_bundle = source_bundle_provenance(
        str(bundle_path), final_revisions, root, False
    )
    if (
        verified_bundle is None
        or verified_bundle["path"] != str(bundle_path)
        or verified_bundle["sha256"] != bundle_sha
        or verified_bundle["verified_revisions"] != final_revisions
        or source_authority_bundle_advertised_tip(bundle_path, bundle_sha)
        != (promoted_revision, "HEAD")
        or source_authority_bundle_advertised_tip(
            root / Path(str(bridge["relative"])), str(bridge["sha256"])
        ) != (F116_BRIDGE_REVISION, "refs/heads/feature/cgl-landau-fluid")
    ):
        raise ValueError("F118 source bundle bytes or advertised tips differ")

    verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "predecessor_current_source_bundle_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": bundle_sha,
        "final_head": promoted_revision,
        "historical_f115_preserved": True,
        "historical_f116_preserved": True,
    }
    reviewers = set()
    candidate_paths = set()
    reviewed_times = []
    for key, kind, decision, label in (
        (
            "provenance_review", "provenance-security", "approved-for-publication",
            "F118 provenance/security review",
        ),
        (
            "plasma_review", "plasma-scientific-continuation", "approved",
            "F118 plasma/scientific review",
        ),
    ):
        review, _ = loaded[key]
        review = require_r17_exact_keys(
            review,
            {
                "schema_version", "record_type", "checkpoint", "execution_epoch",
                "review_kind", "decision", "reviewed_candidate", "published_f118",
                "reviewer", "reviewed_utc", "findings", "limitations", "verified",
            },
            label,
        )
        reviewer = require_r17_exact_keys(
            review["reviewer"],
            {"agent_id", "identity"},
            f"{label} reviewer",
        )
        reviewer_id = require_r17_nonempty_string(
            reviewer["agent_id"], f"{label} reviewer ID"
        )
        # These retained strings declare process separation; they are not
        # treated as cryptographic identities.
        require_r17_nonempty_string(reviewer["identity"], f"{label} reviewer identity")
        reviewed_candidate = require_r17_exact_keys(
            review["reviewed_candidate"], {"path", "sha256"},
            f"{label} reviewed candidate",
        )
        candidate_text = require_r17_nonempty_string(
            reviewed_candidate["path"], f"{label} reviewed candidate path"
        )
        candidate_path = Path(candidate_text)
        reviewed = require_r17_utc(review["reviewed_utc"], f"{label} timestamp")
        if (
            review["schema_version"] != 1
            or review["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or review["checkpoint"] != "F-118"
            or review["execution_epoch"] != EXECUTION_EPOCH
            or review["review_kind"] != kind
            or review["decision"] != decision
            or not candidate_path.is_absolute()
            or ".." in candidate_path.parts
            or str(candidate_path) != candidate_text
            or reviewed_candidate["sha256"] != evidence_sha
            or review["published_f118"]
            != {"path": str(fixed_paths["evidence"]), "sha256": evidence_sha}
            or review["verified"] != verified
            or not isinstance(review["findings"], list)
            or not review["findings"]
            or any(not isinstance(item, str) or not item for item in review["findings"])
            or not isinstance(review["limitations"], list)
            or not review["limitations"]
            or any(not isinstance(item, str) or not item for item in review["limitations"])
            or reviewed < generated
            or reviewed > now + R17_AUTHORIZATION_FUTURE_SKEW
        ):
            raise ValueError(f"{label} identity, independence, or verification differs")
        if reviewer_id in reviewers:
            raise ValueError(
                "F118 independent reviews must declare distinct reviewers as processes"
            )
        reviewers.add(reviewer_id)
        candidate_paths.add(candidate_text)
        reviewed_times.append(reviewed)
    if len(candidate_paths) != 1:
        raise ValueError("F118 independent reviews must bind one exact candidate path")

    def declared(value: object, path: Path, digest: str, mode: str,
                 label: str) -> None:
        binding = require_r17_exact_keys(
            value, {"path", "sha256", "mode", "links"}, label
        )
        if binding != {
            "path": str(path),
            "sha256": digest,
            "mode": mode,
            "links": 1,
        }:
            raise ValueError(f"{label} differs")

    audit, audit_sha = loaded["publication_audit"]
    audit = require_r17_exact_keys(
        audit,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_reviews",
            "historical_f116_authority", "source_archive_catalog",
            "authority_and_enforcement", "publication",
        },
        "F118 publication audit",
    )
    audit_reviews = require_r17_exact_keys(
        audit["independent_reviews"],
        {
            "reviews_bind_exact_published_f118_sha256",
            "provenance_security", "plasma_scientific_continuation",
        },
        "F118 publication audit independent reviews",
    )
    declared(
        audit["artifact"], fixed_paths["evidence"], evidence_sha, "0444",
        "F118 publication audit artifact",
    )
    declared(
        audit_reviews["provenance_security"],
        fixed_paths["provenance_review"], loaded["provenance_review"][1], "0444",
        "F118 publication audit provenance review",
    )
    declared(
        audit_reviews["plasma_scientific_continuation"],
        fixed_paths["plasma_review"], loaded["plasma_review"][1], "0444",
        "F118 publication audit plasma review",
    )
    published = require_r17_utc(audit["published_utc"], "F118 publication time")
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-current-source-authority-supersession-publication-audit"
        or audit["checkpoint"] != "F-118"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit_reviews["reviews_bind_exact_published_f118_sha256"] != evidence_sha
        or audit["historical_f116_authority"]
        != {
            "evidence_sha256": historical_evidence_sha,
            "publication_audit_sha256": historical_loaded["publication_audit"][1],
            "provenance_review_sha256": historical_loaded["provenance_review"][1],
            "plasma_review_sha256": historical_loaded["plasma_review"][1],
        }
        or audit["authority_and_enforcement"] != F118_AUTHORIZATION
        or audit["publication"]
        != "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
        or published < max(reviewed_times)
        or published > now + R17_AUTHORIZATION_FUTURE_SKEW
    ):
        raise ValueError("F118 publication audit identity, chronology, or authority differs")

    source_catalog = require_r17_exact_keys(
        evidence["source_archive_catalog"], {"before", "after"},
        "F118 source-archive catalog evidence",
    )
    before = require_r17_exact_keys(
        source_catalog["before"],
        {
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "predecessor_current_source_bundle_listed_exactly_once",
            "final_bundle_listed", "corrupt_c7_listed", "historical_f115_preserved",
        },
        "F118 predecessor source-archive catalog",
    )
    after = require_r17_exact_keys(
        source_catalog["after"],
        {
            "readme_sha256", "sha256sums_sha256", "bridge_listed_exactly_once",
            "predecessor_current_source_bundle_listed_exactly_once",
            "final_bundle_listed_exactly_once", "corrupt_c7_listed",
            "historical_f115_preserved", "historical_f116_preserved",
            "all_prior_checksum_entries_preserved", "sole_current_source_bundle",
        },
        "F118 published source-archive catalog",
    )
    for key in ("readme_sha256", "sha256sums_sha256"):
        require_r17_sha256(before[key], f"F116 predecessor catalog {key}")
        require_r17_sha256(after[key], f"F116 published catalog {key}")
    if (
        before["bridge_listed_exactly_once"] is not True
        or before["predecessor_current_source_bundle_listed_exactly_once"] is not True
        or before["final_bundle_listed"] is not False
        or before["corrupt_c7_listed"] is not False
        or before["historical_f115_preserved"] is not True
        or after["bridge_listed_exactly_once"] is not True
        or after["predecessor_current_source_bundle_listed_exactly_once"] is not True
        or after["final_bundle_listed_exactly_once"] is not True
        or after["corrupt_c7_listed"] is not False
        or after["historical_f115_preserved"] is not True
        or after["historical_f116_preserved"] is not True
        or after["all_prior_checksum_entries_preserved"] is not True
        or after["sole_current_source_bundle"] != bundle_relative.as_posix()
    ):
        raise ValueError("F118 source-archive catalog claims differ")

    audit_catalog = require_r17_exact_keys(
        audit["source_archive_catalog"],
        {
            "readme", "sha256sums", "bridge_bundle",
            "predecessor_current_source_bundle", "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger",
            "sole_current_source_bundle",
        },
        "F118 publication audit source-archive catalog",
    )
    readme_path = root / "source-archives/README.md"
    sums_path = root / "source-archives/SHA256SUMS"
    readme_payload = read_f116_authority_bytes(
        readme_path, "F116 source-archive README", mode=0o644
    )
    sums_payload = read_f116_authority_bytes(
        sums_path, "F116 source-archive checksum ledger", mode=0o644
    )
    readme_sha = hashlib.sha256(readme_payload).hexdigest()
    sums_sha = hashlib.sha256(sums_payload).hexdigest()
    declared(
        audit_catalog["readme"], readme_path, readme_sha, "0644",
        "F116 publication audit source-archive README",
    )
    declared(
        audit_catalog["sha256sums"], sums_path, sums_sha, "0644",
        "F116 publication audit source-archive checksum ledger",
    )
    expected_bridge_audit = {
        "path": str(root / Path(str(bridge["relative"]))),
        "sha256": bridge["sha256"],
        "mode": "0644",
        "links": 1,
        "head": bridge["head"],
        "role": "retained-non-current-bridge",
        "selected_as_current": False,
    }
    expected_final_audit = {
        "path": str(bundle_path),
        "sha256": bundle_sha,
        "mode": "0644",
        "links": 1,
        "head": promoted_revision,
        "selected_as_current": True,
    }
    predecessor_relative = Path(str(predecessor_bundle["path"]))
    expected_predecessor_audit = {
        "path": str(root / predecessor_relative),
        "sha256": predecessor_bundle["sha256"],
        "mode": "0644",
        "links": 1,
        "head": predecessor_bundle["head"],
        "role": "retained-non-current-predecessor",
        "selected_as_current": False,
    }
    if (
        readme_sha != after["readme_sha256"]
        or sums_sha != after["sha256sums_sha256"]
        or audit_catalog["bridge_bundle"] != expected_bridge_audit
        or audit_catalog["predecessor_current_source_bundle"]
        != expected_predecessor_audit
        or audit_catalog["current_source_bundle"] != expected_final_audit
        or audit_catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
        or audit_catalog["sole_current_source_bundle"] != str(bundle_path)
    ):
        raise ValueError("F118 publication audit source-archive catalog differs")
    try:
        ledger_lines = sums_payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("F116 source-archive checksum ledger is not UTF-8") from error
    ledger = validate_f116_source_archive_checksum_ledger(root, ledger_lines)
    if (
        ledger.get(F116_BRIDGE_NAME) != F116_BRIDGE_SHA256
        or ledger.get(predecessor_relative.name) != predecessor_bundle["sha256"]
        or ledger.get(bundle_relative.name) != bundle_sha
        or ledger.get(R03_F115_SOURCE_BUNDLE_RELATIVE.name)
        != R03_F115_SOURCE_BUNDLE_SHA256
        or F116_CORRUPT_C7_NAME in ledger
    ):
        raise ValueError("F118 source-archive checksum ledger authority differs")
    require_authenticated_f118_source_authority_staging(
        root / "accounting"
        / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F118_source_authority_transactions",
    )
    return {
        "checkpoint": "F-118",
        "evidence": {
            "path": F118_CURRENT_SOURCE_AUTHORITY_RELATIVE.as_posix(),
            "sha256": evidence_sha,
        },
        "provenance_review": {
            "path": F118_PROVENANCE_REVIEW_RELATIVE.as_posix(),
            "sha256": loaded["provenance_review"][1],
        },
        "plasma_review": {
            "path": F118_PLASMA_REVIEW_RELATIVE.as_posix(),
            "sha256": loaded["plasma_review"][1],
        },
        "publication_audit": {
            "path": F118_PUBLICATION_AUDIT_RELATIVE.as_posix(),
            "sha256": audit_sha,
        },
        "final_source_bundle": {
            "path": bundle_relative.as_posix(),
            "sha256": bundle_sha,
            "verified_revisions": final_revisions,
        },
    }


def r17_directory_inventory(path: Path, label: str) -> list[dict[str, str]]:
    """Return the recost-compatible inventory of one flat build manifest."""

    if not path.is_dir():
        raise ValueError(f"{label} is not a directory: {path}")
    entries = []
    for item in sorted(path.iterdir(), key=lambda candidate: candidate.name):
        profile = item.lstat()
        if item.is_symlink() or not stat.S_ISREG(profile.st_mode):
            raise ValueError(f"{label} must contain only direct regular files")
        entries.append({
            "name": item.name,
            "mode": f"{stat.S_IMODE(profile.st_mode):04o}",
            "sha256": sha256(item),
        })
    if not entries:
        raise ValueError(f"{label} must not be empty")
    return entries


def r17_directory_inventory_sha256(path: Path, label: str) -> str:
    """Return the recost-compatible digest of one flat build manifest."""

    return hashlib.sha256(
        (json.dumps(r17_directory_inventory(path, label), sort_keys=True) + "\n").encode()
    ).hexdigest()


def r17_available_storage_bytes(path: Path) -> int:
    """Return live available storage using the shared recost/wave-planner rule."""

    try:
        profile = os.statvfs(path)
    except OSError as error:
        raise ValueError("R17 live storage availability is unavailable") from error
    available = profile.f_bavail * profile.f_frsize
    if available < 0:
        raise ValueError("R17 live storage availability is invalid")
    return available


def current_r17_lineages(paths: dict[str, Path]
                         ) -> dict[str, list[dict[str, object]]]:
    """Reconstruct the recost-compatible current scientific lineages."""

    retained: dict[str, list[dict[str, object]]] = {}
    for case_id in sorted(AUTHORIZED_CASE_IDS):
        segments = accepted_case_segments(paths, case_id)
        if not segments:
            continue
        maximum = max(
            float(item["scientific_inspection"]["final_time"]) for item in segments
        )
        terminals = [
            item for item in segments
            if float(item["scientific_inspection"]["final_time"]) == maximum
        ]
        if len(terminals) != 1:
            raise ValueError(f"{case_id} latest recorded lineage is ambiguous")
        indexed = {
            Path(str(item["_manifest_path"])).resolve(): item for item in segments
        }
        lineage = []
        current = terminals[0]
        visited: set[Path] = set()
        while True:
            current_path = Path(str(current["_manifest_path"])).resolve()
            if current_path in visited:
                raise ValueError(f"{case_id} recorded lineage contains a cycle")
            visited.add(current_path)
            current["_manifest_sha256"] = sha256(current_path)
            lineage.append(current)
            parent = current["command"].get("parent_segment")
            if parent is None:
                break
            if not isinstance(parent, dict) or "manifest" not in parent:
                raise ValueError(f"{case_id} recorded lineage lacks a parent manifest")
            parent_path = Path(str(parent["manifest"])).resolve()
            if parent_path not in indexed:
                raise ValueError(f"{case_id} recorded lineage parent is missing")
            current = indexed[parent_path]
        lineage.reverse()
        if visited != set(indexed):
            raise ValueError(f"{case_id} recorded manifests do not form one lineage")
        retained[case_id] = lineage
    return retained


def r17_lineage_summary_sha256(
    lineages: dict[str, list[dict[str, object]]],
) -> str:
    """Return the exact recost lineage-summary digest."""

    value = {
        case_id: [
            {
                "manifest": item["_manifest_path"],
                "sha256": item["_manifest_sha256"],
                "job_id": item.get("job_id"),
                "segment": item["run"]["segment"],
                "result": item["accounting"]["result"],
                "final_time": float(item["scientific_inspection"]["final_time"]),
            }
            for item in lineage
        ]
        for case_id, lineage in sorted(lineages.items())
    }
    return hashlib.sha256(
        (json.dumps(value, sort_keys=True) + "\n").encode()
    ).hexdigest()


def require_clean_r17_predecessor_state(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    *,
    excluded_manifest_path: Path | None = None,
    source_authority_transactions_authenticated: bool = False,
) -> tuple[str, list[dict[str, str]]]:
    """Require a drained exact-t10 R02-R16 barrier and return its bindings."""

    if active_reservations(reservations):
        raise ValueError("R17 preparation requires zero active Stage I reservations")
    transaction_roots = [
        paths["transactions"],
        paths["accounting"]
        / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions",
    ]
    if not source_authority_transactions_authenticated:
        transaction_roots.append(
            paths["accounting"]
            / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F118_source_authority_transactions"
        )
    for path in transaction_roots:
        require_empty_transaction_directory(path, "R17 preparation transaction store")
    lineages = current_r17_lineages(paths)
    incomplete = []
    for case_id in R17_PREDECESSOR_CASE_IDS:
        lineage = lineages.get(case_id, [])
        if (
            not lineage
            or lineage[-1]["accounting"]["result"] != "accepted"
            or abs(
                float(lineage[-1]["scientific_inspection"]["final_time"])
                - REQUIRED_CASE_FINAL_TIME
            ) > 1.0e-10
        ):
            incomplete.append(case_id)
    if incomplete:
        raise ValueError(
            "R17 readiness requires accepted exact-t10 predecessors: "
            + ", ".join(incomplete)
        )
    manifest_paths = sorted(
        paths["runs"].glob("*/*/manifest/prepared_run.json")
    )
    if excluded_manifest_path is not None:
        excluded = excluded_manifest_path.resolve()
        matches = [path for path in manifest_paths if path.resolve() == excluded]
        if len(matches) != 1:
            raise ValueError(
                "R17 submission predecessor snapshot lacks its exact prepared manifest"
            )
        manifest_paths.remove(matches[0])
    manifest_bindings = [
        {
            "path": path.relative_to(paths["root"]).as_posix(),
            "sha256": sha256(path),
        }
        for path in manifest_paths
    ]
    return r17_lineage_summary_sha256(lineages), manifest_bindings


def compact_json_sha256(value: object) -> str:
    """Return the qualification producer's compact canonical JSON digest."""

    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    ).hexdigest()


def r17_retained_inventory_sha256(value: object) -> str:
    """Return the qualification producer's exact retained-inventory digest."""

    return hashlib.sha256(
        (json.dumps(value, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
    ).hexdigest()


def validate_r17_build_manifest_inventory(
    value: object, inventory_sha256: object,
) -> list[dict[str, str]]:
    """Authenticate the producer's complete sorted build-manifest inventory."""

    if not isinstance(value, list) or not value:
        raise ValueError("R17 build-manifest inventory is empty")
    inventory = []
    names = set()
    for item in value:
        record = require_r17_exact_keys(
            item, {"name", "mode", "sha256"}, "R17 build-manifest inventory record"
        )
        name = require_r17_nonempty_string(
            record["name"], "R17 build-manifest file name"
        )
        if Path(name).name != name or name in names or record["mode"] != "0644":
            raise ValueError("R17 build-manifest inventory is incomplete or ambiguous")
        names.add(name)
        inventory.append({
            "name": name,
            "mode": "0644",
            "sha256": require_r17_sha256(
                record["sha256"], "R17 build-manifest file SHA-256"
            ),
        })
    digest = require_r17_sha256(
        inventory_sha256, "R17 build-manifest inventory SHA-256"
    )
    if (
        inventory != sorted(inventory, key=lambda item: item["name"])
        or r17_retained_inventory_sha256(inventory) != digest
    ):
        raise ValueError("R17 build-manifest inventory digest or ordering differs")
    return inventory


def authenticate_r17_absolute_file_binding(
    root: Path,
    value: object,
    label: str,
    *,
    expected_path: Path | None = None,
    modes: frozenset[int] = frozenset({
        0o440, 0o444, 0o640, 0o644, 0o750, 0o755,
    }),
) -> tuple[Path, str, bytes]:
    """Authenticate one producer absolute retained-file binding."""

    binding = require_r17_exact_keys(value, {"path", "sha256", "size_bytes"}, label)
    path_text = require_r17_nonempty_string(binding["path"], f"{label} path")
    path = Path(path_text)
    if (
        not path.is_absolute()
        or path != path.absolute()
        or ".." in path.parts
        or str(path) != path_text
    ):
        raise ValueError(f"{label} path must be normalized and absolute")
    try:
        if path.resolve().relative_to(root.resolve()) != path.relative_to(root):
            raise ValueError(f"{label} path escapes or aliases the Stage I root")
    except ValueError as error:
        raise ValueError(f"{label} path escapes or aliases the Stage I root") from error
    if expected_path is not None and path != expected_path:
        raise ValueError(f"{label} path differs")
    try:
        mode = stat.S_IMODE(path.lstat().st_mode)
    except OSError as error:
        raise ValueError(f"{label} is unavailable: {path}") from error
    if mode not in modes:
        raise ValueError(f"{label} mode differs: {mode:04o}")
    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        owner_controlled=True,
        symlink_free=True,
    )
    digest = require_r17_sha256(binding["sha256"], f"{label} SHA-256")
    size = require_r17_integer(binding["size_bytes"], f"{label} size", 0)
    if len(retained) != size or hashlib.sha256(retained).hexdigest() != digest:
        raise ValueError(f"{label} retained file binding differs")
    return path, digest, retained


def authenticate_r17_nested_file_bindings(root: Path, value: object,
                                          label: str) -> None:
    """Authenticate every producer retained-file binding nested in evidence."""

    if isinstance(value, dict):
        if {"path", "sha256", "size_bytes"}.issubset(value):
            path = value.get("path")
            if isinstance(path, str) and Path(path).is_absolute():
                authenticate_r17_absolute_file_binding(
                    root,
                    {key: value[key] for key in ("path", "sha256", "size_bytes")},
                    label,
                )
        for key, item in value.items():
            authenticate_r17_nested_file_bindings(root, item, f"{label}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            authenticate_r17_nested_file_bindings(root, item, f"{label}[{index}]")


def validate_r17_frozen_science_build_contract(
    value: object, profile: dict[str, object],
) -> dict[str, object]:
    """Authenticate the producer's full frozen science, source, and build contract."""

    contract = require_r17_exact_keys(
        value,
        {
            "case_id", "case_name", "profile_class", "resolution", "mesh_shape",
            "meshblock_shape", "target_time", "scientific_policy", "run_basename",
            "source_revision", "source_bundle_sha256", "matrix_sha256",
            "input_sha256", "provenance_sha256", "execution_intent_sha256",
            "execution_contract_sha256", "parameter_contract",
            "parameter_contract_sha256", "executable_revision",
            "executable_sha256", "build_manifest_inventory_sha256",
        },
        "R17 frozen science/build contract",
    )
    expected = {
        "case_id": R17_CASE_ID,
        "case_name": "paper_scale_separation_active_alfvenic_beta10_nperp384",
        "profile_class": "scale_separation_384x384x768",
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "target_time": 0.25,
        "scientific_policy": "active_hardwall",
        "run_basename": "qualification_R17_n08",
        "source_revision": QUALIFIED_SOURCE_REVISION,
        "source_bundle_sha256": profile.get("source_bundle_sha256"),
        "matrix_sha256": R17_FROZEN_MATRIX_SHA256,
        "input_sha256": R17_FROZEN_INPUT_SHA256,
        "parameter_contract": R17_FROZEN_PARAMETER_CONTRACT,
        "parameter_contract_sha256": compact_json_sha256(
            R17_FROZEN_PARAMETER_CONTRACT
        ),
        "executable_revision": profile.get("executable_revision"),
        "executable_sha256": profile.get("executable_sha256"),
        "build_manifest_inventory_sha256": profile.get("build_manifest_sha256"),
    }
    for key in (
        "source_bundle_sha256", "provenance_sha256", "execution_intent_sha256",
        "execution_contract_sha256", "parameter_contract_sha256",
        "executable_sha256", "build_manifest_inventory_sha256",
    ):
        require_r17_sha256(contract[key], f"R17 frozen contract {key}")
    if (
        any(contract.get(key) != expected_value for key, expected_value in expected.items())
        or contract["executable_revision"] != contract["source_revision"]
        or (
            contract["executable_revision"],
            contract["executable_sha256"],
        ) not in QUALIFIED_RESTART_BINARY_ABIS
    ):
        raise ValueError("R17 frozen science/build contract differs")
    return contract


def validate_r17_prepared_wave_binding(
    root: Path,
    wave_path: Path,
    payload: bytes,
    frozen: dict[str, object],
    build_inventory: list[dict[str, str]],
) -> tuple[dict[str, object], dict[str, object]]:
    """Authenticate the exact producer wave and its sole R17 execution packet."""

    try:
        wave = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 prepared wave is not valid JSON") from error
    wave = require_r17_exact_keys(
        wave,
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "policy", "provenance", "provenance_sha256", "waves",
        },
        "R17 prepared wave",
    )
    policy = require_r17_exact_keys(
        wave["policy"],
        {
            "canonical_acceptance_eligible", "submission_action_available",
            "submission_performed", "operational_only", "max_concurrent_nodes",
            "selection_minimum_runtime_savings_seconds",
            "selection_maximum_node_hour_ratio", "selection_runtime_tie_fraction",
        },
        "R17 prepared-wave policy",
    )
    provenance = require_r17_exact_keys(
        wave["provenance"],
        {
            "source", "source_bundle", "matrix", "executable",
            "build_manifest", "qualification_helper",
        },
        "R17 prepared-wave provenance",
    )
    source = require_r17_exact_keys(
        provenance["source"], {"directory", "revision"}, "R17 prepared-wave source"
    )
    source_bundle = require_r17_exact_keys(
        provenance["source_bundle"],
        {"path", "sha256", "size_bytes", "verified_revisions"},
        "R17 prepared-wave source bundle",
    )
    matrix = require_r17_exact_keys(
        provenance["matrix"], {"path", "sha256", "size_bytes"},
        "R17 prepared-wave matrix",
    )
    executable = require_r17_exact_keys(
        provenance["executable"], {"path", "sha256", "size_bytes", "revision"},
        "R17 prepared-wave executable",
    )
    build = require_r17_exact_keys(
        provenance["build_manifest"],
        {"path", "athena_sha256", "environment", "inventory", "inventory_sha256"},
        "R17 prepared-wave build manifest",
    )
    helper = require_r17_exact_keys(
        provenance["qualification_helper"],
        {"path", "sha256", "size_bytes", "revision", "committed"},
        "R17 prepared-wave qualification helper",
    )
    helper_revision = require_source_authority_revision(
        helper["revision"], "R17 prepared-wave qualification helper revision"
    )
    require_r17_sha256(helper["sha256"], "R17 prepared-wave qualification helper SHA-256")
    revisions = source_bundle["verified_revisions"]
    if (
        helper["committed"] is not True
        or not isinstance(revisions, list)
        or revisions != sorted({frozen["source_revision"], helper_revision})
    ):
        raise ValueError("R17 prepared-wave source bundle revisions differ")
    source_bundle_path, _, _ = authenticate_r17_absolute_file_binding(
        root,
        {
            key: source_bundle[key] for key in ("path", "sha256", "size_bytes")
        },
        "R17 prepared-wave source bundle",
        modes=frozenset({0o644}),
    )
    executable_path, _, _ = authenticate_r17_absolute_file_binding(
        root,
        {key: executable[key] for key in ("path", "sha256", "size_bytes")},
        "R17 prepared-wave executable",
        modes=frozenset({0o755}),
    )
    _, _, digest_payload = authenticate_r17_absolute_file_binding(
        root,
        build["athena_sha256"],
        "R17 prepared-wave build executable digest",
        modes=frozenset({0o644}),
    )
    authenticate_r17_absolute_file_binding(
        root,
        build["environment"],
        "R17 prepared-wave build environment",
        modes=frozenset({0o644}),
    )
    build_path = Path(require_r17_nonempty_string(
        build["path"], "R17 prepared-wave build-manifest path"
    ))
    if (
        not build_path.is_absolute()
        or build_path != build_path.absolute()
        or ".." in build_path.parts
    ):
        raise ValueError("R17 prepared-wave build-manifest path differs")
    try:
        build_path.resolve().relative_to(root.resolve())
    except ValueError as error:
        raise ValueError("R17 prepared-wave build manifest escapes the Stage I root") from error
    if trusted_directory_entries(
        build_path, "R17 prepared-wave build-manifest directory"
    ) != [record["name"] for record in build_inventory]:
        raise ValueError("R17 prepared-wave build-manifest inventory is incomplete")
    for record in build_inventory:
        retained = read_r17_evidence_bytes(
            build_path / record["name"],
            f"R17 build-manifest file {record['name']}",
            mode=0o644,
            owner_controlled=True,
            symlink_free=True,
        )
        if hashlib.sha256(retained).hexdigest() != record["sha256"]:
            raise ValueError("R17 build-manifest file checksum differs")
    waves = wave["waves"]
    if not isinstance(waves, list) or len(waves) != 1:
        raise ValueError("R17 prepared wave must contain one exclusive wave")
    wave_record = require_r17_exact_keys(
        waves[0], {"wave", "total_nodes", "packets"}, "R17 prepared launch wave"
    )
    packets = wave_record["packets"]
    if not isinstance(packets, list) or len(packets) != 1:
        raise ValueError("R17 prepared wave must contain one exact packet")
    packet = require_r17_exact_keys(
        packets[0], {"execution_intent", "execution_intent_sha256"},
        "R17 prepared packet",
    )
    intent = require_r17_exact_keys(
        packet["execution_intent"],
        {
            "packet_id", "project_root", "qualification_root", "case_id",
            "case_name", "profile_class", "target_time", "run_basename",
            "job_name", "scientific_policy", "input", "provenance_sha256",
            "allocation", "paths", "execution_contract_sha256",
            "authenticated_commands", "execution_intent_sha256",
        },
        "R17 prepared execution intent",
    )
    allocation = require_r17_exact_keys(
        intent["allocation"],
        {
            "nodes", "walltime", "walltime_seconds", "athena_walltime",
            "athena_walltime_seconds", "ranks_per_node", "cpus_per_task",
        },
        "R17 prepared allocation",
    )
    input_record = require_r17_exact_keys(
        intent["input"], {"path", "sha256", "size_bytes"}, "R17 prepared input"
    )
    retained_intent_sha = require_r17_sha256(
        packet["execution_intent_sha256"], "R17 prepared execution-intent SHA-256"
    )
    intent_without_digest = dict(intent)
    embedded_intent_sha = intent_without_digest.pop("execution_intent_sha256")
    intent_core = dict(intent_without_digest)
    intent_core.pop("authenticated_commands")
    execution_contract_sha = intent_core.pop("execution_contract_sha256")
    provenance_sha = require_r17_sha256(
        wave["provenance_sha256"], "R17 prepared provenance SHA-256"
    )
    if (
        wave_path.name != "prepared_wave.json"
        or wave_path.parent != Path(str(wave["qualification_root"]))
        or wave["schema_version"] != 1
        or wave["record_type"] != "cgl_lf_stage_i_qualification_wave"
        or wave["project_root"] != str(root)
        or policy["canonical_acceptance_eligible"] is not False
        or policy["submission_action_available"] is not True
        or policy["submission_performed"] is not False
        or policy["operational_only"] is not True
        or policy["max_concurrent_nodes"] != 8
        or wave_record["wave"] != 1
        or wave_record["total_nodes"] != 8
        or compact_json_sha256(provenance) != provenance_sha
        or source["revision"] != frozen["source_revision"]
        or source_bundle_path != Path(str(source_bundle["path"]))
        or source_bundle["sha256"] != frozen["source_bundle_sha256"]
        or frozen["source_revision"] not in source_bundle["verified_revisions"]
        or matrix["sha256"] != frozen["matrix_sha256"]
        or executable["revision"] != frozen["executable_revision"]
        or executable_path != Path(str(executable["path"]))
        or executable["sha256"] != frozen["executable_sha256"]
        or not digest_payload.split()
        or digest_payload.split()[0].decode("ascii", errors="ignore")
        != frozen["executable_sha256"]
        or build["inventory"] != build_inventory
        or build["inventory_sha256"] != frozen["build_manifest_inventory_sha256"]
        or intent["case_id"] != R17_CASE_ID
        or intent["case_name"] != frozen["case_name"]
        or intent["profile_class"] != frozen["profile_class"]
        or intent["target_time"] != frozen["target_time"]
        or intent["scientific_policy"] != frozen["scientific_policy"]
        or intent["run_basename"] != frozen["run_basename"]
        or intent["project_root"] != str(root)
        or intent["qualification_root"] != str(wave_path.parent)
        or input_record["sha256"] != frozen["input_sha256"]
        or intent["provenance_sha256"] != provenance_sha
        or allocation["nodes"] != 8
        or allocation["ranks_per_node"] != 8
        or embedded_intent_sha != retained_intent_sha
        or compact_json_sha256(intent_without_digest) != retained_intent_sha
        or execution_contract_sha != frozen["execution_contract_sha256"]
        or compact_json_sha256(intent_core) != execution_contract_sha
        or retained_intent_sha != frozen["execution_intent_sha256"]
        or provenance_sha != frozen["provenance_sha256"]
    ):
        raise ValueError("R17 prepared wave or frozen binding differs")
    return wave, packet


def validate_r17_decomposition_evidence(
    value: object, output_inventory_sha256: str,
) -> dict[str, object]:
    """Authenticate the inline exact 1728-block/64-rank R17 decomposition proof."""

    evidence = require_r17_exact_keys(
        value,
        {
            "schema_version", "record_type", "resolution", "mesh_shape",
            "meshblock_shape", "logical_meshblock_grid", "logical_meshblocks",
            "ranks", "meshblocks_per_rank", "complete_block_rank_inventory",
            "complete_block_rank_inventory_sha256",
            "terminal_rank_local_output_inventory_sha256", "checks",
        },
        "R17 decomposition evidence",
    )
    inventory = evidence["complete_block_rank_inventory"]
    if not isinstance(inventory, list) or len(inventory) != 64:
        raise ValueError("R17 decomposition rank inventory is incomplete")
    expected_locations = {
        (lx1, lx2, lx3, 0)
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    }
    observed_locations = []
    for rank, item in enumerate(inventory):
        record = require_r17_exact_keys(
            item, {"rank", "rank_name", "logical_meshblocks"},
            "R17 decomposition rank record",
        )
        locations = record["logical_meshblocks"]
        if (
            record["rank"] != rank
            or record["rank_name"] != f"rank_{rank:08d}"
            or not isinstance(locations, list)
            or len(locations) != 27
            or locations != sorted(locations)
            or len({tuple(location) for location in locations}) != 27
            or any(
                not isinstance(location, list)
                or len(location) != 4
                or any(
                    isinstance(index, bool) or not isinstance(index, int)
                    for index in location
                )
                or tuple(location) not in expected_locations
                for location in locations
            )
        ):
            raise ValueError("R17 decomposition does not retain 27 unique blocks per rank")
        observed_locations.extend(tuple(location) for location in locations)
    checks = require_r17_exact_keys(
        evidence["checks"],
        {
            "exact_resolution", "exact_rank_count",
            "exact_meshblocks_per_rank", "complete_unique_logical_inventory",
        },
        "R17 decomposition checks",
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-r17-decomposition-evidence"
        or evidence["resolution"] != "384x384x768"
        or evidence["mesh_shape"] != [384, 384, 768]
        or evidence["meshblock_shape"] != [32, 32, 64]
        or evidence["logical_meshblock_grid"] != [12, 12, 12]
        or evidence["logical_meshblocks"] != 1728
        or evidence["ranks"] != 64
        or evidence["meshblocks_per_rank"] != 27
        or evidence["complete_block_rank_inventory_sha256"]
        != compact_json_sha256(inventory)
        or evidence["terminal_rank_local_output_inventory_sha256"]
        != output_inventory_sha256
        or len(observed_locations) != 1728
        or set(observed_locations) != expected_locations
        or checks != {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        }
    ):
        raise ValueError("R17 decomposition evidence differs from the reviewed contract")
    return evidence


def parse_r17_account_scheduler_evidence(payload: bytes) -> list[dict[str, object]]:
    """Parse the exact all-users account scheduler evidence retained for R17."""

    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("R17 account scheduler evidence is not UTF-8") from error
    expected_header = (
        "JobIDRaw|JobName|State|ExitCode|NNodes|ElapsedRaw|Submit|Start|End|"
        "Partition|Account|User"
    )
    if len(lines) < 2 or lines[0] != expected_header:
        raise ValueError("R17 account scheduler evidence header differs")
    records = []
    seen = set()
    missing_times = {"", "Unknown", "N/A", "None"}

    def optional_time(value: str, label: str) -> datetime | None:
        if value in missing_times:
            return None
        try:
            parsed = datetime.fromisoformat(value)
        except ValueError as error:
            raise ValueError(f"{label} is not an ISO timestamp") from error
        if parsed.tzinfo is None:
            raise ValueError(f"{label} lacks an explicit offset")
        return parsed.astimezone(timezone.utc)

    for line in lines[1:]:
        fields = line.split("|")
        if len(fields) != 12 or any(field.strip() != field for field in fields):
            raise ValueError("R17 account scheduler evidence row is malformed")
        (
            job_id, job_name, state_value, exit_code, nodes_text, elapsed_text,
            submit, start, end, partition, account, owner,
        ) = fields
        try:
            nodes = int(nodes_text)
            elapsed = int(elapsed_text)
        except ValueError as error:
            raise ValueError("R17 account scheduler evidence numeric fields differ") from error
        state_parts = state_value.split()
        state = state_parts[0].split("+")[0] if state_parts else ""
        submit_time = optional_time(submit, "R17 account scheduler submit time")
        start_time = optional_time(start, "R17 account scheduler start time")
        end_time = optional_time(end, "R17 account scheduler end time")
        if (
            re.fullmatch(r"[1-9][0-9]*(?:_[0-9]+)?", job_id) is None
            or job_id in seen
            or not job_name
            or not state
            or not exit_code
            or nodes < 0
            or elapsed < 0
            or submit_time is None
            or not partition
            or account.casefold() != ACCOUNT.casefold()
            or not owner
            or (
                start_time is None
                and (
                    elapsed != 0
                    or state in {
                        "RUNNING", "COMPLETING", "SUSPENDED", "COMPLETED",
                    }
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
            raise ValueError("R17 account scheduler evidence identity differs")
        seen.add(job_id)
        records.append({
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
        })
    return sorted(records, key=lambda record: str(record["job_id"]))


def validate_r17_account_exclusivity(
    root: Path,
    qualification: dict[str, object],
    job_id: str,
    scheduler_fields: list[str],
) -> dict[str, object]:
    """Authenticate the full-interval account-wide R17 exclusivity proof."""

    _, raw_sha, raw_payload = require_r17_root_binding(
        root,
        qualification["account_scheduler_evidence"],
        "R17 account scheduler evidence",
        mode=0o444,
    )
    records = parse_r17_account_scheduler_evidence(raw_payload)
    _, _, evidence_payload = require_r17_root_binding(
        root,
        qualification["account_exclusivity_evidence"],
        "R17 account exclusivity evidence",
        mode=0o444,
    )
    try:
        evidence = json.loads(evidence_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 account exclusivity evidence is not valid JSON") from error
    evidence = require_r17_exact_keys(
        evidence,
        {
            "schema_version", "record_type", "execution_epoch", "measured_utc",
            "query_contract", "visibility_contract",
            "raw_account_scheduler_sha256", "qualification_job",
            "qualification_job_sha256", "account_jobs", "account_jobs_sha256",
            "overlapping_job_ids", "exclusive_entire_execution_interval",
        },
        "R17 account exclusivity evidence",
    )
    query = require_r17_exact_keys(
        evidence["query_contract"],
        {
            "account", "all_users", "allocations_only", "expanded_arrays",
            "start_utc", "end_utc", "start_argument", "end_argument",
            "start_scheduler_offset", "end_scheduler_offset", "fields",
        },
        "R17 account exclusivity query contract",
    )
    visibility = require_r17_exact_keys(
        evidence["visibility_contract"],
        {"private_data", "all_users_job_visibility"},
        "R17 account exclusivity visibility contract",
    )
    target_jobs = [record for record in records if record["job_id"] == job_id]
    qualification_job = require_r17_exact_keys(
        evidence["qualification_job"],
        {
            "job_id", "job_name", "state", "exit_code", "nodes",
            "elapsed_seconds", "submit_utc", "start_utc", "end_utc",
            "partition", "account",
        },
        "R17 account exclusivity qualification job",
    )
    settings = {
        item.strip().casefold()
        for item in str(visibility["private_data"]).split(",")
        if item.strip()
    }
    scheduler_job = {
        "job_id": scheduler_fields[0],
        "job_name": scheduler_fields[1],
        "state": scheduler_fields[2],
        "exit_code": scheduler_fields[3],
        "nodes": int(scheduler_fields[4]),
        "elapsed_seconds": int(scheduler_fields[5]),
        "submit_utc": scheduler_fields[6],
        "end_utc": scheduler_fields[7],
    }
    target = target_jobs[0] if len(target_jobs) == 1 else {}
    expected_qualification_job = {
        key: value for key, value in target.items() if key != "owner"
    }
    try:
        interval_start = datetime.fromisoformat(str(target["start_utc"])).astimezone(
            timezone.utc
        )
        interval_end = datetime.fromisoformat(str(target["end_utc"])).astimezone(
            timezone.utc
        )
    except (KeyError, ValueError) as error:
        raise ValueError("R17 account exclusivity target interval differs") from error
    query_start = interval_start - timedelta(seconds=1)
    query_end = interval_end + timedelta(seconds=1)
    expected_query = {
        "account": ACCOUNT,
        "all_users": True,
        "allocations_only": True,
        "expanded_arrays": True,
        "start_utc": query_start.isoformat(timespec="seconds"),
        "end_utc": query_end.isoformat(timespec="seconds"),
        "start_argument": query_start.strftime("%Y-%m-%dT%H:%M:%S"),
        "end_argument": query_end.strftime("%Y-%m-%dT%H:%M:%S"),
        "start_scheduler_offset": query_start.strftime("%z"),
        "end_scheduler_offset": query_end.strftime("%z"),
        "fields": [
            "JobIDRaw", "JobName", "State", "ExitCode", "NNodes", "ElapsedRaw",
            "Submit", "Start", "End", "Partition", "Account", "User",
        ],
    }
    overlapping = []
    for record in records:
        try:
            start = datetime.fromisoformat(str(record["start_utc"])).astimezone(
                timezone.utc
            )
        except ValueError:
            continue
        end_text = str(record["end_utc"])
        end = (
            None
            if end_text in {"", "Unknown", "N/A", "None"}
            else datetime.fromisoformat(end_text).astimezone(timezone.utc)
        )
        if start < interval_end and (end is None or end > interval_start):
            overlapping.append(str(record["job_id"]))
    measured = require_r17_utc(
        evidence["measured_utc"], "R17 account exclusivity measurement"
    )
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-r17-account-exclusivity-evidence"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or evidence["raw_account_scheduler_sha256"] != raw_sha
        or evidence["account_jobs"] != records
        or evidence["account_jobs_sha256"] != compact_json_sha256(records)
        or evidence["qualification_job_sha256"] != compact_json_sha256(qualification_job)
        or len(target_jobs) != 1
        or any(target_jobs[0].get(key) != value for key, value in scheduler_job.items())
        or qualification_job != expected_qualification_job
        or any(qualification_job.get(key) != value for key, value in scheduler_job.items())
        or overlapping != [job_id]
        or evidence["overlapping_job_ids"] != overlapping
        or evidence["exclusive_entire_execution_interval"] is not True
        or query != expected_query
        or measured < interval_end
        or visibility["all_users_job_visibility"] is not True
        or not settings
        or settings.intersection({"all", "jobs"})
    ):
        raise ValueError("R17 account exclusivity evidence differs")
    return evidence


def validate_r17_physics_measurements(
    value: object,
    label: str,
    *,
    output_inventory_sha256: str | None = None,
) -> dict[str, object]:
    """Independently enforce every retained R17 physics measurement."""

    keys = {
        "finite_rank_outputs", "mass_relative_drift_max",
        "mhd_user_mass_mismatch_max", "lf_bad_counts_total",
        "normalized_ct_divb_max", "normalized_ct_divb_threshold",
        "normalized_ct_divb_below_threshold",
    }
    if output_inventory_sha256 is not None:
        keys.add("rank_local_output_inventory_sha256")
    measurements = require_r17_exact_keys(value, keys, label)
    mass_drift = require_r17_decimal(
        measurements["mass_relative_drift_max"], f"{label} mass drift"
    )
    mass_mismatch = require_r17_decimal(
        measurements["mhd_user_mass_mismatch_max"], f"{label} mass mismatch"
    )
    normalized_ct_divb = require_r17_decimal(
        measurements["normalized_ct_divb_max"], f"{label} normalized CT divB"
    )
    normalized_ct_threshold = require_r17_decimal(
        measurements["normalized_ct_divb_threshold"],
        f"{label} normalized CT divB threshold",
    )
    if (
        require_r17_integer(
            measurements["finite_rank_outputs"], f"{label} finite rank outputs"
        )
        != 64
        or mass_drift < 0
        or mass_drift > Decimal("1e-12")
        or mass_mismatch < 0
        or mass_mismatch > Decimal("1e-12")
        or require_r17_integer(
            measurements["lf_bad_counts_total"], f"{label} LF bad-count total"
        )
        != 0
        or normalized_ct_divb < 0
        or normalized_ct_threshold != R17_MAX_NORMALIZED_CT_DIVB
        or normalized_ct_divb >= normalized_ct_threshold
        or measurements["normalized_ct_divb_below_threshold"] is not True
        or (
            output_inventory_sha256 is not None
            and require_r17_sha256(
                measurements["rank_local_output_inventory_sha256"],
                f"{label} rank-local output inventory SHA-256",
            )
            != output_inventory_sha256
        )
    ):
        raise ValueError(f"{label} differ")
    return measurements


def validate_r17_operational_qualification(
    root: Path,
    readiness: dict[str, object],
    profile: dict[str, object],
    readiness_reviewer: str,
    authorization_time: datetime | None = None,
) -> dict[str, object]:
    """Reauthenticate the measured and independently reviewed 64-rank proof."""

    authorization_time = (
        datetime.now(timezone.utc)
        if authorization_time is None
        else authorization_time.astimezone(timezone.utc)
    )
    qualification_path, qualification_sha, retained = require_r17_root_binding(
        root, readiness["operational_qualification"],
        "R17 operational qualification", mode=0o444,
    )
    if (
        qualification_path
        != root / "accounting" / (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_"
            "R17_operational_qualification.json"
        )
    ):
        raise ValueError("R17 operational qualification must be retained under accounting")
    try:
        qualification = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 operational qualification is not valid JSON") from error
    qualification = require_r17_exact_keys(
        qualification,
        {
            "schema_version", "record_type", "execution_epoch", "completed_utc",
            "measured_utc", "measured_by", "job_id", "state", "exit_code",
            "nodes", "ranks", "prepared_wave", "qualification_evidence",
            "scientific_evidence", "scheduler_evidence",
            "account_scheduler_evidence", "account_exclusivity_evidence",
            "executable_sha256", "build_manifest_inventory",
            "build_manifest_inventory_sha256", "rank_local_outputs",
            "rank_local_output_inventory_sha256", "rank_local_restarts",
            "rank_local_restart_inventory_sha256", "decomposition_evidence",
            "restart_load_evidence", "physics_validation_evidence",
            "frozen_science_build_contract", "independent_review_contract",
            "authority",
        },
        "R17 operational qualification",
    )
    completed = require_r17_scheduler_utc(
        qualification["completed_utc"], "R17 qualification completion"
    )
    measured = require_r17_utc(
        qualification["measured_utc"], "R17 qualification measurement"
    )
    measured_by = require_r17_nonempty_string(
        qualification["measured_by"], "R17 qualification measurement author"
    )
    build_inventory = validate_r17_build_manifest_inventory(
        qualification["build_manifest_inventory"],
        qualification["build_manifest_inventory_sha256"],
    )
    build_inventory_sha = str(qualification["build_manifest_inventory_sha256"])
    frozen = validate_r17_frozen_science_build_contract(
        qualification["frozen_science_build_contract"], profile
    )
    if (
        qualification["schema_version"] != 2
        or qualification["record_type"] != "stage-i-r17-operational-qualification"
        or qualification["execution_epoch"] != EXECUTION_EPOCH
        or qualification["state"] != "COMPLETED"
        or qualification["exit_code"] != "0:0"
        or require_r17_integer(
            qualification["nodes"], "R17 qualification nodes", 1
        ) != 8
        or require_r17_integer(
            qualification["ranks"], "R17 qualification ranks", 1
        ) != 64
        or qualification["executable_sha256"] != profile["executable_sha256"]
        or build_inventory_sha != profile["build_manifest_sha256"]
        or qualification["executable_sha256"] != frozen["executable_sha256"]
        or build_inventory_sha != frozen["build_manifest_inventory_sha256"]
        or measured < completed
        or measured > authorization_time
        or completed > authorization_time
    ):
        raise ValueError("R17 operational qualification identity or build binding differs")
    job_id = require_r17_nonempty_string(
        qualification["job_id"], "R17 qualification job ID"
    )
    require_numeric_job_id(job_id)
    authority = require_r17_exact_keys(
        qualification["authority"],
        {
            "r17_launch_authorized", "scheduler_mutation_authorized",
            "canonical_mutation_authorized",
        },
        "R17 operational qualification authority",
    )
    independent_contract = require_r17_exact_keys(
        qualification["independent_review_contract"],
        {
            "required", "path", "mode", "schema_version", "record_type",
            "execution_epoch", "decision", "candidate_path",
            "candidate_sha256_required", "reviewer_must_differ_from",
            "reviewed_after_utc",
        },
        "R17 operational qualification independent-review contract",
    )
    expected_review_path = qualification_path.with_name(
        f"{qualification_path.name}.independent_review.json"
    )
    if (
        authority
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or independent_contract
        != {
            "required": True,
            "path": expected_review_path.relative_to(root).as_posix(),
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
    ):
        raise ValueError("R17 independent-review or authority contract differs")

    wave_path, _, wave_payload = require_r17_root_binding(
        root, qualification["prepared_wave"], "R17 prepared wave", mode=0o640
    )
    wave, packet = validate_r17_prepared_wave_binding(
        root, wave_path, wave_payload, frozen, build_inventory
    )
    intent = packet["execution_intent"]
    assert isinstance(intent, dict)
    intent_paths = require_r17_exact_keys(
        intent["paths"],
        {
            "run_dir", "output_dir", "manifest_dir", "environment_log",
            "execution_contract", "prepared_manifest", "archived_input",
            "batch_script", "scheduler_batch_script", "submission_receipt",
            "submission_recovery", "job_binding", "scheduler_raw",
            "scheduler_evidence", "scientific_evidence", "output_binding",
            "restart_smoke_dir", "restart_smoke_log", "restart_smoke_result",
            "slurm_log",
        },
        "R17 prepared execution paths",
    )
    qualification_evidence_path, _, qualification_evidence_payload = (
        require_r17_root_binding(
            root,
            qualification["qualification_evidence"],
            "R17 qualification evidence",
            mode=0o640,
        )
    )
    scientific_path, _, scientific_payload = require_r17_root_binding(
        root,
        qualification["scientific_evidence"],
        "R17 scientific evidence",
        mode=0o640,
    )
    try:
        qualification_evidence = json.loads(qualification_evidence_payload)
        scientific = json.loads(scientific_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 producer-bound evidence is not valid JSON") from error
    qualification_evidence = require_r17_exact_keys(
        qualification_evidence,
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        },
        "R17 qualification evidence",
    )
    prepared_wave_binding = require_r17_exact_keys(
        qualification_evidence["prepared_wave"],
        {"path", "sha256", "size_bytes"},
        "R17 qualification evidence prepared-wave binding",
    )
    results = qualification_evidence["results"]
    if not isinstance(results, list) or len(results) != 1:
        raise ValueError("R17 qualification evidence must contain one exact result")
    result = require_r17_exact_keys(
        results[0],
        {
            "nodes", "elapsed_seconds", "job_id", "scheduler_end_utc",
            "execution_intent_sha256", "job_binding", "scheduler_evidence",
            "scientific_evidence", "output_binding",
            "scientific_agreement_signature",
        },
        "R17 qualification result",
    )
    scientific = require_r17_exact_keys(
        scientific,
        {
            "schema_version", "record_type", "case_id", "nodes", "target_time",
            "scientific_policy", "executable", "execution_intent_sha256",
            "execution_contract_sha256", "forcing_closure_normalized_residual",
            "mhd_history", "user_history", "snapshots", "snapshot_times",
            "restarts", "restart_times", "restart_load_smoke",
            "terminal_rank_local_outputs",
            "terminal_rank_local_output_inventory_sha256", "r17_decomposition",
            "terminal_rank_local_restarts",
            "terminal_rank_local_restart_inventory_sha256",
            "physics_measurements", "scientific_agreement_signature", "checks",
            "accepted_for_profile_selection",
            "accepted_for_operational_qualification",
        },
        "R17 scientific evidence",
    )
    authenticate_r17_nested_file_bindings(root, qualification_evidence, "R17 evidence")
    authenticate_r17_nested_file_bindings(root, scientific, "R17 scientific evidence")
    expected_result_paths = {
        "job_binding": Path(str(intent_paths["job_binding"])),
        "scheduler_evidence": Path(str(intent_paths["scheduler_evidence"])),
        "scientific_evidence": scientific_path,
        "output_binding": Path(str(intent_paths["output_binding"])),
    }
    expected_result_modes = {
        "job_binding": frozenset({0o640}),
        "scheduler_evidence": frozenset({0o640}),
        "scientific_evidence": frozenset({0o640}),
        "output_binding": frozenset({0o440}),
    }
    result_payloads: dict[str, bytes] = {}
    for key in expected_result_paths:
        _, _, result_payloads[key] = authenticate_r17_absolute_file_binding(
            root,
            result[key],
            f"R17 qualification result {key}",
            expected_path=expected_result_paths[key],
            modes=expected_result_modes[key],
        )
    try:
        job_binding = json.loads(result_payloads["job_binding"])
        producer_scheduler = json.loads(result_payloads["scheduler_evidence"])
        output_binding = json.loads(result_payloads["output_binding"])
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 producer result binding is not valid JSON") from error
    job_binding = require_r17_exact_keys(
        job_binding,
        {
            "schema_version", "record_type", "submitted_utc", "wave", "job_id",
            "job_name", "nodes", "prior_wave_job_ids", "execution_intent_sha256",
            "execution_contract_sha256", "submission_argv",
            "submission_stdin_sha256", "submission_stdin_size_bytes",
            "sbatch_stdout", "sbatch_stderr", "submission_receipt",
            "scheduler_batch_script", "post_submission_queue",
            "submission_boundary", "submission_result", "prepared_wave",
            "prepared_manifest", "archived_input", "batch_script",
        },
        "R17 producer job binding",
    )
    producer_scheduler = require_r17_exact_keys(
        producer_scheduler,
        {
            "schema_version", "record_type", "job_id", "job_name", "state",
            "exit_code", "nodes", "elapsed_seconds", "submit_utc", "start_utc",
            "end_utc", "partition", "account", "alloc_tres", "req_tres",
            "timelimit_minutes", "comment", "submit_line",
            "execution_intent_sha256", "execution_contract_sha256",
            "job_binding_sha256", "raw_scheduler",
        },
        "R17 producer scheduler evidence",
    )
    output_binding = require_r17_exact_keys(
        output_binding,
        {
            "schema_version", "record_type", "job_id", "execution_contract_sha256",
            "executable_sha256", "completed_utc", "output_dir", "files",
            "output_tree_sha256",
        },
        "R17 producer output binding",
    )
    output_files = output_binding["files"]
    output_dir = Path(str(output_binding["output_dir"]))
    if not isinstance(output_files, list) or not output_files:
        raise ValueError("R17 producer output-tree inventory is empty")
    retained_output_files = []
    for index, item in enumerate(output_files):
        record = require_r17_exact_keys(
            item, {"path", "sha256", "size_bytes"},
            f"R17 producer output-tree file {index}",
        )
        relative_text = require_r17_nonempty_string(
            record["path"], f"R17 producer output-tree file {index} path"
        )
        relative = Path(relative_text)
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative.as_posix() != relative_text
        ):
            raise ValueError("R17 producer output-tree inventory path differs")
        absolute = output_dir / relative
        authenticate_r17_absolute_file_binding(
            root,
            {
                "path": str(absolute),
                "sha256": record["sha256"],
                "size_bytes": record["size_bytes"],
            },
            f"R17 producer output-tree file {index}",
        )
        retained_output_files.append(record)
    if len({item["path"] for item in retained_output_files}) != len(
        retained_output_files
    ):
        raise ValueError("R17 producer output-tree inventory contains duplicates")
    authenticate_r17_nested_file_bindings(root, job_binding, "R17 job binding")
    authenticate_r17_nested_file_bindings(
        root, producer_scheduler, "R17 producer scheduler evidence"
    )
    expected_job_identity = {
        "job_id": job_id,
        "job_name": intent["job_name"],
        "nodes": 8,
        "execution_intent_sha256": frozen["execution_intent_sha256"],
        "execution_contract_sha256": frozen["execution_contract_sha256"],
    }
    if (
        job_binding["schema_version"] != 4
        or job_binding["record_type"] != "cgl_lf_stage_i_qualification_job_binding"
        or job_binding["wave"] != 1
        or job_binding["prior_wave_job_ids"] != []
        or any(job_binding.get(key) != value for key, value in expected_job_identity.items())
        or producer_scheduler["schema_version"] != 3
        or producer_scheduler["record_type"]
        != "cgl_lf_stage_i_qualification_scheduler_evidence"
        or producer_scheduler["state"] != "COMPLETED"
        or producer_scheduler["exit_code"] != "0:0"
        or producer_scheduler["end_utc"] != qualification["completed_utc"]
        or producer_scheduler["job_binding_sha256"] != compact_json_sha256(job_binding)
        or any(
            producer_scheduler.get(key) != value
            for key, value in expected_job_identity.items()
        )
        or output_binding["schema_version"] != 1
        or output_binding["record_type"]
        != "cgl_lf_stage_i_qualification_output_binding"
        or output_binding["job_id"] != job_id
        or output_binding["execution_contract_sha256"]
        != frozen["execution_contract_sha256"]
        or output_binding["executable_sha256"] != frozen["executable_sha256"]
        or output_binding["output_dir"] != intent_paths["output_dir"]
        or output_binding["output_tree_sha256"]
        != compact_json_sha256(output_binding["files"])
        or result["scientific_evidence"]["sha256"]
        != qualification["scientific_evidence"]["sha256"]
    ):
        raise ValueError("R17 producer scheduler/job/output binding differs")
    executable = require_r17_exact_keys(
        scientific["executable"], {"revision", "sha256"}, "R17 scientific executable"
    )
    if (
        qualification_evidence_path
        != Path(str(wave["qualification_root"])) / "R17.qualification_evidence.json"
        or qualification_evidence["schema_version"] != 2
        or qualification_evidence["record_type"]
        != "cgl_lf_stage_i_qualification_evidence"
        or qualification_evidence["project_root"] != str(root)
        or qualification_evidence["qualification_root"] != wave["qualification_root"]
        or prepared_wave_binding
        != {
            "path": str(wave_path),
            "sha256": qualification["prepared_wave"]["sha256"],
            "size_bytes": len(wave_payload),
        }
        or qualification_evidence["case_id"] != R17_CASE_ID
        or qualification_evidence["target_time"] != frozen["target_time"]
        or result["nodes"] != 8
        or result["job_id"] != job_id
        or result["scheduler_end_utc"] != qualification["completed_utc"]
        or result["execution_intent_sha256"] != frozen["execution_intent_sha256"]
        or result["scientific_agreement_signature"]
        != scientific["scientific_agreement_signature"]
        or scientific_path != Path(str(intent_paths["scientific_evidence"]))
        or scientific["schema_version"] != 6
        or scientific["record_type"]
        != "cgl_lf_stage_i_qualification_scientific_evidence"
        or scientific["case_id"] != R17_CASE_ID
        or scientific["nodes"] != 8
        or scientific["target_time"] != frozen["target_time"]
        or scientific["scientific_policy"] != frozen["scientific_policy"]
        or scientific["execution_intent_sha256"] != frozen["execution_intent_sha256"]
        or scientific["execution_contract_sha256"] != frozen["execution_contract_sha256"]
        or executable
        != {
            "revision": frozen["executable_revision"],
            "sha256": frozen["executable_sha256"],
        }
        or scientific["accepted_for_profile_selection"] is not False
        or scientific["accepted_for_operational_qualification"] is not True
    ):
        raise ValueError("R17 producer evidence binding differs")

    def inventory(value: object, label: str, *,
                  exactly: int | None) -> tuple[list[dict[str, str]], str]:
        if not isinstance(value, list) or len(value) < 64:
            raise ValueError(f"{label} must retain at least one file per rank")
        if exactly is not None and len(value) != exactly:
            raise ValueError(f"{label} must retain exactly {exactly} files")
        ranks = set()
        paths = set()
        retained_inventory = []
        for index, binding in enumerate(value):
            path, digest, payload = require_r17_root_binding(
                root, binding, f"{label} {index}", mode=0o644
            )
            matches = [
                part for part in path.parts if re.fullmatch(r"rank_[0-9]{8}", part)
            ]
            if len(matches) != 1 or path in paths or not payload:
                raise ValueError(f"{label} inventory is invalid")
            ranks.add(matches[0])
            paths.add(path)
            retained_inventory.append({
                "path": path.relative_to(root).as_posix(),
                "sha256": digest,
            })
        if ranks != {f"rank_{rank:08d}" for rank in range(64)}:
            raise ValueError(f"{label} rank inventory is incomplete")
        if retained_inventory != sorted(
            retained_inventory, key=lambda item: item["path"]
        ):
            raise ValueError(f"{label} rank inventory ordering differs")
        digest = r17_retained_inventory_sha256(retained_inventory)
        return retained_inventory, digest

    outputs, output_inventory = inventory(
        qualification["rank_local_outputs"],
        "R17 qualification rank-local outputs",
        exactly=64,
    )
    restarts, restart_inventory = inventory(
        qualification["rank_local_restarts"],
        "R17 qualification rank-local restarts",
        exactly=64,
    )
    decomposition = validate_r17_decomposition_evidence(
        qualification["decomposition_evidence"], output_inventory
    )
    if (
        qualification["rank_local_output_inventory_sha256"] != output_inventory
        or qualification["rank_local_restart_inventory_sha256"] != restart_inventory
        or scientific["terminal_rank_local_outputs"] != outputs
        or scientific["terminal_rank_local_output_inventory_sha256"] != output_inventory
        or scientific["terminal_rank_local_restarts"] != restarts
        or scientific["terminal_rank_local_restart_inventory_sha256"] != restart_inventory
        or scientific["r17_decomposition"] != decomposition
        or scientific["physics_measurements"] is None
    ):
        raise ValueError("R17 producer rank inventory or scientific binding differs")
    scientific_physics = validate_r17_physics_measurements(
        scientific["physics_measurements"],
        "R17 producer scientific physics measurements",
    )
    scientific_checks = require_r17_exact_keys(
        scientific["checks"],
        {
            "exact_endpoint", "complete_rank_inventory",
            "finite_synchronized_histories", "mass_conserved",
            "strict_lf_failure_counters_zero", "hard_volume_zero",
            "snapshot_hard_bounds_independently_verified",
            "normalized_ct_divb_below_threshold", "interval_cap_counts_valid",
            "nontrivial_forcing_and_pressure_work", "case_aware_policy_passed",
            "snapshot_cadence_complete", "terminal_snapshot_unique",
            "snapshots_structurally_complete", "restart_headers_authenticated",
            "terminal_restart_unique", "restart_load_smoke_passed",
        },
        "R17 producer scientific checks",
    )
    if any(value is not True for value in scientific_checks.values()):
        raise ValueError("R17 producer scientific checks differ")
    measured_times = [measured]
    evidence_bindings: dict[str, dict[str, object]] = {}
    for key, record_type, expected_measurements in (
        (
            "restart_load_evidence",
            "stage-i-r17-restart-load-evidence",
            {
                "rank_local_restart_inventory_sha256": restart_inventory,
                "loaded_rank_count": 64,
                "load_state": "COMPLETED",
                "load_exit_code": "0:0",
            },
        ),
        (
            "physics_validation_evidence",
            "stage-i-r17-physics-validation-evidence",
            None,
        ),
    ):
        _, _, payload = require_r17_root_binding(
            root, qualification[key], f"R17 qualification {key}", mode=0o444
        )
        try:
            evidence = json.loads(payload)
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError(f"R17 qualification {key} is not valid JSON") from error
        evidence = require_r17_exact_keys(
            evidence,
            {
                "schema_version", "record_type", "execution_epoch", "measured_utc",
                "measured_by", "job_id", "executable_sha256",
                "build_manifest_inventory_sha256", "passed", "measurements",
            },
            f"R17 qualification {key}",
        )
        author = require_r17_nonempty_string(
            evidence["measured_by"], f"R17 qualification {key} author"
        )
        measured = require_r17_utc(
            evidence["measured_utc"], f"R17 qualification {key} timestamp"
        )
        measured_times.append(measured)
        if (
            evidence["schema_version"] != 1
            or evidence["record_type"] != record_type
            or evidence["execution_epoch"] != EXECUTION_EPOCH
            or evidence["job_id"] != job_id
            or evidence["executable_sha256"] != profile["executable_sha256"]
            or evidence["build_manifest_inventory_sha256"] != build_inventory_sha
            or evidence["passed"] is not True
            or not isinstance(evidence["measurements"], dict)
            or evidence["measured_utc"] != qualification["measured_utc"]
            or evidence["measured_by"] != measured_by
            or measured < completed
            or measured > authorization_time
        ):
            raise ValueError(f"R17 qualification {key} differs")
        measurements = evidence["measurements"]
        if expected_measurements is not None:
            require_r17_exact_keys(
                measurements, set(expected_measurements),
                "R17 qualification restart-load measurements",
            )
            if measurements != expected_measurements:
                raise ValueError(f"R17 qualification {key} measurements differ")
        else:
            measurements = validate_r17_physics_measurements(
                measurements,
                "R17 qualification physics measurements",
                output_inventory_sha256=output_inventory,
            )
            expected_scientific_measurements = {
                "rank_local_output_inventory_sha256": output_inventory,
                **scientific_physics,
            }
            if measurements != expected_scientific_measurements:
                raise ValueError("R17 qualification physics measurements differ")
        binding = qualification[key]
        assert isinstance(binding, dict)
        evidence_bindings[key] = {
            "path": binding["path"],
            "sha256": binding["sha256"],
            "measured_utc": measured.isoformat(),
            "measured_by": author,
            "measurements": measurements,
        }

    review_path, review_sha, review_payload = require_r17_root_binding(
        root, readiness["operational_qualification_review"],
        "R17 operational qualification independent review", mode=0o444,
    )
    if review_path != expected_review_path:
        raise ValueError("R17 operational qualification review path differs")
    try:
        review = json.loads(review_payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("R17 operational qualification review is not valid JSON") from error
    review = require_r17_exact_keys(
        review,
        {
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate",
        },
        "R17 operational qualification review",
    )
    reviewer = require_r17_nonempty_string(
        review["reviewer"], "R17 operational qualification reviewer"
    )
    reviewed = require_r17_utc(
        review["reviewed_utc"], "R17 operational qualification review time"
    )
    if (
        review["schema_version"] != 1
        or review["record_type"]
        != "stage-i-r17-operational-qualification-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved"
        or review["candidate"]
        != {"path": str(qualification_path), "sha256": qualification_sha}
        # This is the producer's declared process-separation contract.  A
        # reviewer label is not treated as a cryptographic identity.
        or reviewer in independent_contract["reviewer_must_differ_from"]
        or reviewed < max([completed, *measured_times])
        or reviewed > authorization_time + R17_AUTHORIZATION_FUTURE_SKEW
    ):
        raise ValueError(
            "R17 operational qualification review differs or lacks declared "
            "process independence"
        )

    scheduler_path, _, scheduler_payload = require_r17_root_binding(
        root, qualification["scheduler_evidence"],
        "R17 qualification scheduler evidence", mode=0o444,
    )
    if (
        scheduler_path.parent != root / "accounting"
        or scheduler_path.name != f"{job_id}.r17_qualification.sacct.txt"
    ):
        raise ValueError("R17 qualification scheduler evidence path differs")
    try:
        lines = scheduler_payload.decode().splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("R17 qualification scheduler evidence is not UTF-8") from error
    if len(lines) != 1:
        raise ValueError("R17 qualification scheduler evidence must contain one row")
    fields = lines[0].split("|")
    try:
        elapsed = int(fields[5])
    except (IndexError, ValueError) as error:
        raise ValueError("R17 qualification scheduler elapsed time is invalid") from error
    submitted = (
        require_r17_scheduler_utc(fields[6], "R17 qualification scheduler submit time")
        if len(fields) == 8 else None
    )
    scheduler_completed = (
        require_r17_scheduler_utc(fields[7], "R17 qualification scheduler completion time")
        if len(fields) == 8 else None
    )
    if (
        len(fields) != 8
        or fields[0] != job_id
        or re.fullmatch(r"[A-Za-z0-9_.-]+", fields[1]) is None
        or fields[2:5] != ["COMPLETED", "0:0", "8"]
        or any(not field or field.strip() != field for field in fields)
        or elapsed <= 0
        or scheduler_completed != completed
        or submitted is None
        or scheduler_completed <= submitted
    ):
        raise ValueError("R17 qualification scheduler evidence differs")
    if root.resolve() == DEFAULT_ROOT.expanduser().resolve():
        try:
            live = subprocess.run(
                [
                    str(SACCT), "-X", "-j", job_id,
                    "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                    "ElapsedRaw,Submit,End",
                    "-n", "-P",
                ],
                check=True,
                capture_output=True,
                text=True,
                env=scheduler_environment(),
            ).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            raise ValueError(
                "live R17 qualification scheduler evidence is unavailable"
            ) from error
        live_rows = [
            line.rstrip("|").split("|")
            for line in live.splitlines()
            if line.rstrip("|").split("|")[0] == job_id
        ]
        if live_rows != [fields]:
            raise ValueError("live R17 qualification scheduler evidence differs")
    account_exclusivity = validate_r17_account_exclusivity(
        root, qualification, job_id, fields
    )
    if account_exclusivity["measured_utc"] != qualification["measured_utc"]:
        raise ValueError("R17 account-exclusivity measurement binding differs")
    return {
        **qualification,
        "rank_local_output_inventory_sha256": output_inventory,
        "rank_local_restart_inventory_sha256": restart_inventory,
        "build_manifest_inventory_sha256": build_inventory_sha,
        "authenticated_rank_local_outputs": outputs,
        "authenticated_rank_local_restarts": restarts,
        "authenticated_validation_evidence": evidence_bindings,
        "reviewed_by": reviewer,
    }


def validate_fixed_r17_readiness_chain(
    paths: dict[str, Path],
    recost: dict[str, object],
    profile: dict[str, object],
    lineage_sha256: str,
    current_source_authority: dict[str, object],
    now: datetime,
) -> tuple[dict[str, object], dict[str, str]]:
    """Authenticate the fixed reviewed/published R17 readiness chain."""

    readiness_path = paths["root"] / R17_READINESS_RELATIVE
    readiness, readiness_sha = read_controller_publication_json(
        readiness_path, "R17 readiness evidence", mode=0o444
    )
    readiness = require_r17_exact_keys(
        readiness,
        {
            "schema_version", "record_type", "execution_epoch", "root",
            "generated_utc", "expires_utc", "reviewed_by",
            "predecessor_lineages_sha256", "storage_evidence_sha256",
            "computed_projection_sha256", "executable_sha256",
            "build_manifest_sha256", "required_retained_bytes", "nodes", "ranks",
            "storage_ready", "node_hour_ready", "rank_64_ready",
            "operational_qualification", "operational_qualification_review",
            "current_source_authority",
        },
        "R17 readiness evidence",
    )
    generated = require_r17_utc(readiness["generated_utc"], "R17 readiness generation")
    expires = require_r17_utc(readiness["expires_utc"], "R17 readiness expiry")
    recost_generated = require_r17_utc(
        recost["generated_utc"], "latest recost generation"
    )
    recost_expires = require_r17_utc(recost["expires_utc"], "latest recost expiry")
    if (
        generated > now + R17_AUTHORIZATION_FUTURE_SKEW
        or now - generated > R17_AUTHORIZATION_MAX_LIFETIME
        or expires <= now
        or expires <= generated
        or expires - generated > R17_AUTHORIZATION_MAX_LIFETIME
        or generated > recost_generated
        or expires < recost_expires
    ):
        raise ValueError("R17 readiness evidence is stale or has an invalid lifetime")
    readiness_reviewer = require_r17_nonempty_string(
        readiness["reviewed_by"], "R17 readiness reviewer"
    )
    provenance = recost["provenance"]
    assert isinstance(provenance, dict)
    expected = {
        "schema_version": 1,
        "record_type": "stage-i-r17-readiness",
        "execution_epoch": EXECUTION_EPOCH,
        "root": str(paths["root"]),
        "predecessor_lineages_sha256": lineage_sha256,
        "storage_evidence_sha256": provenance.get("storage_evidence_sha256"),
        "computed_projection_sha256": provenance.get("computed_projection_sha256"),
        "executable_sha256": profile["executable_sha256"],
        "build_manifest_sha256": profile["build_manifest_sha256"],
        "required_retained_bytes": R17_MINIMUM_RETAINED_BYTES,
        "nodes": 8,
        "ranks": 64,
        "storage_ready": True,
        "node_hour_ready": True,
        "rank_64_ready": True,
        "current_source_authority": current_source_authority,
    }
    if (
        any(readiness.get(key) != value for key, value in expected.items())
        or r17_available_storage_bytes(paths["root"]) < R17_MINIMUM_RETAINED_BYTES
    ):
        raise ValueError("R17 readiness evidence current-state bindings differ")
    if provenance.get("r17_readiness_evidence_sha256") != readiness_sha:
        raise ValueError("latest recost does not bind the fixed R17 readiness bytes")

    review_path = readiness_path.with_name(
        f"{readiness_path.name}.independent_review.json"
    )
    review, review_sha = read_controller_publication_json(
        review_path, "R17 readiness independent review", mode=0o444
    )
    review = require_r17_exact_keys(
        review,
        {
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate",
        },
        "R17 readiness independent review",
    )
    reviewed = require_r17_utc(review["reviewed_utc"], "R17 readiness review time")
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-r17-readiness-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or review["reviewer"] != readiness_reviewer
        or review["candidate"] != {"path": str(readiness_path), "sha256": readiness_sha}
        or reviewed < generated
        or reviewed > recost_generated
        or reviewed > now + R17_AUTHORIZATION_FUTURE_SKEW
    ):
        raise ValueError("R17 readiness independent review differs")

    audit_path = readiness_path.with_name(f"{readiness_path.name}.publication_audit.json")
    audit, audit_sha = read_controller_publication_json(
        audit_path, "R17 readiness publication audit", mode=0o444
    )
    audit = require_r17_exact_keys(
        audit,
        {
            "schema_version", "record_type", "execution_epoch", "published_utc",
            "artifact", "independent_review", "authority",
        },
        "R17 readiness publication audit",
    )
    published = require_r17_utc(audit["published_utc"], "R17 readiness publication time")
    if (
        audit["schema_version"] != 1
        or audit["record_type"] != "stage-i-r17-readiness-publication-audit"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"] != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or published < reviewed
        or published > recost_generated
        or published > now + R17_AUTHORIZATION_FUTURE_SKEW
    ):
        raise ValueError("R17 readiness publication audit differs or broadens authority")
    require_r17_declared_publication(
        audit["artifact"], readiness_path, readiness_sha,
        "R17 readiness publication artifact",
    )
    require_r17_declared_publication(
        audit["independent_review"], review_path, review_sha,
        "R17 readiness publication review",
    )
    qualification_expanded = validate_r17_operational_qualification(
        paths["root"],
        readiness,
        profile,
        readiness_reviewer,
        authorization_time=recost_generated,
    )
    embedded = require_r17_exact_keys(
        recost.get("r17_readiness"),
        {*readiness, "operational_qualification_evidence", "publication_chain"},
        "latest recost embedded R17 readiness",
    )
    if (
        embedded["operational_qualification_evidence"] != qualification_expanded
        or any(embedded.get(key) != value for key, value in readiness.items())
    ):
        raise ValueError("latest recost embedded R17 readiness differs")
    publication_chain = embedded.get("publication_chain")
    if not isinstance(publication_chain, dict) or publication_chain != {
        "independent_review_sha256": review_sha,
        "publication_audit_sha256": audit_sha,
        "reviewed_by": readiness_reviewer,
        "published_utc": published.isoformat(),
    }:
        raise ValueError("latest recost R17 readiness publication chain differs")
    return readiness, {
        "path": str(readiness_path),
        "sha256": readiness_sha,
        "independent_review_path": str(review_path),
        "independent_review_sha256": review_sha,
        "publication_audit_path": str(audit_path),
        "publication_audit_sha256": audit_sha,
    }


def expected_f119_predecessor_recost(root: Path) -> dict[str, object]:
    """Return the exact F119-only failed-F117 supersession identity."""

    root = root.expanduser().resolve()
    legacy = root / "accounting" / LEGACY_F114_RECOST_NAME
    failed = {
        "checkpoint": "F-117",
        **{
            key: {
                "path": relative.as_posix(),
                "sha256": F117_FAILED_ATTEMPT_SHA256[key],
            }
            for key, relative in F117_FAILED_ATTEMPT_RELATIVES.items()
        },
        "status": "authenticated-unpromoted-failed-attempt",
    }
    return {
        "path": str(legacy),
        "artifact_name": LEGACY_F114_RECOST_NAME,
        "sha256": LEGACY_F114_RECOST_SHA256,
        "publication_audit_path": str(
            legacy.with_name(f"{LEGACY_F114_RECOST_NAME}.publication_audit.json")
        ),
        "publication_audit_sha256": LEGACY_F114_PUBLICATION_AUDIT_SHA256,
        "independent_review_path": None,
        "independent_review_sha256": None,
        "checkpoint": "F-114",
        "generated_utc": LEGACY_F114_GENERATED_UTC,
        "published_utc": LEGACY_F114_PUBLISHED_UTC,
        "bootstrap": F119_LEGACY_F114_BOOTSTRAP,
        "superseded_failed_attempt": failed,
    }


def require_f119_recost_supersession(
    root: Path, recost_path: Path, recost: dict[str, object]
) -> None:
    """Authenticate the exact F119-only recovery from the failed F117 attempt."""

    root = root.expanduser().resolve()
    recost_path = recost_path.expanduser().resolve()
    predecessor = recost.get("predecessor_recost")
    path_is_f117 = recost_path.name == F117_RECOST_ARTIFACT_NAME
    declares_f117 = (
        recost.get("checkpoint") == "F-117"
        or recost.get("artifact_name") == F117_RECOST_ARTIFACT_NAME
    )
    if path_is_f117 or declares_f117:
        raise ValueError("the retained failed F117 attempt must never be promoted or consumed")

    path_is_f119 = recost_path.name == F119_RECOST_ARTIFACT_NAME
    declares_f119 = (
        recost.get("checkpoint") == "F-119"
        or recost.get("artifact_name") == F119_RECOST_ARTIFACT_NAME
    )
    if path_is_f119 != (
        recost.get("checkpoint") == "F-119"
        and recost.get("artifact_name") == F119_RECOST_ARTIFACT_NAME
    ):
        raise ValueError("F119 recost path, checkpoint, and artifact identity differ")

    legacy_claim = isinstance(predecessor, dict) and (
        predecessor.get("checkpoint") == "F-114"
        or predecessor.get("sha256") == LEGACY_F114_RECOST_SHA256
        or predecessor.get("publication_audit_sha256")
        == LEGACY_F114_PUBLICATION_AUDIT_SHA256
        or predecessor.get("bootstrap") == F119_LEGACY_F114_BOOTSTRAP
        or "superseded_failed_attempt" in predecessor
    )
    if not path_is_f119:
        if declares_f119 or legacy_claim:
            raise ValueError("the failed-F117 legacy-F114 supersession is reserved for F119")
        return

    expected = expected_f119_predecessor_recost(root)
    if predecessor != expected:
        raise ValueError("F119 predecessor recost or failed-F117 quartet differs")
    for relative in F117_FORBIDDEN_PROMOTION_RELATIVES:
        if os.path.lexists(root / relative):
            raise ValueError("F119 requires the retained failed F117 attempt to remain unpromoted")
    for key, relative in F117_FAILED_ATTEMPT_RELATIVES.items():
        _, digest = read_controller_publication_json(
            root / relative, f"retained failed F117 {key}", mode=0o644
        )
        if digest != F117_FAILED_ATTEMPT_SHA256[key]:
            raise ValueError(f"retained failed F117 {key} bytes differ")


def latest_published_r17_recost(paths: dict[str, Path], now: datetime
                               ) -> tuple[Path, str, dict[str, object], Path, str]:
    """Authenticate the latest independently reviewed recost publication."""

    publications = []
    accounting_entries = trusted_directory_entries(
        paths["accounting"], "F117-era accounting publication store"
    )
    for name in accounting_entries:
        audit_path = paths["accounting"] / name
        match = R17_RECOST_PUBLICATION_AUDIT_PATTERN.fullmatch(audit_path.name)
        if match is None:
            continue
        if match.group(1) == "117":
            raise ValueError("the retained failed F117 attempt must never be promoted")
        audit, audit_sha = read_controller_publication_json(
            audit_path, f"recost publication audit {audit_path.name}", mode=0o444
        )
        if audit.get("record_type") != "stage-i-recost-recommendation-publication-audit":
            continue
        publications.append((
            require_r17_utc(
                audit.get("published_utc"),
                f"recost publication audit {audit_path.name} time",
            ),
            audit_path,
            audit_sha,
            audit,
        ))
    if not publications:
        raise ValueError("R17 preparation requires a promoted recost publication")
    latest_time = max(item[0] for item in publications)
    latest = [item for item in publications if item[0] == latest_time]
    if len(latest) != 1:
        raise ValueError("latest promoted recost publication is ambiguous")
    published, audit_path, audit_sha, audit = latest[0]
    if published > now + R17_AUTHORIZATION_FUTURE_SKEW:
        raise ValueError("latest promoted recost publication is implausibly future")
    artifact_path = audit_path.with_name(
        audit_path.name.removesuffix(".publication_audit.json")
    )
    review_path = artifact_path.with_name(f"{artifact_path.name}.independent_review.json")
    recost, recost_sha = read_controller_publication_json(
        artifact_path, "latest promoted recost", mode=0o444
    )
    require_f119_recost_supersession(paths["root"], artifact_path, recost)
    review, review_sha = read_controller_publication_json(
        review_path, "latest promoted recost independent review", mode=0o444
    )
    base_audit_keys = {
        "schema_version", "record_type", "execution_epoch", "published_utc",
        "artifact", "independent_review", "authority",
    }
    full_audit_keys = base_audit_keys | {
        "transaction_id", "recost_recommendations", "counts", "generator",
        "scheduler_evidence", "source_bundle", "stage_i_helper", "utility",
        "forensic_copy", "publication", "generalized_publication_context",
    }
    if (
        not isinstance(audit, dict)
        or (set(audit) != base_audit_keys and set(audit) != full_audit_keys)
    ):
        raise ValueError("latest promoted recost publication audit schema differs")
    if (
        audit["schema_version"] != 1
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"] != {
            "action_authority": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("latest promoted recost publication audit broadens authority")
    require_r17_declared_publication(
        audit["artifact"], artifact_path, recost_sha,
        "latest promoted recost publication artifact",
    )
    require_r17_declared_publication(
        audit["independent_review"], review_path, review_sha,
        "latest promoted recost publication review",
    )
    if set(audit) == full_audit_keys and (
        audit["recost_recommendations"] != recost.get("recommendations")
        or audit["publication"] != R17_RECOST_PUBLICATION_METHOD
        or not require_r17_nonempty_string(
            audit["transaction_id"], "latest promoted recost transaction ID"
        )
    ):
        raise ValueError("latest promoted recost full publication audit differs")
    review = require_r17_exact_keys(
        review,
        {
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate", "scope",
        },
        "latest promoted recost independent review",
    )
    reviewer = require_r17_exact_keys(
        review["reviewer"], {"agent_id", "independent_from_generator"},
        "latest promoted recost reviewer",
    )
    review_scope = review["scope"]
    if (
        review_scope != {"non_authorizing": True}
    ):
        raise ValueError("latest promoted recost review scope differs")
    reviewed = require_r17_utc(review["reviewed_utc"], "latest promoted recost review time")
    generated = require_r17_utc(recost.get("generated_utc"), "latest recost generation")
    expires = require_r17_utc(recost.get("expires_utc"), "latest recost expiry")
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-recost-recommendation-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or review["candidate"] != {"path": str(artifact_path), "sha256": recost_sha}
        or reviewer["independent_from_generator"] is not True
        or not require_r17_nonempty_string(
            reviewer["agent_id"], "latest promoted recost reviewer agent ID"
        )
        or generated > reviewed
        or reviewed > published
        or generated > now + R17_AUTHORIZATION_FUTURE_SKEW
        or now - generated > R17_AUTHORIZATION_MAX_LIFETIME
        or expires <= now
        or expires <= generated
        or expires - generated > R17_AUTHORIZATION_MAX_LIFETIME
    ):
        raise ValueError("latest promoted recost review, chronology, or freshness differs")
    return artifact_path, recost_sha, recost, audit_path, audit_sha


def require_fresh_r12_profile_contract(
    profile: dict[str, object],
    parent_segment: dict[str, object] | None,
    restart: Path | None,
    time_tlim_target: float | None,
) -> None:
    """Require the exact reviewed fresh R12 rerun before later continuations."""

    if profile.get("case_id") != "R12" or parent_segment is not None:
        return
    expected_null = {
        "parent_job_id": None,
        "parent_result": None,
        "parent_segment": None,
        "restart_file": None,
        "restart_file_sha256": None,
        "restart_time": None,
    }
    if (
        profile.get("segment") != R12_FRESH_RERUN_SEGMENT
        or profile.get("nodes") != R12_FRESH_RERUN_NODES
        or profile.get("ranks_per_node") != R12_FRESH_RERUN_RANKS_PER_NODE
        or (
            int(profile.get("nodes", 0))
            * int(profile.get("ranks_per_node", 0))
            != R12_FRESH_RERUN_TOTAL_RANKS
        )
        or profile.get("walltime") != R12_FRESH_RERUN_WALLTIME
        or profile.get("athena_walltime") != R12_FRESH_RERUN_ATHENA_WALLTIME
        or restart is not None
        or time_tlim_target is None
        or abs(float(time_tlim_target) - R12_FRESH_RERUN_TARGET) > 1.0e-12
        or any(profile.get(key) != value for key, value in expected_null.items())
    ):
        raise ValueError(
            "R12 must restart fresh as s01_rankio_t0_t0p12 on 4 nodes / "
            "32 ranks with Slurm 02:00:00 and Athena 01:50:00 from t=0 "
            "to t=0.12 without a parent or restart"
        )


def require_promoted_profile_for_prepare(
    paths: dict[str, Path],
    args: argparse.Namespace,
    *,
    source_dir: Path,
    matrix_path: Path,
    input_path: Path,
    input_revision: str,
    utility_provenance: dict[str, object],
    build_provenance: dict[str, str],
    build_manifest: Path,
    qualification_approval: dict[str, object] | None,
    bundle_provenance: dict[str, object] | None,
    current_source_authority: dict[str, object] | None,
    restart: Path | None,
    parent_segment: dict[str, object] | None,
    time_tlim_target: float | None,
    now: datetime | None = None,
) -> dict[str, object] | None:
    """Bind every canonical R03-R16 prepare to its reviewed F117-era profile."""

    if (
        args.case_id not in CONCURRENT_CASE_IDS
        or paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve()
    ):
        return None
    if (
        qualification_approval is None
        or bundle_provenance is None
        or current_source_authority is None
    ):
        raise ValueError("canonical prepare requires promoted-profile provenance")
    now = datetime.now(timezone.utc) if now is None else now.astimezone(timezone.utc)
    recost_path, recost_sha, recost, audit_path, audit_sha = latest_published_r17_recost(
        paths, now
    )
    match = R17_RECOST_PUBLICATION_AUDIT_PATTERN.fullmatch(audit_path.name)
    if (
        match is None
        or int(match.group(1)) < 117
        or recost.get("schema_version") != 2
        or recost.get("record_type") != "stage-i-recost-recommendation-evidence"
        or recost.get("checkpoint") != f"F-{match.group(1)}"
        or recost.get("artifact_name") != recost_path.name
        or recost.get("execution_epoch") != EXECUTION_EPOCH
        or recost.get("authority")
        != {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or recost.get("publication_requirements")
        != {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        }
    ):
        raise ValueError("latest promoted F117-era recost identity or authority differs")
    recommendations = recost.get("recommendations")
    required = {
        "mode", "authorizing", "recommended_next_profiles",
        "bounded_concurrency", "controller_consumption_state",
        "non_authorizing_reason",
    }
    if (
        not isinstance(recommendations, dict)
        or not required.issubset(recommendations)
        or not set(recommendations).issubset(
            required | {"sole_next_segment_recommendation"}
        )
    ):
        raise ValueError("latest promoted F117-era profile schema differs")
    profiles = recommendations["recommended_next_profiles"]
    bounded = recommendations["bounded_concurrency"]
    if (
        recommendations["mode"] not in {"sole-next-profile", "bounded-wave"}
        or recommendations["authorizing"] is not False
        or not isinstance(profiles, list)
        or not profiles
        or len(profiles) > MAX_ACTIVE_STAGE_I_SEGMENTS
        or not isinstance(bounded, dict)
        or bounded
        != {
            "max_active_segments": MAX_ACTIVE_STAGE_I_SEGMENTS,
            "max_wave_nodes": sum(
                int(profile.get("nodes", 0))
                for profile in profiles if isinstance(profile, dict)
            ),
            "r17_exclusive_and_last": True,
        }
        or int(bounded["max_wave_nodes"]) > MAX_ACTIVE_STAGE_I_NODES
        or not require_r17_nonempty_string(
            recommendations["controller_consumption_state"],
            "latest promoted F117-era controller state",
        )
        or not require_r17_nonempty_string(
            recommendations["non_authorizing_reason"],
            "latest promoted F117-era non-authorizing reason",
        )
    ):
        raise ValueError("latest promoted F117-era bounded profile differs")
    if (
        recommendations["mode"] == "sole-next-profile" and len(profiles) != 1
    ) or (
        recommendations["mode"] == "bounded-wave" and len(profiles) < 2
    ):
        raise ValueError("latest promoted F117-era recommendation mode differs")
    validated_profiles = [
        require_r17_exact_keys(profile, R17_PROFILE_KEYS, "promoted Stage I profile")
        for profile in profiles
    ]
    profile_cases = [profile["case_id"] for profile in validated_profiles]
    if (
        len(profile_cases) != len(set(profile_cases))
        or any(
            case_id not in CONCURRENT_CASE_IDS | {R17_CASE_ID}
            for case_id in profile_cases
        )
    ):
        raise ValueError("latest promoted F117-era profiles duplicate or broaden cases")
    if (
        recommendations["mode"] == "sole-next-profile"
        and recommendations.get("sole_next_segment_recommendation")
        != {
            key: validated_profiles[0][key]
            for key in R17_SOLE_PROFILE_KEYS
        }
    ):
        raise ValueError("latest promoted F117-era sole-profile compatibility differs")
    selected = [
        profile for profile in validated_profiles
        if profile.get("case_id") == args.case_id
        and profile.get("segment") == args.segment
    ]
    if len(selected) != 1:
        raise ValueError("latest promoted F117-era recost lacks one exact requested profile")
    profile = selected[0]
    parent_job_id = None
    parent_result = None
    parent_segment_id = None
    restart_time = None
    restart_sha = None
    if parent_segment is not None:
        parent_job_id = read_manifest(
            Path(str(parent_segment["manifest"]))
        ).get("job_id")
        parent_result = parent_segment["result"]
        parent_segment_id = parent_segment["segment"]
        restart_time = parent_segment["final_time"]
        if restart is None:
            raise ValueError("promoted continuation profile lacks a restart")
        restart_sha = sha256(restart)
    require_fresh_r12_profile_contract(
        profile, parent_segment, restart, time_tlim_target
    )
    try:
        input_relative = input_path.relative_to(source_dir).as_posix()
    except ValueError as error:
        raise ValueError("promoted profile input is outside its source directory") from error
    expected = {
        "acceptance_criterion": args.acceptance_criterion,
        "athena_walltime": args.athena_walltime,
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": r17_directory_inventory_sha256(
            build_manifest, "promoted profile build manifest"
        ),
        "case_id": args.case_id,
        "controller_walltime_max_seconds": MAX_SEGMENT_SECONDS,
        "cpus_per_task": args.cpus_per_task,
        "executable": str(Path(args.executable).expanduser().resolve()),
        "executable_revision": build_provenance["revision"],
        "executable_sha256": build_provenance["sha256"],
        "input_file": input_relative,
        "input_revision": input_revision,
        "input_sha256": sha256(input_path),
        "nodes": args.nodes,
        "output_layout": "rank-local",
        "segment": args.segment,
        "parent_job_id": parent_job_id,
        "parent_result": parent_result,
        "parent_segment": parent_segment_id,
        "restart_file": str(restart) if restart else None,
        "restart_file_sha256": restart_sha,
        "restart_time": restart_time,
        "ranks_per_node": args.ranks_per_node,
        "time_tlim_target": time_tlim_target,
        "walltime": args.walltime,
        "source_bundle": bundle_provenance["path"],
        "source_bundle_sha256": bundle_provenance["sha256"],
    }
    if any(profile.get(key) != value for key, value in expected.items()):
        raise ValueError("promoted F117-era profile differs from prepare arguments")
    if (
        not isinstance(profile["acceptance_policy"], str)
        or not profile["acceptance_policy"]
        or not isinstance(profile["recommendation_basis"], dict)
        or not profile["recommendation_basis"]
        or require_r17_integer(
            profile["estimated_storage_bytes"], "promoted profile storage", 1
        ) < 1
    ):
        raise ValueError("promoted F117-era profile lacks reviewed scientific policy")
    tooling_revision = require_source_authority_revision(
        utility_provenance.get("revision"), "current F117 tooling revision"
    )
    provenance = recost.get("provenance")
    if (
        not isinstance(provenance, dict)
        or provenance.get("source_authority") != current_source_authority
        or provenance.get("stage_i_helper_sha256") != utility_provenance.get("sha256")
        or provenance.get("stage_i_helper_revision") != tooling_revision
        or provenance.get("matrix_sha256") != sha256(matrix_path)
        or provenance.get("matrix_revision") != tooling_revision
        or provenance.get("source_bundle_sha256") != bundle_provenance.get("sha256")
        or provenance.get("qualification_approval_sha256")
        != qualification_approval.get("sha256")
    ):
        raise ValueError("promoted F117-era profile provenance is stale")
    authority = {
        "recost_path": str(recost_path),
        "recost_sha256": recost_sha,
        "recost_publication_audit_path": str(audit_path),
        "recost_publication_audit_sha256": audit_sha,
        "profile_sha256": stable_json_sha256(profile),
        "profile": profile,
        "current_source_authority": current_source_authority,
    }
    return authority


def require_r17_readiness_for_prepare(
    paths: dict[str, Path],
    args: argparse.Namespace,
    reservations: list[dict[str, object]],
    *,
    source_dir: Path,
    matrix_path: Path,
    input_path: Path,
    input_revision: str,
    utility_provenance: dict[str, object],
    build_provenance: dict[str, str],
    build_manifest: Path,
    qualification_approval: dict[str, object] | None,
    bundle_provenance: dict[str, object] | None,
    current_source_authority: dict[str, object] | None,
    restart: Path | None,
    parent_segment: dict[str, object] | None,
    time_tlim_target: float | None,
    now: datetime | None = None,
    reservation_snapshot_sha256: str | None = None,
    excluded_manifest_path: Path | None = None,
) -> dict[str, object] | None:
    """Require the strong published readiness chain before canonical R17 prepare."""

    if args.case_id != R17_CASE_ID:
        return None
    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return None
    now = datetime.now(timezone.utc) if now is None else now.astimezone(timezone.utc)
    if excluded_manifest_path is None:
        lineage_sha, manifest_bindings = require_clean_r17_predecessor_state(
            paths,
            reservations,
            source_authority_transactions_authenticated=(
                current_source_authority is not None
            ),
        )
    else:
        lineage_sha, manifest_bindings = require_clean_r17_predecessor_state(
            paths,
            reservations,
            excluded_manifest_path=excluded_manifest_path,
            source_authority_transactions_authenticated=(
                current_source_authority is not None
            ),
        )
    if reservation_snapshot_sha256 is None:
        reservation_snapshot_sha256 = sha256(paths["reservations"])
    elif (
        require_r17_sha256(
            reservation_snapshot_sha256, "R17 predecessor reservation snapshot SHA-256"
        )
        != stable_json_sha256(reservations)
    ):
        raise ValueError("R17 predecessor reservation snapshot checksum differs")
    recost_path, recost_sha, recost, audit_path, audit_sha = latest_published_r17_recost(
        paths, now
    )
    recost = require_r17_exact_keys(
        recost,
        {
            "schema_version", "record_type", "checkpoint", "artifact_name",
            "execution_epoch", "generated_utc", "expires_utc", "requested_by",
            "scope", "predecessor_recost", "authority", "publication_requirements",
            "recommendations", "barrier", "budget", "storage", "ledger",
            "reservations", "manifests", "r17_readiness", "promoted_f113",
            "reconcile", "provenance",
        },
        "latest promoted recost",
    )
    match = R17_RECOST_PUBLICATION_AUDIT_PATTERN.fullmatch(audit_path.name)
    if (
        match is None
        or recost["schema_version"] != 2
        or recost["record_type"] != "stage-i-recost-recommendation-evidence"
        or recost["checkpoint"] != f"F-{match.group(1)}"
        or recost["artifact_name"] != recost_path.name
        or recost["execution_epoch"] != EXECUTION_EPOCH
        or recost["authority"] != {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
        or recost["publication_requirements"] != {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        }
    ):
        raise ValueError("latest promoted recost identity or authority differs")
    recommendations = require_r17_exact_keys(
        recost["recommendations"],
        {
            "mode", "authorizing", "recommended_next_profiles",
            "bounded_concurrency", "controller_consumption_state",
            "sole_next_segment_recommendation", "non_authorizing_reason",
        },
        "latest promoted recost recommendations",
    )
    concurrency = require_r17_exact_keys(
        recommendations["bounded_concurrency"],
        {"max_active_segments", "max_wave_nodes", "r17_exclusive_and_last"},
        "latest promoted recost bounded concurrency",
    )
    profiles = recommendations["recommended_next_profiles"]
    if (
        recommendations["mode"] != "sole-next-profile"
        or recommendations["authorizing"] is not False
        or recommendations["controller_consumption_state"]
        != "non-authorizing until exact independent review and publication audit"
        or recommendations["non_authorizing_reason"]
        != (
            "A recost request and generated candidate are evidence only; neither may "
            "authorize prepare, submit, scheduler, controller, or canonical mutations."
        )
        or concurrency["max_active_segments"] != MAX_ACTIVE_STAGE_I_SEGMENTS
        or concurrency["max_wave_nodes"] != 8
        or concurrency["r17_exclusive_and_last"] is not True
        or not isinstance(profiles, list)
        or len(profiles) != 1
    ):
        raise ValueError("latest promoted recost is not an exact non-broadening R17 profile")
    profile = require_r17_exact_keys(profiles[0], R17_PROFILE_KEYS, "sole R17 profile")
    if (
        recommendations["sole_next_segment_recommendation"]
        != {key: profile[key] for key in R17_SOLE_PROFILE_KEYS}
    ):
        raise ValueError("latest promoted recost sole R17 compatibility profile differs")

    build_manifest_sha = r17_directory_inventory_sha256(
        build_manifest, "R17 build manifest"
    )
    if bundle_provenance is None or qualification_approval is None:
        raise ValueError("R17 preparation requires retained source and qualification provenance")
    if current_source_authority is None:
        raise ValueError("R17 preparation requires authenticated F118 source authority")
    current_source_authority = require_r17_exact_keys(
        current_source_authority,
        {
            "checkpoint", "evidence", "provenance_review", "plasma_review",
            "publication_audit", "final_source_bundle",
        },
        "authenticated F118 current-source authority",
    )
    authority_bundle = require_r17_exact_keys(
        current_source_authority["final_source_bundle"],
        {"path", "sha256", "verified_revisions"},
        "authenticated F118 final source bundle",
    )
    authority_bundle_relative_text = require_r17_nonempty_string(
        authority_bundle["path"], "authenticated F118 final source bundle path"
    )
    authority_bundle_relative = Path(authority_bundle_relative_text)
    authority_revisions = authority_bundle["verified_revisions"]
    if (
        current_source_authority["checkpoint"] != "F-118"
        or authority_bundle_relative.is_absolute()
        or ".." in authority_bundle_relative.parts
        or authority_bundle_relative.as_posix() != authority_bundle_relative_text
        or (paths["root"] / authority_bundle_relative).resolve()
        != Path(str(bundle_provenance.get("path", ""))).resolve()
        or authority_bundle["sha256"] != bundle_provenance.get("sha256")
        or not isinstance(authority_revisions, list)
        or any(
            not isinstance(revision, str)
            or GIT_REVISION_PATTERN.fullmatch(revision) is None
            for revision in authority_revisions
        )
    ):
        raise ValueError("authenticated current-source authority identity differs")
    parent_job_id = None
    parent_result = None
    parent_segment_id = None
    restart_time = None
    restart_sha = None
    if parent_segment is not None:
        parent_manifest = read_manifest(Path(str(parent_segment["manifest"])))
        parent_job_id = parent_manifest.get("job_id")
        parent_result = parent_segment["result"]
        parent_segment_id = parent_segment["segment"]
        restart_time = parent_segment["final_time"]
        if restart is None:
            raise ValueError("R17 continuation profile lacks a restart")
        restart_sha = sha256(restart)
    try:
        input_relative = input_path.relative_to(source_dir).as_posix()
    except ValueError as error:
        raise ValueError("R17 input is outside its source directory") from error
    expected_profile = {
        "acceptance_criterion": args.acceptance_criterion,
        "athena_walltime": args.athena_walltime,
        "build_manifest": str(build_manifest),
        "build_manifest_sha256": build_manifest_sha,
        "case_id": R17_CASE_ID,
        "controller_walltime_max_seconds": MAX_SEGMENT_SECONDS,
        "cpus_per_task": args.cpus_per_task,
        "executable": str(Path(args.executable).expanduser().resolve()),
        "executable_revision": build_provenance["revision"],
        "executable_sha256": build_provenance["sha256"],
        "input_file": input_relative,
        "input_revision": input_revision,
        "input_sha256": sha256(input_path),
        "nodes": args.nodes,
        "output_layout": "rank-local",
        "segment": args.segment,
        "parent_job_id": parent_job_id,
        "parent_result": parent_result,
        "parent_segment": parent_segment_id,
        "restart_file": str(restart) if restart else None,
        "restart_file_sha256": restart_sha,
        "restart_time": restart_time,
        "ranks_per_node": args.ranks_per_node,
        "time_tlim_target": time_tlim_target,
        "walltime": args.walltime,
        "source_bundle": bundle_provenance["path"],
        "source_bundle_sha256": bundle_provenance["sha256"],
    }
    if any(profile.get(key) != value for key, value in expected_profile.items()):
        raise ValueError("sole published R17 profile differs from prepare arguments or provenance")
    if (
        not isinstance(profile["acceptance_policy"], str)
        or not profile["acceptance_policy"]
        or not isinstance(profile["recommendation_basis"], dict)
        or not profile["recommendation_basis"]
        or require_r17_integer(
            profile["estimated_storage_bytes"], "R17 estimated storage bytes", 1
        ) < R17_MINIMUM_RETAINED_BYTES
    ):
        raise ValueError("sole published R17 profile lacks reviewed policy or storage")

    provenance = require_r17_exact_keys(
        recost["provenance"],
        {
            "request_sha256", "request_independent_review_sha256", "generator_sha256",
            "generator_revision", "stage_i_helper_sha256", "stage_i_helper_revision",
            "matrix_sha256", "matrix_revision", "source_bundle_sha256",
            "source_bundle_verified_revisions", "source_authority",
            "qualification_approval_sha256", "ceiling_evidence_sha256",
            "ceiling_publication_audit_sha256", "f113_historical_helper_revision",
            "f113_historical_helper_sha256", "storage_evidence_sha256",
            "reconciliation_sha256", "ledger_sha256", "reservations_sha256",
            "scheduler_evidence", "predecessor_recost_sha256",
            "predecessor_recost_independent_review_sha256",
            "predecessor_recost_publication_audit_sha256",
            "authenticated_lineages_sha256", "computed_projection_sha256",
            "r17_readiness_evidence_sha256", "scheduler_sha256",
        },
        "latest promoted recost provenance",
    )
    tooling_revision = require_source_authority_revision(
        utility_provenance["revision"], "current R17 tooling revision"
    )
    current_bindings = {
        "generator_sha256": source_authority_committed_sha256(
            tooling_revision,
            "scripts/frontier/cgl_lf_stage_i_recost.py",
            "current R17 recost generator",
        ),
        "generator_revision": tooling_revision,
        "stage_i_helper_sha256": utility_provenance["sha256"],
        "stage_i_helper_revision": tooling_revision,
        "matrix_sha256": sha256(matrix_path),
        "matrix_revision": tooling_revision,
        "source_bundle_sha256": bundle_provenance["sha256"],
        "qualification_approval_sha256": qualification_approval["sha256"],
        "ledger_sha256": sha256(paths["ledger"]),
        "reservations_sha256": reservation_snapshot_sha256,
        "authenticated_lineages_sha256": lineage_sha,
    }
    if any(provenance.get(key) != value for key, value in current_bindings.items()):
        raise ValueError("latest promoted recost current controller/state bindings differ")
    recost_revisions = provenance["source_bundle_verified_revisions"]
    required_revisions = bundle_provenance["verified_revisions"]
    if (
        provenance["source_authority"] != current_source_authority
        or not isinstance(recost_revisions, list)
        or len(recost_revisions) != len(set(recost_revisions))
        or any(
            not isinstance(revision, str)
            or GIT_REVISION_PATTERN.fullmatch(revision) is None
            for revision in recost_revisions
        )
        or not isinstance(required_revisions, list)
        or recost_revisions != authority_revisions
        or any(revision not in recost_revisions for revision in required_revisions)
    ):
        raise ValueError("latest promoted recost F118 source authority or bundle differs")
    current_rows = read_ledger(paths)
    ledger = recost["ledger"]
    recost_reservations = recost["reservations"]
    manifests = recost["manifests"]
    if (
        not isinstance(ledger, dict)
        or ledger.get("sha256") != current_bindings["ledger_sha256"]
        or ledger.get("rows") != len(current_rows)
        or not isinstance(recost_reservations, dict)
        or recost_reservations.get("sha256") != current_bindings["reservations_sha256"]
        or recost_reservations.get("rows") != len(reservations)
        or recost_reservations.get("active") != 0
        or not isinstance(manifests, dict)
        or manifests.get("rows") != len(manifest_bindings)
        or manifests.get("bindings") != manifest_bindings
        or manifests.get("authenticated_lineages_sha256") != lineage_sha
    ):
        raise ValueError("latest promoted recost does not bind the exact drained canonical state")

    budget = require_r17_exact_keys(
        recost["budget"],
        {
            "method", "measurement_basis", "actual_stage_i_node_hours",
            "authorized_wave_reserved_node_hours",
            "actual_plus_authorized_wave_node_hours",
            "computed_remaining_stage_i_node_hours",
            "computed_stage_i_total_node_hours",
            "promoted_stage_i_envelope_node_hours", "project_ceiling_node_hours",
            "computed_stage_i_margin_node_hours", "case_breakdown",
            "storage_projection_evidence",
        },
        "latest promoted recost budget",
    )
    projection_sha = hashlib.sha256(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    ).hexdigest()
    actual = sum(
        (
            require_r17_decimal(row["actual_node_hours"], "ledger actual node-hours")
            for row in current_rows
        ),
        Decimal("0"),
    )
    authorized = Decimal(args.nodes * parse_walltime(args.walltime)) / Decimal(3600)
    actual_retained = require_r17_decimal(
        budget["actual_stage_i_node_hours"], "recost actual Stage I node-hours"
    )
    authorized_retained = require_r17_decimal(
        budget["authorized_wave_reserved_node_hours"],
        "recost authorized R17 node-hours",
    )
    committed = require_r17_decimal(
        budget["actual_plus_authorized_wave_node_hours"],
        "recost actual plus authorized node-hours",
    )
    remaining = require_r17_decimal(
        budget["computed_remaining_stage_i_node_hours"],
        "recost remaining Stage I node-hours",
    )
    projected = require_r17_decimal(
        budget["computed_stage_i_total_node_hours"],
        "recost computed Stage I total node-hours",
    )
    envelope = require_r17_decimal(
        budget["promoted_stage_i_envelope_node_hours"],
        "recost promoted Stage I envelope",
    )
    project = require_r17_decimal(
        budget["project_ceiling_node_hours"], "recost project ceiling"
    )
    margin = require_r17_decimal(
        budget["computed_stage_i_margin_node_hours"], "recost Stage I margin"
    )
    validate_r17_scoped_budget_projection(budget, current_rows)
    if (
        provenance["computed_projection_sha256"] != projection_sha
        or actual_retained != actual
        or authorized_retained != authorized
        or committed != actual + authorized
        or projected != actual + remaining
        or projected < committed
        or envelope != Decimal(str(CURRENT_STAGE_I_RESERVED_NODE_HOURS))
        or project != Decimal(str(PROJECT_BUDGET_NODE_HOURS))
        or projected > envelope
        or projected > project
        or margin != envelope - projected
    ):
        raise ValueError("latest promoted recost budget is stale or arithmetically invalid")
    storage = require_r17_exact_keys(
        recost["storage"],
        {
            "available_bytes", "retained_stage_i_bytes", "required_safety_bytes",
            "projected_authorized_wave_growth_bytes",
            "headroom_after_authorized_wave_and_safety_bytes",
        },
        "latest promoted recost storage",
    )
    available = require_r17_integer(storage["available_bytes"], "R17 available bytes")
    require_r17_integer(storage["retained_stage_i_bytes"], "R17 retained bytes")
    safety = require_r17_integer(storage["required_safety_bytes"], "R17 safety bytes", 1)
    projected_growth = require_r17_integer(
        storage["projected_authorized_wave_growth_bytes"],
        "R17 projected growth bytes",
    )
    headroom = require_r17_integer(
        storage["headroom_after_authorized_wave_and_safety_bytes"],
        "R17 storage headroom bytes",
    )
    if (
        projected_growth != profile["estimated_storage_bytes"]
        or headroom != available - safety - projected_growth
        or headroom < 0
    ):
        raise ValueError("latest promoted recost storage readiness differs")
    _, readiness_binding = validate_fixed_r17_readiness_chain(
        paths, recost, profile, lineage_sha, current_source_authority, now
    )
    return {
        "recost_path": str(recost_path),
        "recost_sha256": recost_sha,
        "recost_publication_audit_path": str(audit_path),
        "recost_publication_audit_sha256": audit_sha,
        "current_source_authority": current_source_authority,
        **readiness_binding,
    }


def require_positive_finite_float(value: object, label: str) -> float:
    """Require one positive finite command or API threshold."""

    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{label} must be positive and finite")
    return result


def positive_finite_float_arg(value: str) -> float:
    """Parse a positive finite command-line floating-point value."""

    try:
        return require_positive_finite_float(value, "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def layout(root: Path) -> dict[str, Path]:
    """Return retained current-epoch Stage I production locations."""

    accounting = root / "accounting"
    return {
        "root": root,
        "accounting": accounting,
        "ledger": accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv",
        "reservations": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json"
        ),
        "summary": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_budget_summary.md"
        ),
        "transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions"
        ),
        "qualification": (
            accounting
            / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_qualification_approval.json"
        ),
        "runs": root / "runs" / "mks24-stage-i" / EXECUTION_EPOCH,
        "logs_slurm": root / "logs" / "slurm",
    }


def require_safe_manifest_path(paths: dict[str, Path],
                               value: object) -> Path:
    """Require the exact E03 manifest location for one authorized segment."""

    declared = Path(str(value)).expanduser()
    if not declared.is_absolute():
        raise ValueError(f"manifest target must be absolute: {declared}")
    manifest_path = declared.resolve()
    try:
        relative = manifest_path.relative_to(paths["runs"].resolve())
    except ValueError as error:
        raise ValueError(
            f"manifest target is outside the E03 run store: {manifest_path}"
        ) from error
    if len(relative.parts) != 4 or relative.parts[2:] != (
        "manifest", "prepared_run.json"
    ):
        raise ValueError(f"manifest target is not an exact E03 segment path: {manifest_path}")
    case_id, segment = relative.parts[:2]
    require_authorized_case(case_id)
    require_safe_segment(segment)
    expected = paths["runs"] / case_id / segment / "manifest" / "prepared_run.json"
    if declared != expected or manifest_path != expected.resolve():
        raise ValueError(f"manifest target does not resolve exactly beneath E03: {manifest_path}")
    return manifest_path


def orphaned_segment_run_directories(paths: dict[str, Path]) -> list[Path]:
    """Return segment directories left behind before a manifest was retained."""

    orphans = []
    if not paths["runs"].is_dir():
        return orphans
    for case_dir in sorted(paths["runs"].iterdir()):
        if not case_dir.is_dir() or case_dir.name not in AUTHORIZED_CASE_IDS:
            continue
        for run_dir in sorted(case_dir.iterdir()):
            if run_dir.name == "analysis":
                continue
            if (
                run_dir.is_dir()
                and not (run_dir / "manifest" / "prepared_run.json").is_file()
            ):
                orphans.append(run_dir)
    return orphans


def require_no_orphaned_segment_runs(paths: dict[str, Path]) -> None:
    """Fail closed after an interrupted prepare leaves unaudited run content."""

    orphans = orphaned_segment_run_directories(paths)
    if orphans:
        raise ValueError(
            "Stage I interrupted-prepare cleanup is required before mutation: "
            + ", ".join(str(path) for path in orphans)
        )


def read_qualification_approval(path: Path) -> dict[str, object]:
    """Read and validate one corrected-build Frontier qualification token."""

    if not path.is_file():
        raise ValueError(f"E03 qualification approval token is absent: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"E03 qualification approval token is invalid: {path}")
    if value.get("schema_version") != 1:
        raise ValueError(f"E03 qualification approval token has wrong schema: {path}")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(f"E03 qualification approval token has wrong epoch: {path}")
    executable_sha256 = value.get("approved_executable_sha256")
    revision = value.get("approved_executable_revision")
    if (
        not isinstance(executable_sha256, str)
        or SHA256_PATTERN.fullmatch(executable_sha256) is None
    ):
        raise ValueError(
            f"E03 qualification approval token has invalid executable digest: {path}"
        )
    if (
        not isinstance(revision, str)
        or GIT_REVISION_PATTERN.fullmatch(revision) is None
    ):
        raise ValueError(
            f"E03 qualification approval token has invalid git revision: {path}"
        )
    for key in ("approved_utc", "approved_by", "review_notes"):
        if not isinstance(value.get(key), str) or not str(value[key]).strip():
            raise ValueError(
                f"E03 qualification approval token lacks nonempty {key}: {path}"
            )
    return value


def qualification_approval_status(paths: dict[str, Path]) -> dict[str, object]:
    """Return a summary-safe view of the current E03 qualification token."""

    path = paths["qualification"]
    if not path.is_file():
        return {
            "state": "pending",
            "path": str(path),
            "reason": "approval token is absent",
        }
    try:
        approval = read_qualification_approval(path)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return {
            "state": "invalid",
            "path": str(path),
            "reason": str(error),
        }
    return {
        "state": "approved",
        "path": str(path),
        "sha256": sha256(path),
        "approved_executable_sha256": approval["approved_executable_sha256"],
        "approved_executable_revision": approval["approved_executable_revision"],
        "approved_utc": approval["approved_utc"],
        "approved_by": approval["approved_by"],
    }


def require_qualification_approval(
    paths: dict[str, Path],
    executable_sha256: str,
    executable_revision: str,
    offline_local_root: bool,
) -> dict[str, object] | None:
    """Require the fixed-location token to approve this exact corrected build."""

    path = paths["qualification"]
    if offline_local_root and not path.is_file():
        return None
    approval = read_qualification_approval(path)
    if approval["approved_executable_sha256"] != executable_sha256:
        raise ValueError("E03 qualification token does not approve this executable")
    if approval["approved_executable_revision"] != executable_revision:
        raise ValueError("E03 qualification token does not approve this git revision")
    return {
        "path": str(path),
        "sha256": sha256(path),
        "execution_epoch": EXECUTION_EPOCH,
        "approved_executable_sha256": approval["approved_executable_sha256"],
        "approved_executable_revision": approval["approved_executable_revision"],
        "token": approval,
    }


def initialize(root: Path) -> dict[str, Path]:
    """Create the production layout and accounting stores."""

    with canonical_root_lock(root):
        paths = layout(root)
        for key in ("accounting", "runs", "logs_slurm", "transactions"):
            mkdir_durable(paths[key])
        if not paths["ledger"].exists():
            write_text(paths["ledger"], ledger_csv_text([]), mode=0o644)
        if not paths["reservations"].exists():
            write_json(paths["reservations"], [])
        read_ledger(paths)
        read_reservations(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        refresh_summary(paths)
        return paths


def read_ledger(paths: dict[str, Path]) -> list[dict[str, str]]:
    """Read retained Stage I allocation records."""

    with paths["ledger"].open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Stage I ledger is empty: {paths['ledger']}") from error
        if tuple(header) != LEDGER_COLUMNS:
            raise ValueError(
                f"Stage I ledger header is invalid: {paths['ledger']}: {header!r}"
            )
        rows = []
        cumulative = 0.0
        for index, row in enumerate(reader, start=2):
            if len(row) != len(LEDGER_COLUMNS):
                raise ValueError(
                    f"Stage I ledger row {index} has {len(row)} columns; "
                    f"expected {len(LEDGER_COLUMNS)}"
                )
            retained = dict(zip(LEDGER_COLUMNS, row))
            cumulative = validate_ledger_cumulative_fields(
                retained, cumulative, f"Stage I ledger row {index}"
            )
            rows.append(retained)
        return rows


def read_reservations(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Read all Stage I preparation records."""

    value = json.loads(paths["reservations"].read_text(encoding="utf-8"))
    if not isinstance(value, list):
        raise ValueError("Stage I reservation store must contain a list")
    return value


def transaction_store_entries(paths: dict[str, Path]) -> list[str]:
    """Return a final repeated scan from one bound transaction-store descriptor.

    The descriptor stays open through the last content-identity check.  POSIX
    cannot freeze a writable directory after this list-returning interface
    releases its descriptor, so every returned journal is authenticated again
    when opened.
    """

    directory = paths["transactions"].absolute()
    label = "Stage I transaction store"
    require_owner_symlink_free_path(directory.parent, label)
    try:
        directory.lstat()
    except FileNotFoundError:
        return []
    except OSError as error:
        raise ValueError(f"{label} is unavailable") from error
    require_owner_symlink_free_path(directory, label)
    try:
        directory_fd = os.open(
            directory,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as error:
        raise ValueError(f"{label} is unavailable") from error
    try:
        directory_profile = require_directory_descriptor_binding(
            directory, directory_fd, label
        )
        if (
            directory_profile.st_uid != os.geteuid()
            or stat.S_IMODE(directory_profile.st_mode) & 0o022
        ):
            raise ValueError(f"{label} is not owner controlled")
        entries, _scan_profile = stable_bound_directory_entries(
            directory, directory_fd, directory_profile, label
        )
        temporaries = {
            name for name in entries
            if transaction_metadata_temporary_target(name) is not None
        }
        if temporaries:
            require_canonical_mutation_lock(directory)
        for name in sorted(temporaries):
            target = transaction_metadata_temporary_target(name)
            if target is None:
                raise ValueError(f"{label} temporary identity changed")
            retire_metadata_temporary(
                directory,
                directory_fd,
                name,
                directory_profile,
                f"Stage I transaction metadata temporary for {target}",
            )
        entries, final_profile = stable_bound_directory_entries(
            directory, directory_fd, directory_profile, label
        )
        confirmed_entries, confirmed_profile = stable_bound_directory_entries(
            directory, directory_fd, directory_profile, label
        )
        if (
            entries != confirmed_entries
            or directory_content_identity(final_profile)
            != directory_content_identity(confirmed_profile)
        ):
            raise ValueError(f"{label} final contents changed while returning result")
        require_directory_descriptor_binding(
            directory, directory_fd, label, directory_profile
        )
        returned_profile = os.fstat(directory_fd)
        if (
            directory_content_identity(confirmed_profile)
            != directory_content_identity(returned_profile)
        ):
            raise ValueError(f"{label} final contents changed while returning result")
        return confirmed_entries
    finally:
        os.close(directory_fd)


def pending_transaction_paths(paths: dict[str, Path]) -> list[Path]:
    """Return durable metadata transitions awaiting completion."""

    directory = paths["transactions"]
    entries = transaction_store_entries(paths)
    if any(
        Path(name).name != name
        or name in {".", ".."}
        or not name.endswith(".json")
        for name in entries
    ):
        raise ValueError("Stage I transaction store contains an unexpected entry")
    return [directory / name for name in entries]


def require_no_pending_transactions(paths: dict[str, Path]) -> None:
    """Fail closed while a prior cross-file transition awaits recovery."""

    pending = pending_transaction_paths(paths)
    if pending:
        raise ValueError(
            "Stage I metadata recovery is required before mutation: "
            + ", ".join(str(path) for path in pending)
        )


def validate_reservation_record(paths: dict[str, Path],
                                reservation: object) -> dict[str, object]:
    """Require one journal-controlled reservation to use the retained schema."""

    if not isinstance(reservation, dict):
        raise ValueError("transaction reservation record is not an object")
    columns = frozenset(reservation)
    if not RESERVATION_REQUIRED_COLUMNS.issubset(columns):
        missing = sorted(RESERVATION_REQUIRED_COLUMNS - columns)
        raise ValueError(f"transaction reservation lacks columns: {missing}")
    supported = RESERVATION_REQUIRED_COLUMNS | RESERVATION_OPTIONAL_COLUMNS
    if not columns.issubset(supported):
        raise ValueError(
            "transaction reservation has unsupported columns: "
            f"{sorted(columns - supported)}"
        )
    if reservation.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("transaction reservation has wrong execution epoch")
    manifest_path = require_safe_manifest_path(paths, reservation["manifest"])
    case_id = str(reservation["case_id"])
    segment = require_safe_segment(str(reservation["segment"]))
    if (
        case_id != manifest_path.parents[2].name
        or segment != manifest_path.parents[1].name
    ):
        raise ValueError("transaction reservation identity differs from manifest path")
    require_authorized_case(case_id)
    try:
        nodes = int(reservation["nodes"])
        requested_seconds = parse_walltime(str(reservation["requested_walltime"]))
        reserved = float(reservation["reserved_node_hours"])
    except (TypeError, ValueError) as error:
        raise ValueError("transaction reservation allocation is invalid") from error
    if (
        nodes < 1
        or requested_seconds <= 0
        or not math.isfinite(reserved)
        or reserved <= 0.0
    ):
        raise ValueError("transaction reservation allocation must be positive")
    if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
        require_case_node_count(case_id, nodes)
    if abs(reserved - node_hours(nodes, requested_seconds)) > 5.0e-12:
        raise ValueError("transaction reservation node-hours differ from allocation")
    parse_utc_timestamp(reservation["prepared_utc"], "reservation prepared_utc")
    state = reservation.get("state")
    if state not in {"prepared", "submitted", "recorded", "cancelled"}:
        raise ValueError(f"transaction reservation has invalid state: {state!r}")
    if state in {"submitted", "recorded"}:
        require_numeric_job_id(str(reservation.get("job_id", "")))
    if state == "recorded":
        try:
            actual = float(reservation["actual_node_hours"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("recorded transaction reservation lacks actual use") from error
        if not math.isfinite(actual) or actual < 0.0:
            raise ValueError("recorded transaction reservation has invalid actual use")
        if reservation.get("result") not in {
            "accepted", "clean_partial", "rejected", "failed", "aborted"
        }:
            raise ValueError("recorded transaction reservation has invalid result")
    if state == "cancelled" and not isinstance(reservation.get("notes"), str):
        raise ValueError("cancelled transaction reservation lacks notes")
    intent = reservation.get("execution_intent_sha256")
    if (
        paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve()
        and intent is None
    ):
        raise ValueError(
            "canonical transaction reservation lacks execution intent digest"
        )
    if intent is not None and (
        not isinstance(intent, str) or SHA256_PATTERN.fullmatch(intent) is None
    ):
        raise ValueError("transaction reservation execution intent is invalid")
    return reservation


def validate_transaction_ledger_row(row: object,
                                    manifest: dict[str, object]) -> None:
    """Require one journal-controlled ledger row to match its target manifest."""

    if not isinstance(row, dict) or frozenset(row) != frozenset(LEDGER_COLUMNS):
        raise ValueError("recorded transaction ledger row has invalid columns")
    validate_ledger_numeric_fields(row, "recorded transaction ledger row")
    require_numeric_job_id(str(row["job_id"]))
    if row["execution_epoch"] != EXECUTION_EPOCH:
        raise ValueError("recorded transaction ledger row has wrong execution epoch")
    run = manifest.get("run")
    command = manifest.get("command")
    manifest_paths = manifest.get("paths")
    allocation = manifest.get("allocation")
    if (
        not isinstance(run, dict)
        or not isinstance(command, dict)
        or not isinstance(manifest_paths, dict)
        or not isinstance(allocation, dict)
    ):
        raise ValueError("recorded transaction manifest lacks ledger provenance")
    expected = {
        "job_id": manifest.get("job_id"),
        "case_id": run.get("case_id"),
        "case_name": run.get("case_name"),
        "segment": run.get("segment"),
        "nodes": str(allocation.get("nodes")),
        "requested_walltime": allocation.get("requested_walltime"),
        "executable_revision": command.get("executable_revision"),
        "executable_sha256": command.get("executable_sha256"),
        "input_revision": command.get("input_revision"),
        "input_file": command.get("input_file"),
        "output_dir": manifest_paths.get("output_dir"),
    }
    for key, value in expected.items():
        if row.get(key) != value:
            raise ValueError(f"recorded transaction ledger {key} differs from manifest")
    if manifest.get("accounting") != row:
        raise ValueError("recorded transaction manifest accounting differs from ledger")


def validate_submission_audit(paths: dict[str, Path], value: object) -> None:
    """Require fixed machine-readable submit-policy evidence."""

    if not isinstance(value, dict):
        raise ValueError("submission journal lacks submission audit")
    required = {
        "created_utc",
        "offline_local_root",
        "skip_slurm_test",
        "slurm_test_only",
        "acknowledged_shared_root_campaigns",
    }
    optional = {
        "legacy_mark_submitted",
        "initial_queue_authentication",
        "final_queue_authentication",
        "batch_script",
        "shared_root_isolation_clearance",
    }
    if not required.issubset(value) or not frozenset(value).issubset(required | optional):
        raise ValueError("submission journal has invalid submission audit columns")
    parse_utc_timestamp(value["created_utc"], "submission audit created_utc")
    if (
        not isinstance(value["offline_local_root"], bool)
        or not isinstance(value["skip_slurm_test"], bool)
        or not isinstance(value["slurm_test_only"], str)
        or not isinstance(value["acknowledged_shared_root_campaigns"], list)
        or not all(
            isinstance(item, str)
            for item in value["acknowledged_shared_root_campaigns"]
        )
        or (
            "legacy_mark_submitted" in value
            and value["legacy_mark_submitted"] is not True
        )
    ):
        raise ValueError("submission journal has invalid submission audit values")
    offline_local_root = (
        paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve()
    )
    if value["offline_local_root"] is not offline_local_root:
        raise ValueError("submission journal root mode differs from submission audit")
    clearance = value.get("shared_root_isolation_clearance")
    if clearance is not None:
        clearance = require_r17_exact_keys(
            clearance,
            {
                "checkpoint", "artifact", "independent_review",
                "publication_audit", "publication_audit_review",
                "source_authority_publications", "supersedes",
                "retained_terminal_sacct", "live_queue_absence", "consumption",
            },
            "submission audit shared-root clearance",
        )
        require_r17_nonempty_string(
            clearance["checkpoint"], "submission audit shared-root clearance checkpoint"
        )
        for key in (
            "artifact", "independent_review", "publication_audit",
            "publication_audit_review",
        ):
            binding = require_r17_exact_keys(
                clearance[key],
                {"path", "sha256", "mode", "links"},
                f"submission audit shared-root clearance {key}",
            )
            require_r17_nonempty_string(
                binding["path"], f"submission audit shared-root clearance {key} path"
            )
            require_r17_sha256(
                binding["sha256"],
                f"submission audit shared-root clearance {key} SHA-256",
            )
            if binding["mode"] != "0444" or binding["links"] != 1:
                raise ValueError(
                    f"submission audit shared-root clearance {key} profile differs"
                )
        for key in ("source_authority_publications", "supersedes"):
            if not isinstance(clearance[key], list) or not clearance[key]:
                raise ValueError(
                    f"submission audit shared-root clearance {key} must not be empty"
                )
            for index, binding in enumerate(clearance[key]):
                binding = require_r17_exact_keys(
                    binding,
                    {"path", "sha256", "mode", "links"},
                    f"submission audit shared-root clearance {key} {index}",
                )
                require_r17_nonempty_string(
                    binding["path"],
                    f"submission audit shared-root clearance {key} {index} path",
                )
                require_r17_sha256(
                    binding["sha256"],
                    f"submission audit shared-root clearance {key} {index} SHA-256",
                )
                if binding["mode"] != "0444" or binding["links"] != 1:
                    raise ValueError(
                        f"submission audit shared-root clearance {key} profile differs"
                    )
        terminal_sacct = require_r17_exact_keys(
            clearance["retained_terminal_sacct"],
            {
                "source_publication", "json_path", "sacct_sha256", "job_id",
                "job_name", "state", "exit_code", "completed_utc",
            },
            "submission audit shared-root clearance retained sacct",
        )
        source_publication = require_r17_exact_keys(
            terminal_sacct["source_publication"],
            {"path", "sha256", "mode", "links"},
            "submission audit shared-root clearance retained sacct publication",
        )
        if (
            source_publication["path"]
            != str(paths["root"] / SHARED_ROOT_CLEARANCE_SUPERSESSION_RELATIVE)
            or source_publication["mode"] != "0444"
            or source_publication["links"] != 1
            or terminal_sacct["json_path"] != "stale_campaign.scheduler.sacct"
            or terminal_sacct["job_id"] != SHARED_ROOT_STALE_JOB_ID
            or terminal_sacct["job_name"] != SHARED_ROOT_STALE_JOB_NAME
            or terminal_sacct["state"] != "COMPLETED"
            or terminal_sacct["exit_code"] != "0:0"
        ):
            raise ValueError(
                "submission audit shared-root clearance retained sacct differs"
            )
        require_r17_sha256(
            source_publication["sha256"],
            "submission audit shared-root clearance retained publication SHA-256",
        )
        require_r17_sha256(
            terminal_sacct["sacct_sha256"],
            "submission audit shared-root clearance retained sacct SHA-256",
        )
        require_r17_scheduler_utc(
            terminal_sacct["completed_utc"],
            "submission audit shared-root clearance retained completion time",
        )
        live_absence = require_r17_exact_keys(
            clearance["live_queue_absence"],
            {"job_id", "absent", "checked_utc", "rows_sha256"},
            "submission audit shared-root clearance live queue absence",
        )
        initial_queue = value.get("initial_queue_authentication")
        if (
            not isinstance(initial_queue, dict)
            or live_absence["job_id"] != SHARED_ROOT_STALE_JOB_ID
            or live_absence["absent"] is not True
            or live_absence["checked_utc"] != initial_queue.get("checked_utc")
            or live_absence["rows_sha256"] != initial_queue.get("rows_sha256")
        ):
            raise ValueError(
                "submission audit shared-root clearance live queue absence differs"
            )
        validate_queue_authentication_evidence(
            initial_queue, "submission audit shared-root clearance complete queue"
        )
        consumption = require_r17_exact_keys(
            clearance["consumption"],
            {"action", "case_id", "authorized_stale_manifest"},
            "submission audit shared-root clearance consumption",
        )
        if (
            consumption["action"] not in {"check-submit", "submit"}
            or consumption["case_id"]
            not in {f"R{number:02d}" for number in range(3, 18)}
            or consumption["authorized_stale_manifest"]
            != str(
                paths["root"] / "runs" / SHARED_ROOT_STALE_CAMPAIGN_ID
                / "manifest/prepared_run.json"
            )
        ):
            raise ValueError(
                "submission audit shared-root clearance consumption differs"
            )
    acknowledged = value["acknowledged_shared_root_campaigns"]
    if (
        not offline_local_root
        and bool(acknowledged) != (clearance is not None)
    ):
        raise ValueError(
            "production shared-root acknowledgement and managed clearance differ"
        )
    for key in ("initial_queue_authentication", "final_queue_authentication"):
        if key in value:
            validate_queue_authentication_evidence(value[key], key)
    binding = value.get("batch_script")
    if binding is None:
        if not (
            offline_local_root
            and value.get("legacy_mark_submitted") is True
        ):
            raise ValueError("submission journal lacks descriptor-bound batch script")
    else:
        validate_batch_script_binding(binding, "submission-journal batch script")


def open_batch_script_from_binding(value: object) -> dict[str, object]:
    """Reopen the exact script retained by a submission journal."""

    binding = validate_batch_script_binding(value, "submission-journal batch script")
    path = Path(str(binding["path"]))
    require_symlink_free_path(path, "submission-journal batch script")
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ValueError("submission-journal batch script is unavailable") from error
    authenticated = {
        "fd": descriptor,
        "descriptor_path": f"/proc/self/fd/{descriptor}",
        "path": path,
        "binding": binding,
    }
    try:
        reauthenticate_open_batch_script(authenticated, binding)
    except BaseException:
        close_authenticated_batch_script(authenticated)
        raise
    return authenticated


def require_transaction_batch_script(
    transaction: dict[str, object],
    authenticated: dict[str, object] | None = None,
) -> dict[str, object] | None:
    """Authenticate a journal's exact script before recovery or commit."""

    audit = transaction.get("submission_audit")
    if not isinstance(audit, dict):
        raise ValueError("submission transition lacks retained policy audit")
    binding = audit.get("batch_script")
    if binding is None:
        if (
            audit.get("offline_local_root") is True
            and audit.get("legacy_mark_submitted") is True
        ):
            return None
        raise ValueError("submission transition lacks descriptor-bound batch script")
    if authenticated is not None:
        reauthenticate_open_batch_script(authenticated, binding)
        return authenticated
    reopened = open_batch_script_from_binding(binding)
    close_authenticated_batch_script(reopened)
    return None


def scheduler_output_contains_job(output: str, expected_job_name: str) -> bool:
    """Return whether retained scheduler rows name the prepared allocation."""

    for row in csv.reader(output.splitlines(), delimiter="|"):
        if not row:
            continue
        if len(row) < 2:
            raise ValueError("retained scheduler evidence contains a malformed row")
        if row[1] == expected_job_name:
            return True
    return False


def validate_scheduler_absence_evidence(paths: dict[str, Path],
                                        value: object,
                                        transaction: dict[str, object]) -> None:
    """Require retained scheduler-clear evidence to use one fixed schema."""

    if not isinstance(value, dict):
        raise ValueError("cleared submission transaction lacks scheduler evidence")
    common = {
        "mode", "checked_utc", "expected_job_name", "ambiguity_created_utc",
    }
    mode = value.get("mode")
    optional = {
        "live scheduler absence query": {
            "squeue_command", "squeue_output", "sacct_command", "sacct_output",
        },
        "offline-local fixture": {"fixture"},
        "offline-local operator confirmation": set(),
        "break-glass after scheduler query failure": {
            "operator_evidence", "query_error",
        },
    }
    if mode not in optional or frozenset(value) != frozenset(common | optional[mode]):
        raise ValueError("cleared submission scheduler evidence has invalid schema")
    parse_utc_timestamp(value["checked_utc"], "scheduler absence checked_utc")
    manifest = transaction.get("manifest")
    if not isinstance(manifest, dict):
        raise ValueError("cleared submission transaction lacks manifest evidence")
    expected_name = expected_job_name(manifest)
    if (
        value["ambiguity_created_utc"] != transaction.get("created_utc")
        or value["expected_job_name"] != expected_name
    ):
        raise ValueError("cleared submission scheduler evidence is inconsistent")
    offline_local_root = (
        paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve()
    )
    if mode in {
        "offline-local fixture", "offline-local operator confirmation",
    } and not offline_local_root:
        raise ValueError(
            "canonical cleared submission may not retain offline scheduler evidence"
        )
    if mode == "live scheduler absence query" and (
        not isinstance(value["squeue_command"], list)
        or not isinstance(value["sacct_command"], list)
        or not isinstance(value["squeue_output"], str)
        or not isinstance(value["sacct_output"], str)
    ):
        raise ValueError("cleared submission live scheduler evidence is invalid")
    if mode == "live scheduler absence query" and (
        scheduler_output_contains_job(value["squeue_output"], expected_name)
        or scheduler_output_contains_job(value["sacct_output"], expected_name)
    ):
        raise ValueError(
            "retained scheduler evidence still reports the prepared allocation"
        )
    if mode == "offline-local fixture" and (
        not isinstance(value["fixture"], dict)
        or value["fixture"].get("absent") is not True
    ):
        raise ValueError("cleared submission offline fixture is invalid")
    if mode == "break-glass after scheduler query failure" and (
        not isinstance(value["operator_evidence"], str)
        or not value["operator_evidence"].strip()
        or not isinstance(value["query_error"], str)
        or not value["query_error"].strip()
    ):
        raise ValueError("cleared submission break-glass evidence is invalid")


def transaction_expected_columns(kind: str) -> frozenset[str]:
    """Return the fixed journal schema for one transition kind."""

    if kind == "submit_pending":
        return TRANSACTION_COMMON_COLUMNS | frozenset({
            "prepared_manifest_sha256", "submission_audit",
        })
    columns = TRANSACTION_COMMON_COLUMNS | TRANSACTION_PAYLOAD_COLUMNS
    if kind == "submitted":
        return columns | frozenset({
            "prepared_manifest_sha256", "submission_audit",
            "job_id", "submitted_recorded_utc",
        })
    if kind == "submit_cleared":
        return columns | frozenset({
            "prepared_manifest_sha256", "submission_audit",
            "recovery_notes", "scheduler_absence_evidence",
        })
    return columns


def read_transaction(paths: dict[str, Path], path: Path) -> dict[str, object]:
    """Read and fully validate one durable Stage I transition journal."""

    declared = path.expanduser().absolute()
    transaction_root = paths["transactions"].expanduser().absolute()
    if (
        declared.parent != transaction_root
        or declared.suffix != ".json"
        or declared.name in {"", ".", ".."}
    ):
        raise ValueError(f"transaction journal is outside the E03 store: {declared}")
    retained = read_r17_evidence_bytes(
        declared,
        "Stage I transaction journal",
        mode=0o644,
        owner_controlled=True,
        symlink_free=True,
    )
    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid Stage I transaction journal: {declared}") from error
    if not isinstance(value, dict):
        raise ValueError(f"invalid Stage I transaction journal: {declared}")
    kind = value.get("kind")
    if kind not in TRANSACTION_KINDS:
        raise ValueError(f"transaction journal has invalid kind: {declared}")
    if frozenset(value) != transaction_expected_columns(str(kind)):
        raise ValueError(f"transaction journal has invalid schema for {kind}: {declared}")
    if value.get("schema_version") != 1:
        raise ValueError(f"transaction journal has wrong schema version: {declared}")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(f"transaction journal has wrong execution epoch: {declared}")
    if value.get("transaction_id") != declared.stem:
        raise ValueError(f"transaction journal ID differs from filename: {declared}")
    parse_utc_timestamp(value.get("created_utc"), "transaction created_utc")
    prior_reservations_sha256 = value.get("prior_reservations_sha256")
    if (
        not isinstance(prior_reservations_sha256, str)
        or SHA256_PATTERN.fullmatch(prior_reservations_sha256) is None
    ):
        raise ValueError(f"transaction journal has invalid reservation baseline: {declared}")
    prior_reservations = value.get("prior_reservations")
    if (
        not isinstance(prior_reservations, list)
        or stable_json_sha256(prior_reservations) != prior_reservations_sha256
    ):
        raise ValueError(f"transaction journal reservation baseline is invalid: {declared}")
    reservation_records_by_manifest(
        paths, prior_reservations, "transaction prior reservation snapshot"
    )
    require_active_reservation_policy(prior_reservations)
    manifest_path = require_safe_manifest_path(paths, value["manifest_path"])
    if kind == "submit_pending":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("pending submission journal has invalid manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        authenticate_canonical_reservation_snapshot(paths, prior_reservations)
        require_retained_r17_policy(paths, prior_reservations)
        require_reservation_budget(
            read_ledger(paths), prior_reservations,
            "pending submission reservation snapshot",
        )
        require_transaction_batch_script(value)
        return value
    manifest = value.get("manifest")
    reservations = value.get("reservations")
    if not isinstance(manifest, dict) or not isinstance(reservations, list):
        raise ValueError(f"transaction payload is incomplete: {declared}")
    require_current_epoch(manifest, "transaction manifest")
    if Path(str(manifest.get("project_root", ""))).resolve() != paths["root"].resolve():
        raise ValueError("transaction manifest project root differs from E03 root")
    run = manifest.get("run")
    if (
        not isinstance(run, dict)
        or run.get("case_id") != manifest_path.parents[2].name
        or run.get("segment") != manifest_path.parents[1].name
    ):
        raise ValueError("transaction manifest identity differs from target path")
    validated = [
        validate_reservation_record(paths, reservation)
        for reservation in reservations
    ]
    validate_reservation_snapshot_transition(
        paths, str(kind), manifest_path, prior_reservations, validated
    )
    matches = [
        reservation for reservation in validated
        if Path(str(reservation["manifest"])).resolve() == manifest_path
    ]
    if len(matches) != 1 or matches[0].get("state") != manifest.get("state"):
        raise ValueError("transaction target reservation state differs from manifest")
    expected_state = {
        "prepared": "prepared",
        "submitted": "submitted",
        "submit_cleared": "prepared",
        "recorded": "recorded",
        "cancelled": "cancelled",
        "cancelled_submitted": "cancelled",
    }[str(kind)]
    if manifest.get("state") != expected_state:
        raise ValueError(f"transaction manifest state is invalid for {kind}")
    if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
        validate_prepared_resources(manifest, canonical_production=True)
        require_reservation_matches_manifest(matches[0], manifest)
        require_reserved_execution_intent(
            matches[0], manifest, allow_legacy_local=False
        )
        if kind == "prepared":
            require_prepare_case_policy(
                paths, str(run["case_id"]), int(matches[0]["nodes"]),
                offline_local_root=False,
            )
    row = value.get("ledger_row")
    if kind == "recorded":
        validate_transaction_ledger_row(row, manifest)
        if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
            require_recorded_scheduler_evidence(paths, row, manifest)
        difference = abs(
            float(matches[0].get("actual_node_hours", -1.0))
            - float(row.get("actual_node_hours", -2.0))
        )
        if (
            matches[0].get("result") != row.get("result")
            or not math.isfinite(difference)
            or difference > 5.0e-7
        ):
            raise ValueError("recorded transaction reservation differs from ledger")
    elif row is not None:
        raise ValueError(f"{kind} transaction unexpectedly controls a ledger row")
    if kind == "submitted":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("submitted journal has invalid prepared manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        require_numeric_job_id(str(value.get("job_id", "")))
        if manifest.get("job_id") != value.get("job_id"):
            raise ValueError("submitted journal job ID differs from manifest")
        parse_utc_timestamp(
            value.get("submitted_recorded_utc"), "transaction submitted_recorded_utc"
        )
        require_transaction_batch_script(value)
    if kind == "submit_cleared":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("cleared journal has invalid prepared manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        if not isinstance(value.get("recovery_notes"), str):
            raise ValueError("cleared submission transaction lacks recovery notes")
        validate_scheduler_absence_evidence(
            paths, value.get("scheduler_absence_evidence"), value
        )
        require_transaction_batch_script(value)
    if kind == "cancelled_submitted":
        validate_submitted_cancellation_metadata(paths, manifest_path, manifest)
    baseline_ledger, prospective_ledger = transaction_ledger_views(paths, row)
    for label, snapshot, ledger in (
        ("transaction prior reservation snapshot", prior_reservations, baseline_ledger),
        ("transaction payload reservation snapshot", validated, prospective_ledger),
    ):
        require_retained_r17_policy(paths, snapshot)
        require_reservation_budget(ledger, snapshot, label)
    return value


def unlink_trusted_transaction(paths: dict[str, Path], transaction_path: Path,
                               expected: dict[str, object]) -> None:
    """Retire only the exact journal read through the trusted transaction store."""

    directory = paths["transactions"].absolute()
    declared = transaction_path.absolute()
    if declared.parent != directory:
        raise ValueError("transaction journal removal target is outside the store")
    require_owner_symlink_free_path(directory, "Stage I transaction store")
    try:
        directory_fd = os.open(
            directory,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as error:
        raise ValueError("Stage I transaction store is unavailable for removal") from error
    journal = None
    try:
        directory_profile = require_directory_descriptor_binding(
            directory, directory_fd, "Stage I transaction store"
        )
        if (
            directory_profile.st_uid != os.geteuid()
            or stat.S_IMODE(directory_profile.st_mode) & 0o022
        ):
            raise ValueError("Stage I transaction store is not owner controlled")
        expected_payload = (
            json.dumps(expected, indent=2, sort_keys=True) + "\n"
        ).encode("utf-8")
        journal = open_regular_file_binding(
            directory_fd,
            declared.name,
            "Stage I transaction journal",
            expected_payload=expected_payload,
        )
        if journal is None:
            raise ValueError("Stage I transaction journal changed before removal")
        profile = journal.profile
        if (
            not stat.S_ISREG(profile.st_mode)
            or profile.st_uid != os.geteuid()
            or stat.S_IMODE(profile.st_mode) != 0o644
            or profile.st_nlink != 1
        ):
            raise ValueError("Stage I transaction journal changed before removal")
        quarantine_bound_predecessor(
            directory_fd,
            declared.name,
            journal,
            directory,
            directory_profile,
            declared,
            "Stage I transaction journal retirement",
        )
    finally:
        if journal is not None:
            journal.close()
        os.close(directory_fd)


def expected_pretransition_reservation(kind: str,
                                       reservation: dict[str, object],
                                       ) -> dict[str, object] | None:
    """Return the only permitted pre-transition form of one target reservation."""

    if kind == "prepared":
        return None
    result = dict(reservation)
    if kind == "submitted":
        result["state"] = "prepared"
        result.pop("job_id", None)
    elif kind == "recorded":
        result["state"] = "submitted"
        result.pop("actual_node_hours", None)
        result.pop("result", None)
    elif kind == "cancelled":
        result["state"] = "prepared"
        result.pop("notes", None)
    elif kind == "cancelled_submitted":
        result["state"] = "submitted"
        result.pop("notes", None)
    elif kind != "submit_cleared":
        raise ValueError(f"transaction kind has no reservation transition: {kind}")
    return result


def reservation_records_by_manifest(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    label: str,
) -> dict[Path, dict[str, object]]:
    """Index one validated reservation snapshot without duplicate targets."""

    indexed: dict[Path, dict[str, object]] = {}
    for reservation in reservations:
        validated = validate_reservation_record(paths, reservation)
        path = Path(str(validated["manifest"])).resolve()
        if path in indexed:
            raise ValueError(f"{label} contains duplicate manifest records")
        indexed[path] = validated
    return indexed


def authenticate_canonical_reservation_snapshot(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    exempt_paths: frozenset[Path] = frozenset(),
) -> None:
    """Bind canonical baseline reservations to their retained launch intents."""

    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    indexed = reservation_records_by_manifest(
        paths, reservations, "canonical reservation snapshot"
    )
    for path, reservation in indexed.items():
        if path in exempt_paths:
            continue
        if not path.is_file():
            raise ValueError(
                f"canonical reservation lacks retained manifest during replay: {path}"
            )
        manifest = read_manifest(path)
        require_current_epoch(manifest, "canonical replay reservation manifest")
        validate_prepared_resources(manifest, canonical_production=True)
        require_reservation_matches_manifest(reservation, manifest)
        require_reserved_execution_intent(
            reservation, manifest, allow_legacy_local=False
        )


def validate_reservation_snapshot_transition(
    paths: dict[str, Path],
    kind: str,
    target: Path,
    prior: list[dict[str, object]],
    payload: list[dict[str, object]],
) -> None:
    """Require a journal snapshot to preserve every unrelated reservation."""

    prior_by_manifest = reservation_records_by_manifest(
        paths, prior, "transaction prior reservation snapshot"
    )
    payload_by_manifest = reservation_records_by_manifest(
        paths, payload, "transaction reservation snapshot"
    )
    require_active_reservation_policy(prior)
    require_active_reservation_policy(payload)
    authenticate_canonical_reservation_snapshot(
        paths, prior, exempt_paths=frozenset({target})
    )
    authenticate_canonical_reservation_snapshot(
        paths, payload, exempt_paths=frozenset({target})
    )
    require_retained_r17_policy(paths, prior)
    require_retained_r17_policy(paths, payload)
    if (
        {path: record for path, record in prior_by_manifest.items() if path != target}
        != {path: record for path, record in payload_by_manifest.items() if path != target}
    ):
        raise ValueError("transaction snapshot would alter unrelated reservations")
    payload_target = payload_by_manifest.get(target)
    if payload_target is None:
        raise ValueError("transaction snapshot lacks target reservation")
    if prior_by_manifest.get(target) != expected_pretransition_reservation(
        kind, payload_target
    ):
        raise ValueError("transaction snapshot target reservation transition is invalid")


def validate_transaction_reservation_baseline(
    paths: dict[str, Path],
    transaction: dict[str, object],
) -> None:
    """Authenticate reservation replay against its retained pre-transition store."""

    reservations = transaction.get("reservations")
    if not isinstance(reservations, list):
        raise ValueError("transaction reservation snapshot is invalid")
    current = read_reservations(paths)
    require_active_reservation_policy(current)
    current_sha256 = sha256(paths["reservations"])
    payload_sha256 = stable_json_sha256(reservations)
    if current_sha256 not in {
        transaction.get("prior_reservations_sha256"), payload_sha256,
    }:
        raise ValueError("reservation store differs from transaction replay baseline")
    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    target = require_safe_manifest_path(paths, transaction["manifest_path"])
    reservation_records_by_manifest(paths, current, "reservation store")
    authenticate_canonical_reservation_snapshot(
        paths, current, exempt_paths=frozenset({target})
    )
    require_retained_r17_policy(paths, current)
    baseline_ledger, prospective_ledger = transaction_ledger_views(
        paths, transaction.get("ledger_row")
    )
    ledger = (
        baseline_ledger
        if current_sha256 == transaction.get("prior_reservations_sha256")
        else prospective_ledger
    )
    require_reservation_budget(
        ledger, current,
        "reservation replay baseline",
    )


def apply_transaction(paths: dict[str, Path], transaction_path: Path) -> None:
    """Idempotently complete one journaled metadata transition."""

    transaction = read_transaction(paths, transaction_path)
    kind = transaction.get("kind")
    if kind == "submit_pending":
        raise ValueError(
            "submission outcome is ambiguous; use recover-submit with the "
            f"scheduler job ID: {transaction_path}"
        )
    if kind not in {
        "prepared", "submitted", "submit_cleared", "recorded", "cancelled",
        "cancelled_submitted",
    }:
        raise ValueError(f"transaction journal has invalid kind: {transaction_path}")
    manifest_path = Path(str(transaction["manifest_path"])).resolve()
    manifest = transaction.get("manifest")
    reservations = transaction.get("reservations")
    if not isinstance(manifest, dict) or not isinstance(reservations, list):
        raise ValueError(f"transaction payload is incomplete: {transaction_path}")
    validate_transaction_reservation_baseline(paths, transaction)
    row = transaction.get("ledger_row")
    _, prospective_ledger = transaction_ledger_views(paths, row)
    require_reservation_budget(
        prospective_ledger, reservations, "transaction payload reservation snapshot"
    )
    require_retained_r17_policy(paths, reservations)
    if row is not None:
        if not isinstance(row, dict):
            raise ValueError(f"transaction ledger row is invalid: {transaction_path}")
        ledger = read_ledger(paths)
        matches = [item for item in ledger if item.get("job_id") == row.get("job_id")]
        if len(matches) > 1:
            raise ValueError(f"transaction ledger job is duplicated: {transaction_path}")
        if matches and matches[0] != row:
            raise ValueError(f"transaction ledger row conflicts: {transaction_path}")
        if not matches:
            append_ledger_row(paths["ledger"], row, ledger)
    write_json(paths["reservations"], reservations)
    write_json(manifest_path, manifest)
    refresh_summary(paths)
    unlink_trusted_transaction(paths, transaction_path, transaction)


def durable_transition(paths: dict[str, Path], kind: str, manifest_path: Path,
                       manifest: dict[str, object],
                       reservations: list[dict[str, object]],
                       ledger_row: dict[str, object] | None = None) -> None:
    """Journal and apply one recoverable cross-file metadata transition."""

    require_no_pending_transactions(paths)
    require_safe_manifest_path(paths, manifest_path)
    prior_reservations = read_reservations(paths)
    transaction_path = (
        paths["transactions"] / f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}.json"
    )
    write_json(transaction_path, {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_path.stem,
        "kind": kind,
        "created_utc": utc_now(),
        "manifest_path": str(manifest_path),
        "prior_reservations": prior_reservations,
        "prior_reservations_sha256": stable_json_sha256(prior_reservations),
        "manifest": manifest,
        "reservations": reservations,
        "ledger_row": ledger_row,
    }, mode=0o644)
    apply_transaction(paths, transaction_path)


def write_submit_pending_transaction(paths: dict[str, Path],
                                     manifest_path: Path,
                                     submission_audit: dict[str, object],
                                     ) -> Path:
    """Persist an ambiguity barrier immediately before invoking sbatch."""

    require_no_pending_transactions(paths)
    require_safe_manifest_path(paths, manifest_path)
    prior_reservations = read_reservations(paths)
    transaction_path = (
        paths["transactions"] / f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}.json"
    )
    write_json(transaction_path, {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_path.stem,
        "kind": "submit_pending",
        "created_utc": utc_now(),
        "manifest_path": str(manifest_path),
        "prior_reservations": prior_reservations,
        "prior_reservations_sha256": stable_json_sha256(prior_reservations),
        "prepared_manifest_sha256": sha256(manifest_path),
        "submission_audit": submission_audit,
    }, mode=0o644)
    return transaction_path


def submit_pending_transaction(paths: dict[str, Path],
                               manifest_path: Path) -> Path:
    """Return the one ambiguous sbatch journal for a prepared manifest."""

    matches = []
    for path in pending_transaction_paths(paths):
        transaction = read_transaction(paths, path)
        if (
            transaction.get("kind") == "submit_pending"
            and Path(str(transaction.get("manifest_path", ""))).resolve()
            == manifest_path.resolve()
        ):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(
            f"expected one pending submission journal for {manifest_path}, "
            f"found {len(matches)}"
        )
    return matches[0]


def finish_submit_transaction(paths: dict[str, Path], transaction_path: Path,
                              manifest_path: Path,
                              manifest: dict[str, object],
                              reservations: list[dict[str, object]],
                              job_id: str,
                              authenticated_batch: dict[str, object] | None = None,
                              ) -> None:
    """Attach a scheduler ID to an ambiguity journal and commit submission."""

    transaction = read_transaction(paths, transaction_path)
    if transaction.get("kind") != "submit_pending":
        raise ValueError(f"transaction is not an ambiguous submission: {transaction_path}")
    if transaction.get("prepared_manifest_sha256") != sha256(manifest_path):
        raise ValueError("prepared manifest changed after the sbatch boundary")
    if transaction.get("prior_reservations_sha256") != sha256(paths["reservations"]):
        raise ValueError("reservation store changed after the sbatch boundary")
    require_transaction_batch_script(transaction, authenticated_batch)
    require_numeric_job_id(job_id)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if manifest.get("state") != "prepared" or reservation.get("state") != "prepared":
        raise ValueError("submission recovery requires matching prepared state")
    require_reserved_execution_intent(
        reservation, manifest,
        allow_legacy_local=paths["root"].resolve() != DEFAULT_ROOT.resolve(),
    )
    submitted_utc = utc_now()
    audit = transaction.get("submission_audit")
    if not isinstance(audit, dict):
        raise ValueError("submission transition lacks retained policy audit")
    audits = manifest.setdefault("submission_audits", [])
    if not isinstance(audits, list):
        raise ValueError("prepared manifest has invalid submission audits")
    audits.append(audit)
    manifest["state"] = "submitted"
    manifest["job_id"] = job_id
    manifest["submitted_recorded_utc"] = submitted_utc
    reservation["state"] = "submitted"
    reservation["job_id"] = job_id
    transaction.update({
        "kind": "submitted",
        "job_id": job_id,
        "submitted_recorded_utc": submitted_utc,
        "manifest": manifest,
        "reservations": reservations,
        "ledger_row": None,
    })
    write_json(transaction_path, transaction, mode=0o644)
    apply_transaction(paths, transaction_path)


def active_reservations(reservations: list[dict[str, object]]
                        ) -> list[dict[str, object]]:
    """Return unaccounted prepared or submitted segment reservations."""

    return [
        item for item in reservations
        if isinstance(item, dict)
        and item.get("state") in {"prepared", "submitted"}
    ]


def require_reservation_budget(ledger: list[dict[str, str]],
                               reservations: list[dict[str, object]],
                               label: str) -> tuple[float, float]:
    """Require one prospective reservation snapshot to fit both budget ceilings."""

    try:
        actual = sum(float(row["actual_node_hours"]) for row in ledger)
        reserved = sum(
            float(item["reserved_node_hours"])
            for item in active_reservations(reservations)
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label} usage is invalid") from error
    if not math.isfinite(actual) or not math.isfinite(reserved):
        raise ValueError(f"{label} usage is invalid")
    if actual + reserved > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError(f"{label} exceeds the Stage I reservation ceiling")
    if actual + reserved > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError(f"{label} exceeds the incremental project ceiling")
    return actual, reserved


def transaction_ledger_views(paths: dict[str, Path],
                             row: object
                             ) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Return the modeled ledgers immediately before and after one replay."""

    ledger = read_ledger(paths)
    if row is None:
        return ledger, ledger
    if not isinstance(row, dict):
        raise ValueError("transaction ledger row is invalid")
    matches = [item for item in ledger if item.get("job_id") == row.get("job_id")]
    if len(matches) > 1:
        raise ValueError("transaction ledger job is duplicated")
    if matches and matches[0] != row:
        raise ValueError("transaction ledger row conflicts")
    baseline = [item for item in ledger if item.get("job_id") != row.get("job_id")]
    return baseline, [*baseline, row]


def require_active_reservation_policy(
    reservations: list[dict[str, object]],
    candidate_case_id: str | None = None,
    candidate_nodes: int | None = None,
) -> list[dict[str, object]]:
    """Enforce bounded distinct-case/node overlap while keeping preparation serial."""

    active = active_reservations(reservations)
    case_ids = []
    active_nodes = 0
    for reservation in active:
        case_id = str(reservation.get("case_id", ""))
        require_authorized_case(case_id)
        if case_id in case_ids:
            raise ValueError(f"active Stage I reservations duplicate case {case_id}")
        case_ids.append(case_id)
        nodes = reservation.get("nodes")
        if isinstance(nodes, bool) or not isinstance(nodes, int) or nodes < 1:
            raise ValueError(
                f"active Stage I reservation {case_id} has invalid nodes"
            )
        active_nodes += nodes
    prepared = [
        reservation for reservation in active
        if reservation.get("state") == "prepared"
    ]
    if len(active) > MAX_ACTIVE_STAGE_I_SEGMENTS:
        raise ValueError(
            f"active Stage I reservations exceed the "
            f"{MAX_ACTIVE_STAGE_I_SEGMENTS}-segment concurrency limit"
        )
    if active_nodes > MAX_ACTIVE_STAGE_I_NODES:
        raise ValueError(
            f"active Stage I reservations exceed the "
            f"{MAX_ACTIVE_STAGE_I_NODES}-node concurrency limit"
        )
    if len(prepared) > 1:
        raise ValueError("only one Stage I segment may be prepared at a time")
    if R17_CASE_ID in case_ids and len(active) != 1:
        raise ValueError(f"{R17_CASE_ID} requires exclusive Stage I execution")
    if len(active) > 1 and any(
        case_id not in CONCURRENT_CASE_IDS for case_id in case_ids
    ):
        raise ValueError(
            "overlapping Stage I reservations are restricted to distinct R03-R16 cases"
        )
    if candidate_case_id is None:
        if candidate_nodes is not None:
            raise ValueError("candidate nodes require a candidate Stage I case")
        return active
    require_authorized_case(candidate_case_id)
    if (
        isinstance(candidate_nodes, bool)
        or not isinstance(candidate_nodes, int)
        or candidate_nodes < 1
    ):
        raise ValueError(
            f"candidate Stage I reservation {candidate_case_id} has invalid nodes"
        )
    if candidate_case_id in case_ids:
        raise ValueError(
            f"{candidate_case_id} already has an active Stage I reservation"
        )
    if active and (
        candidate_case_id == R17_CASE_ID or R17_CASE_ID in case_ids
    ):
        raise ValueError(f"{R17_CASE_ID} requires exclusive Stage I execution")
    if active and (
        candidate_case_id not in CONCURRENT_CASE_IDS
        or any(case_id not in CONCURRENT_CASE_IDS for case_id in case_ids)
    ):
        raise ValueError(
            "overlapping Stage I reservations are restricted to distinct R03-R16 cases"
        )
    if prepared:
        raise ValueError(
            "another Stage I segment is prepared; submit or cancel it before "
            "preparing a new segment"
        )
    if len(active) >= MAX_ACTIVE_STAGE_I_SEGMENTS:
        raise ValueError(
            f"active Stage I reservations reached the "
            f"{MAX_ACTIVE_STAGE_I_SEGMENTS}-segment concurrency limit"
        )
    if active_nodes + candidate_nodes > MAX_ACTIVE_STAGE_I_NODES:
        raise ValueError(
            f"candidate Stage I reservation would exceed the "
            f"{MAX_ACTIVE_STAGE_I_NODES}-node concurrency limit"
        )
    return active


def reservation_usage(paths: dict[str, Path]) -> tuple[float, float]:
    """Return actual and actively reserved current-epoch Stage I node-hours."""

    return require_reservation_budget(
        read_ledger(paths), read_reservations(paths), "active Stage I reservation"
    )


def refresh_summary(paths: dict[str, Path]) -> None:
    """Regenerate the human-readable current-epoch production summary."""

    with canonical_root_lock(paths["root"]):
        ledger = read_ledger(paths)
        reservations = read_reservations(paths)
        actual = sum(float(row["actual_node_hours"]) for row in ledger)
        active = active_reservations(reservations)
        reserved = sum(float(item["reserved_node_hours"]) for item in active)
        stage_remaining = max(
            0.0, CURRENT_STAGE_I_RESERVED_NODE_HOURS - actual - reserved
        )
        project_remaining = PROJECT_BUDGET_NODE_HOURS - actual - reserved
        qualification = qualification_approval_status(paths)
        if qualification["state"] == "approved":
            qualification_line = (
                "- E03 corrected-build Frontier qualification: approved by "
                f"`{qualification['approved_by']}` at "
                f"`{qualification['approved_utc']}` for executable "
                f"`{qualification['approved_executable_sha256']}`."
            )
        else:
            qualification_line = (
                "- E03 corrected-build Frontier qualification: pending until "
                f"a valid approval token exists at `{qualification['path']}` "
                f"({qualification['reason']})."
            )
        lines = [
            f"# MKS24 Stage I Frontier {EXECUTION_EPOCH} Budget",
            "",
            f"- Updated UTC: `{utc_now()}`",
            f"- Execution epoch: `{EXECUTION_EPOCH}`",
            qualification_line,
            f"- Fresh incremental project ceiling: "
            f"`{PROJECT_BUDGET_NODE_HOURS:.6f}` node-hours",
            f"- Historical debug qualification use, reported but not charged to E03: "
            f"`{HISTORICAL_DEBUG_NODE_HOURS:.6f}` node-hours",
            f"- Historical E01 Stage I use, reported but not charged to E03: "
            f"`{HISTORICAL_E01_STAGE_I_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 pipeline evidence, reported but not charged to E03: "
            f"`{HISTORICAL_E02_PIPELINE_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R16 use: `{COMPLETED_R16_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R02 standard-layout timing-pilot use: "
            f"`{COMPLETED_R02_STANDARD_LAYOUT_PILOT_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R17 high-resolution timing-pilot use: "
            f"`{COMPLETED_R17_HIGH_RESOLUTION_PILOT_NODE_HOURS:.6f}` node-hours",
            f"- Current E03 mapped-matrix planning envelope: "
            f"`{CURRENT_STAGE_I_RESERVED_NODE_HOURS:.6f}` node-hours",
            f"- {EXECUTION_EPOCH} Stage I actual use: `{actual:.6f}` node-hours",
            f"- Active segment reservations: `{reserved:.6f}` node-hours",
            f"- Unreserved E03 mapped-matrix remainder: "
            f"`{stage_remaining:.6f}` node-hours",
            f"- Incremental project remainder after active E03 Stage I use: "
            f"`{project_remaining:.6f}` node-hours",
            "",
            "## Recorded Segments",
            "",
        ]
        if ledger:
            lines.extend([
                "| Job | Case/segment | State | Node-hours | Result |",
                "| --- | --- | --- | ---: | --- |",
            ])
            for row in ledger:
                lines.append(
                    "| `{job_id}` | `{case_id}/{segment}` | {state} | "
                    "`{actual_node_hours}` | {result} |".format(**row)
                )
        else:
            lines.append(
                f"No {EXECUTION_EPOCH} Stage I production allocation has been recorded."
            )
        lines.extend(["", "## Active Reservations", ""])
        if active:
            lines.extend([
                "| Case/segment | Nodes | Walltime | Node-hours | State |",
                "| --- | ---: | --- | ---: | --- |",
            ])
            for item in active:
                lines.append(
                    "| `{case_id}/{segment}` | `{nodes}` | `{requested_walltime}` | "
                    "`{reserved_node_hours:.6f}` | {state} |".format(**item)
                )
        else:
            lines.append(
                f"No prepared or submitted {EXECUTION_EPOCH} Stage I segment "
                "is reserved."
            )
        lines.extend([
            "",
            f"Up to {MAX_ACTIVE_STAGE_I_SEGMENTS} distinct R03-R16 cases may "
            f"be active concurrently within {MAX_ACTIVE_STAGE_I_NODES} total "
            "nodes, with at most one prepared submission packet. "
            f"{R17_CASE_ID} remains exclusive and last. "
            "E02 is retained only as pipeline and cost evidence after the Phase A "
            "forcing-policy audit. "
            "Jobs use the `batch` partition with Frontier's default production "
            "`normal` QOS; the `debug` QOS is not used for paper production.",
            "",
        ])
        write_text(paths["summary"], "\n".join(lines))


def load_matrix(path: Path) -> dict[str, object]:
    """Read the mapped Stage I matrix."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid matrix manifest: {path}")
    return value


def parse_input_parameters(path: Path) -> dict[str, str]:
    """Read block-qualified input values while ignoring comments."""

    block = ""
    values: dict[str, str] = {}
    for original in path.read_text(encoding="utf-8").splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if "=" in line and block:
            parameter, value = line.split("=", 1)
            values[f"{block}/{parameter.strip()}"] = value.strip()
    return values


def time_tlim_override_target(overrides: list[str],
                              allow_missing: bool = False) -> float | None:
    """Return the single positive finite ``time/tlim`` execution target."""

    targets = [
        override.split("=", 1)[1]
        for override in overrides
        if "=" in override and override.split("=", 1)[0] == "time/tlim"
    ]
    if not targets and allow_missing:
        return None
    if len(targets) != 1:
        raise ValueError(
            "prepare requires exactly one --override time/tlim=<numeric-target>"
        )
    try:
        target = float(targets[0])
    except ValueError as error:
        raise ValueError("time/tlim override target must be numeric") from error
    if not math.isfinite(target) or target <= 0.0:
        raise ValueError("time/tlim override target must be positive and finite")
    return target


def validate_prepare_overrides(path: Path, overrides: list[str],
                               allow_missing_time_target: bool = False,
                               canonical_production: bool = False,
                               ) -> float | None:
    """Require declared input keys and one explicit segment time target."""

    keys = set(parse_input_parameters(path))
    if canonical_production and (
        len(overrides) != 1 or not overrides[0].startswith("time/tlim=")
    ):
        raise ValueError(
            "canonical production requires exactly one --override "
            "time/tlim=<positive-target>"
        )
    for override in overrides:
        if "=" not in override:
            raise ValueError(f"override must use block/name=value: {override}")
        key, _ = override.split("=", 1)
        if key not in keys:
            raise ValueError(
                f"override targets parameter absent from input deck: {key}"
            )
    return time_tlim_override_target(
        overrides, allow_missing=allow_missing_time_target
    )


def prepared_time_tlim_target(manifest: dict[str, object]) -> float:
    """Return and verify the prepared segment's retained execution target."""

    command = manifest.get("command")
    if not isinstance(command, dict):
        raise ValueError("prepared manifest lacks command metadata")
    overrides = command.get("overrides")
    if not isinstance(overrides, list) or not all(
        isinstance(item, str) for item in overrides
    ):
        raise ValueError("prepared manifest lacks command overrides")
    override_target = time_tlim_override_target(overrides)
    retained_target = command.get("time_tlim_target", override_target)
    try:
        target = float(retained_target)
    except (TypeError, ValueError) as error:
        raise ValueError("prepared manifest has invalid time/tlim target") from error
    if (
        not math.isfinite(target)
        or override_target is None
        or abs(target - override_target) > 1.0e-12
    ):
        raise ValueError("prepared manifest time/tlim target is inconsistent")
    return target


def validate_prepared_continuation_target(manifest: dict[str, object]) -> None:
    """Require a retained continuation target to advance beyond its parent."""

    command = manifest.get("command")
    if not isinstance(command, dict):
        raise ValueError("prepared manifest lacks command metadata")
    parent = command.get("parent_segment")
    if parent is None:
        return
    if not isinstance(parent, dict):
        raise ValueError("prepared continuation parent metadata is invalid")
    try:
        final_time = float(parent["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("prepared continuation lacks parent final time") from error
    target = prepared_time_tlim_target(manifest)
    if not math.isfinite(final_time) or target <= final_time + 1.0e-12:
        raise ValueError(
            "prepared continuation time/tlim target does not advance beyond "
            f"parent inspection time {final_time:.12g}"
        )


def validate_matrix(matrix_path: Path, source_dir: Path) -> dict[str, object]:
    """Validate unique case mapping and admitted file-level aliases."""

    matrix = load_matrix(matrix_path)
    cases = matrix.get("cases")
    if not isinstance(cases, list):
        raise ValueError("matrix cases must be a list")
    expected_ids = [f"R{number:02d}" for number in range(2, 18)]
    ids = [str(case.get("id")) for case in cases]
    if ids != expected_ids:
        raise ValueError(f"Stage I case identifiers must be {expected_ids}: {ids}")
    inputs = [str(case.get("input")) for case in cases]
    if len(set(inputs)) != len(inputs):
        raise ValueError("canonical Stage I inputs are not unique")
    if matrix.get("authorization", {}).get("mapped_unique_runs") != len(cases):
        raise ValueError("mapped run count disagrees with matrix cases")
    for path_text in inputs:
        if not (source_dir / path_text).is_file():
            raise ValueError(f"matrix input does not exist: {path_text}")
    excluded = {
        str(item.get("input"))
        for item in matrix.get("excluded_definitions", [])
    }
    if excluded.intersection(inputs):
        raise ValueError("an excluded input appears in the canonical run matrix")
    for alias in matrix.get("aliases", []):
        equivalent = alias.get("equivalent_input")
        if equivalent is None:
            continue
        canonical = source_dir / str(alias["canonical_input"])
        comparison = source_dir / str(equivalent)
        if not canonical.is_file() or not comparison.is_file():
            raise ValueError(f"alias inputs are missing for {alias['canonical_case']}")
        ignored = set(alias.get("comparison_ignore_parameters", []))
        canonical_values = parse_input_parameters(canonical)
        comparison_values = parse_input_parameters(comparison)
        for key in ignored:
            canonical_values.pop(key, None)
            comparison_values.pop(key, None)
        if canonical_values != comparison_values:
            raise ValueError(
                f"alias is not physically identical for {alias['canonical_case']}"
            )
    return matrix


def case_for_id(matrix: dict[str, object], case_id: str) -> dict[str, object]:
    """Select one unique canonical Stage I case."""

    matches = [
        case for case in matrix["cases"]
        if str(case.get("id")) == case_id
    ]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate Stage I case ID: {case_id}")
    return matches[0]


def git_revision_for_input(source_dir: Path, input_path: Path,
                           matrix_path: Path) -> str:
    """Require submitted input and matrix to match committed source content."""

    revision = subprocess.run(
        [
            authenticated_git_binary(), "--no-replace-objects",
            *controller_git_worktree_prefix(source_dir), "rev-parse", "HEAD",
        ],
        check=True, capture_output=True, text=True,
        env=hardened_git_environment(),
    ).stdout.strip()
    relative_paths = [
        str(input_path.relative_to(source_dir)),
        str(matrix_path.relative_to(source_dir)),
    ]
    tracked = subprocess.run(
        [
            authenticated_git_binary(), "--no-replace-objects",
            *controller_git_worktree_prefix(source_dir),
            "ls-files", "--error-unmatch", "--", *relative_paths,
        ],
        check=False,
        env=hardened_git_environment(),
    )
    if tracked.returncode != 0:
        raise ValueError("submitted input and matrix must be tracked by Git")
    for diff_args in (
        ["diff", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
        ["diff", "--cached", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
    ):
        result = subprocess.run(
            [
                authenticated_git_binary(), "--no-replace-objects",
                *controller_git_worktree_prefix(source_dir),
                *diff_args, *relative_paths,
            ],
            check=False,
            env=hardened_git_environment(),
        )
        if result.returncode != 0:
            raise ValueError(
                "submitted input and matrix must be committed before production"
            )
    return revision


def source_bundle_provenance(source_bundle_value: str | None,
                             revisions: list[str], root: Path,
                             allow_local_root: bool
                             ) -> dict[str, object] | None:
    """Require retained bundle provenance for real Stage I preparation."""

    if source_bundle_value is None:
        if allow_local_root:
            return None
        raise ValueError(
            "--source-bundle is required for retained Stage I source provenance"
        )
    source_bundle = Path(source_bundle_value).expanduser().resolve()
    require_beneath_root(
        source_bundle, root, "source bundle", allow_local_root
    )
    if not source_bundle.is_file():
        raise ValueError(f"source bundle is missing: {source_bundle}")
    with tempfile.TemporaryDirectory(prefix="cgl_lf_bundle_verify_") as directory:
        repository = Path(directory) / "source.git"
        try:
            subprocess.run(
                [authenticated_git_binary(), "--no-replace-objects",
                 "clone", "--bare", "--quiet",
                 str(source_bundle), str(repository)],
                check=True, capture_output=True, text=True,
                env=hardened_git_environment(),
            )
        except subprocess.CalledProcessError as error:
            raise ValueError(
                f"source bundle cannot be cloned: {source_bundle}"
            ) from error
        for revision in revisions:
            present = subprocess.run(
                [
                    authenticated_git_binary(), "--no-replace-objects",
                    "-C", str(repository),
                    "cat-file", "-e", f"{revision}^{{commit}}",
                ],
                check=False, capture_output=True, text=True,
                env=hardened_git_environment(),
            )
            if present.returncode != 0:
                raise ValueError(
                    f"source bundle does not contain revision {revision}: "
                    f"{source_bundle}"
                )
    return {
        "path": str(source_bundle),
        "sha256": sha256(source_bundle),
        "verified_revisions": revisions,
    }


def production_utility_provenance(allow_uncommitted: bool = False
                                  ) -> dict[str, object]:
    """Require the production-control script itself to be committed."""

    script_path = Path(__file__).resolve()
    relative = str(script_path.relative_to(ROOT_DIR))
    if Path(relative) != PRODUCTION_UTILITY_RELATIVE:
        raise ValueError("production utility path is inconsistent")
    revision = subprocess.run(
        [
            authenticated_git_binary(), "--no-replace-objects",
            *controller_git_worktree_prefix(ROOT_DIR),
            "rev-parse", "--verify", "HEAD",
        ],
        check=True, capture_output=True, text=True,
        env=hardened_git_environment(),
    ).stdout.strip()
    if GIT_REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError("production utility revision is invalid")
    tracked = subprocess.run(
        [
            authenticated_git_binary(), "--no-replace-objects",
            *controller_git_worktree_prefix(ROOT_DIR),
            "ls-files", "--error-unmatch", "--", relative,
        ],
        check=False, capture_output=True, text=True,
        env=hardened_git_environment(),
    )
    if tracked.returncode != 0:
        raise ValueError("production utility is not tracked by Git")
    committed = True
    for diff_args in (
        ["diff", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
        ["diff", "--cached", "--no-ext-diff", "--no-textconv", "--quiet", "--"],
    ):
        result = subprocess.run(
            [
                authenticated_git_binary(), "--no-replace-objects",
                *controller_git_worktree_prefix(ROOT_DIR),
                *diff_args, relative,
            ],
            check=False,
            env=hardened_git_environment(),
        )
        if result.returncode == 1:
            committed = False
        elif result.returncode != 0:
            raise ValueError("cannot determine production utility worktree status")
    if not committed and not allow_uncommitted:
        raise ValueError(
            "production utility must be committed before preparing a segment"
        )
    return {
        "path": str(script_path),
        "revision": revision,
        "sha256": sha256(script_path),
        "committed": committed,
    }


def authenticate_production_utility(
    record: object,
    *,
    source_bundle: object = None,
    allow_uncommitted: bool = False,
    allow_historical: bool = False,
) -> None:
    """Revalidate the retained production helper used to prepare a segment."""

    if not isinstance(record, dict):
        raise ValueError("prepared manifest lacks production utility provenance")
    path = Path(str(record.get("path", ""))).resolve()
    if path != Path(__file__).resolve():
        raise ValueError("prepared production utility path is inconsistent")
    try:
        relative = path.relative_to(ROOT_DIR)
    except ValueError as error:
        raise ValueError("prepared production utility path is inconsistent") from error
    if relative != PRODUCTION_UTILITY_RELATIVE:
        raise ValueError("prepared production utility path is inconsistent")
    revision = record.get("revision")
    if not isinstance(revision, str) or GIT_REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError("prepared production utility revision is invalid")
    if record.get("committed") is not True and not (
        allow_uncommitted and record.get("committed") is False
    ):
        raise ValueError("prepared production utility is not committed")
    expected = record.get("sha256")
    if not isinstance(expected, str) or SHA256_PATTERN.fullmatch(expected) is None:
        raise ValueError("prepared production utility checksum is invalid")
    production_utility_provenance(allow_uncommitted=allow_uncommitted)
    if not allow_historical:
        require_file_sha256(path, expected, "production utility")
        return
    if source_bundle is None and allow_uncommitted:
        require_file_sha256(path, expected, "production utility")
        return
    if not isinstance(source_bundle, dict):
        raise ValueError(
            "recorded manifest lacks historical production utility bundle provenance"
        )
    bundle_path = Path(str(source_bundle.get("path", ""))).resolve()
    require_file_sha256(bundle_path, source_bundle.get("sha256"), "source bundle")
    revisions = source_bundle.get("verified_revisions")
    if (
        not isinstance(revisions, list)
        or revision not in revisions
    ):
        raise ValueError(
            "recorded source bundle lacks the historical production utility revision"
        )
    with tempfile.TemporaryDirectory(prefix="cgl_lf_utility_verify_") as directory:
        repository = Path(directory) / "source.git"
        try:
            subprocess.run(
                [authenticated_git_binary(), "--no-replace-objects",
                 "clone", "--bare", "--quiet",
                 str(bundle_path), str(repository)],
                check=True, capture_output=True, text=True,
                env=hardened_git_environment(),
            )
        except subprocess.CalledProcessError as error:
            raise ValueError(
                f"recorded source bundle cannot be cloned: {bundle_path}"
            ) from error
        historical = subprocess.run(
            [
                authenticated_git_binary(), "--no-replace-objects",
                "-C", str(repository), "show",
                f"{revision}:{PRODUCTION_UTILITY_RELATIVE}",
            ],
            check=False, capture_output=True,
            env=hardened_git_environment(),
        )
    if (
        historical.returncode != 0
        or hashlib.sha256(historical.stdout).hexdigest() != expected
    ):
        raise ValueError(f"prepared production utility checksum has changed: {path}")


def read_build_provenance(executable: Path,
                          build_manifest: Path) -> dict[str, str]:
    """Verify an archived executable against its immutable build manifest."""

    digest_file = build_manifest / "athena.sha256"
    environment_file = build_manifest / "environment.txt"
    if not digest_file.is_file() or not environment_file.is_file():
        raise ValueError(f"incomplete build manifest directory: {build_manifest}")
    recorded_sha = digest_file.read_text(encoding="utf-8").split()[0]
    actual_sha = sha256(executable)
    if recorded_sha != actual_sha:
        raise ValueError("executable digest does not match build manifest")
    revision_match = re.search(
        r"(?m)^git_revision=([0-9a-f]{40})\s*$",
        environment_file.read_text(encoding="utf-8"),
    )
    if revision_match is None:
        raise ValueError("build manifest does not record a full git revision")
    return {
        "revision": revision_match.group(1),
        "sha256": actual_sha,
        "manifest_dir": str(build_manifest),
    }


@locked_root_action
def approve_qualification(args: argparse.Namespace) -> int:
    """Atomically retain reviewed corrected-build qualification approval."""

    if not args.confirm_corrected_build_frontier_qualified:
        raise ValueError("--confirm-corrected-build-frontier-qualified is required")
    if not args.approved_by.strip() or not args.review_notes.strip():
        raise ValueError("--approved-by and --review-notes must be nonempty")
    root = require_root(Path(args.root), args.allow_local_root)
    paths = initialize(root)
    require_reconciled_store_consistency(
        paths,
        allow_absent_qualification=True,
        allow_invalid_qualification=args.replace_existing_approval,
    )
    executable = Path(args.executable).expanduser().resolve()
    build_manifest = Path(args.build_manifest).expanduser().resolve()
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(
        build_manifest, root, "build manifest", args.allow_local_root
    )
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is unavailable: {executable}")
    provenance = read_build_provenance(executable, build_manifest)
    retained_manifests = sorted(
        paths["runs"].glob("*/*/manifest/prepared_run.json")
    )
    if retained_manifests:
        raise ValueError(
            "E03 qualification approval cannot change after segment preparation: "
            + ", ".join(str(path) for path in retained_manifests)
        )
    if paths["qualification"].exists() and not args.replace_existing_approval:
        raise ValueError(
            "E03 qualification approval token already exists; pass "
            "--replace-existing-approval only after reviewing the replacement build"
        )
    approval = {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "approval_scope": "corrected-build Frontier qualification for E03 prepare",
        "approved_utc": utc_now(),
        "approved_by": args.approved_by,
        "review_notes": args.review_notes,
        "approved_executable": str(executable),
        "approved_executable_sha256": provenance["sha256"],
        "approved_executable_revision": provenance["revision"],
        "build_manifest": provenance["manifest_dir"],
    }
    write_json(paths["qualification"], approval)
    refresh_summary(paths)
    print(f"Wrote E03 qualification approval token: {paths['qualification']}")
    print(f"Approved executable sha256: {provenance['sha256']}")
    print(f"Approved git revision: {provenance['revision']}")
    return 0


def quote(value: str | Path) -> str:
    """Quote a shell literal in the generated Slurm script."""

    return shlex.quote(str(value))


def expected_job_name(manifest: dict[str, object]) -> str:
    """Return the exact safe Slurm job name for one prepared segment."""

    run = manifest["run"]
    segment = require_safe_segment(str(run["segment"]))
    job_name = (
        f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_{run['case_id']}_{segment}"
    )
    job_name = re.sub(r"[^A-Za-z0-9_]+", "_", job_name)[:60]
    if re.fullmatch(r"[A-Za-z0-9_]+", job_name) is None:
        raise ValueError(f"generated Stage I job name is unsafe: {job_name}")
    return job_name


def normalized_batch_script_text(value: str) -> str:
    """Normalize the embedded self-digest before hashing a batch script."""

    return BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}", value
    )


def normalized_batch_script_sha256(path: Path) -> str:
    """Return the digest authenticated by a generated batch script itself."""

    return hashlib.sha256(
        normalized_batch_script_text(path.read_text(encoding="utf-8")).encode("utf-8")
    ).hexdigest()


def batch_script_binding(path: Path, profile: os.stat_result,
                         payload: bytes) -> dict[str, object]:
    """Return the exact stable identity retained across the sbatch boundary."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"prepared batch script is not UTF-8: {path}") from error
    return {
        "path": str(path),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "normalized_sha256": hashlib.sha256(
            normalized_batch_script_text(text).encode("utf-8")
        ).hexdigest(),
        "size_bytes": len(payload),
        "mode": f"{stat.S_IMODE(profile.st_mode):04o}",
        "links": profile.st_nlink,
        "device": profile.st_dev,
        "inode": profile.st_ino,
        "owner_uid": profile.st_uid,
    }


def validate_batch_script_binding(value: object, label: str) -> dict[str, object]:
    """Require one exact descriptor-bound batch-script identity."""

    value = require_r17_exact_keys(
        value,
        {
            "path", "sha256", "normalized_sha256", "size_bytes", "mode",
            "links", "device", "inode", "owner_uid",
        },
        label,
    )
    path = Path(require_r17_nonempty_string(value["path"], f"{label} path"))
    if (
        not path.is_absolute()
        or path != path.absolute()
        or ".." in path.parts
        or require_r17_sha256(value["sha256"], f"{label} SHA-256")
        != value["sha256"]
        or require_r17_sha256(
            value["normalized_sha256"], f"{label} normalized SHA-256"
        )
        != value["normalized_sha256"]
        or value["mode"] != "0750"
        or require_r17_integer(value["links"], f"{label} links", 1) != 1
        or any(
            require_r17_integer(value[key], f"{label} {key}") != value[key]
            for key in ("size_bytes", "device", "inode", "owner_uid")
        )
    ):
        raise ValueError(f"{label} differs")
    return value


def read_batch_script_descriptor(descriptor: int, path: Path,
                                 label: str) -> tuple[bytes, os.stat_result]:
    """Read and authenticate one already-open batch script."""

    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or stat.S_IMODE(before.st_mode) != 0o750
            or before.st_uid != os.geteuid()
            or before.st_nlink != 1
        ):
            raise ValueError(
                f"{label} must be an owner-controlled 0750 single-link file: {path}"
            )
        os.lseek(descriptor, 0, os.SEEK_SET)
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        after = os.fstat(descriptor)
    except OSError as error:
        raise ValueError(f"{label} cannot be read: {path}") from error
    stable = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if any(getattr(before, key) != getattr(after, key) for key in stable):
        raise ValueError(f"{label} changed while reading: {path}")
    payload = b"".join(blocks)
    if len(payload) != before.st_size:
        raise ValueError(f"{label} size changed while reading: {path}")
    return payload, before


def open_authenticated_batch_script(manifest: dict[str, object],
                                    manifest_path: Path) -> dict[str, object]:
    """Open the exact prepared script once for test-only and real submission."""

    command = manifest.get("command")
    paths = manifest.get("paths")
    if not isinstance(command, dict) or not isinstance(paths, dict):
        raise ValueError("prepared manifest lacks batch-script metadata")
    path = Path(str(paths.get("batch_script", ""))).expanduser().absolute()
    expected = (manifest_path.parent / "cgl_lf_stage_i.sbatch").absolute()
    if path != expected:
        raise ValueError("prepared batch script path is inconsistent")
    require_symlink_free_path(path, "prepared batch script")
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise ValueError(f"prepared batch script is unavailable: {path}") from error
    try:
        payload, profile = read_batch_script_descriptor(
            descriptor, path, "prepared batch script"
        )
        binding = batch_script_binding(path, profile, payload)
        text = payload.decode("utf-8")
        embedded = BATCH_SCRIPT_DIGEST_PATTERN.findall(text)
        if (
            embedded != [command.get("batch_script_sha256")]
            or binding["normalized_sha256"] != command.get("batch_script_sha256")
            or normalized_batch_script_text(text)
            != normalized_batch_script_text(generated_batch_script(manifest, manifest_path))
        ):
            raise ValueError("prepared batch script differs from retained launch intent")
        named = path.lstat()
        if (named.st_dev, named.st_ino) != (profile.st_dev, profile.st_ino):
            raise ValueError("prepared batch script pathname changed while opening")
    except BaseException:
        os.close(descriptor)
        raise
    return {
        "fd": descriptor,
        "descriptor_path": f"/proc/self/fd/{descriptor}",
        "path": path,
        "binding": binding,
    }


def reauthenticate_open_batch_script(authenticated: dict[str, object],
                                     expected_binding: object | None = None,
                                     ) -> dict[str, object]:
    """Recheck exact bytes and pathname identity for an open script."""

    try:
        descriptor = int(authenticated["fd"])
        path = Path(str(authenticated["path"]))
        retained = authenticated["binding"]
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("authenticated batch-script descriptor is invalid") from error
    retained = validate_batch_script_binding(retained, "authenticated batch script")
    payload, profile = read_batch_script_descriptor(
        descriptor, path, "authenticated batch script"
    )
    current = batch_script_binding(path, profile, payload)
    require_symlink_free_path(path, "authenticated batch script")
    try:
        named = path.lstat()
    except OSError as error:
        raise ValueError("authenticated batch-script pathname is unavailable") from error
    if (
        current != retained
        or (named.st_dev, named.st_ino) != (profile.st_dev, profile.st_ino)
        or (
            expected_binding is not None
            and current
            != validate_batch_script_binding(
                expected_binding, "submission-journal batch script"
            )
        )
    ):
        raise ValueError("authenticated batch script changed across submission")
    return current


def close_authenticated_batch_script(authenticated: object) -> None:
    """Close one controller-owned batch-script descriptor."""

    if isinstance(authenticated, dict) and isinstance(authenticated.get("fd"), int):
        os.close(authenticated["fd"])
        authenticated["fd"] = None


def finalize_batch_script(value: str) -> tuple[str, str]:
    """Embed a normalized self-digest in one generated batch script."""

    matches = BATCH_SCRIPT_DIGEST_PATTERN.findall(value)
    if matches != [BATCH_SCRIPT_DIGEST_PLACEHOLDER]:
        raise ValueError("generated batch script lacks one self-digest placeholder")
    digest = hashlib.sha256(normalized_batch_script_text(value).encode("utf-8")).hexdigest()
    return value.replace(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}",
        f"BATCH_SCRIPT_SHA256={digest}",
        1,
    ), digest


def require_file_sha256(path: Path, expected: object, label: str) -> None:
    """Require one retained prepared artifact to preserve its digest."""

    if not path.is_file():
        raise ValueError(f"prepared {label} is missing: {path}")
    if not isinstance(expected, str) or sha256(path) != expected:
        raise ValueError(f"prepared {label} checksum has changed: {path}")


def prepared_restart_inventory(manifest_path: Path,
                               command: dict[str, object]) -> list[Path]:
    """Return the exact retained restart archive inventory."""

    records = command.get("restart_files")
    if not isinstance(records, list):
        raise ValueError("prepared manifest lacks restart sibling metadata")
    paths = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("prepared restart sibling metadata is invalid")
        paths.append(Path(str(record.get("path", ""))).resolve())
    restart = command.get("restart_file")
    if restart is None:
        if paths:
            raise ValueError("prepared manifest retains restart siblings without a restart")
        return []
    restart_path = Path(str(restart)).resolve()
    if restart_path not in paths:
        raise ValueError("prepared primary restart is absent from retained siblings")
    archive_root = manifest_path.parent / "submitted_restart"
    if archive_root.is_dir():
        actual = sorted(path.resolve() for path in archive_root.rglob("*") if path.is_file())
    else:
        actual = [restart_path] if restart_path.is_file() else []
    if sorted(paths) != actual:
        raise ValueError("prepared restart archive inventory has changed")
    return paths


def restart_parameter_dump(path: Path) -> tuple[str, int]:
    """Read restart parameter text and return its exact binary payload offset."""

    marker = b"<par_end>\n"
    header = b""
    with path.open("rb") as stream:
        while marker not in header and len(header) < MAX_RESTART_PARAMETER_DUMP_BYTES:
            block = stream.read(
                min(4096, MAX_RESTART_PARAMETER_DUMP_BYTES - len(header))
            )
            if not block:
                break
            header += block
    end = header.find(marker)
    if end < 0:
        raise ValueError(
            f"restart parameter dump lacks loadable <par_end> terminator: {path}"
        )
    payload_offset = end + len(marker)
    try:
        text = header[:end].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"restart parameter dump is not UTF-8 text: {path}") from error
    return text, payload_offset


def restart_time_marker_text(path: Path) -> str:
    """Read the explicit physical-time marker text from a restart parameter dump."""

    text, _ = restart_parameter_dump(path)
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
        raise ValueError(
            f"restart parameter dump must contain one time/restart_time marker: {path}"
        )
    return markers[0]


def restart_time_marker(path: Path) -> float:
    """Read the explicit physical-time marker from a restart parameter dump."""

    try:
        result = float(restart_time_marker_text(path))
    except ValueError as error:
        raise ValueError(f"restart time marker is not numeric: {path}") from error
    if not math.isfinite(result):
        raise ValueError(f"restart time marker is not finite: {path}")
    return result


def restart_product_time(paths: list[Path],
                         allow_missing_marker: bool = False) -> float | None:
    """Require every selected restart sibling to retain one physical time."""

    try:
        times = [restart_time_marker(path) for path in paths]
    except (OSError, ValueError):
        if allow_missing_marker:
            return None
        raise
    if not times:
        raise ValueError("restart product has no selected siblings")
    if any(abs(value - times[0]) > 1.0e-12 for value in times[1:]):
        raise ValueError("restart sibling physical-time markers disagree")
    return times[0]


def restart_binary_abi(command: object) -> dict[str, object]:
    """Return the qualified binary-restart ABI for one prepared executable."""

    if not isinstance(command, dict):
        raise ValueError("prepared manifest lacks command metadata")
    key = (
        command.get("executable_revision"),
        command.get("executable_sha256"),
    )
    abi = QUALIFIED_RESTART_BINARY_ABIS.get(key)
    if abi is None:
        raise ValueError("prepared executable has no qualified binary-restart ABI")
    return abi


def restart_binary_time_encoding(path: Path, command: object) -> tuple[bytes, str]:
    """Read encoded Mesh::time from a restart file using its qualified ABI."""

    abi = restart_binary_abi(command)
    _, payload_offset = restart_parameter_dump(path)
    offset = payload_offset + int(abi["mesh_time_offset_after_parameter_dump"])
    value_format = str(abi["mesh_time_format"])
    value_size = struct.calcsize(value_format)
    with path.open("rb") as stream:
        stream.seek(offset)
        encoded = stream.read(value_size)
    if len(encoded) != value_size:
        raise ValueError(f"restart binary header is truncated: {path}")
    return encoded, value_format


def restart_binary_time(path: Path, command: object) -> float:
    """Read Mesh::time from a restart file using its qualified executable ABI."""

    encoded, value_format = restart_binary_time_encoding(path, command)
    value = float(struct.unpack(value_format, encoded)[0])
    if not math.isfinite(value):
        raise ValueError(f"restart binary physical time is not finite: {path}")
    return value


def authenticated_restart_product_time(
    paths: list[Path], command: object,
) -> dict[str, object]:
    """Bind restart marker text to synchronized binary Mesh::time evidence."""

    if not paths:
        raise ValueError("restart product has no selected siblings")
    abi = restart_binary_abi(command)
    encodings = [restart_binary_time_encoding(path, command) for path in paths]
    if any(encoded != encodings[0] for encoded in encodings[1:]):
        raise ValueError("restart sibling binary physical times disagree")
    binary_times = []
    for encoded, value_format in encodings:
        value = float(struct.unpack(value_format, encoded)[0])
        if not math.isfinite(value):
            raise ValueError("restart binary physical time is not finite")
        binary_times.append(value)
    marker_modes = []
    for path, binary_time in zip(paths, binary_times):
        marker_text = restart_time_marker_text(path)
        try:
            marker_time = float(marker_text)
        except ValueError as error:
            raise ValueError(f"restart time marker is not numeric: {path}") from error
        if not math.isfinite(marker_time):
            raise ValueError(f"restart time marker is not finite: {path}")
        if marker_text == format(binary_time, ".17g"):
            marker_modes.append("full_precision")
        elif marker_text == format(binary_time, ".6g"):
            marker_modes.append("legacy_default_precision")
        else:
            raise ValueError(
                f"restart time marker does not authenticate binary physical time: {path}"
            )
    if any(mode not in abi["allowed_marker_modes"] for mode in marker_modes):
        raise ValueError("restart time marker mode is not qualified for this executable")
    if any(mode != marker_modes[0] for mode in marker_modes[1:]):
        raise ValueError("restart sibling physical-time marker modes disagree")
    return {
        "binary_time": binary_times[0],
        "marker_modes": marker_modes,
    }


def validate_prepared_resources(manifest: dict[str, object],
                                canonical_production: bool) -> None:
    """Require retained Slurm and Athena resource policy to remain valid."""

    allocation = manifest.get("allocation")
    command = manifest.get("command")
    if not isinstance(allocation, dict) or not isinstance(command, dict):
        raise ValueError("prepared manifest lacks resource metadata")
    try:
        nodes = int(allocation["nodes"])
        requested_seconds = parse_walltime(str(allocation["requested_walltime"]))
        athena_seconds = parse_walltime(str(command["athena_walltime"]))
        ranks_per_node = int(allocation["ranks_per_node"])
        cpus_per_task = int(allocation["cpus_per_task"])
        reserved_node_hours = float(allocation["reserved_node_hours"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("prepared manifest has invalid resource metadata") from error
    if (
        nodes < 1
        or requested_seconds <= 0
        or athena_seconds <= 0
        or not math.isfinite(reserved_node_hours)
        or reserved_node_hours <= 0.0
    ):
        raise ValueError("prepared manifest resource values must be positive")
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError("prepared Slurm walltime exceeds the Stage I limit")
    if athena_seconds >= requested_seconds:
        raise ValueError("prepared Athena walltime is not shorter than Slurm walltime")
    if athena_seconds > requested_seconds - 600:
        raise ValueError("prepared Athena walltime lacks the shutdown margin")
    if (
        allocation.get("requested_seconds") != requested_seconds
        or abs(reserved_node_hours - node_hours(nodes, requested_seconds)) > 5.0e-12
    ):
        raise ValueError("prepared resource accounting is inconsistent")
    if canonical_production and (
        ranks_per_node != EXPECTED_RANKS_PER_NODE
        or cpus_per_task != EXPECTED_CPUS_PER_TASK
    ):
        raise ValueError("canonical Frontier resource shape is inconsistent")
    if canonical_production:
        run = manifest.get("run")
        if not isinstance(run, dict):
            raise ValueError("prepared manifest lacks run metadata")
        require_case_node_count(str(run.get("case_id")), nodes)


def require_reservation_matches_manifest(reservation: dict[str, object],
                                         manifest: dict[str, object]) -> None:
    """Bind one canonical reservation snapshot to its manifest allocation."""

    allocation = manifest.get("allocation")
    run = manifest.get("run")
    if not isinstance(allocation, dict) or not isinstance(run, dict):
        raise ValueError("transaction manifest lacks reservation metadata")
    try:
        matches = (
            reservation.get("case_id") == run["case_id"]
            and reservation.get("case_name") == run["case_name"]
            and reservation.get("segment") == run["segment"]
            and int(reservation["nodes"]) == int(allocation["nodes"])
            and reservation.get("requested_walltime")
            == allocation["requested_walltime"]
            and abs(
                float(reservation["reserved_node_hours"])
                - float(allocation["reserved_node_hours"])
            ) <= 5.0e-12
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "transaction manifest has invalid reservation metadata"
        ) from error
    if not matches:
        raise ValueError("transaction reservation allocation differs from manifest")


def historical_submitted_utility_transition_authorized(
    manifest: dict[str, object],
) -> bool:
    """Admit only the reviewed pre-promotion R03 submission."""

    command = manifest.get("command")
    run = manifest.get("run")
    if not isinstance(command, dict) or not isinstance(run, dict):
        return False
    utility = command.get("production_utility")
    bundle = command.get("source_bundle")
    if not isinstance(utility, dict) or not isinstance(bundle, dict):
        return False
    actual = {
        "project_root": manifest.get("project_root"),
        "state": manifest.get("state"),
        "job_id": manifest.get("job_id"),
        "case_id": run.get("case_id"),
        "segment": run.get("segment"),
        "production_utility_revision": utility.get("revision"),
        "production_utility_sha256": utility.get("sha256"),
        "source_bundle_sha256": bundle.get("sha256"),
        "executable_revision": command.get("executable_revision"),
        "executable_sha256": command.get("executable_sha256"),
    }
    return actual == HISTORICAL_SUBMITTED_R03_UTILITY_TRANSITION


def authenticate_prepared_execution(manifest: dict[str, object],
                                    manifest_path: Path,
                                    allow_legacy_local: bool = False) -> None:
    """Authenticate every prepared artifact used by one production launch."""

    command = manifest.get("command")
    paths = manifest.get("paths")
    if allow_legacy_local and (
        not isinstance(command, dict) or not isinstance(paths, dict)
    ):
        return
    if not isinstance(command, dict) or not isinstance(paths, dict):
        raise ValueError("prepared manifest lacks execution metadata")
    required = (
        "batch_script_sha256",
        "input_file",
        "input_sha256",
        "matrix_file",
        "matrix_sha256",
        "executable",
        "executable_sha256",
        "restart_files",
        "production_utility",
    )
    if allow_legacy_local and any(key not in command for key in required):
        return
    overrides = command.get("overrides")
    if not allow_legacy_local and (
        not isinstance(overrides, list)
        or len(overrides) != 1
        or not isinstance(overrides[0], str)
        or not overrides[0].startswith("time/tlim=")
    ):
        raise ValueError(
            "canonical prepared manifest requires exactly one time/tlim override"
        )
    if not allow_legacy_local or isinstance(overrides, list):
        prepared_time_tlim_target(manifest)
    validate_prepared_continuation_target(manifest)
    validate_prepared_resources(manifest, canonical_production=not allow_legacy_local)
    bundle = command.get("source_bundle")
    recorded = manifest.get("state") == "recorded"
    reviewed_historical_submission = (
        historical_submitted_utility_transition_authorized(manifest)
    )
    authenticate_production_utility(
        command.get("production_utility"),
        source_bundle=bundle,
        allow_uncommitted=allow_legacy_local,
        allow_historical=recorded or reviewed_historical_submission,
    )
    batch_script = Path(str(paths.get("batch_script", ""))).resolve()
    if batch_script != (manifest_path.parent / "cgl_lf_stage_i.sbatch").resolve():
        raise ValueError("prepared batch script path is inconsistent")
    if not batch_script.is_file():
        raise ValueError(f"prepared batch script is missing: {batch_script}")
    script_text = batch_script.read_text(encoding="utf-8")
    embedded = BATCH_SCRIPT_DIGEST_PATTERN.findall(script_text)
    if embedded != [command.get("batch_script_sha256")]:
        raise ValueError("prepared batch script self-digest is inconsistent")
    if normalized_batch_script_sha256(batch_script) != command.get(
        "batch_script_sha256"
    ):
        raise ValueError("prepared batch script normalized checksum has changed")
    if not recorded and not reviewed_historical_submission:
        expected_script = generated_batch_script(manifest, manifest_path)
        if normalized_batch_script_text(script_text) != normalized_batch_script_text(
            expected_script
        ):
            raise ValueError("prepared batch script differs from retained launch intent")
    require_file_sha256(
        Path(str(command.get("input_file", ""))).resolve(),
        command.get("input_sha256"),
        "input",
    )
    require_file_sha256(
        Path(str(command.get("matrix_file", ""))).resolve(),
        command.get("matrix_sha256"),
        "matrix",
    )
    executable = Path(str(command.get("executable", ""))).resolve()
    require_file_sha256(executable, command.get("executable_sha256"), "executable")
    if not os.access(executable, os.X_OK):
        raise ValueError(f"prepared executable is not executable: {executable}")
    qualification = command.get("qualification_approval")
    if qualification is None:
        if not allow_legacy_local:
            raise ValueError("prepared manifest lacks E03 qualification approval")
    elif not isinstance(qualification, dict):
        raise ValueError("prepared E03 qualification approval metadata is invalid")
    else:
        approval_path = Path(str(qualification.get("path", ""))).resolve()
        expected_path = layout(
            Path(str(manifest.get("project_root", ""))).resolve()
        )["qualification"].resolve()
        if approval_path != expected_path:
            raise ValueError("prepared E03 qualification token path is inconsistent")
        require_file_sha256(
            approval_path, qualification.get("sha256"), "E03 qualification token"
        )
        approval = read_qualification_approval(approval_path)
        if (
            qualification.get("execution_epoch") != EXECUTION_EPOCH
            or approval["approved_executable_sha256"]
            != command.get("executable_sha256")
            or approval["approved_executable_revision"]
            != command.get("executable_revision")
            or qualification.get("approved_executable_sha256")
            != approval["approved_executable_sha256"]
            or qualification.get("approved_executable_revision")
            != approval["approved_executable_revision"]
            or qualification.get("token") != approval
        ):
            raise ValueError(
                "prepared E03 qualification approval does not match the executable"
            )
    if bundle is None:
        if not allow_legacy_local:
            raise ValueError("prepared manifest lacks source bundle provenance")
    elif not isinstance(bundle, dict):
        raise ValueError("prepared source bundle metadata is invalid")
    else:
        revisions = bundle.get("verified_revisions")
        expected_revisions = {
            command.get("input_revision"),
            command.get("executable_revision"),
            command["production_utility"].get("revision"),
        }
        if (
            not isinstance(revisions, list)
            or not all(isinstance(item, str) for item in revisions)
            or not expected_revisions.issubset(set(revisions))
        ):
            raise ValueError(
                "prepared source bundle lacks launch provenance revisions"
            )
        require_file_sha256(
            Path(str(bundle.get("path", ""))).resolve(),
            bundle.get("sha256"),
            "source bundle",
        )
    records = command.get("restart_files")
    if not isinstance(records, list):
        raise ValueError("prepared manifest lacks restart sibling metadata")
    for record in records:
        revalidate_retained_file(record, label="prepared restart sibling")
    restart_paths = prepared_restart_inventory(manifest_path, command)
    if restart_paths:
        marker_bypass = (
            allow_legacy_local
            and command.get("allow_missing_restart_time_marker") is True
        )
        if (
            command.get("allow_missing_restart_time_marker") is True
            and not allow_legacy_local
        ):
            raise ValueError(
                "canonical prepared restart may not bypass physical-time markers"
            )
        parent = command.get("parent_segment")
        if not isinstance(parent, dict):
            if not allow_legacy_local:
                raise ValueError("prepared restart lacks inspected parent metadata")
        else:
            marker = (
                restart_product_time(restart_paths, allow_missing_marker=True)
                if marker_bypass or allow_legacy_local
                else float(
                    authenticated_restart_product_time(
                        restart_paths, command
                    )["binary_time"]
                )
            )
            try:
                parent_final_time = float(parent["final_time"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError("prepared restart parent lacks final time") from error
            if marker is not None and abs(marker - parent_final_time) > 1.0e-12:
                raise ValueError(
                    "prepared restart physical time differs from parent inspection"
                )


def generated_batch_script(manifest: dict[str, object],
                           manifest_path: Path) -> str:
    """Generate an authenticated normal-QOS production segment."""

    run = manifest["run"]
    allocation = manifest["allocation"]
    command = manifest["command"]
    paths = manifest["paths"]
    overrides = " ".join(quote(value) for value in command["overrides"])
    if overrides:
        overrides = " " + overrides
    restart = command.get("restart_file")
    restart_literal = quote(str(restart)) if restart else "''"
    job_name = expected_job_name(manifest)
    bundle = command.get("source_bundle")
    runtime_checks = [
        f"require_sha {quote(command['production_utility']['sha256'])} "
        f"{quote(command['production_utility']['path'])} production_utility",
        f"require_sha {quote(command['input_sha256'])} \"${{INPUT}}\" input",
        f"require_sha {quote(command['matrix_sha256'])} {quote(command['matrix_file'])} matrix",
        f"require_sha {quote(command['executable_sha256'])} \"${{ATHENA}}\" executable",
    ]
    if isinstance(bundle, dict):
        runtime_checks.append(
            f"require_sha {quote(bundle['sha256'])} {quote(bundle['path'])} source_bundle"
        )
    qualification = command.get("qualification_approval")
    if isinstance(qualification, dict):
        runtime_checks.append(
            f"require_sha {quote(qualification['sha256'])} "
            f"{quote(qualification['path'])} qualification_approval"
        )
    restart_records = command.get("restart_files", [])
    if not isinstance(restart_records, list):
        raise ValueError("prepared restart sibling metadata is invalid")
    for index, record in enumerate(restart_records):
        runtime_checks.append(
            f"require_sha {quote(record['sha256'])} {quote(record['path'])} "
            f"restart_{index:04d}"
        )
    archive_root = manifest_path.parent / "submitted_restart"
    if restart_records and archive_root.is_dir():
        runtime_checks.append(
            f'test "$(find {quote(archive_root)} -type f -print | wc -l)" '
            f'-eq {len(restart_records)} || {{ echo "restart inventory changed" >&2; '
            "exit 1; }"
        )
    runtime_checks_text = "\n".join(runtime_checks)
    return f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
# Frontier default normal QOS is required for production; do not add -q debug.
#SBATCH -t {allocation["requested_walltime"]}
#SBATCH -N {allocation["nodes"]}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail

RUN_MANIFEST={quote(manifest_path)}
BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}
ATHENA={quote(command["executable"])}
INPUT={quote(command["input_file"])}
RESTART={restart_literal}
OUT_DIR={quote(paths["output_dir"])}
ENV_LOG={quote(paths["environment_log"])}
RANKS_PER_NODE={allocation["ranks_per_node"]}
CPUS_PER_TASK={allocation["cpus_per_task"]}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * RANKS_PER_NODE))"

require_sha() {{
  local expected="$1"
  local path="$2"
  local label="$3"
  local actual
  test -f "${{path}}" || {{ echo "missing ${{label}}: ${{path}}" >&2; exit 1; }}
  actual="$(sha256sum "${{path}}" | awk '{{print $1}}')"
  test "${{actual}}" = "${{expected}}" || {{
    echo "checksum mismatch for ${{label}}: ${{path}}" >&2
    exit 1
  }}
}}

normalized_script_sha256() {{
  sed -E 's/^BATCH_SCRIPT_SHA256=[0-9a-f]{{64}}$/BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}/' "$0" |
    sha256sum | awk '{{print $1}}'
}}

test -f "${{RUN_MANIFEST}}"
test -x "${{ATHENA}}"
test -f "${{INPUT}}"
if [[ -n "${{RESTART}}" ]]; then
  test -f "${{RESTART}}"
fi
test "$(normalized_script_sha256)" = "${{BATCH_SCRIPT_SHA256}}" || {{
  echo "batch script checksum mismatch" >&2
  exit 1
}}
{runtime_checks_text}
mkdir -p "${{OUT_DIR}}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1

{{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${{SLURM_JOB_ID:?Missing SLURM_JOB_ID}}"
  echo "prepared_manifest=${{RUN_MANIFEST}}"
  echo "nodes=${{NNODES}}"
  echo "ranks=${{NRANKS}}"
  module -t list 2>&1
  env | LC_ALL=C sort | grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)' || true
}} > "${{ENV_LOG}}"

RUN_ARGS=(-i "${{INPUT}}")
if [[ -n "${{RESTART}}" ]]; then
  RUN_ARGS=(-r "${{RESTART}}")
fi

srun -N "${{NNODES}}" -n "${{NRANKS}}" --ntasks-per-node="${{RANKS_PER_NODE}}" \
  -c "${{CPUS_PER_TASK}}" --threads-per-core=1 --cpu-bind=threads \
  --gpus-per-task=1 --gpu-bind=closest \
  "${{ATHENA}}" "${{RUN_ARGS[@]}}" -d "${{OUT_DIR}}" \
  -t {quote(command["athena_walltime"])} \
  job/basename={quote(run["run_basename"])}{overrides}

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${{ENV_LOG}}"
"""


def execution_intent_sha256(manifest: dict[str, object]) -> str:
    """Digest the immutable prepared execution intent across lifecycle updates."""

    try:
        intent = {
            key: manifest[key]
            for key in (
                "schema_version", "execution_epoch", "project_root", "policy",
                "run", "allocation", "command", "paths",
            )
        }
    except KeyError as error:
        raise ValueError(
            f"prepared manifest lacks immutable execution intent field {error.args[0]}"
        ) from error
    return hashlib.sha256(
        json.dumps(intent, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def require_reserved_execution_intent(reservation: dict[str, object],
                                      manifest: dict[str, object],
                                      allow_legacy_local: bool) -> None:
    """Require reservation metadata to bind the immutable prepared launch."""

    expected = reservation.get("execution_intent_sha256")
    if expected is None and allow_legacy_local:
        return
    if (
        not isinstance(expected, str)
        or expected != execution_intent_sha256(manifest)
    ):
        raise ValueError("prepared execution intent differs from its reservation")


@locked_root_action
def prepare(args: argparse.Namespace) -> Path:
    """Create one retained production segment under the bounded overlap policy."""

    root = require_root(Path(args.root), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    if offline_local_root:
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    require_authorized_case(args.case_id)
    require_safe_segment(args.segment)
    reservations = read_reservations(paths)
    if isinstance(args.nodes, bool) or args.nodes < 1:
        raise ValueError("--nodes must be positive")
    require_active_reservation_policy(reservations, args.case_id, args.nodes)
    require_prepare_case_policy(
        paths, args.case_id, args.nodes, offline_local_root
    )
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    executable = Path(args.executable).expanduser().resolve()
    build_manifest = Path(args.build_manifest).expanduser().resolve()
    restart = (
        Path(args.restart_file).expanduser().resolve()
        if args.restart_file else None
    )
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(
        build_manifest, root, "build manifest", args.allow_local_root
    )
    if restart is not None:
        require_beneath_root(restart, root, "restart file", args.allow_local_root)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is unavailable: {executable}")
    matrix = validate_matrix(matrix_path, source_dir)
    case = case_for_id(matrix, args.case_id)
    input_path = source_dir / str(case["input"])
    allow_missing_time_target = getattr(
        args, "allow_missing_time_target", False
    )
    if allow_missing_time_target and not offline_local_root:
        raise ValueError(
            "--allow-missing-time-target is restricted to offline validation"
        )
    time_tlim_target = validate_prepare_overrides(
        input_path, args.override,
        allow_missing_time_target=allow_missing_time_target,
        canonical_production=not offline_local_root,
    )
    input_revision = git_revision_for_input(source_dir, input_path, matrix_path)
    utility_provenance = production_utility_provenance(
        allow_uncommitted=offline_local_root
    )
    provenance = read_build_provenance(executable, build_manifest)
    qualification_approval = require_qualification_approval(
        paths, provenance["sha256"], provenance["revision"], offline_local_root
    )
    bundle_provenance = source_bundle_provenance(
        args.source_bundle,
        list(dict.fromkeys([
            input_revision,
            provenance["revision"],
            utility_provenance["revision"],
        ])),
        root,
        offline_local_root,
    )
    current_source_authority = require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=bundle_provenance,
        utility_provenance=utility_provenance,
        matrix_path=matrix_path,
        input_revision=input_revision,
        offline_local_root=offline_local_root,
    )
    if restart is not None and not restart.is_file():
        raise ValueError(f"restart file is unavailable: {restart}")
    allow_missing_restart_time_marker = getattr(
        args, "allow_missing_restart_time_marker", False
    )
    if allow_missing_restart_time_marker and not offline_local_root:
        raise ValueError(
            "--allow-missing-restart-time-marker is restricted to offline validation"
        )
    parent_segment = (
        verify_continuation_restart(
            restart,
            allow_missing_restart_time_marker=allow_missing_restart_time_marker,
        )
        if restart else None
    )
    if parent_segment is not None:
        if parent_segment["case_id"] != args.case_id:
            raise ValueError("continuation restart belongs to a different case")
        if parent_segment["input_sha256"] != sha256(input_path):
            raise ValueError("continuation restart input differs from its parent")
        if parent_segment["executable_sha256"] != provenance["sha256"]:
            raise ValueError("continuation executable differs from its parent")
        if (
            time_tlim_target is not None
            and time_tlim_target <= float(parent_segment["final_time"]) + 1.0e-12
        ):
            raise ValueError(
                "continuation time/tlim target must advance beyond its parent "
                f"inspection time {float(parent_segment['final_time']):.12g}"
            )
    requested_seconds = parse_walltime(args.walltime)
    if requested_seconds <= 0:
        raise ValueError("--walltime must be positive")
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError(
            "a Stage I normal-QOS segment may not request more than two hours"
        )
    athena_seconds = parse_walltime(args.athena_walltime)
    if athena_seconds <= 0:
        raise ValueError("--athena-walltime must be positive")
    if athena_seconds >= requested_seconds:
        raise ValueError("Athena walltime must be shorter than the Slurm walltime")
    if athena_seconds > requested_seconds - 600:
        raise ValueError("Athena walltime must leave ten minutes for shutdown")
    if not offline_local_root and (
        args.ranks_per_node != EXPECTED_RANKS_PER_NODE
        or args.cpus_per_task != EXPECTED_CPUS_PER_TASK
    ):
        raise ValueError(
            "canonical Frontier production requires --ranks-per-node=8 "
            "and --cpus-per-task=7"
        )
    promoted_profile_authority = require_promoted_profile_for_prepare(
        paths,
        args,
        source_dir=source_dir,
        matrix_path=matrix_path,
        input_path=input_path,
        input_revision=input_revision,
        utility_provenance=utility_provenance,
        build_provenance=provenance,
        build_manifest=build_manifest,
        qualification_approval=qualification_approval,
        bundle_provenance=bundle_provenance,
        current_source_authority=current_source_authority,
        restart=restart,
        parent_segment=parent_segment,
        time_tlim_target=time_tlim_target,
    )
    r17_readiness_evidence_chain = require_r17_readiness_for_prepare(
        paths,
        args,
        reservations,
        source_dir=source_dir,
        matrix_path=matrix_path,
        input_path=input_path,
        input_revision=input_revision,
        utility_provenance=utility_provenance,
        build_provenance=provenance,
        build_manifest=build_manifest,
        qualification_approval=qualification_approval,
        bundle_provenance=bundle_provenance,
        current_source_authority=current_source_authority,
        restart=restart,
        parent_segment=parent_segment,
        time_tlim_target=time_tlim_target,
    )
    segment_hours = node_hours(args.nodes, requested_seconds)
    actual, reserved = reservation_usage(paths)
    if actual + reserved + segment_hours > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("proposed segment exceeds the Stage I reservation")
    if (
        actual + reserved + segment_hours > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError("proposed segment exceeds the project ceiling")
    run_dir = paths["runs"] / args.case_id / args.segment
    manifest_dir = run_dir / "manifest"
    manifest_path = manifest_dir / "prepared_run.json"
    if run_dir.exists():
        raise ValueError(f"segment run directory already exists: {run_dir}")
    output_dir = run_dir / "output"
    for path in (manifest_dir, output_dir):
        mkdir_durable(path)
    archived_input = manifest_dir / "submitted_input.athinput"
    archived_matrix = manifest_dir / "mks24_stage_i_manifest.json"
    copy_file(input_path, archived_input)
    copy_file(matrix_path, archived_matrix)
    archived_restart = None
    archived_restart_files: list[Path] = []
    if restart is not None:
        source_restart_files = (
            [Path(str(path)) for path in parent_segment["restart_files"]]
            if parent_segment is not None
            else [restart]
        )
        if restart.parent.name.startswith("rank_"):
            archive_root = manifest_dir / "submitted_restart"
            for source in source_restart_files:
                target = archive_root / source.parent.name / source.name
                mkdir_durable(target.parent)
                copy_file(source, target)
                archived_restart_files.append(target)
            archived_restart = (
                archive_root / restart.parent.name / restart.name
            )
        else:
            archived_restart = manifest_dir / "submitted_restart.rst"
            copy_file(restart, archived_restart)
            archived_restart_files.append(archived_restart)
    batch_script = manifest_dir / "cgl_lf_stage_i.sbatch"
    manifest: dict[str, object] = {
        "schema_version": 3,
        "execution_epoch": EXECUTION_EPOCH,
        "state": "prepared",
        "prepared_utc": utc_now(),
        "project_root": str(root),
        "policy": {
            "account": ACCOUNT,
            "partition": PARTITION,
            "qos": PRODUCTION_QOS,
            "project_budget_node_hours": PROJECT_BUDGET_NODE_HOURS,
            "historical_debug_node_hours": HISTORICAL_DEBUG_NODE_HOURS,
            "historical_e01_stage_i_node_hours": HISTORICAL_E01_STAGE_I_NODE_HOURS,
            "stage_i_reserved_node_hours": CURRENT_STAGE_I_RESERVED_NODE_HOURS,
            "stage_i_authorization": (
                "frozen mapped Stage I matrix R02-R17 with bounded distinct-case "
                "overlap and exclusive R17"
            ),
            "atomic_submission_required": True,
        },
        "run": {
            "case_id": args.case_id,
            "case_name": case["name"],
            "segment": args.segment,
            "run_basename": f"{EXECUTION_EPOCH_SLUG}_{case['name']}_{args.segment}",
            "resolution": case["resolution"],
            "figure_roles": case["figure_roles"],
            "acceptance_criterion": args.acceptance_criterion,
        },
        "allocation": {
            "nodes": args.nodes,
            "requested_walltime": args.walltime,
            "requested_seconds": requested_seconds,
            "reserved_node_hours": segment_hours,
            "ranks_per_node": args.ranks_per_node,
            "cpus_per_task": args.cpus_per_task,
        },
        "command": {
            "production_utility": utility_provenance,
            "qualification_approval": qualification_approval,
            "source_dir": str(source_dir),
            "source_bundle": bundle_provenance,
            "input_revision": input_revision,
            "source_input_file": str(input_path),
            "input_file": str(archived_input),
            "input_sha256": sha256(archived_input),
            "matrix_file": str(archived_matrix),
            "matrix_sha256": sha256(archived_matrix),
            "executable": str(executable),
            "executable_revision": provenance["revision"],
            "executable_sha256": provenance["sha256"],
            "build_manifest": provenance["manifest_dir"],
            "source_restart_file": str(restart) if restart else None,
            "restart_file": str(archived_restart) if archived_restart else None,
            "restart_sha256": sha256(archived_restart) if archived_restart else None,
            "restart_files": [
                retained_file(path) for path in archived_restart_files
            ],
            "parent_segment": parent_segment,
            "allow_missing_restart_time_marker": (
                allow_missing_restart_time_marker
            ),
            "athena_walltime": args.athena_walltime,
            "overrides": args.override,
            "time_tlim_target": time_tlim_target,
        },
        "paths": {
            "run_dir": str(run_dir),
            "output_dir": str(output_dir),
            "environment_log": str(manifest_dir / "run_environment.txt"),
            "batch_script": str(batch_script),
            "slurm_log": str(paths["logs_slurm"] / "%x.%j.log"),
        },
    }
    if r17_readiness_evidence_chain is not None:
        manifest["command"]["r17_readiness_evidence_chain"] = (
            r17_readiness_evidence_chain
        )
    if promoted_profile_authority is not None:
        manifest["command"]["promoted_profile_authority"] = (
            promoted_profile_authority
        )
    if current_source_authority is not None:
        manifest["command"]["current_source_authority"] = current_source_authority
    script, script_sha256 = finalize_batch_script(
        generated_batch_script(manifest, manifest_path)
    )
    manifest["command"]["batch_script_sha256"] = script_sha256
    write_text(batch_script, script, mode=0o750)
    authenticate_prepared_execution(
        manifest, manifest_path, allow_legacy_local=offline_local_root
    )
    reservations = read_reservations(paths)
    reservations.append({
        "execution_epoch": EXECUTION_EPOCH,
        "manifest": str(manifest_path),
        "case_id": args.case_id,
        "case_name": case["name"],
        "segment": args.segment,
        "nodes": args.nodes,
        "requested_walltime": args.walltime,
        "reserved_node_hours": segment_hours,
        "execution_intent_sha256": execution_intent_sha256(manifest),
        "state": "prepared",
        "prepared_utc": manifest["prepared_utc"],
    })
    durable_transition(paths, "prepared", manifest_path, manifest, reservations)
    print(f"Prepared Stage I segment: {manifest_path}")
    print(f"Reserved node-hours: {segment_hours:.6f}")
    return manifest_path


def read_manifest(path: Path) -> dict[str, object]:
    """Read one segment preparation manifest."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid run manifest: {path}")
    return value


def verify_continuation_restart(
    restart: Path,
    allow_missing_restart_time_marker: bool = False,
) -> dict[str, object]:
    """Require an inspected production parent for a continuation restart."""

    parent_manifest_path = None
    for ancestor in restart.parents:
        candidate = ancestor / "manifest" / "prepared_run.json"
        if candidate.is_file():
            parent_manifest_path = candidate
            break
    if parent_manifest_path is None:
        raise ValueError("restart does not belong to a retained Stage I segment")
    parent = read_manifest(parent_manifest_path)
    require_current_epoch(parent, "continuation parent")
    parent_root = Path(str(parent.get("project_root", ""))).resolve()
    if (
        allow_missing_restart_time_marker
        and parent_root == DEFAULT_ROOT.expanduser().resolve()
    ):
        raise ValueError(
            "canonical continuation may not bypass restart physical-time markers"
        )
    authenticate_prepared_execution(
        parent, parent_manifest_path,
        allow_legacy_local=parent_root != DEFAULT_ROOT.expanduser().resolve(),
    )
    if parent_root == DEFAULT_ROOT.expanduser().resolve():
        parent_paths = layout(parent_root)
        require_existing_layout(parent_paths)
        require_reserved_execution_intent(
            reservation_for_manifest(
                read_reservations(parent_paths), parent_manifest_path
            ),
            parent,
            allow_legacy_local=False,
        )
    accounting = parent.get("accounting", {})
    inspection = parent.get("scientific_inspection", {})
    if (
        parent.get("state") != "recorded"
        or not isinstance(accounting, dict)
        or accounting.get("result") not in {"accepted", "clean_partial"}
        or not isinstance(inspection, dict)
        or inspection.get("clean_for_continuation") is not True
    ):
        raise ValueError("restart parent has not passed continuation inspection")
    revalidate_inspection_files(inspection, parent)
    if accounting.get("result") == "clean_partial":
        inspection = revalidate_continuation_plasma_evidence(inspection, parent)
        plasma_evidence = inspection.get("plasma_continuation_evidence")
        if (
            isinstance(plasma_evidence, dict)
            and plasma_evidence.get("eligibility_only") is True
        ):
            parent_key = (
                str(parent.get("run", {}).get("case_id", "")),
                str(parent.get("job_id", "")),
            )
            if parent_key == ("R03", "4762472"):
                raise ValueError(
                    "R03 frozen-E03 clean partial is inventory-compatible only; "
                    "future continuation must use accepted job "
                    f"{R03_ACCEPTED_CONTINUATION_JOB_ID}"
                )
            if parent_key == ("R12", R12_HISTORICAL_CLEAN_PARTIAL_JOB_ID):
                raise ValueError(
                    "R12 frozen-E03 clean partial is inventory-compatible only; "
                    "R12 must restart fresh as "
                    f"{R12_FRESH_RERUN_SEGMENT} from t=0"
                )
            raise ValueError(
                "frozen-E03 clean partial is inventory-compatible only"
            )
    terminal = inspection.get("terminal_restart")
    if (
        not isinstance(terminal, dict)
        or Path(str(terminal.get("path", ""))).resolve() != restart
        or terminal.get("sha256") != sha256(restart)
    ):
        raise ValueError("continuation must use the inspected terminal restart")
    revalidate_retained_product(terminal)
    restart_files = retained_product_paths(terminal)
    restart_time = (
        restart_product_time(
            restart_files,
            allow_missing_marker=allow_missing_restart_time_marker,
        )
        if allow_missing_restart_time_marker or parent_root != DEFAULT_ROOT.resolve()
        else float(
            authenticated_restart_product_time(
                restart_files, parent.get("command")
            )["binary_time"]
        )
    )
    try:
        final_time = float(inspection["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("continuation parent inspection lacks final time") from error
    if (
        restart_time is not None
        and abs(restart_time - final_time) > 1.0e-12
    ):
        raise ValueError(
            "continuation terminal restart physical time differs from inspection"
        )
    result = {
        "execution_epoch": EXECUTION_EPOCH,
        "manifest": str(parent_manifest_path),
        "case_id": parent["run"]["case_id"],
        "segment": parent["run"]["segment"],
        "result": accounting["result"],
        "restart_sha256": terminal["sha256"],
        "restart_files": [str(path) for path in restart_files],
        "final_time": final_time,
        "restart_time": restart_time,
        "input_sha256": parent["command"]["input_sha256"],
        "executable_sha256": parent["command"]["executable_sha256"],
    }
    if accounting.get("result") == "clean_partial":
        result["plasma_continuation_policy"] = inspection[
            "plasma_continuation_policy"
        ]
        result["plasma_continuation_evidence"] = inspection[
            "plasma_continuation_evidence"
        ]
    return result


def reservation_for_manifest(reservations: list[dict[str, object]],
                             manifest_path: Path) -> dict[str, object]:
    """Return exactly one reservation corresponding to a run manifest."""

    matches = [
        item for item in reservations
        if Path(str(item.get("manifest", ""))).resolve() == manifest_path.resolve()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one reservation for {manifest_path}, found {len(matches)}"
        )
    if matches[0].get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("reservation does not belong to the current execution epoch")
    return matches[0]


def scheduler_environment() -> dict[str, str]:
    """Return scheduler state without caller routing, interpreter, or loader hooks."""

    try:
        account = pwd.getpwuid(os.geteuid())
    except KeyError as error:
        raise ValueError(
            "effective UID has no account record; cannot construct scheduler environment"
        ) from error
    return {
        "HOME": account.pw_dir,
        "LC_ALL": "C",
        "LOGNAME": account.pw_name,
        "PATH": "/usr/bin:/bin",
        "USER": account.pw_name,
    }


def scheduler_account_name() -> str:
    """Return the account name bound to this process effective UID."""

    try:
        return pwd.getpwuid(os.geteuid()).pw_name
    except KeyError as error:
        raise ValueError(
            "effective UID has no account name; cannot query scheduler"
        ) from error


def production_queue_output(args: argparse.Namespace,
                            offline_local_root: bool) -> str:
    """Query all of this user's queued jobs before a root-writing submission."""

    fixture = getattr(args, "squeue_file", None)
    if fixture:
        if not offline_local_root:
            raise ValueError("--squeue-file is restricted to offline local roots")
        return Path(fixture).read_text(encoding="utf-8")
    if offline_local_root:
        raise ValueError("offline local-root submission requires --squeue-file")
    user = scheduler_account_name()
    try:
        return subprocess.run(
            [str(SQUEUE), "-h", "-u", user, "-o", "%i|%P|%T|%j"],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        raise ValueError(
            "squeue is unavailable; refusing Stage I submission"
        ) from error


def cancelled_job_queue_output(args: argparse.Namespace,
                               offline_local_root: bool) -> str:
    """Query one cancelled job ID and require it to be absent from squeue."""

    fixture = getattr(args, "squeue_file", None)
    if fixture:
        if not offline_local_root:
            raise ValueError("--squeue-file is restricted to offline local roots")
        return Path(fixture).read_text(encoding="utf-8")
    if offline_local_root:
        raise ValueError(
            "offline local-root submitted cancellation requires --squeue-file"
        )
    return live_cancelled_job_queue_output(args.job_id)


def live_cancelled_job_queue_output(job_id: str) -> str:
    """Query one exact job ID from trusted live squeue."""

    require_numeric_job_id(job_id)
    try:
        return subprocess.run(
            [
                str(SQUEUE), "-h", "-j", job_id,
                "-o", "%i|%P|%T|%j",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        raise ValueError(
            "squeue is unavailable; refusing submitted Stage I cancellation"
        ) from error


def expected_submitted_queue_jobs(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
) -> dict[str, str]:
    """Map retained submitted reservations to their exact scheduler names."""

    expected = {}
    for reservation in active_reservations(reservations):
        if reservation.get("state") != "submitted":
            continue
        job_id = require_numeric_job_id(str(reservation.get("job_id", "")))
        if job_id in expected:
            raise ValueError(f"submitted Stage I job ID is duplicated: {job_id}")
        manifest_path = require_safe_manifest_path(paths, reservation["manifest"])
        if not manifest_path.is_file():
            raise ValueError(
                f"submitted Stage I reservation lacks retained manifest: {manifest_path}"
            )
        manifest = read_manifest(manifest_path)
        require_current_epoch(manifest, "submitted queue reservation")
        if (
            manifest.get("state") != "submitted"
            or str(manifest.get("job_id", "")) != job_id
        ):
            raise ValueError(
                f"submitted Stage I reservation differs from manifest: {manifest_path}"
            )
        require_reservation_matches_manifest(reservation, manifest)
        expected[job_id] = expected_job_name(manifest)
    return expected


def authenticate_production_queue(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    lines: list[str],
) -> None:
    """Authenticate every CGL row while retaining unrelated user jobs."""

    expected = expected_submitted_queue_jobs(paths, reservations)
    seen = set()
    conflicts = []
    for line in lines:
        row = next(csv.reader([line], delimiter="|"))
        if len(row) != 4:
            conflicts.append(line)
            continue
        if any(not value or value.strip() != value for value in row):
            conflicts.append(line)
            continue
        job_id, partition, state, job_name = row
        if (
            not job_name.startswith(CGL_JOB_NAME_PREFIX)
            and job_id not in expected
        ):
            continue
        if (
            JOB_ID_PATTERN.fullmatch(job_id) is None
            or partition != PARTITION
            or state not in NONTERMINAL_STATES
            or expected.get(job_id) != job_name
            or job_id in seen
        ):
            conflicts.append(line)
            continue
        seen.add(job_id)
    if conflicts:
        raise ValueError(
            "a queued CGL job or malformed scheduler row lacks an authenticated submitted "
            f"{EXECUTION_EPOCH} Stage I reservation: " + "; ".join(conflicts)
        )


def validate_queue_authentication_evidence(value: object, label: str) -> None:
    """Require retained queue rows and their deterministic digest."""

    if not isinstance(value, dict) or frozenset(value) != {
        "checked_utc", "rows", "rows_sha256",
    }:
        raise ValueError(f"{label} has invalid queue authentication columns")
    parse_utc_timestamp(value["checked_utc"], f"{label} checked_utc")
    rows = value["rows"]
    if (
        not isinstance(rows, list)
        or not all(isinstance(row, str) and row.strip() == row for row in rows)
        or stable_json_sha256(rows) != value["rows_sha256"]
    ):
        raise ValueError(f"{label} has invalid queue authentication evidence")


def authenticated_production_queue_evidence(
    args: argparse.Namespace,
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    offline_local_root: bool,
) -> dict[str, object]:
    """Query, authenticate, and retain one complete scheduler queue snapshot."""

    rows = [
        line for line in production_queue_output(args, offline_local_root).splitlines()
        if line.strip()
    ]
    authenticate_production_queue(paths, reservations, rows)
    return {
        "checked_utc": utc_now(),
        "rows": rows,
        "rows_sha256": stable_json_sha256(rows),
    }


def validate_submission_fixture_options(args: argparse.Namespace,
                                        offline_local_root: bool) -> None:
    """Keep scheduler fixture injection and test bypass out of production."""

    if getattr(args, "squeue_file", None) and not offline_local_root:
        raise ValueError("--squeue-file is restricted to offline local roots")
    if getattr(args, "skip_slurm_test", False) and not offline_local_root:
        raise ValueError("--skip-slurm-test is restricted to offline local roots")


def scheduler_test_only_output(script: dict[str, object] | Path) -> str:
    """Run the live Slurm submission validator with trusted routing."""

    if isinstance(script, dict):
        reauthenticate_open_batch_script(script)
        command_path = str(script["descriptor_path"])
        pass_fds = (int(script["fd"]),)
    else:
        command_path = str(script)
        pass_fds = ()
    return subprocess.run(
        [str(SBATCH), "--test-only", command_path],
        check=True,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
        pass_fds=pass_fds,
    ).stdout


def require_clearance_file_binding(value: object, path: Path, label: str, *,
                                   mode: int) -> str:
    """Authenticate one exact path, digest, and inode bound by a clearance."""

    binding = require_r17_exact_keys(
        value, {"path", "sha256", "mode", "links", "device", "inode"}, label
    )
    expected_digest = require_r17_sha256(binding["sha256"], f"{label} SHA-256")
    if (
        binding["path"] != str(path)
        or binding["mode"] != f"{mode:04o}"
        or binding["links"] != 1
    ):
        raise ValueError(f"{label} differs from its exact retained path profile")
    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        owner_controlled=True,
        symlink_free=True,
    )
    profile = path.lstat()
    if (
        binding["device"] != profile.st_dev
        or binding["inode"] != profile.st_ino
        or hashlib.sha256(retained).hexdigest() != expected_digest
    ):
        raise ValueError(f"{label} differs from its exact retained inode or bytes")
    return expected_digest


def require_clearance_directory_binding(value: object, path: Path,
                                        label: str) -> None:
    """Authenticate one exact non-writable shared-root namespace inode."""

    binding = require_r17_exact_keys(
        value, {"path", "mode", "device", "inode"}, label
    )
    require_owner_symlink_free_path(path, label)
    profile = path.lstat()
    if (
        not stat.S_ISDIR(profile.st_mode)
        or stat.S_IMODE(profile.st_mode) & 0o022
        or binding != {
            "path": str(path),
            "mode": f"{stat.S_IMODE(profile.st_mode):04o}",
            "device": profile.st_dev,
            "inode": profile.st_ino,
        }
    ):
        raise ValueError(f"{label} differs from its exact retained directory inode")


def require_clearance_publication_bindings(value: object, expected: list[Path],
                                           label: str
                                           ) -> list[dict[str, object]]:
    """Authenticate an ordered list of exact immutable publications."""

    if not isinstance(value, list) or len(value) != len(expected):
        raise ValueError(f"{label} must bind every exact retained publication")
    retained = []
    for index, path in enumerate(expected):
        _, digest = read_controller_publication_json(
            path, f"{label} publication {index}", mode=0o444
        )
        require_r17_declared_publication(
            value[index], path, digest, f"{label} publication {index}"
        )
        retained.append({
            "path": str(path),
            "sha256": digest,
            "mode": "0444",
            "links": 1,
        })
    return retained


def require_clearance_review(
    path: Path,
    *,
    checkpoint: str,
    candidate_path: Path,
    candidate_sha: str,
    record_type: str,
    decision: str,
    label: str,
) -> tuple[str, str, datetime]:
    """Authenticate one exact independent clearance or publication review."""

    review, review_sha = read_controller_publication_json(path, label, mode=0o444)
    review = require_r17_exact_keys(
        review,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "reviewed_utc", "decision", "reviewer",
            "independent_of_implementation", "candidate",
        },
        label,
    )
    reviewer = require_r17_exact_keys(
        review["reviewer"], {"role", "reviewer_id"}, f"{label} reviewer"
    )
    reviewer_id = require_r17_nonempty_string(
        reviewer["reviewer_id"], f"{label} reviewer ID"
    )
    require_r17_nonempty_string(reviewer["role"], f"{label} reviewer role")
    if (
        review["schema_version"] != 1
        or review["record_type"] != record_type
        or review["checkpoint"] != checkpoint
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != decision
        or review["independent_of_implementation"] is not True
    ):
        raise ValueError(f"{label} decision or authority differs")
    require_r17_declared_publication(
        review["candidate"], candidate_path, candidate_sha, f"{label} candidate"
    )
    return reviewer_id, review_sha, require_r17_utc(
        review["reviewed_utc"], f"{label} review time"
    )


def managed_shared_root_clearance_chain_paths(
    accounting: Path, checkpoint: int
) -> tuple[Path, Path, Path, Path]:
    """Return the four canonical paths for one managed clearance checkpoint."""

    artifact = accounting / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F{checkpoint}_"
        "shared_root_isolation_clearance_refresh.json"
    )
    review = Path(f"{artifact}.independent_review.json")
    audit = Path(f"{artifact}.publication_audit.json")
    audit_review = Path(f"{audit}.independent_review.json")
    return artifact, review, audit, audit_review


def managed_shared_root_clearance_scope(paths: dict[str, Path]) -> dict[str, object]:
    """Return the fixed non-broadening authority of a managed refresh."""

    return {
        "authorized_actions": ["check-submit", "submit"],
        "authorized_case_ids": [f"R{number:02d}" for number in range(3, 18)],
        "authorized_execution_epoch": EXECUTION_EPOCH,
        "authorized_stage_i_namespace": str(
            paths["root"] / "runs/mks24-stage-i" / EXECUTION_EPOCH
        ),
        "all_controller_preflights_required": True,
        "all_user_queue_authentication_required": True,
        "root_locked_submit_required": True,
        "exact_controller_authenticated_manifest_required": True,
        "prepare_authorized": False,
        "direct_sbatch_authorized": False,
        "other_campaign_acknowledgement_authorized": False,
        "outside_namespace_authorized": False,
        "r17_last_policy_unchanged": True,
        "binding_refresh_only": True,
    }


def managed_shared_root_clearance_publication_requirements() -> dict[str, bool]:
    """Return the fixed independent-review requirements of a refresh."""

    return {
        "independent_review_required": True,
        "publication_audit_required": True,
        "independent_publication_review_required": True,
        "distinct_reviewers_required": True,
    }


def current_clearance_file_binding(path: Path, label: str, *,
                                   mode: int) -> dict[str, object]:
    """Return one exact current path, digest, and inode binding."""

    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        owner_controlled=True,
        symlink_free=True,
    )
    profile = path.lstat()
    return {
        "path": str(path),
        "sha256": hashlib.sha256(retained).hexdigest(),
        "mode": f"{mode:04o}",
        "links": 1,
        "device": profile.st_dev,
        "inode": profile.st_ino,
    }


def current_clearance_publication_binding(path: Path,
                                          label: str) -> dict[str, object]:
    """Return one exact current immutable-publication binding."""

    _, digest = read_controller_publication_json(path, label, mode=0o444)
    return {"path": str(path), "sha256": digest, "mode": "0444", "links": 1}


def current_clearance_directory_binding(path: Path,
                                        label: str) -> dict[str, object]:
    """Return one exact current non-writable directory-inode binding."""

    require_owner_symlink_free_path(path, label)
    profile = path.lstat()
    if not stat.S_ISDIR(profile.st_mode) or stat.S_IMODE(profile.st_mode) & 0o022:
        raise ValueError(f"{label} is not an owner-controlled non-writable directory")
    return {
        "path": str(path),
        "mode": f"{stat.S_IMODE(profile.st_mode):04o}",
        "device": profile.st_dev,
        "inode": profile.st_ino,
    }


def retained_stale_terminal_sacct_binding(
    paths: dict[str, Path],
) -> dict[str, object]:
    """Authenticate the immutable terminal sacct evidence for the stale job."""

    source = paths["root"] / SHARED_ROOT_CLEARANCE_SUPERSESSION_RELATIVE
    publication, digest = read_controller_publication_json(
        source, "managed shared-root clearance retained scheduler evidence",
        mode=0o444,
    )
    stale = publication.get("stale_campaign")
    if not isinstance(stale, dict):
        raise ValueError(
            "managed shared-root clearance retained scheduler evidence is absent"
        )
    scheduler = require_r17_exact_keys(
        stale.get("scheduler"),
        {"observed_utc", "sacct", "squeue", "terminal_completed_0_0_absent"},
        "managed shared-root clearance retained scheduler evidence",
    )
    sacct = require_r17_exact_keys(
        scheduler["sacct"], {"argv", "returncode", "stderr", "stdout"},
        "managed shared-root clearance retained sacct evidence",
    )
    expected_argv = [
        str(SACCT), "-X", "-n", "-P", "-j", SHARED_ROOT_STALE_JOB_ID,
        "-o", "JobIDRaw,JobName,State,Elapsed,NNodes,ExitCode,Submit,Start,End",
    ]
    if (
        stale.get("campaign_id") != SHARED_ROOT_STALE_CAMPAIGN_ID
        or stale.get("job_id") != SHARED_ROOT_STALE_JOB_ID
        or stale.get("retained_state") != "running"
        or scheduler["terminal_completed_0_0_absent"] is not True
        or sacct["argv"] != expected_argv
        or sacct["returncode"] != 0
        or sacct["stderr"] != ""
        or not isinstance(sacct["stdout"], str)
    ):
        raise ValueError("managed shared-root clearance retained sacct evidence differs")
    rows = [
        row for row in csv.reader(sacct["stdout"].splitlines(), delimiter="|")
        if row
    ]
    if len(rows) != 1 or len(rows[0]) != 9:
        raise ValueError(
            "managed shared-root clearance retained sacct evidence must contain "
            "one exact top-level row"
        )
    row = rows[0]
    observed = require_r17_utc(
        scheduler["observed_utc"],
        "managed shared-root clearance retained scheduler observation",
    )
    submitted = require_r17_scheduler_utc(
        row[6], "managed shared-root clearance retained sacct submit time"
    )
    started = require_r17_scheduler_utc(
        row[7], "managed shared-root clearance retained sacct start time"
    )
    completed = require_r17_scheduler_utc(
        row[8], "managed shared-root clearance retained sacct completion time"
    )
    try:
        elapsed = parse_walltime(row[3])
        nodes = int(row[4])
    except ValueError as error:
        raise ValueError(
            "managed shared-root clearance retained sacct allocation differs"
        ) from error
    if (
        row[0] != SHARED_ROOT_STALE_JOB_ID
        or row[1] != SHARED_ROOT_STALE_JOB_NAME
        or row[2] != "COMPLETED"
        or row[5] != "0:0"
        or elapsed <= 0
        or nodes <= 0
        or not submitted <= started <= completed <= observed
    ):
        raise ValueError("managed shared-root clearance retained sacct evidence differs")
    return {
        "source_publication": {
            "path": str(source),
            "sha256": digest,
            "mode": "0444",
            "links": 1,
        },
        "json_path": "stale_campaign.scheduler.sacct",
        "sacct_sha256": stable_json_sha256(sacct),
        "job_id": SHARED_ROOT_STALE_JOB_ID,
        "job_name": SHARED_ROOT_STALE_JOB_NAME,
        "state": "COMPLETED",
        "exit_code": "0:0",
        "completed_utc": row[8],
    }


def require_stale_job_absent_from_complete_queue(
    value: object, *, now: datetime,
) -> dict[str, object]:
    """Require the stale job ID absent from a fresh complete-user squeue snapshot."""

    validate_queue_authentication_evidence(value, "managed clearance live squeue")
    if not isinstance(value, dict):
        raise ValueError("managed clearance live squeue evidence is invalid")
    checked = parse_utc_timestamp(
        value["checked_utc"], "managed clearance live squeue checked_utc"
    )
    if (
        checked > now + timedelta(minutes=5)
        or now - checked > timedelta(minutes=5)
    ):
        raise ValueError("managed clearance live squeue evidence is not fresh")
    rows = value["rows"]
    if not isinstance(rows, list):
        raise ValueError("managed clearance live squeue rows are invalid")
    for line in rows:
        row = next(csv.reader([line], delimiter="|"))
        if len(row) != 4:
            raise ValueError("managed clearance live squeue row is malformed")
        if row[0] == SHARED_ROOT_STALE_JOB_ID:
            raise ValueError(
                "managed clearance stale job remains present in complete squeue evidence"
            )
    return {
        "job_id": SHARED_ROOT_STALE_JOB_ID,
        "absent": True,
        "checked_utc": value["checked_utc"],
        "rows_sha256": value["rows_sha256"],
    }


def authenticate_managed_shared_root_clearance_chain_commit(
    accounting: Path, checkpoint: int,
) -> None:
    """Authenticate one complete audit-review-last publication commit."""

    artifact_path, review_path, audit_path, audit_review_path = (
        managed_shared_root_clearance_chain_paths(accounting, checkpoint)
    )
    artifact, artifact_sha = read_controller_publication_json(
        artifact_path, f"managed shared-root clearance F-{checkpoint} artifact",
        mode=0o444,
    )
    review, review_sha = read_controller_publication_json(
        review_path, f"managed shared-root clearance F-{checkpoint} review",
        mode=0o444,
    )
    audit, audit_sha = read_controller_publication_json(
        audit_path, f"managed shared-root clearance F-{checkpoint} audit",
        mode=0o444,
    )
    audit_review, _ = read_controller_publication_json(
        audit_review_path,
        f"managed shared-root clearance F-{checkpoint} audit review",
        mode=0o444,
    )
    if not all(isinstance(value, dict) for value in (
        artifact, review, audit, audit_review,
    )):
        raise ValueError("managed shared-root clearance committed chain differs")
    require_r17_declared_publication(
        review.get("candidate"), artifact_path, artifact_sha,
        f"managed shared-root clearance F-{checkpoint} review candidate",
    )
    require_r17_declared_publication(
        audit.get("artifact"), artifact_path, artifact_sha,
        f"managed shared-root clearance F-{checkpoint} audit artifact",
    )
    require_r17_declared_publication(
        audit.get("independent_review"), review_path, review_sha,
        f"managed shared-root clearance F-{checkpoint} audit review",
    )
    require_r17_declared_publication(
        audit_review.get("candidate"), audit_path, audit_sha,
        f"managed shared-root clearance F-{checkpoint} audit-review candidate",
    )


def managed_shared_root_clearance_checkpoints(accounting: Path) -> list[int]:
    """Return audit-review-last complete managed-clearance checkpoints."""

    entries = set(trusted_directory_entries(
        accounting, "managed shared-root clearance accounting directory"
    ))
    checkpoints = []
    for entry in sorted(entries):
        match = SHARED_ROOT_CLEARANCE_REFRESH_AUDIT_REVIEW_PATTERN.fullmatch(entry)
        if match is not None:
            checkpoint = int(match.group("checkpoint"))
            chain = managed_shared_root_clearance_chain_paths(accounting, checkpoint)
            if not all(path.name in entries for path in chain):
                continue
            authenticate_managed_shared_root_clearance_chain_commit(
                accounting, checkpoint
            )
            checkpoints.append(checkpoint)
    return sorted(set(checkpoints))


def build_managed_shared_root_clearance_refresh(
    paths: dict[str, Path], checkpoint: int, *,
    now: datetime | None = None,
) -> dict[str, object]:
    """Build a read-only exact-binding candidate for independent review."""

    checkpoint = require_r17_integer(
        checkpoint, "managed shared-root clearance checkpoint", minimum=1
    )
    accounting = paths["root"] / "accounting"
    existing = managed_shared_root_clearance_checkpoints(accounting)
    if existing and checkpoint <= max(existing):
        raise ValueError(
            "managed shared-root clearance checkpoint must exceed every publication"
        )
    chain_paths = managed_shared_root_clearance_chain_paths(accounting, checkpoint)
    if any(path.exists() for path in chain_paths):
        raise ValueError("managed shared-root clearance checkpoint already exists")

    source_paths = [
        paths["root"] / relative
        for relative in SHARED_ROOT_CLEARANCE_SOURCE_AUTHORITY_RELATIVES
    ]
    prior_paths = [
        paths["root"] / relative for relative in SHARED_ROOT_CLEARANCE_LEGACY_RELATIVES
    ]
    for prior_checkpoint in existing:
        prior_paths.extend(
            managed_shared_root_clearance_chain_paths(accounting, prior_checkpoint)
        )
    stale_namespace = paths["root"] / "runs" / SHARED_ROOT_STALE_CAMPAIGN_ID
    stale_manifest = stale_namespace / "manifest/prepared_run.json"
    stale = read_manifest(stale_manifest)
    if (
        stale.get("campaign_id") != SHARED_ROOT_STALE_CAMPAIGN_ID
        or stale.get("slurm_job_id") != SHARED_ROOT_STALE_JOB_ID
        or stale.get("state") != "running"
    ):
        raise ValueError("managed shared-root clearance stale campaign differs")
    stage_namespace = paths["root"] / "runs/mks24-stage-i" / EXECUTION_EPOCH
    stage_binding = current_clearance_directory_binding(
        stage_namespace, "managed shared-root clearance Stage I namespace"
    )
    stale_binding = current_clearance_directory_binding(
        stale_namespace, "managed shared-root clearance stale namespace"
    )
    if (
        stage_namespace == stale_namespace
        or stage_namespace in stale_namespace.parents
        or stale_namespace in stage_namespace.parents
    ):
        raise ValueError("managed shared-root clearance namespaces are not disjoint")
    generated = datetime.now(timezone.utc) if now is None else now.astimezone(timezone.utc)
    return {
        "schema_version": 1,
        "record_type": "stage-i-managed-shared-root-isolation-clearance-refresh",
        "checkpoint": f"F-{checkpoint}",
        "execution_epoch": EXECUTION_EPOCH,
        "generated_utc": generated.isoformat(),
        "requested_acknowledgement": SHARED_ROOT_STALE_CAMPAIGN_ID,
        "scope": managed_shared_root_clearance_scope(paths),
        "controller_binding": current_clearance_file_binding(
            retained_controller_source_path(),
            "managed shared-root clearance controller",
            mode=0o644,
        ),
        "source_authority_publications": [
            current_clearance_publication_binding(
                path, f"managed shared-root clearance source authority {index}"
            )
            for index, path in enumerate(source_paths)
        ],
        "supersedes": [
            current_clearance_publication_binding(
                path, f"managed shared-root clearance superseded publication {index}"
            )
            for index, path in enumerate(prior_paths)
        ],
        "stale_campaign": {
            "campaign_id": SHARED_ROOT_STALE_CAMPAIGN_ID,
            "job_id": SHARED_ROOT_STALE_JOB_ID,
            "retained_state": "running",
            "manifest": current_clearance_file_binding(
                stale_manifest,
                "managed shared-root clearance stale manifest",
                mode=0o644,
            ),
            "terminal_sacct": retained_stale_terminal_sacct_binding(paths),
        },
        "isolation": {
            "stage_i_namespace": stage_binding,
            "stale_namespace": stale_binding,
            "lexically_disjoint": True,
        },
        "publication_requirements": (
            managed_shared_root_clearance_publication_requirements()
        ),
    }


def render_managed_shared_root_clearance_refresh(args: argparse.Namespace) -> int:
    """Print one canonical-JSON managed clearance candidate without mutation."""

    root = require_root(Path(args.root), args.allow_local_root)
    candidate = build_managed_shared_root_clearance_refresh(
        layout(root), args.checkpoint
    )
    print(json.dumps(candidate, indent=2, sort_keys=True))
    return 0


def require_managed_shared_root_clearance(
    paths: dict[str, Path], requested_campaigns: set[str], *,
    action: str,
    case_id: str,
    queue_evidence: dict[str, object],
    now: datetime | None = None,
) -> dict[str, object] | None:
    """Require the latest fully reviewed clearance for a production overlap."""

    if not requested_campaigns:
        return None
    if requested_campaigns != {SHARED_ROOT_STALE_CAMPAIGN_ID}:
        raise ValueError(
            "managed shared-root clearance authorizes exactly the retained campaign"
        )
    if action not in {"check-submit", "submit"}:
        raise ValueError("managed shared-root clearance action is not authorized")
    authorized_cases = {f"R{number:02d}" for number in range(3, 18)}
    if case_id not in authorized_cases:
        raise ValueError("managed shared-root clearance case is not authorized")
    current = datetime.now(timezone.utc) if now is None else now.astimezone(timezone.utc)
    accounting = paths["root"] / "accounting"
    checkpoints = managed_shared_root_clearance_checkpoints(accounting)
    if not checkpoints:
        raise ValueError("managed shared-root clearance is not canonically published")
    checkpoint_number = max(checkpoints)
    checkpoint = f"F-{checkpoint_number}"
    artifact_path, review_path, audit_path, audit_review_path = (
        managed_shared_root_clearance_chain_paths(accounting, checkpoint_number)
    )
    artifact, artifact_sha = read_controller_publication_json(
        artifact_path, "managed shared-root clearance", mode=0o444
    )
    artifact = require_r17_exact_keys(
        artifact,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "generated_utc", "requested_acknowledgement", "scope",
            "controller_binding", "source_authority_publications", "supersedes",
            "stale_campaign", "isolation", "publication_requirements",
        },
        "managed shared-root clearance",
    )
    scope = managed_shared_root_clearance_scope(paths)
    publication_requirements = managed_shared_root_clearance_publication_requirements()
    generated = require_r17_utc(
        artifact["generated_utc"], "managed shared-root clearance generation time"
    )
    if (
        artifact["schema_version"] != 1
        or artifact["record_type"]
        != "stage-i-managed-shared-root-isolation-clearance-refresh"
        or artifact["checkpoint"] != checkpoint
        or artifact["execution_epoch"] != EXECUTION_EPOCH
        or artifact["requested_acknowledgement"] != SHARED_ROOT_STALE_CAMPAIGN_ID
        or artifact["scope"] != scope
        or artifact["publication_requirements"] != publication_requirements
    ):
        raise ValueError("managed shared-root clearance scope or identity differs")

    controller_path = retained_controller_source_path()
    require_clearance_file_binding(
        artifact["controller_binding"],
        controller_path,
        "managed shared-root clearance controller",
        mode=0o644,
    )
    source_authority = require_clearance_publication_bindings(
        artifact["source_authority_publications"],
        [
            paths["root"] / relative
            for relative in SHARED_ROOT_CLEARANCE_SOURCE_AUTHORITY_RELATIVES
        ],
        "managed shared-root clearance source authority",
    )
    prior_paths = [paths["root"] / relative for relative in SHARED_ROOT_CLEARANCE_LEGACY_RELATIVES]
    for prior_checkpoint in sorted(set(checkpoints)):
        if prior_checkpoint >= checkpoint_number:
            continue
        prior_paths.extend(
            managed_shared_root_clearance_chain_paths(accounting, prior_checkpoint)
        )
    supersedes = require_clearance_publication_bindings(
        artifact["supersedes"], prior_paths, "managed shared-root clearance supersedes"
    )

    stale = require_r17_exact_keys(
        artifact["stale_campaign"],
        {"campaign_id", "job_id", "retained_state", "manifest", "terminal_sacct"},
        "managed shared-root clearance stale campaign",
    )
    stale_manifest = (
        paths["root"] / "runs" / SHARED_ROOT_STALE_CAMPAIGN_ID
        / "manifest/prepared_run.json"
    )
    require_clearance_file_binding(
        stale["manifest"],
        stale_manifest,
        "managed shared-root clearance stale manifest",
        mode=0o644,
    )
    retained_stale = read_manifest(stale_manifest)
    if (
        stale["campaign_id"] != SHARED_ROOT_STALE_CAMPAIGN_ID
        or stale["job_id"] != SHARED_ROOT_STALE_JOB_ID
        or stale["retained_state"] != "running"
        or retained_stale.get("campaign_id") != SHARED_ROOT_STALE_CAMPAIGN_ID
        or retained_stale.get("slurm_job_id") != SHARED_ROOT_STALE_JOB_ID
        or retained_stale.get("state") != "running"
    ):
        raise ValueError("managed shared-root clearance stale campaign differs")
    terminal_sacct = retained_stale_terminal_sacct_binding(paths)
    if stale["terminal_sacct"] != terminal_sacct:
        raise ValueError("managed shared-root clearance retained sacct binding differs")

    isolation = require_r17_exact_keys(
        artifact["isolation"],
        {"stage_i_namespace", "stale_namespace", "lexically_disjoint"},
        "managed shared-root clearance isolation",
    )
    stage_namespace = paths["root"] / "runs/mks24-stage-i" / EXECUTION_EPOCH
    stale_namespace = paths["root"] / "runs" / SHARED_ROOT_STALE_CAMPAIGN_ID
    require_clearance_directory_binding(
        isolation["stage_i_namespace"], stage_namespace,
        "managed shared-root clearance Stage I namespace",
    )
    require_clearance_directory_binding(
        isolation["stale_namespace"], stale_namespace,
        "managed shared-root clearance stale namespace",
    )
    if (
        isolation["lexically_disjoint"] is not True
        or stage_namespace == stale_namespace
        or stage_namespace in stale_namespace.parents
        or stale_namespace in stage_namespace.parents
    ):
        raise ValueError("managed shared-root clearance namespaces are not disjoint")

    reviewer_id, review_sha, reviewed = require_clearance_review(
        review_path,
        checkpoint=checkpoint,
        candidate_path=artifact_path,
        candidate_sha=artifact_sha,
        record_type="stage-i-managed-shared-root-isolation-clearance-review",
        decision="approved",
        label="managed shared-root clearance independent review",
    )
    audit, audit_sha = read_controller_publication_json(
        audit_path, "managed shared-root clearance publication audit", mode=0o444
    )
    audit = require_r17_exact_keys(
        audit,
        {
            "schema_version", "record_type", "checkpoint", "execution_epoch",
            "published_utc", "artifact", "independent_review", "authority",
        },
        "managed shared-root clearance publication audit",
    )
    authority = {
        "binding_refresh_only": True,
        "scope_expanded": False,
        "bypasses_controller_preflights": False,
    }
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-managed-shared-root-isolation-clearance-publication-audit"
        or audit["checkpoint"] != checkpoint
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"] != authority
    ):
        raise ValueError("managed shared-root clearance publication audit differs")
    require_r17_declared_publication(
        audit["artifact"], artifact_path, artifact_sha,
        "managed shared-root clearance publication audit artifact",
    )
    require_r17_declared_publication(
        audit["independent_review"], review_path, review_sha,
        "managed shared-root clearance publication audit review",
    )
    published = require_r17_utc(
        audit["published_utc"], "managed shared-root clearance publication time"
    )
    audit_reviewer_id, audit_review_sha, audit_reviewed = require_clearance_review(
        audit_review_path,
        checkpoint=checkpoint,
        candidate_path=audit_path,
        candidate_sha=audit_sha,
        record_type=(
            "stage-i-managed-shared-root-isolation-clearance-"
            "publication-audit-review"
        ),
        decision="approved-for-publication",
        label="managed shared-root clearance publication-audit review",
    )
    if (
        reviewer_id == audit_reviewer_id
        or not generated < reviewed < published < audit_reviewed
        or audit_reviewed > current + timedelta(minutes=5)
    ):
        raise ValueError("managed shared-root clearance review chain differs")
    live_queue_absence = require_stale_job_absent_from_complete_queue(
        queue_evidence, now=current
    )
    return {
        "checkpoint": checkpoint,
        "artifact": {
            "path": str(artifact_path), "sha256": artifact_sha,
            "mode": "0444", "links": 1,
        },
        "independent_review": {
            "path": str(review_path), "sha256": review_sha,
            "mode": "0444", "links": 1,
        },
        "publication_audit": {
            "path": str(audit_path), "sha256": audit_sha,
            "mode": "0444", "links": 1,
        },
        "publication_audit_review": {
            "path": str(audit_review_path), "sha256": audit_review_sha,
            "mode": "0444", "links": 1,
        },
        "source_authority_publications": source_authority,
        "supersedes": supersedes,
        "retained_terminal_sacct": terminal_sacct,
        "live_queue_absence": live_queue_absence,
        "consumption": {
            "action": action,
            "case_id": case_id,
            "authorized_stale_manifest": str(stale_manifest),
        },
    }


def shared_root_campaign_conflicts(root: Path,
                                   allowed_manifests: set[Path]) -> list[str]:
    """Return top-level CGL campaigns requiring an explicit overlap review."""

    conflicts = []
    for manifest_path in sorted(
        (root / "runs").glob("*/manifest/prepared_run.json")
    ):
        manifest = read_manifest(manifest_path)
        state = str(manifest.get("state", "")).lower()
        if state not in SHARED_ROOT_ACTIVE_STATES:
            continue
        campaign_id = str(
            manifest.get("campaign_id", manifest_path.parents[1].name)
        )
        if manifest_path not in allowed_manifests:
            conflicts.append(f"{campaign_id}|{state}|{manifest_path}")
    return conflicts


def require_existing_layout(paths: dict[str, Path]) -> None:
    """Require initialized stores without modifying them."""

    for key in ("ledger", "reservations"):
        if not paths[key].is_file():
            raise ValueError(f"Stage I store is missing: {paths[key]}")
    for key in ("runs", "transactions"):
        if not paths[key].is_dir():
            raise ValueError(f"Stage I directory is missing: {paths[key]}")
    read_ledger(paths)
    read_reservations(paths)


def require_reconciled_store_consistency(
    paths: dict[str, Path],
    *,
    allow_absent_qualification: bool = False,
    allow_invalid_qualification: bool = False,
    ignored_issue_prefixes: tuple[str, ...] = (),
) -> None:
    """Fail closed before canonical metadata or bundle mutations."""

    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    report = reconcile_report(paths["root"])
    issues = []
    for issue in report["issues"]:
        if (
            (
                allow_absent_qualification
                and issue.startswith(
                    "E03 corrected-build Frontier qualification is pending:"
                )
            )
            or (
                allow_invalid_qualification
                and issue.startswith(
                    "E03 corrected-build Frontier qualification is invalid:"
                )
            )
        ):
            continue
        if any(issue.startswith(prefix) for prefix in ignored_issue_prefixes):
            continue
        issues.append(issue)
    if issues:
        raise ValueError(
            "canonical E03 store reconciliation failed before mutation: "
            + "; ".join(issues)
        )


def reauthenticate_submission_authority(
    paths: dict[str, Path],
    manifest_path: Path,
    manifest: dict[str, object],
    reservations: list[dict[str, object]],
    reservation: dict[str, object],
    offline_local_root: bool,
) -> None:
    """Reauthenticate F118 and exact R17 launch evidence before scheduler access."""

    if offline_local_root:
        return
    command = manifest.get("command")
    run = manifest.get("run")
    allocation = manifest.get("allocation")
    if (
        not isinstance(command, dict)
        or not isinstance(run, dict)
        or not isinstance(allocation, dict)
    ):
        raise ValueError("canonical prepared manifest lacks submission authority metadata")
    bundle = command.get("source_bundle")
    utility = command.get("production_utility")
    if not isinstance(bundle, dict) or not isinstance(utility, dict):
        raise ValueError("canonical prepared manifest lacks F118 provenance")
    matrix_path = Path(str(command.get("matrix_file", ""))).expanduser().resolve()
    current_source_authority = require_current_source_authority_for_prepare(
        paths,
        bundle_provenance=bundle,
        utility_provenance=utility,
        matrix_path=matrix_path,
        input_revision=str(command.get("input_revision", "")),
        offline_local_root=False,
    )
    if (
        current_source_authority is None
        or command.get("current_source_authority") != current_source_authority
    ):
        raise ValueError("prepared manifest F118 current-source authority is stale")

    source_dir = Path(str(command.get("source_dir", ""))).expanduser().resolve()
    source_input = Path(
        str(command.get("source_input_file", ""))
    ).expanduser().resolve()
    try:
        source_input.relative_to(source_dir)
    except ValueError as error:
        raise ValueError("prepared source input is outside its source directory") from error
    require_file_sha256(
        source_input, command.get("input_sha256"), "source input"
    )
    restart_value = command.get("source_restart_file")
    restart = (
        Path(str(restart_value)).expanduser().resolve()
        if restart_value is not None else None
    )
    build_manifest = Path(
        str(command.get("build_manifest", ""))
    ).expanduser().resolve()
    qualification = command.get("qualification_approval")
    parent_segment = command.get("parent_segment")
    if qualification is not None and not isinstance(qualification, dict):
        raise ValueError("prepared qualification provenance is invalid")
    if parent_segment is not None and not isinstance(parent_segment, dict):
        raise ValueError("prepared parent-segment provenance is invalid")
    profile_args = argparse.Namespace(
        case_id=run.get("case_id"),
        segment=run.get("segment"),
        acceptance_criterion=run.get("acceptance_criterion"),
        athena_walltime=command.get("athena_walltime"),
        executable=command.get("executable"),
        nodes=allocation.get("nodes"),
        ranks_per_node=allocation.get("ranks_per_node"),
        cpus_per_task=allocation.get("cpus_per_task"),
        walltime=allocation.get("requested_walltime"),
    )
    build_provenance = {
        "revision": str(command.get("executable_revision", "")),
        "sha256": str(command.get("executable_sha256", "")),
        "manifest_dir": str(build_manifest),
    }
    if run.get("case_id") != R17_CASE_ID:
        if command.get("r17_readiness_evidence_chain") is not None:
            raise ValueError("non-R17 prepared manifest retains R17 readiness authority")
        if run.get("case_id") not in CONCURRENT_CASE_IDS:
            if command.get("promoted_profile_authority") is not None:
                raise ValueError("R02 prepared manifest retains unsupported profile authority")
            return
        fresh_profile = require_promoted_profile_for_prepare(
            paths,
            profile_args,
            source_dir=source_dir,
            matrix_path=matrix_path,
            input_path=source_input,
            input_revision=str(command.get("input_revision", "")),
            utility_provenance=utility,
            build_provenance=build_provenance,
            build_manifest=build_manifest,
            qualification_approval=qualification,
            bundle_provenance=bundle,
            current_source_authority=current_source_authority,
            restart=restart,
            parent_segment=parent_segment,
            time_tlim_target=command.get("time_tlim_target"),
        )
        if (
            fresh_profile is None
            or command.get("promoted_profile_authority") != fresh_profile
        ):
            raise ValueError("prepared F117-era reviewed profile authority is stale")
        return

    matching_indices = [
        index for index, item in enumerate(reservations)
        if isinstance(item, dict)
        and item.get("manifest") == str(manifest_path)
    ]
    if (
        len(matching_indices) != 1
        or reservations[matching_indices[0]] != reservation
        or reservation.get("case_id") != R17_CASE_ID
        or reservation.get("state") != "prepared"
    ):
        raise ValueError("R17 submission reservation snapshot differs")
    predecessor_reservations = [
        item for index, item in enumerate(reservations)
        if index != matching_indices[0]
    ]
    fresh_readiness = require_r17_readiness_for_prepare(
        paths,
        profile_args,
        predecessor_reservations,
        source_dir=source_dir,
        matrix_path=matrix_path,
        input_path=source_input,
        input_revision=str(command.get("input_revision", "")),
        utility_provenance=utility,
        build_provenance=build_provenance,
        build_manifest=build_manifest,
        qualification_approval=qualification,
        bundle_provenance=bundle,
        current_source_authority=current_source_authority,
        restart=restart,
        parent_segment=parent_segment,
        time_tlim_target=command.get("time_tlim_target"),
        reservation_snapshot_sha256=stable_json_sha256(predecessor_reservations),
        excluded_manifest_path=manifest_path,
    )
    if (
        fresh_readiness is None
        or command.get("r17_readiness_evidence_chain") != fresh_readiness
    ):
        raise ValueError("prepared R17 readiness/recost authority is stale")


def submission_preflight(args: argparse.Namespace, manifest_path: Path,
                         manifest: dict[str, object],
                         run_slurm_test: bool,
                         ) -> tuple[
                             dict[str, Path], dict[str, object], dict[str, object]
                         ]:
    """Authenticate and check one prepared segment immediately before submission."""

    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    validate_submission_fixture_options(args, offline_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    if manifest.get("state") != "prepared":
        raise ValueError("only a prepared segment can be submission-checked")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is unavailable for submission")
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    active = require_active_reservation_policy(reservations)
    prepared = [
        item for item in active if item.get("state") == "prepared"
    ]
    if prepared != [reservation]:
        raise ValueError("submission requires one matching prepared reservation")
    authenticate_prepared_execution(
        manifest, manifest_path, allow_legacy_local=offline_local_root
    )
    reauthenticate_submission_authority(
        paths,
        manifest_path,
        manifest,
        reservations,
        reservation,
        offline_local_root,
    )
    if (
        not offline_local_root
        or isinstance(manifest.get("command", {}).get("overrides"), list)
    ):
        prepared_time_tlim_target(manifest)
    actual, reserved = reservation_usage(paths)
    if actual + reserved > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("active Stage I reservation exceeds its ceiling")
    initial_queue_authentication = authenticated_production_queue_evidence(
        args, paths, reservations, offline_local_root
    )
    requested_campaigns = set(getattr(args, "allow_shared_root_campaign", []))
    shared_root_clearance = None
    allowed_shared_root_manifests: set[Path] = set()
    if not offline_local_root and requested_campaigns:
        run = manifest.get("run")
        if not isinstance(run, dict):
            raise ValueError("managed shared-root clearance requires manifest case")
        shared_root_clearance = require_managed_shared_root_clearance(
            paths,
            requested_campaigns,
            action=str(getattr(args, "action", "")),
            case_id=str(run.get("case_id", "")),
            queue_evidence=initial_queue_authentication,
        )
        allowed_shared_root_manifests.add(Path(
            str(shared_root_clearance["consumption"]["authorized_stale_manifest"])
        ))
    elif (
        offline_local_root
        and requested_campaigns == {SHARED_ROOT_STALE_CAMPAIGN_ID}
    ):
        allowed_shared_root_manifests.add(
            root / "runs" / SHARED_ROOT_STALE_CAMPAIGN_ID
            / "manifest/prepared_run.json"
        )
    conflicts = shared_root_campaign_conflicts(root, allowed_shared_root_manifests)
    if conflicts:
        raise ValueError(
            "shared-root campaign records require explicit review; pass "
            "--allow-shared-root-campaign only after confirming isolation: "
            + "; ".join(conflicts)
        )
    script = open_authenticated_batch_script(manifest, manifest_path)
    try:
        slurm_test_outcome = "not requested"
        if (
            run_slurm_test
            and not offline_local_root
            and not getattr(args, "skip_slurm_test", False)
        ):
            print(scheduler_test_only_output(script).strip())
            slurm_test_outcome = "passed"
        elif run_slurm_test and offline_local_root:
            slurm_test_outcome = "offline local-root fixture"
        elif run_slurm_test and getattr(args, "skip_slurm_test", False):
            slurm_test_outcome = "operator skipped with --skip-slurm-test"
        reauthenticate_open_batch_script(script)
        audit = {
            "created_utc": utc_now(),
            "offline_local_root": offline_local_root,
            "skip_slurm_test": bool(getattr(args, "skip_slurm_test", False)),
            "slurm_test_only": slurm_test_outcome,
            "acknowledged_shared_root_campaigns": sorted(requested_campaigns),
            "initial_queue_authentication": initial_queue_authentication,
            "batch_script": script["binding"],
        }
        if shared_root_clearance is not None:
            audit["shared_root_isolation_clearance"] = shared_root_clearance
        return paths, script, audit
    except BaseException:
        close_authenticated_batch_script(script)
        raise


def check_submit(args: argparse.Namespace) -> int:
    """Print a read-only, explicitly non-authoritative submission preview."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    _, script, _ = submission_preflight(
        args, manifest_path, manifest, run_slurm_test=True
    )
    try:
        print("Non-authoritative preview passed. Re-run the checks atomically with:")
        acknowledgements = "".join(
            " " + shlex.quote(f"--allow-shared-root-campaign={campaign}")
            for campaign in sorted(set(getattr(args, "allow_shared_root_campaign", [])))
        )
        print(
            f"  {shlex.quote(authenticated_python_binary())} -I -S -B "
            f"{shlex.quote(str(Path(__file__).resolve()))} submit "
            f"--manifest {shlex.quote(str(manifest_path))}{acknowledgements}"
        )
        print(f"Prepared script: {script['path']}")
        return 0
    finally:
        close_authenticated_batch_script(script)


def parse_sbatch_job_id(output: str) -> str:
    """Return the top-level numeric ID from ``sbatch --parsable`` output."""

    match = re.fullmatch(r"\s*([1-9][0-9]*)(?:;[^\s;]+)?\s*", output)
    if match is None:
        raise ValueError(f"sbatch did not return one numeric job ID: {output!r}")
    return require_numeric_job_id(match.group(1))


def scheduler_submit_time_evidence(value: object,
                                   transaction: dict[str, object]) -> str:
    """Require scheduler submit time within five minutes of the sbatch barrier."""

    submitted = parse_utc_timestamp(value, "scheduler submit time")
    barrier = parse_utc_timestamp(
        transaction.get("created_utc"), "submission ambiguity barrier time"
    )
    if abs((submitted - barrier).total_seconds()) > (
        SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS
    ):
        raise ValueError(
            "scheduler recovery submit time is outside the symmetric "
            f"{SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS}-second ambiguity "
            "barrier tolerance"
        )
    return submitted.isoformat()


def verify_recovered_scheduler_job(manifest: dict[str, object], job_id: str,
                                   offline_local_root: bool,
                                   transaction: dict[str, object],
                                   ) -> dict[str, object]:
    """Require an operator-supplied recovery ID to name the prepared job."""

    require_numeric_job_id(job_id)
    if offline_local_root:
        return {
            "mode": "offline-local fixture",
            "checked_utc": utc_now(),
            "job_id": job_id,
        }
    expected_name = expected_job_name(manifest)
    expected_script = Path(str(manifest["paths"]["batch_script"])).resolve()
    control = subprocess.run(
        [str(SCONTROL), "show", "job", "-o", job_id],
        check=False,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
    )
    if control.returncode == 0:
        fields = {
            key: value
            for key, value in (
                token.split("=", 1)
                for token in shlex.split(control.stdout)
                if "=" in token
            )
        }
        command = fields.get("Command")
        if (
            fields.get("JobId") == job_id
            and fields.get("JobName") == expected_name
            and str(fields.get("Account", "")).casefold() == ACCOUNT.casefold()
            and fields.get("Partition") == PARTITION
            and (
                command in {None, "", "(null)", "N/A"}
                or Path(command).expanduser().resolve() == expected_script
            )
        ):
            submitted = scheduler_submit_time_evidence(
                fields.get("SubmitTime"), transaction
            )
            return {
                "mode": "scontrol",
                "checked_utc": utc_now(),
                "job_id": job_id,
                "job_name": expected_name,
                "account": fields["Account"],
                "partition": fields["Partition"],
                "command": command,
                "submit_time": submitted,
            }
    accounting = subprocess.run(
        [
            str(SACCT), "-X", "-j", job_id,
            "--format=JobIDRaw,JobName,Account,Partition,Submit", "-n", "-P",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
    ).stdout
    rows = [
        row for row in csv.reader(accounting.splitlines(), delimiter="|")
        if row and row[0] == job_id
    ]
    if (
        len(rows) != 1
        or len(rows[0]) < 5
        or rows[0][1] != expected_name
        or rows[0][2].casefold() != ACCOUNT.casefold()
        or rows[0][3] != PARTITION
    ):
        raise ValueError("scheduler recovery job does not match the prepared segment")
    submitted = scheduler_submit_time_evidence(rows[0][4], transaction)
    return {
        "mode": "sacct",
        "checked_utc": utc_now(),
        "job_id": job_id,
        "job_name": expected_name,
        "account": rows[0][2],
        "partition": rows[0][3],
        "submit_time": submitted,
    }


def scheduler_absence_evidence(args: argparse.Namespace,
                               manifest: dict[str, object],
                               transaction: dict[str, object],
                               offline_local_root: bool,
                               ) -> dict[str, object]:
    """Retain scheduler-side evidence before clearing an ambiguous submission."""

    expected_name = expected_job_name(manifest)
    checked_utc = utc_now()
    fixture = getattr(args, "scheduler_absence_evidence_file", None)
    if fixture:
        if not offline_local_root:
            raise ValueError(
                "--scheduler-absence-evidence-file is restricted to offline validation"
            )
        value = json.loads(Path(fixture).read_text(encoding="utf-8"))
        if not isinstance(value, dict) or value.get("absent") is not True:
            raise ValueError("offline scheduler absence fixture must assert absent=true")
        return {
            "mode": "offline-local fixture",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
            "fixture": value,
        }
    if offline_local_root:
        return {
            "mode": "offline-local operator confirmation",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
        }
    user = scheduler_account_name()
    barrier = parse_utc_timestamp(
        transaction.get("created_utc"), "submission ambiguity barrier time"
    )
    scheduler_start = (
        barrier - timedelta(seconds=SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS)
    ).strftime(
        "%Y-%m-%dT%H:%M:%S"
    )
    squeue_command = [
        str(SQUEUE), "-h", "-u", user, "-o", "%i|%j|%a|%P|%V|%o",
    ]
    sacct_command = [
        str(SACCT), "-X", "-S", scheduler_start,
        "--format=JobIDRaw,JobName,Account,Partition,Submit", "-n", "-P",
    ]
    try:
        queued = subprocess.run(
            squeue_command,
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
        accounted = subprocess.run(
            sacct_command,
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        break_glass = getattr(args, "break_glass_clear_evidence", "")
        if (
            not getattr(args, "confirm_break_glass_clear", False)
            or not str(break_glass).strip()
        ):
            raise ValueError(
                "scheduler absence query failed; explicit break-glass evidence "
                "and --confirm-break-glass-clear are required"
            ) from error
        return {
            "mode": "break-glass after scheduler query failure",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
            "operator_evidence": break_glass,
            "query_error": str(error),
        }
    if (
        scheduler_output_contains_job(queued, expected_name)
        or scheduler_output_contains_job(accounted, expected_name)
    ):
        raise ValueError("scheduler still reports a matching ambiguous submission")
    return {
        "mode": "live scheduler absence query",
        "checked_utc": checked_utc,
        "expected_job_name": expected_name,
        "ambiguity_created_utc": transaction["created_utc"],
        "squeue_command": squeue_command,
        "squeue_output": queued,
        "sacct_command": sacct_command,
        "sacct_output": accounted,
    }


@locked_manifest_action
def submit(args: argparse.Namespace) -> int:
    """Atomically authenticate, submit, and retain one scheduler job ID."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    output_file = getattr(args, "sbatch_output_file", None)
    if output_file and not offline_local_root:
        raise ValueError("--sbatch-output-file is restricted to offline validation")
    paths, script, audit = submission_preflight(
        args, manifest_path, manifest, run_slurm_test=True
    )
    try:
        audit["final_queue_authentication"] = authenticated_production_queue_evidence(
            args, paths, read_reservations(paths), offline_local_root
        )
        reauthenticate_open_batch_script(script, audit["batch_script"])
        transaction_path = write_submit_pending_transaction(
            paths, manifest_path, audit
        )
        reauthenticate_open_batch_script(script, audit["batch_script"])
        if output_file:
            output = Path(output_file).read_text(encoding="utf-8")
        else:
            output = subprocess.run(
                [str(SBATCH), "--parsable", str(script["descriptor_path"])],
                check=True,
                capture_output=True,
                text=True,
                env=scheduler_environment(),
                pass_fds=(int(script["fd"]),),
            ).stdout
        reauthenticate_open_batch_script(script, audit["batch_script"])
        job_id = parse_sbatch_job_id(output)
        finish_submit_transaction(
            paths, transaction_path, manifest_path, manifest,
            read_reservations(paths), job_id, authenticated_batch=script,
        )
        print(f"Submitted Stage I job {job_id}: {script['path']}")
        return 0
    finally:
        close_authenticated_batch_script(script)


@locked_manifest_action
def mark_submitted(args: argparse.Namespace) -> int:
    """Attach a fake scheduler ID only for legacy offline validation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    if not is_offline_local_root(root, args.allow_local_root):
        raise ValueError("production submissions must use the atomic submit action")
    paths = initialize(root)
    require_numeric_job_id(args.job_id)
    if manifest.get("state") != "prepared":
        raise ValueError("manifest is not in prepared state")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is not in prepared state")
    transaction_path = write_submit_pending_transaction(
        paths, manifest_path, {
            "created_utc": utc_now(),
            "offline_local_root": True,
            "legacy_mark_submitted": True,
            "skip_slurm_test": True,
            "slurm_test_only": "offline local-root fixture",
            "acknowledged_shared_root_campaigns": [],
        }
    )
    finish_submit_transaction(
        paths, transaction_path, manifest_path, manifest, reservations, args.job_id
    )
    print(f"Marked Stage I job {args.job_id} submitted.")
    return 0


@locked_manifest_action
def recover_submit(args: argparse.Namespace) -> int:
    """Resolve an ambiguous sbatch boundary with its retained scheduler ID."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    transaction_path = submit_pending_transaction(paths, manifest_path)
    transaction = read_transaction(paths, transaction_path)
    evidence = verify_recovered_scheduler_job(
        manifest, args.job_id, offline_local_root, transaction
    )
    retained_evidence = manifest.setdefault("scheduler_recovery_evidence", [])
    if not isinstance(retained_evidence, list):
        raise ValueError("prepared manifest has invalid scheduler recovery evidence")
    retained_evidence.append(evidence)
    finish_submit_transaction(
        paths, transaction_path, manifest_path, manifest,
        read_reservations(paths), args.job_id,
    )
    print(f"Recovered submitted Stage I job {args.job_id}.")
    return 0


@locked_manifest_action
def clear_submit_pending(args: argparse.Namespace) -> int:
    """Clear an ambiguous submit barrier after confirming no job was launched."""

    if not args.confirm_no_job_submitted:
        raise ValueError("--confirm-no-job-submitted is required")
    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    transaction_path = submit_pending_transaction(paths, manifest_path)
    transaction = read_transaction(paths, transaction_path)
    if transaction.get("prepared_manifest_sha256") != sha256(manifest_path):
        raise ValueError("prepared manifest changed after the sbatch boundary")
    if transaction.get("prior_reservations_sha256") != sha256(paths["reservations"]):
        raise ValueError("reservation store changed after the sbatch boundary")
    absence_evidence = scheduler_absence_evidence(
        args, manifest, transaction, offline_local_root
    )
    notes = manifest.setdefault("submission_recovery_notes", [])
    if not isinstance(notes, list):
        raise ValueError("prepared manifest has invalid submission recovery notes")
    notes.append({
        "cleared_utc": utc_now(),
        "notes": args.notes,
        "outcome": "operator confirmed no scheduler job was submitted",
        "scheduler_absence_evidence": absence_evidence,
    })
    transaction.update({
        "kind": "submit_cleared",
        "manifest": manifest,
        "reservations": read_reservations(paths),
        "ledger_row": None,
        "recovery_notes": args.notes,
        "scheduler_absence_evidence": absence_evidence,
    })
    write_json(transaction_path, transaction, mode=0o644)
    apply_transaction(paths, transaction_path)
    print(f"Cleared ambiguous submission barrier: {manifest_path}")
    return 0


@locked_root_action
def recover_transactions(args: argparse.Namespace) -> int:
    """Replay deterministic journals while preserving ambiguous submissions."""

    root = require_root(Path(args.root), args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    recovered = 0
    for transaction_path in pending_transaction_paths(paths):
        apply_transaction(paths, transaction_path)
        recovered += 1
    print(f"Recovered {recovered} deterministic Stage I transaction(s).")
    return 0


def parse_history_bytes(retained: bytes, label: str) -> dict[str, list[float]]:
    """Parse exact retained AthenaK history bytes keyed by labeled columns."""

    try:
        text = retained.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"history file is not valid UTF-8: {label}") from error
    labels: list[str] = []
    rows: list[list[float]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            found = re.findall(r"\[(\d+)\]=([^\s]+)", line)
            if found:
                indices = [int(index) for index, _ in found]
                candidate = [label for _, label in found]
                if (
                    indices != list(range(1, len(candidate) + 1))
                    or len(candidate) != len(set(candidate))
                    or (labels and candidate != labels)
                ):
                    raise ValueError(
                        f"history labels are ambiguous or invalid: {label}"
                    )
                labels = candidate
            continue
        if line.strip():
            try:
                rows.append([float(value) for value in line.split()])
            except ValueError as error:
                raise ValueError(
                    f"history file contains invalid numeric data: {label}"
                ) from error
    if not labels or not rows:
        raise ValueError(f"history file is missing labels or data: {label}")
    if any(len(row) != len(labels) for row in rows):
        raise ValueError(f"history row width does not match labels: {label}")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise ValueError(f"history file contains a non-finite value: {label}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def parse_history(path: Path) -> dict[str, list[float]]:
    """Read one AthenaK history file keyed by its labeled columns."""

    return parse_history_bytes(path.read_bytes(), str(path))


def _continuation_plasma_evidence(
    case_id: str,
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    *,
    maximum_divb: float | None,
    normalized_ct_divb_evidence: dict[str, object],
) -> dict[str, object]:
    """Evaluate clean-partial histories against accepted plasma policy."""

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
    if (
        mhd_mass_drift > CONTINUATION_MASS_TOLERANCE
        or user_mass_drift > CONTINUATION_MASS_TOLERANCE
        or mass_mismatch > CONTINUATION_MASS_TOLERANCE
    ):
        raise ValueError("continuation mass conservation exceeds accepted policy")
    if (
        maximum_divb is not None
        and (
            maximum_divb < 0.0
            or maximum_divb >= CONTINUATION_MAX_NORMALIZED_CT_DIVB
        )
    ):
        raise ValueError("continuation normalized CT divB exceeds accepted policy")
    if any(value != 0.0 for value in user["hard_vol"]):
        raise ValueError("continuation hard_vol is nonzero")
    strict_maxima = {
        name: max(mhd[name]) for name in STRICT_LF_FAILURE_COLUMNS
    }
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
    passive = case_id in {"R06", "R07", "R08", "R09"}
    finite_limiter = case_id in {"R14", "R15"}
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
    return {
        "schema_version": 2,
        "policy": CONTINUATION_PLASMA_POLICY,
        "case_id": case_id,
        "continuation_authorized": True,
        "checks": {
            "finite_synchronized_histories": True,
            "mass_conservation": True,
            "normalized_ct_divb": maximum_divb is not None,
            "strict_lf_policy": True,
            "forcing_policy": True,
            "pressure_feedback_policy": True,
            "limiter_policy": True,
        },
        "normalized_ct_divb_evidence": normalized_ct_divb_evidence,
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


def continuation_plasma_evidence(
    case_id: str,
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
) -> dict[str, object]:
    """Require ordinary clean-partial histories to satisfy accepted plasma policy."""

    if "max_ndiv" not in user:
        raise ValueError(
            "continuation normalized CT divB is irreducibly unavailable for the "
            "qualified historical E03 executable: user history lacks max_ndiv, "
            "retained mhd_w_bcc snapshots are cell-centered, and retained restart "
            "face fields have no authenticated CT-divergence interpretation"
        )
    maximum_divb = max(user["max_ndiv"])
    return _continuation_plasma_evidence(
        case_id,
        mhd,
        user,
        maximum_divb=maximum_divb,
        normalized_ct_divb_evidence={
            "source": "authenticated user history max_ndiv",
            "comparison": "strictly-less-than",
            "threshold": CONTINUATION_MAX_NORMALIZED_CT_DIVB,
            "maximum": maximum_divb,
            "passed": True,
        },
    )


def read_frozen_e03_migration_binding(
    binding: dict[str, object],
    label: str,
) -> tuple[Path, bytes, dict[str, object]]:
    """Read one exact owner-controlled frozen-E03 migration artifact."""

    relative_text = str(binding.get("path", ""))
    relative = Path(relative_text)
    if (
        not relative_text
        or relative.is_absolute()
        or ".." in relative.parts
        or relative.as_posix() != relative_text
    ):
        raise ValueError(f"{label} path is not normalized and root-relative")
    path = DEFAULT_ROOT.expanduser().absolute() / relative
    try:
        mode = int(binding["mode"])
        size_bytes = int(binding["size_bytes"])
        expected_sha256 = str(binding["sha256"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label} binding is invalid") from error
    retained = read_r17_evidence_bytes(
        path,
        label,
        mode=mode,
        owner_controlled=True,
        symlink_free=True,
    )
    if len(retained) != size_bytes:
        raise ValueError(f"{label} size differs from frozen-E03 migration binding")
    if hashlib.sha256(retained).hexdigest() != expected_sha256:
        raise ValueError(f"{label} checksum differs from frozen-E03 migration binding")
    return path, retained, {
        "path": str(path),
        "sha256": expected_sha256,
        "size_bytes": size_bytes,
        "mode": f"{mode:04o}",
        "links": 1,
    }


def read_frozen_e03_migration_json(
    binding: dict[str, object],
    label: str,
) -> tuple[dict[str, object], dict[str, object]]:
    """Read and decode one exact frozen-E03 migration JSON artifact."""

    _, retained, public_binding = read_frozen_e03_migration_binding(binding, label)
    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be a JSON object")
    return value, public_binding


def validate_frozen_e03_independent_validation(
    contract: dict[str, object],
    validation: dict[str, object],
    review: dict[str, object] | None,
    bindings: dict[str, dict[str, object] | None],
) -> None:
    """Require exact independent evidence to bind one migrated clean partial."""

    manifest = contract["manifest"]
    inspection = contract["inspection"]
    validation_contract = contract["independent_validation"]
    if not all(
        isinstance(item, dict)
        for item in (manifest, inspection, validation_contract)
    ):
        raise ValueError("frozen-E03 migration contract is invalid")
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
        raise ValueError(
            "frozen-E03 independent validation differs from eligible segment"
        )
    product_sha256 = validation.get("product_sha256")
    product_sizes = validation.get("product_sizes")
    if not isinstance(product_sha256, dict) or not isinstance(product_sizes, dict):
        raise ValueError("frozen-E03 independent validation lacks product inventory")
    for key in ("mhd_history", "user_history"):
        history_binding = bindings.get(key)
        if (
            not isinstance(history_binding, dict)
            or product_sha256.get(history_binding["path"])
            != history_binding["sha256"]
            or product_sizes.get(history_binding["path"])
            != history_binding["size_bytes"]
        ):
            raise ValueError(
                "frozen-E03 independent validation history binding differs"
            )
    zero_strict = {name: 0 for name in STRICT_LF_FAILURE_COLUMNS}
    if validation.get("strict_failure_maxima") != zero_strict:
        raise ValueError(
            "frozen-E03 independent validation strict LF evidence differs"
        )
    if case_id == "R03":
        if (
            validation.get("clean_for_continuation") is not True
            or validation.get("formal_inspection_accepted") is not False
            or validation.get("product_count")
            != validation.get("expected_product_count")
            or review is not None
        ):
            raise ValueError("R03 frozen-E03 migration validation differs")
        return
    expected_ct = {
        "authorizing": False,
        "executable_revision": FROZEN_E03_EXECUTABLE["revision"],
        "executable_sha256": FROZEN_E03_EXECUTABLE["sha256"],
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


def frozen_e03_migrated_inspection(
    inspection: dict[str, object],
    manifest: dict[str, object],
) -> dict[str, object] | None:
    """Reproduce full schema-4 plasma evidence for exact retained E03 partials."""

    run = manifest.get("run")
    if not isinstance(run, dict):
        return None
    case_id = str(run.get("case_id", ""))
    job_id = str(manifest.get("job_id", ""))
    retained_contract = FROZEN_E03_CONTINUATION_MIGRATIONS.get((case_id, job_id))
    if retained_contract is None:
        return None
    contract = {
        **retained_contract,
        "case_id": case_id,
        "job_id": job_id,
    }
    accounting = manifest.get("accounting")
    command = manifest.get("command")
    if (
        manifest.get("project_root") != str(DEFAULT_ROOT)
        or manifest.get("execution_epoch") != EXECUTION_EPOCH
        or run.get("segment") != contract["segment"]
        or manifest.get("state") != "recorded"
        or not isinstance(accounting, dict)
        or accounting.get("result") != "clean_partial"
        or accounting.get("job_id") != job_id
        or not isinstance(command, dict)
        or command.get("executable_revision")
        != FROZEN_E03_EXECUTABLE["revision"]
        or command.get("executable_sha256")
        != FROZEN_E03_EXECUTABLE["sha256"]
    ):
        raise ValueError("frozen-E03 migration manifest identity differs")
    exact_manifest, manifest_binding = read_frozen_e03_migration_json(
        contract["manifest"], "frozen-E03 migration manifest"
    )
    exact_inspection, inspection_binding = read_frozen_e03_migration_json(
        contract["inspection"], "frozen-E03 migration inspection"
    )
    if exact_manifest != manifest or exact_inspection != inspection:
        raise ValueError("frozen-E03 migration manifest or inspection bytes differ")
    if (
        exact_manifest.get("scientific_inspection") != exact_inspection
        or inspection.get("schema_version") != 4
        or inspection.get("case_id") != case_id
        or inspection.get("job_id") != job_id
        or inspection.get("segment") != contract["segment"]
        or inspection.get("clean_for_continuation") is not True
        or inspection.get("accepted") is not False
        or inspection.get("manifest") != manifest_binding["path"]
    ):
        raise ValueError("frozen-E03 migration inspection identity differs")
    bindings: dict[str, dict[str, object] | None] = {
        "manifest": manifest_binding,
        "inspection": inspection_binding,
    }
    histories: dict[str, dict[str, list[float]]] = {}
    for key in ("mhd_history", "user_history"):
        path, retained, public_binding = read_frozen_e03_migration_binding(
            contract[key], f"frozen-E03 migration {key}"
        )
        bindings[key] = public_binding
        expected_record = {
            "path": str(path),
            "size_bytes": public_binding["size_bytes"],
            "sha256": public_binding["sha256"],
        }
        if inspection.get(key) != expected_record:
            raise ValueError(f"frozen-E03 migration {key} inspection binding differs")
        histories[key] = parse_history_bytes(retained, str(path))
    validation, validation_binding = read_frozen_e03_migration_json(
        contract["independent_validation"],
        "frozen-E03 migration independent validation",
    )
    bindings["independent_validation"] = validation_binding
    review_contract = contract.get("independent_review")
    if review_contract is None:
        review = None
        bindings["independent_review"] = None
    else:
        if not isinstance(review_contract, dict):
            raise ValueError("frozen-E03 migration independent review is invalid")
        review, review_binding = read_frozen_e03_migration_json(
            review_contract, "frozen-E03 migration independent review"
        )
        bindings["independent_review"] = review_binding
    validate_frozen_e03_independent_validation(
        contract, validation, review, bindings
    )
    normalized_ct_divb_evidence = {
        "ct_divergence_claimed": False,
        "status": "historically_unavailable",
        "reason": FROZEN_E03_CT_DIVERGENCE_REASON,
        "authorizing": False,
        "source": "authenticated frozen-E03 migration contract",
        "qualified_executable": dict(FROZEN_E03_EXECUTABLE),
        "independent_validation": validation_binding,
        "independent_review": bindings["independent_review"],
    }
    plasma_evidence = _continuation_plasma_evidence(
        case_id,
        histories["mhd_history"],
        histories["user_history"],
        maximum_divb=None,
        normalized_ct_divb_evidence=normalized_ct_divb_evidence,
    )
    plasma_evidence["policy"] = FROZEN_E03_CONTINUATION_MIGRATION_POLICY
    plasma_evidence["continuation_authorized"] = False
    plasma_evidence["continuation_eligible"] = False
    plasma_evidence["eligibility_only"] = True
    plasma_evidence["ct_divergence_claimed"] = False
    plasma_evidence["ct_divergence_reason"] = FROZEN_E03_CT_DIVERGENCE_REASON
    plasma_evidence["authorization_basis"] = (
        "none; historical frozen-E03 evidence is inventory-only and cannot "
        "authorize continuation"
    )
    plasma_evidence["checks"]["frozen_e03_exact_migration"] = True
    plasma_evidence["migration_contract"] = {
        "schema_version": 1,
        "policy": FROZEN_E03_CONTINUATION_MIGRATION_POLICY,
        "case_id": case_id,
        "job_id": job_id,
        "segment": contract["segment"],
        "frozen_science_executable": dict(FROZEN_E03_EXECUTABLE),
        "bindings": bindings,
        "authority": {
            "continuation_authorized": False,
            "submission_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
    }
    migrated = dict(inspection)
    migrated["checks"] = dict(inspection["checks"])
    migrated["checks"]["plasma_continuation_policy"] = False
    migrated["plasma_continuation_policy"] = (
        FROZEN_E03_CONTINUATION_MIGRATION_POLICY
    )
    migrated["plasma_continuation_authorized"] = False
    migrated["plasma_continuation_eligible"] = False
    migrated["plasma_continuation_evidence"] = plasma_evidence
    return migrated


def revalidate_continuation_plasma_evidence(
    inspection: dict[str, object],
    manifest: dict[str, object],
) -> dict[str, object]:
    """Reparse retained histories before propagating a clean-partial state."""

    migrated = frozen_e03_migrated_inspection(inspection, manifest)
    if migrated is not None:
        return migrated
    mhd_record = inspection.get("mhd_history")
    user_record = inspection.get("user_history")
    if not isinstance(mhd_record, dict) or not isinstance(user_record, dict):
        raise ValueError("continuation inspection lacks retained histories")
    case_id = str(manifest.get("run", {}).get("case_id", ""))
    current = continuation_plasma_evidence(
        case_id,
        parse_history(Path(str(mhd_record.get("path", ""))).resolve()),
        parse_history(Path(str(user_record.get("path", ""))).resolve()),
    )
    if (
        inspection.get("schema_version") == 4
        and inspection.get("plasma_continuation_policy") != CONTINUATION_PLASMA_POLICY
    ):
        raise ValueError("schema-4 continuation plasma policy is absent or changed")
    if (
        inspection.get("plasma_continuation_evidence") != current
        or current.get("continuation_authorized") is not True
    ):
        raise ValueError("continuation plasma evidence is absent, stale, or changed")
    return inspection


def retained_file(path: Path) -> dict[str, object]:
    """Describe an output retained for acceptance review."""

    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def output_product_groups(directory: Path, pattern: str,
                          expected_ranks: int | None = None) -> list[list[Path]]:
    """Collect shared outputs or complete rank-local product sets."""

    groups = [[path] for path in sorted(directory.glob(pattern))]
    rank0_dir = directory / "rank_00000000"
    if rank0_dir.is_dir():
        for rank0 in sorted(rank0_dir.glob(pattern)):
            rank_files = sorted(directory.glob(f"rank_*/{rank0.name}"))
            expected_names = [
                f"rank_{rank:08d}" for rank in range(len(rank_files))
            ]
            if [path.parent.name for path in rank_files] != expected_names:
                raise ValueError(f"rank-local output set is not contiguous: {rank0}")
            if expected_ranks is not None and len(rank_files) != expected_ranks:
                raise ValueError(
                    f"rank-local output set has {len(rank_files)} files; "
                    f"expected {expected_ranks}: {rank0}"
                )
            groups.append(rank_files)
    return sorted(groups, key=lambda group: str(group[0]))


def retained_product(group: list[Path]) -> dict[str, object]:
    """Describe one shared or rank-local retained output product."""

    representative = retained_file(group[0])
    if group[0].parent.name.startswith("rank_"):
        representative["storage"] = "per_rank"
        representative["rank_files"] = [retained_file(path) for path in group]
    else:
        representative["storage"] = "shared_mpiio"
    return representative


def retained_product_paths(record: dict[str, object]) -> list[Path]:
    """Return every file belonging to an inspected retained product."""

    rank_files = record.get("rank_files")
    if isinstance(rank_files, list):
        if not rank_files or not all(isinstance(item, dict) for item in rank_files):
            raise ValueError("inspection-retained rank-local product is invalid")
        return [Path(str(item["path"])).resolve() for item in rank_files]
    return [Path(str(record["path"])).resolve()]


def revalidate_retained_file(record: object,
                             label: str = "inspection-retained file") -> None:
    """Require a retained file to preserve size and digest."""

    if not isinstance(record, dict):
        raise ValueError(f"{label} record is invalid")
    path = Path(str(record.get("path", ""))).resolve()
    if not path.is_file():
        raise ValueError(f"{label} is missing: {path}")
    if record.get("size_bytes") != path.stat().st_size:
        raise ValueError(f"{label} size has changed: {path}")
    if record.get("sha256") != sha256(path):
        raise ValueError(f"{label} checksum has changed: {path}")


def revalidate_retained_product(record: object) -> None:
    """Require every member of one inspected output product to be unchanged."""

    revalidate_retained_file(record)
    if not isinstance(record, dict):
        raise ValueError("inspection-retained product record is invalid")
    rank_files = record.get("rank_files")
    if rank_files is None:
        return
    if not isinstance(rank_files, list) or not rank_files:
        raise ValueError("inspection-retained rank-local product is invalid")
    for rank_file in rank_files:
        revalidate_retained_file(rank_file)


def recorded_product_groups(records: object, label: str) -> list[list[Path]]:
    """Return the exact inspected path groups for one retained output class."""

    if not isinstance(records, list):
        raise ValueError(f"segment inspection lacks retained {label}")
    groups = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError(f"segment inspection has invalid retained {label}")
        groups.append(retained_product_paths(record))
    return groups


def output_group_signature(groups: list[list[Path]]) -> list[tuple[str, ...]]:
    """Return a stable path-only signature for grouped output inventory."""

    return sorted(tuple(str(path.resolve()) for path in group) for group in groups)


def require_complete_output_inventory(directory: Path, pattern: str,
                                      groups: list[list[Path]],
                                      label: str) -> None:
    """Reject product files omitted by shared or rank-local grouping."""

    actual = sorted(
        path.resolve() for path in directory.rglob(pattern) if path.is_file()
    )
    grouped = sorted(path.resolve() for group in groups for path in group)
    if actual != grouped:
        raise ValueError(f"{label} output inventory contains ungrouped files")


def revalidate_inspection_inventory(manifest: dict[str, object],
                                    inspection: dict[str, object]) -> None:
    """Reject additions, removals, or regrouping after formal inspection."""

    output_dir = Path(str(manifest["paths"]["output_dir"])).resolve()
    expected_ranks = (
        int(manifest["allocation"]["nodes"])
        * int(manifest["allocation"].get("ranks_per_node", 1))
    )
    histories = {
        "mhd_history": sorted(path.resolve() for path in output_dir.glob("*.mhd.hst")),
        "user_history": sorted(path.resolve() for path in output_dir.glob("*.user.hst")),
    }
    for key, actual in histories.items():
        record = inspection.get(key)
        if not isinstance(record, dict):
            raise ValueError(f"segment inspection lacks retained {key}")
        expected = [Path(str(record.get("path", ""))).resolve()]
        if actual != expected:
            raise ValueError(f"inspection-retained {key} inventory has changed")
    for key, directory, pattern in (
        ("snapshots", output_dir / "bin", "*.bin"),
        ("restarts", output_dir / "rst", "*.rst"),
    ):
        expected = recorded_product_groups(inspection.get(key), key)
        actual = output_product_groups(directory, pattern, expected_ranks)
        require_complete_output_inventory(directory, pattern, actual, key)
        if output_group_signature(actual) != output_group_signature(expected):
            raise ValueError(f"inspection-retained {key} inventory has changed")


def revalidate_inspection_restart_times(inspection: dict[str, object],
                                        manifest: dict[str, object] | None,
                                        allow_legacy_local: bool) -> None:
    """Reparse retained restart markers and bind the terminal product to final time."""

    if allow_legacy_local and "restart_times" not in inspection:
        return
    groups = recorded_product_groups(inspection.get("restarts"), "restarts")
    bypass = (
        allow_legacy_local
        and inspection.get("restart_time_marker_bypass") is True
    )
    schema_version = int(inspection.get("schema_version", 0))
    if schema_version >= 4 and not bypass:
        if not isinstance(manifest, dict):
            raise ValueError(
                "binary-authenticated inspection lacks prepared manifest"
            )
        authenticated = [
            authenticated_restart_product_time(group, manifest.get("command"))
            for group in groups
        ]
        parsed = [float(item["binary_time"]) for item in authenticated]
        marker_modes = [item["marker_modes"] for item in authenticated]
        if inspection.get("restart_time_marker_modes") != marker_modes:
            raise ValueError(
                "segment inspection restart marker authentication has changed"
            )
    else:
        parsed = [
            restart_product_time(group, allow_missing_marker=bypass)
            for group in groups
        ]
    retained = inspection.get("restart_times")
    if not isinstance(retained, list) or len(retained) != len(parsed):
        raise ValueError("segment inspection restart-time evidence is incomplete")
    for expected, actual in zip(retained, parsed):
        if expected is None or actual is None:
            if not bypass or expected is not None or actual is not None:
                raise ValueError("segment inspection restart-time bypass is inconsistent")
        else:
            try:
                difference = abs(float(expected) - actual)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "segment inspection restart-time evidence is invalid"
                ) from error
            if difference > 1.0e-12:
                raise ValueError("segment inspection restart-time evidence has changed")
    if bypass:
        return
    try:
        final_time = float(inspection["final_time"])
        terminal_time = float(inspection["terminal_restart_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("segment inspection lacks terminal restart time") from error
    matches = [
        index for index, value in enumerate(parsed)
        if value is not None and abs(value - final_time) <= 1.0e-10
    ]
    if len(matches) != 1 or abs(terminal_time - final_time) > 1.0e-10:
        raise ValueError("segment inspection terminal restart time differs from final time")
    terminal = inspection.get("terminal_restart")
    if not isinstance(terminal, dict):
        raise ValueError("segment inspection lacks terminal restart product")
    if (
        output_group_signature([retained_product_paths(terminal)])
        != output_group_signature([groups[matches[0]]])
    ):
        raise ValueError("segment inspection terminal restart product is inconsistent")


def revalidate_inspection_files(inspection: dict[str, object],
                                manifest: dict[str, object] | None = None) -> None:
    """Recheck every file retained by an accepted or clean-partial inspection."""

    for key in ("mhd_history", "user_history"):
        revalidate_retained_file(inspection.get(key))
    for key in ("snapshots", "restarts"):
        records = inspection.get(key)
        if not isinstance(records, list):
            raise ValueError(f"segment inspection lacks retained {key}")
        for record in records:
            revalidate_retained_product(record)
    allow_legacy_local = (
        manifest is not None
        and Path(str(manifest.get("project_root", ""))).resolve()
        != DEFAULT_ROOT.expanduser().resolve()
    )
    revalidate_inspection_restart_times(inspection, manifest, allow_legacy_local)
    if manifest is not None:
        revalidate_inspection_inventory(manifest, inspection)


def merge_history_files(sources: list[Path], destination: Path) -> None:
    """Merge restart-segment histories while removing repeated boundary rows."""

    reference_labels: list[str] | None = None
    header: list[str] = []
    retained_rows: list[str] = []
    last_time = float("-inf")
    for source_index, source in enumerate(sources):
        source_header: list[str] = []
        source_labels: list[str] = []
        source_rows: list[str] = []
        for line in source.read_text(encoding="utf-8").splitlines():
            if line.startswith("#"):
                source_header.append(line)
                found = re.findall(r"\[\d+\]=([^\s]+)", line)
                if found:
                    source_labels = found
            elif line.strip():
                source_rows.append(line)
        if not source_labels:
            raise ValueError(f"history file has no labeled header: {source}")
        if reference_labels is None:
            reference_labels = source_labels
            header = source_header
        elif source_labels != reference_labels:
            raise ValueError("restart-segment history columns do not match")
        for line in source_rows:
            time = float(line.split()[0])
            if time > last_time + 1.0e-12:
                retained_rows.append(line)
                last_time = time
        if source_index == 0 and not retained_rows:
            raise ValueError(f"history file has no retained rows: {source}")
    write_text(destination, "\n".join([*header, *retained_rows]) + "\n")


def binary_snapshot_time(path: Path) -> float:
    """Read a retained Athena binary snapshot time without loading field data."""

    with path.open("rb") as stream:
        code_header = stream.readline().split()
        if not code_header or code_header[0] != b"Athena":
            raise ValueError(f"invalid Athena binary snapshot: {path}")
        pheader_count = int(stream.readline().split(b"=")[-1])
        values: dict[str, str] = {}
        for _ in range(pheader_count - 1):
            key, value = [
                token.strip()
                for token in stream.readline().decode("utf-8").split("=", 1)
            ]
            values[key] = value
    if "time" not in values:
        raise ValueError(f"binary snapshot has no time field: {path}")
    return float(values["time"])


def binary_product_time(group: list[Path]) -> float:
    """Require every rank-local member of a snapshot to record one time."""

    times = [binary_snapshot_time(path) for path in group]
    if any(abs(time - times[0]) > 1.0e-12 for time in times[1:]):
        raise ValueError(f"rank-local snapshot times disagree: {group[0]}")
    return times[0]


def analysis_model_choices(input_path: Path) -> dict[str, str]:
    """Use the workflow's interpretation mapping for a production bundle."""

    workflow_path = ROOT_DIR / "scripts/cgl_lf_workflow.py"
    spec = importlib.util.spec_from_file_location(
        "_cgl_lf_workflow_stage_i", workflow_path
    )
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load workflow model mapping: {workflow_path}")
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
    return workflow.model_choices(input_path.read_text(encoding="utf-8"), [])


@locked_manifest_action
def inspect_segment(args: argparse.Namespace) -> int:
    """Inspect terminal-time output evidence before scientific acceptance."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "submitted segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    allow_missing_restart_time_marker = getattr(
        args, "allow_missing_restart_time_marker", False
    )
    if allow_missing_restart_time_marker and not offline_local_root:
        raise ValueError(
            "--allow-missing-restart-time-marker is restricted to offline validation"
        )
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    reservation = reservation_for_manifest(read_reservations(paths), manifest_path)
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    authenticate_prepared_execution(
        manifest, manifest_path,
        allow_legacy_local=offline_local_root,
    )
    if manifest.get("state") not in {"submitted", "recorded"}:
        raise ValueError("only a submitted or recorded segment may be inspected")
    required_time = float(args.required_time)
    prepared_target = prepared_time_tlim_target(manifest)
    if (
        not math.isfinite(required_time)
        or abs(required_time - prepared_target) > 1.0e-12
    ):
        raise ValueError(
            "--required-time must match the prepared time/tlim target "
            f"{prepared_target:.12g}"
        )
    output_dir = Path(str(manifest["paths"]["output_dir"])).resolve()
    mhd_histories = sorted(output_dir.glob("*.mhd.hst"))
    user_histories = sorted(output_dir.glob("*.user.hst"))
    expected_ranks = (
        int(manifest["allocation"]["nodes"])
        * int(manifest["allocation"].get("ranks_per_node", 1))
    )
    snapshots = output_product_groups(
        output_dir / "bin", "*.bin", expected_ranks
    )
    restarts = output_product_groups(
        output_dir / "rst", "*.rst", expected_ranks
    )
    require_complete_output_inventory(output_dir / "bin", "*.bin", snapshots, "snapshot")
    require_complete_output_inventory(output_dir / "rst", "*.rst", restarts, "restart")
    if len(mhd_histories) != 1 or len(user_histories) != 1:
        raise ValueError("segment must retain exactly one MHD and one user history")
    history = parse_history(mhd_histories[0])
    user_history = parse_history(user_histories[0])
    if "time" not in history:
        raise ValueError("MHD history does not contain a time column")
    missing = [
        label for label in STRICT_LF_FAILURE_COLUMNS
        if label not in history
    ]
    if missing:
        raise ValueError(f"MHD history lacks strict safety counters: {missing}")
    final_time = history["time"][-1]
    maximum_failure_counts = {
        label: max(history[label])
        for label in STRICT_LF_FAILURE_COLUMNS
    }
    plasma_evidence = continuation_plasma_evidence(
        str(manifest["run"]["case_id"]), history, user_history
    )
    snapshot_times = [binary_product_time(group) for group in snapshots]
    if allow_missing_restart_time_marker or offline_local_root:
        restart_times = [
            restart_product_time(
                group, allow_missing_marker=allow_missing_restart_time_marker
            )
            for group in restarts
        ]
        restart_time_marker_modes: list[list[str]] = []
    else:
        authenticated_restarts = [
            authenticated_restart_product_time(group, manifest.get("command"))
            for group in restarts
        ]
        restart_times = [
            float(item["binary_time"]) for item in authenticated_restarts
        ]
        restart_time_marker_modes = [
            list(item["marker_modes"]) for item in authenticated_restarts
        ]
    restart_records = [retained_product(group) for group in restarts]
    restart_markers_verified = bool(restart_times) and all(
        value is not None for value in restart_times
    )
    terminal_restart_matches = [
        index for index, value in enumerate(restart_times)
        if value is not None and abs(value - final_time) <= 1.0e-10
    ]
    if restart_markers_verified:
        if len(terminal_restart_matches) != 1:
            raise ValueError(
                "segment must retain exactly one restart product whose explicit "
                "physical time matches the inspected final history time"
            )
        terminal_restart = restart_records[terminal_restart_matches[0]]
        terminal_restart_time = restart_times[terminal_restart_matches[0]]
    elif allow_missing_restart_time_marker and restart_records:
        terminal_restart = restart_records[-1]
        terminal_restart_time = None
    else:
        terminal_restart = None
        terminal_restart_time = None
    checks = {
        "required_time_reached": final_time >= required_time - 1.0e-10,
        "strict_lf_failure_counters_zero": all(
            value == 0.0 for value in maximum_failure_counts.values()
        ),
        "snapshots_retained": bool(snapshots),
        "terminal_snapshot_retained": bool(snapshot_times)
        and max(snapshot_times) >= final_time - 1.0e-10,
        "restart_retained": bool(restarts),
        "terminal_restart_physical_time_matches_final": (
            restart_markers_verified and len(terminal_restart_matches) == 1
        ) or (
            offline_local_root and allow_missing_restart_time_marker
            and bool(restart_records)
        ),
        "plasma_continuation_policy": True,
    }
    accepted = all(checks.values())
    clean_for_continuation = all(
        checks[key]
        for key in (
            "strict_lf_failure_counters_zero",
            "snapshots_retained",
            "terminal_snapshot_retained",
            "restart_retained",
            "terminal_restart_physical_time_matches_final",
            "plasma_continuation_policy",
        )
    )
    inspection: dict[str, object] = {
        "schema_version": 3 if offline_local_root else 4,
        "execution_epoch": EXECUTION_EPOCH,
        "inspected_utc": utc_now(),
        "manifest": str(manifest_path),
        "job_id": manifest.get("job_id"),
        "case_id": manifest["run"]["case_id"],
        "segment": manifest["run"]["segment"],
        "required_time": required_time,
        "final_time": final_time,
        "maximum_strict_failure_counts": maximum_failure_counts,
        "checks": checks,
        "accepted": accepted,
        "clean_for_continuation": clean_for_continuation,
        "mhd_history": retained_file(mhd_histories[0]),
        "user_history": retained_file(user_histories[0]),
        "plasma_continuation_policy": CONTINUATION_PLASMA_POLICY,
        "plasma_continuation_evidence": plasma_evidence,
        "snapshots": [retained_product(group) for group in snapshots],
        "snapshot_times": snapshot_times,
        "restarts": restart_records,
        "restart_times": restart_times,
        "restart_time_marker_modes": restart_time_marker_modes,
        "terminal_restart": terminal_restart,
        "terminal_restart_time": terminal_restart_time,
        "restart_time_marker_bypass": (
            offline_local_root
            and allow_missing_restart_time_marker
            and not restart_markers_verified
        ),
    }
    if "lf_hwproj" in history:
        inspection["final_hardwall_projection_count"] = history["lf_hwproj"][-1]
    inspection_path = manifest_path.parent / "segment_inspection.json"
    write_json(inspection_path, inspection)
    status = (
        "accepted"
        if accepted
        else "clean partial" if clean_for_continuation else "not accepted"
    )
    print(
        f"Segment inspection {status}: final_time={final_time:.12g}, "
        f"required_time={required_time:.12g}"
    )
    print(f"Wrote {inspection_path}")
    return 0 if accepted else 1


def sacct_output(args: argparse.Namespace, paths: dict[str, Path],
                 allow_fixture: bool) -> str:
    """Read or query a top-level Slurm allocation record."""

    if args.sacct_file:
        if not allow_fixture:
            raise ValueError("--sacct-file is allowed only for offline local roots")
        output = Path(args.sacct_file).read_text(encoding="utf-8")
    else:
        output = subprocess.run(
            [
                str(SACCT), "-X", "-j", args.job_id,
                "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                "ElapsedRaw,Submit,End", "-n", "-P",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    write_text(paths["accounting"] / f"{args.job_id}.stage_i.sacct.txt", output)
    return output


def cancellation_sacct_output(args: argparse.Namespace,
                              allow_fixture: bool) -> str:
    """Read no-start cancellation evidence without publishing it prematurely."""

    if args.sacct_file:
        if not allow_fixture:
            raise ValueError("--sacct-file is allowed only for offline local roots")
        return Path(args.sacct_file).read_text(encoding="utf-8")
    return live_cancellation_sacct_output(args.job_id)


def live_cancellation_sacct_output(job_id: str) -> str:
    """Query one exact job ID from trusted live sacct."""

    require_numeric_job_id(job_id)
    try:
        return subprocess.run(
            [
                str(SACCT), "-X", "-j", job_id,
                "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                "ElapsedRaw,Submit,End", "-n", "-P",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        raise ValueError(
            "sacct is unavailable; refusing submitted Stage I cancellation"
        ) from error


def parse_sacct(output: str, job_id: str) -> dict[str, str]:
    """Select a completed top-level allocation row."""

    require_numeric_job_id(job_id)
    records: list[list[str]] = []
    for row in csv.reader(output.splitlines(), delimiter="|"):
        if row and row[0] == job_id:
            while row and not row[-1]:
                row.pop()
            records.append(row)
    if len(records) != 1 or len(records[0]) != 8:
        raise ValueError(f"expected one eight-column allocation row for {job_id}")
    keys = (
        "job_id", "job_name", "state", "exit_code", "nodes",
        "elapsed_seconds", "submitted_utc", "completed_utc",
    )
    result = dict(zip(keys, records[0]))
    state = result["state"].split()[0].split("+")[0]
    if state in NONTERMINAL_STATES:
        raise ValueError(f"job {job_id} is not complete: {result['state']}")
    result["state"] = state
    return result


def require_recorded_scheduler_evidence(paths: dict[str, Path],
                                        row: dict[str, object],
                                        manifest: dict[str, object]) -> None:
    """Bind canonical recorded accounting to its archived Slurm allocation row."""

    job_id = str(row["job_id"])
    evidence = paths["accounting"] / f"{job_id}.stage_i.sacct.txt"
    if not evidence.is_file():
        raise ValueError(f"recorded transaction lacks scheduler evidence: {evidence}")
    sacct = parse_sacct(evidence.read_text(encoding="utf-8"), job_id)
    if sacct["job_name"] != expected_job_name(manifest):
        raise ValueError("recorded scheduler evidence job name differs from manifest")
    expected = {
        "job_id": row["job_id"],
        "state": row["state"],
        "exit_code": row["exit_code"],
        "nodes": row["nodes"],
        "elapsed_seconds": row["elapsed_seconds"],
        "submitted_utc": row["submitted_utc"],
        "completed_utc": row["completed_utc"],
    }
    for key, value in expected.items():
        if sacct[key] != str(value):
            raise ValueError(
                f"recorded scheduler evidence {key} differs from ledger"
            )


@locked_manifest_action
def record(args: argparse.Namespace) -> int:
    """Account a completed segment and release its reservation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "submitted segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    require_numeric_job_id(args.job_id)
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    authenticate_prepared_execution(
        manifest, manifest_path,
        allow_legacy_local=offline_local_root,
    )
    if manifest.get("state") != "submitted":
        raise ValueError("only a submitted segment can be accounted")
    if str(manifest.get("job_id")) != args.job_id:
        raise ValueError("job ID does not match the submitted segment")
    ledger = read_ledger(paths)
    if any(row["job_id"] == args.job_id for row in ledger):
        raise ValueError(f"job {args.job_id} is already accounted")
    sacct = parse_sacct(
        sacct_output(args, paths, allow_fixture=offline_local_root), args.job_id
    )
    if sacct["job_name"] != expected_job_name(manifest):
        raise ValueError(
            "sacct job name does not match the prepared segment: "
            f"{sacct['job_name']!r}"
        )
    nodes = int(sacct["nodes"])
    if nodes != int(manifest["allocation"]["nodes"]):
        raise ValueError("allocated nodes differ from the prepared reservation")
    inspection = None
    if args.result in {"accepted", "clean_partial"}:
        if sacct["state"] != "COMPLETED" or sacct["exit_code"] != "0:0":
            raise ValueError(
                "scientific continuation output requires a clean COMPLETED job"
            )
        inspection_path = manifest_path.parent / "segment_inspection.json"
        if not inspection_path.is_file():
            raise ValueError(
                "scientific continuation output requires inspect-segment evidence"
            )
        inspection = json.loads(inspection_path.read_text(encoding="utf-8"))
        if (
            not isinstance(inspection, dict)
            or inspection.get("execution_epoch") != EXECUTION_EPOCH
            or inspection.get("job_id") != args.job_id
            or Path(str(inspection.get("manifest", ""))).resolve()
            != manifest_path
        ):
            raise ValueError("segment inspection does not match this submitted job")
        try:
            inspected_target = float(inspection["required_time"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("segment inspection lacks its prepared target") from error
        if (
            not math.isfinite(inspected_target)
            or abs(inspected_target - prepared_time_tlim_target(manifest)) > 1.0e-12
        ):
            raise ValueError("segment inspection target differs from preparation")
        if args.result == "accepted" and inspection.get("accepted") is not True:
            raise ValueError("segment inspection does not accept this submitted job")
        if (
            args.result == "clean_partial"
            and (
                inspection.get("accepted") is True
                or inspection.get("clean_for_continuation") is not True
            )
        ):
            raise ValueError("clean_partial requires clean output short of its target")
        revalidate_inspection_files(inspection, manifest)
        if args.result == "clean_partial":
            revalidate_continuation_plasma_evidence(inspection, manifest)
    actual = node_hours(nodes, int(sacct["elapsed_seconds"]))
    cumulative = sum(float(row["actual_node_hours"]) for row in ledger) + actual
    remaining_reserved = sum(
        float(item["reserved_node_hours"])
        for item in active_reservations(reservations)
        if item is not reservation
    )
    if cumulative > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("actual use exceeds the Stage I reservation")
    if cumulative > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError("actual use exceeds the incremental project ceiling")
    if cumulative + remaining_reserved > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError(
            "accounted use plus active Stage I reservations exceeds the ceiling"
        )
    if cumulative + remaining_reserved > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError(
            "accounted use plus active Stage I reservations exceeds the "
            "incremental project ceiling"
        )
    command = manifest["command"]
    run = manifest["run"]
    allocation = manifest["allocation"]
    row = {
        "execution_epoch": EXECUTION_EPOCH,
        "job_id": args.job_id,
        "submitted_utc": sacct["submitted_utc"],
        "completed_utc": sacct["completed_utc"],
        "case_id": run["case_id"],
        "case_name": run["case_name"],
        "segment": run["segment"],
        "state": sacct["state"],
        "exit_code": sacct["exit_code"],
        "nodes": str(nodes),
        "requested_walltime": allocation["requested_walltime"],
        "elapsed_seconds": sacct["elapsed_seconds"],
        "reserved_node_hours": f"{float(allocation['reserved_node_hours']):.6f}",
        "actual_node_hours": f"{actual:.6f}",
        "cumulative_stage_i_node_hours": f"{cumulative:.6f}",
        "executable_revision": command["executable_revision"],
        "executable_sha256": command["executable_sha256"],
        "input_revision": command["input_revision"],
        "input_file": command["input_file"],
        "output_dir": manifest["paths"]["output_dir"],
        "result": args.result,
        "notes": args.notes,
    }
    reservation["state"] = "recorded"
    reservation["actual_node_hours"] = actual
    reservation["result"] = args.result
    manifest["state"] = "recorded"
    manifest["accounting"] = row
    if inspection is not None:
        manifest["scientific_inspection"] = inspection
    durable_transition(
        paths, "recorded", manifest_path, manifest, reservations, ledger_row=row
    )
    print(
        f"Recorded {actual:.6f} node-hours for {run['case_id']}/{run['segment']}; "
        f"Stage I cumulative={cumulative:.6f}."
    )
    return 0


def accepted_case_segments(paths: dict[str, Path],
                           case_id: str) -> list[dict[str, object]]:
    """Load scientifically qualifying, accounted segments for one mapped case."""

    segments: list[dict[str, object]] = []
    manifests = sorted(
        (paths["runs"] / case_id).glob("*/manifest/prepared_run.json")
    )
    allow_legacy_local = paths["root"].resolve() != DEFAULT_ROOT.resolve()
    reservations = read_reservations(paths)
    for manifest_path in manifests:
        manifest = read_manifest(manifest_path)
        require_current_epoch(manifest, "retained segment")
        accounting = manifest.get("accounting", {})
        if manifest.get("state") != "recorded" or not isinstance(accounting, dict):
            continue
        result = accounting.get("result")
        if result not in {"accepted", "clean_partial"}:
            continue
        inspection = manifest.get("scientific_inspection")
        if (
            not isinstance(inspection, dict)
            or (
                result == "accepted"
                and inspection.get("accepted") is not True
            )
            or (
                result == "clean_partial"
                and inspection.get("clean_for_continuation") is not True
            )
        ):
            raise ValueError(
                f"retained segment lacks qualifying inspection: {manifest_path}"
            )
        authenticate_prepared_execution(
            manifest, manifest_path, allow_legacy_local=allow_legacy_local
        )
        require_reserved_execution_intent(
            reservation_for_manifest(reservations, manifest_path),
            manifest,
            allow_legacy_local=allow_legacy_local,
        )
        if not allow_legacy_local or "mhd_history" in inspection:
            revalidate_inspection_files(inspection, manifest)
        if result == "clean_partial":
            manifest["scientific_inspection"] = (
                revalidate_continuation_plasma_evidence(inspection, manifest)
            )
        manifest["_manifest_path"] = str(manifest_path)
        segments.append(manifest)
    segments.sort(
        key=lambda item: float(item["scientific_inspection"]["final_time"])
    )
    return segments


def accepted_case_lineage(paths: dict[str, Path],
                          case_id: str) -> list[dict[str, object]]:
    """Select the recorded restart lineage ending at the latest accepted segment."""

    segments = accepted_case_segments(paths, case_id)
    if not segments:
        return []
    accepted = [
        segment for segment in segments
        if segment["accounting"]["result"] == "accepted"
    ]
    if not accepted:
        raise ValueError(f"{case_id} has no accepted terminal segment")
    final_time = max(
        float(segment["scientific_inspection"]["final_time"])
        for segment in accepted
    )
    terminals = [
        segment for segment in accepted
        if abs(
            float(segment["scientific_inspection"]["final_time"]) - final_time
        ) <= 1.0e-10
    ]
    if len(terminals) != 1:
        raise ValueError(
            f"{case_id} has multiple latest accepted terminal segments"
        )
    indexed = {
        Path(str(segment["_manifest_path"])).resolve(): segment
        for segment in segments
    }
    lineage: list[dict[str, object]] = []
    current = terminals[0]
    visited: set[Path] = set()
    while True:
        path = Path(str(current["_manifest_path"])).resolve()
        if path in visited:
            raise ValueError(f"{case_id} accepted restart lineage contains a cycle")
        visited.add(path)
        lineage.append(current)
        parent = current["command"].get("parent_segment")
        if parent is None:
            break
        if not isinstance(parent, dict) or "manifest" not in parent:
            raise ValueError(
                f"{case_id} accepted restart lineage lacks a parent manifest"
            )
        parent_path = Path(str(parent["manifest"])).resolve()
        if parent_path not in indexed:
            raise ValueError(
                f"{case_id} accepted restart lineage parent is not a "
                f"qualifying recorded segment: {parent_path}"
            )
        current = indexed[parent_path]
    lineage.reverse()
    return lineage


def one_segment_output(manifest: dict[str, object], pattern: str) -> Path:
    """Select the unique retained output matching one segment product."""

    output_dir = Path(str(manifest["paths"]["output_dir"]))
    matches = sorted(output_dir.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"segment output must have one {pattern} product: {output_dir}"
        )
    return matches[0]


def inspected_history_path(manifest: dict[str, object], key: str,
                           pattern: str, allow_legacy_local: bool) -> Path:
    """Return one authenticated inspection-retained history path."""

    inspection = manifest["scientific_inspection"]
    record = inspection.get(key)
    if isinstance(record, dict):
        revalidate_retained_file(record)
        return Path(str(record["path"])).resolve()
    if allow_legacy_local:
        return one_segment_output(manifest, pattern)
    raise ValueError(f"retained segment inspection lacks {key}")


def inspected_snapshot_groups(manifest: dict[str, object],
                              allow_legacy_local: bool) -> list[list[Path]]:
    """Return only snapshot groups retained by formal inspection."""

    inspection = manifest["scientific_inspection"]
    records = inspection.get("snapshots")
    if isinstance(records, list):
        groups = recorded_product_groups(records, "snapshots")
        for record in records:
            revalidate_retained_product(record)
        return groups
    if allow_legacy_local:
        output_dir = Path(str(manifest["paths"]["output_dir"]))
        expected_ranks = (
            int(manifest["allocation"]["nodes"])
            * int(manifest["allocation"].get("ranks_per_node", 1))
        )
        return output_product_groups(output_dir / "bin", "*.bin", expected_ranks)
    raise ValueError("retained segment inspection lacks snapshots")


def link_distinct_snapshots(sources: list[list[Path]], destination: Path,
                            case_name: str) -> list[Path]:
    """Link one retained binary per physical time into an analysis bundle."""

    timed = sorted((binary_product_time(group), group) for group in sources)
    linked: list[Path] = []
    last_time = float("-inf")
    for time, group in timed:
        if time <= last_time + 1.0e-12:
            continue
        name = f"{case_name}.{len(linked):05d}.bin"
        if group[0].parent.name.startswith("rank_"):
            rank0_target = None
            for source in group:
                rank_dir = destination / source.parent.name
                mkdir_durable(rank_dir)
                target = rank_dir / name
                target.symlink_to(source)
                fsync_directory(target.parent)
                if source.parent.name == "rank_00000000":
                    rank0_target = target
            if rank0_target is None:
                raise ValueError("rank-local snapshot set has no rank zero file")
            linked.append(rank0_target)
        else:
            target = destination / name
            target.symlink_to(group[0])
            fsync_directory(target.parent)
            linked.append(target)
        last_time = time
    if not linked:
        raise ValueError("accepted segments retain no distinct snapshots")
    return linked


@locked_root_action
def bundle_case(args: argparse.Namespace) -> int:
    """Assemble accepted restart segments as one analyzer-compatible bundle."""

    required_final_time = require_positive_finite_float(
        args.required_final_time, "--required-final-time"
    )
    root = require_root(Path(args.root), args.allow_local_root)
    if is_offline_local_root(root, args.allow_local_root):
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    matrix = validate_matrix(matrix_path, source_dir)
    case = case_for_id(matrix, args.case_id)
    segments = accepted_case_lineage(paths, args.case_id)
    if not segments:
        raise ValueError(f"no accepted segments are recorded for {args.case_id}")
    final_time = float(segments[-1]["scientific_inspection"]["final_time"])
    if final_time < required_final_time - 1.0e-10:
        raise ValueError(
            f"{args.case_id} reaches only t={final_time}; "
            f"required final time is {required_final_time}"
        )
    first_command = segments[0]["command"]
    expected_digests = (
        first_command["input_sha256"],
        first_command["executable_sha256"],
    )
    submitted_input = Path(str(first_command["input_file"]))
    allow_legacy_local = root.resolve() != DEFAULT_ROOT.resolve()
    if not allow_legacy_local:
        require_file_sha256(
            submitted_input, first_command["input_sha256"], "accepted input"
        )
    model_choices = analysis_model_choices(submitted_input)
    if model_choices.get("output1_file_type") != "hst":
        raise ValueError("accepted production input does not retain output1 history")
    try:
        history_interval = float(model_choices["output1_dt"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "accepted production input has no numeric history cadence"
        ) from error
    if not math.isfinite(history_interval) or history_interval <= 0.0:
        raise ValueError("accepted production input has invalid history cadence")
    previous_final = None
    previous_history: dict[str, list[float]] | None = None
    for segment in segments:
        command = segment["command"]
        if (
            command["input_sha256"],
            command["executable_sha256"],
        ) != expected_digests:
            raise ValueError("accepted segments do not share input/executable digests")
        history = parse_history(
            inspected_history_path(
                segment, "mhd_history", "*.mhd.hst", allow_legacy_local
            )
        )
        first_time = history["time"][0]
        segment_final = history["time"][-1]
        if previous_final is None:
            if first_time > 1.0e-10:
                raise ValueError("accepted case history does not begin at t=0")
        else:
            if (
                previous_history is None
                or "dt" not in previous_history
                or "dt" not in history
            ):
                raise ValueError(
                    "accepted histories lack timesteps needed to validate restarts"
                )
            boundary_step = max(previous_history["dt"][-1], history["dt"][0])
            if not math.isfinite(boundary_step) or boundary_step < 0.0:
                raise ValueError("accepted histories have invalid boundary timestep")
            # A walltime terminal row can replace one scheduled output at restart.
            permitted_gap = 2.0 * history_interval + boundary_step + 1.0e-10
            if first_time > previous_final + permitted_gap:
                raise ValueError(
                    "accepted restart histories exceed configured sampling cadence"
                )
            if segment_final <= previous_final + 1.0e-10:
                raise ValueError("accepted restart segment does not advance time")
        previous_final = segment_final
        previous_history = history
    bundle = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else paths["runs"] / "bundles" / args.case_id
    )
    require_beneath_root(bundle, root, "analysis bundle", args.allow_local_root)
    if bundle.exists():
        if not args.replace:
            raise ValueError(f"analysis bundle already exists: {bundle}")
        shutil.rmtree(bundle)
    history_dir = bundle / "history"
    input_dir = bundle / "inputs"
    snapshot_dir = bundle / "cases" / str(case["name"]) / "bin"
    for directory in (history_dir, input_dir, snapshot_dir):
        mkdir_durable(directory)
    mhd_history = history_dir / f"{case['name']}.mhd.hst"
    user_history = history_dir / f"{case['name']}.user.hst"
    merge_history_files(
        [
            inspected_history_path(
                segment, "mhd_history", "*.mhd.hst", allow_legacy_local
            )
            for segment in segments
        ],
        mhd_history,
    )
    merge_history_files(
        [
            inspected_history_path(
                segment, "user_history", "*.user.hst", allow_legacy_local
            )
            for segment in segments
        ],
        user_history,
    )
    snapshot_sources: list[list[Path]] = []
    for segment in segments:
        snapshot_sources.extend(inspected_snapshot_groups(segment, allow_legacy_local))
    snapshots = link_distinct_snapshots(
        snapshot_sources, snapshot_dir, str(case["name"])
    )
    archived_input = input_dir / submitted_input.name
    copy_file(submitted_input, archived_input)
    case_entry = {
        "name": case["name"],
        "input": case["input"],
        "execution_input": str(archived_input.relative_to(bundle)),
        "command": ["archived production segments"],
        "overrides": [],
        "log": "retained in production segment archive",
        "status": "passed",
        "outputs": {
            "mhd_history": str(mhd_history.relative_to(bundle)),
            "user_history": str(user_history.relative_to(bundle)),
            "snapshot_paths": [
                str(path.relative_to(bundle)) for path in snapshots
            ],
        },
        "lf_active": True,
        "amr": False,
        "paper_smoke": False,
        "model_choices": model_choices,
    }
    manifest = {
        "workflow": "paper-mks24-stage-i-production",
        "execution_epoch": EXECUTION_EPOCH,
        "created_utc": utc_now(),
        "status": "accepted_for_analysis",
        "git_revision": first_command["input_revision"],
        "git_worktree_dirty": (
            "not recorded; submitted input and matrix verified committed"
        ),
        "executable": first_command["executable"],
        "production_case_id": args.case_id,
        "required_final_time": required_final_time,
        "accepted_final_time": final_time,
        "production_segment_manifests": [
            str(segment["_manifest_path"]) for segment in segments
        ],
        "cases": [case_entry],
    }
    write_json(bundle / "manifest.json", manifest)
    print(f"Wrote accepted Stage I analysis bundle: {bundle}")
    print(f"  final_time={final_time:.12g}, snapshots={len(snapshots)}")
    return 0


def prefix_bundle_case_paths(case: dict[str, object], prefix: Path
                             ) -> dict[str, object]:
    """Rewrite nested per-case paths for a campaign-level bundle manifest."""

    copied = json.loads(json.dumps(case))
    copied["execution_input"] = str(prefix / str(copied["execution_input"]))
    outputs = copied["outputs"]
    for key in ("mhd_history", "user_history"):
        outputs[key] = str(prefix / str(outputs[key]))
    outputs["snapshot_paths"] = [
        str(prefix / str(path)) for path in outputs["snapshot_paths"]
    ]
    return copied


@locked_root_action
def bundle_campaign(args: argparse.Namespace) -> int:
    """Assemble all accepted mapped cases into one paper-analysis bundle."""

    required_final_time = require_positive_finite_float(
        args.required_final_time, "--required-final-time"
    )
    root = require_root(Path(args.root), args.allow_local_root)
    if is_offline_local_root(root, args.allow_local_root):
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    matrix = validate_matrix(matrix_path, source_dir)
    for case in matrix["cases"]:
        segments = accepted_case_lineage(paths, str(case["id"]))
        if not segments:
            raise ValueError(f"no accepted segments are recorded for {case['id']}")
        final_time = float(segments[-1]["scientific_inspection"]["final_time"])
        if final_time < required_final_time - 1.0e-10:
            raise ValueError(
                f"{case['id']} reaches only t={final_time}; "
                f"required final time is {required_final_time}"
            )
    bundle = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else paths["runs"] / "bundles" / "mks24-stage-i-campaign"
    )
    require_beneath_root(bundle, root, "analysis bundle", args.allow_local_root)
    if bundle.exists():
        if not args.replace:
            raise ValueError(f"analysis bundle already exists: {bundle}")
        shutil.rmtree(bundle)
    mkdir_durable(bundle)
    cases: list[dict[str, object]] = []
    segment_manifests: list[str] = []
    case_times: dict[str, float] = {}
    for case in matrix["cases"]:
        case_id = str(case["id"])
        nested = bundle / "case_bundles" / case_id
        case_args = argparse.Namespace(**vars(args))
        case_args.case_id = case_id
        case_args.output_dir = str(nested)
        case_args.replace = False
        bundle_case(case_args)
        nested_manifest = read_manifest(nested / "manifest.json")
        prefix = nested.relative_to(bundle)
        cases.append(prefix_bundle_case_paths(nested_manifest["cases"][0], prefix))
        segment_manifests.extend(nested_manifest["production_segment_manifests"])
        case_times[case_id] = float(nested_manifest["accepted_final_time"])
    campaign_manifest = {
        "workflow": "paper-mks24-stage-i-production",
        "execution_epoch": EXECUTION_EPOCH,
        "created_utc": utc_now(),
        "status": "accepted_for_analysis",
        "git_revision": subprocess.run(
            [
                authenticated_git_binary(), "--no-replace-objects",
                *controller_git_worktree_prefix(source_dir),
                "rev-parse", "HEAD",
            ],
            check=True, capture_output=True, text=True,
            env=hardened_git_environment(),
        ).stdout.strip(),
        "git_worktree_dirty": (
            "not recorded; submitted inputs and matrices verified committed"
        ),
        "executable": "recorded per accepted production segment",
        "required_final_time": required_final_time,
        "accepted_case_final_times": case_times,
        "production_segment_manifests": segment_manifests,
        "cases": cases,
    }
    write_json(bundle / "manifest.json", campaign_manifest)
    print(f"Wrote accepted Stage I campaign bundle: {bundle}")
    print(f"  accepted cases={len(cases)}")
    return 0


def require_accounting_evidence_path(paths: dict[str, Path], value: object,
                                     label: str) -> Path:
    """Require retained cancellation evidence to remain beneath accounting."""

    declared = Path(str(value)).expanduser()
    if not declared.is_absolute():
        raise ValueError(f"{label} path must be absolute")
    resolved = declared.resolve()
    try:
        resolved.relative_to(paths["accounting"].resolve())
    except ValueError as error:
        raise ValueError(f"{label} must remain beneath Stage I accounting") from error
    if declared != resolved:
        raise ValueError(f"{label} path must resolve exactly")
    return resolved


def submitted_cancellation_evidence_paths(
    paths: dict[str, Path], job_id: str
) -> dict[str, Path]:
    """Return the fixed evidence namespace for one no-start cancellation."""

    require_numeric_job_id(job_id)
    prefix = paths["accounting"] / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_{job_id}_"
        "source_bundle_recovery_cancellation"
    )
    return {
        "authorization": prefix.with_suffix(".authorization.json"),
        "publication_audit": prefix.with_suffix(".authorization.json.publication_audit.json"),
        "submitted_manifest": prefix.with_suffix(".submitted_manifest.json"),
        "scheduler_batch_script": prefix.with_suffix(".scheduler_batch_script.sbatch"),
        "pre_cancel_hold": prefix.with_suffix(".pre_cancel_hold.squeue.txt"),
        "post_cancel_sacct": prefix.with_suffix(".post_cancel.sacct.txt"),
        "post_cancel_queue": prefix.with_suffix(".post_cancel.squeue.txt"),
    }


def require_immutable_recovery_evidence(
    path: Path, expected_sha256: object, label: str
) -> None:
    """Require one owner-controlled, immutable, single-link evidence file."""

    try:
        profile = path.stat(follow_symlinks=False)
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or stat.S_IMODE(profile.st_mode) != 0o444
    ):
        raise ValueError(f"{label} must be an owner-controlled immutable file: {path}")
    require_file_sha256(path, expected_sha256, label)


def require_canonical_source_bundle_recovery_identity(
    paths: dict[str, Path],
    manifest_path: Path,
    submitted_manifest: dict[str, object],
) -> None:
    """Pin the sole canonical corruption-recovery cancellation identity."""

    if paths["root"].resolve() != DEFAULT_ROOT.resolve():
        return
    policy = CANONICAL_SOURCE_BUNDLE_RECOVERY
    job_id = str(submitted_manifest.get("job_id", ""))
    if (
        job_id != policy["job_id"]
        or manifest_path
        != (paths["root"] / str(policy["manifest_relative"])).resolve()
    ):
        raise ValueError("submitted cancellation is not the promoted recovery identity")
    command = submitted_manifest.get("command")
    bundle = command.get("source_bundle") if isinstance(command, dict) else None
    if (
        not isinstance(bundle, dict)
        or Path(str(bundle.get("path", ""))).resolve()
        != (paths["root"] / str(policy["source_bundle_relative"])).resolve()
        or bundle.get("sha256") != policy["source_bundle_expected_sha256"]
    ):
        raise ValueError("submitted cancellation source-bundle identity is not promoted")


def validate_pre_cancel_hold_evidence(
    path: Path,
    job_id: str,
    expected_job_name_value: str,
) -> None:
    """Require one exact held, pending, never-started scheduler row."""

    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if len(lines) != 1:
        raise ValueError("pre-cancel hold evidence must contain one scheduler row")
    row = next(csv.reader(lines, delimiter="|"))
    if len(row) != 6:
        raise ValueError("pre-cancel hold evidence scheduler row is malformed")
    retained_job_id, job_name, state, elapsed, nodes, reason = row
    if (
        retained_job_id != job_id
        or job_name != expected_job_name_value
        or state != "PENDING"
        or elapsed != "0:00"
        or nodes != "1"
        or reason not in {"JobHeldUser", "(JobHeldUser)"}
    ):
        raise ValueError("pre-cancel hold evidence does not prove a held no-start job")


def validate_recovery_publication_audit(
    evidence: dict[str, Path],
    authorization_sha256: str,
    expected_audit_sha256: object,
) -> None:
    """Require an independently reviewed fixed-path recovery authorization."""

    audit_path = evidence["publication_audit"]
    require_immutable_recovery_evidence(
        audit_path,
        expected_audit_sha256,
        "submitted cancellation publication audit",
    )
    try:
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("submitted cancellation publication audit is invalid") from error
    required = {
        "schema_version",
        "record_type",
        "execution_epoch",
        "published_utc",
        "authorization",
        "review",
    }
    if not isinstance(audit, dict) or frozenset(audit) != required:
        raise ValueError("submitted cancellation publication audit schema is invalid")
    binding = audit["authorization"]
    review = audit["review"]
    if (
        audit["schema_version"] != 1
        or audit["record_type"]
        != "stage-i-source-bundle-recovery-cancellation-publication-audit"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or not isinstance(binding, dict)
        or frozenset(binding) != {"path", "sha256", "mode"}
        or Path(str(binding["path"])).resolve() != evidence["authorization"]
        or binding["sha256"] != authorization_sha256
        or binding["mode"] != "0444"
        or not isinstance(review, dict)
        or frozenset(review) != {"status", "authority", "reviewed_by"}
        or review["status"] != "approved"
        or review["authority"] != "independent-source-bundle-recovery-review"
        or not isinstance(review["reviewed_by"], str)
        or not review["reviewed_by"].strip()
    ):
        raise ValueError("submitted cancellation publication audit is not approved")
    parse_utc_timestamp(audit["published_utc"], "recovery authorization publication UTC")


def source_bundle_cancellation_authorization(
    paths: dict[str, Path],
    manifest_path: Path,
    submitted_manifest: dict[str, object],
) -> tuple[Path, str]:
    """Authenticate a one-use held-job cancellation caused by bundle corruption."""

    require_canonical_source_bundle_recovery_identity(
        paths, manifest_path, submitted_manifest
    )
    job_id = require_numeric_job_id(str(submitted_manifest.get("job_id", "")))
    evidence = submitted_cancellation_evidence_paths(paths, job_id)
    authorization = evidence["authorization"]
    expected_sha256: object = (
        sha256(authorization) if authorization.is_file() else None
    )
    if paths["root"].resolve() == DEFAULT_ROOT.resolve():
        expected_sha256 = CANONICAL_SOURCE_BUNDLE_RECOVERY["authorization_sha256"]
    require_immutable_recovery_evidence(
        authorization, expected_sha256, "submitted cancellation authorization"
    )
    try:
        retained = json.loads(authorization.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ValueError("submitted cancellation authorization is invalid JSON") from error
    required = {
        "schema_version",
        "record_type",
        "execution_epoch",
        "job_id",
        "decision",
        "reason",
        "manifest",
        "batch_script",
        "pre_cancel_hold",
        "source_bundle",
    }
    if not isinstance(retained, dict) or frozenset(retained) != required:
        raise ValueError("submitted cancellation authorization schema is invalid")
    if (
        retained["schema_version"] != 1
        or retained["record_type"]
        != "stage-i-source-bundle-recovery-cancellation-authorization"
        or retained["execution_epoch"] != EXECUTION_EPOCH
        or retained["job_id"] != submitted_manifest.get("job_id")
        or retained["decision"] != "cancel-held-job-without-start"
        or not isinstance(retained["reason"], str)
        or not retained["reason"].strip()
    ):
        raise ValueError("submitted cancellation authorization identity is invalid")
    manifest_evidence = retained["manifest"]
    if (
        not isinstance(manifest_evidence, dict)
        or frozenset(manifest_evidence)
        != {"live_path", "snapshot_path", "sha256"}
        or Path(str(manifest_evidence["live_path"])).resolve() != manifest_path
        or Path(str(manifest_evidence["snapshot_path"])).resolve()
        != evidence["submitted_manifest"]
    ):
        raise ValueError("submitted cancellation authorization manifest differs")
    snapshot_path = evidence["submitted_manifest"]
    require_immutable_recovery_evidence(
        snapshot_path, manifest_evidence["sha256"], "submitted manifest snapshot"
    )
    if read_manifest(snapshot_path) != submitted_manifest:
        raise ValueError("submitted manifest snapshot differs from launch manifest")
    manifest_paths = submitted_manifest.get("paths")
    command = submitted_manifest.get("command")
    if not isinstance(manifest_paths, dict) or not isinstance(command, dict):
        raise ValueError("submitted manifest lacks retained launch metadata")
    batch_script = Path(str(manifest_paths.get("batch_script", ""))).resolve()
    batch_evidence = retained["batch_script"]
    if (
        not isinstance(batch_evidence, dict)
        or frozenset(batch_evidence)
        != {
            "path",
            "sha256",
            "normalized_sha256",
            "scheduler_snapshot_path",
            "scheduler_snapshot_sha256",
        }
        or Path(str(batch_evidence["path"])).resolve() != batch_script
        or Path(str(batch_evidence["scheduler_snapshot_path"])).resolve()
        != evidence["scheduler_batch_script"]
    ):
        raise ValueError("submitted cancellation authorization batch script differs")
    require_file_sha256(
        batch_script, batch_evidence["sha256"], "submitted cancellation batch script"
    )
    if (
        batch_evidence["normalized_sha256"] != command.get("batch_script_sha256")
        or normalized_batch_script_sha256(batch_script)
        != batch_evidence["normalized_sha256"]
    ):
        raise ValueError("submitted cancellation batch script launch digest differs")
    scheduler_script = evidence["scheduler_batch_script"]
    require_immutable_recovery_evidence(
        scheduler_script,
        batch_evidence["scheduler_snapshot_sha256"],
        "scheduler-retained submitted batch script",
    )
    if (
        batch_evidence["scheduler_snapshot_sha256"] != batch_evidence["sha256"]
        or normalized_batch_script_sha256(scheduler_script)
        != batch_evidence["normalized_sha256"]
    ):
        raise ValueError("scheduler-retained submitted batch script differs")
    hold = retained["pre_cancel_hold"]
    if (
        not isinstance(hold, dict)
        or frozenset(hold) != {"path", "sha256"}
        or Path(str(hold["path"])).resolve() != evidence["pre_cancel_hold"]
    ):
        raise ValueError("pre-cancel hold evidence binding differs")
    require_immutable_recovery_evidence(
        evidence["pre_cancel_hold"], hold["sha256"], "pre-cancel hold evidence"
    )
    validate_pre_cancel_hold_evidence(
        evidence["pre_cancel_hold"], job_id, expected_job_name(submitted_manifest)
    )
    if not isinstance(command.get("source_bundle"), dict):
        raise ValueError("submitted manifest lacks source-bundle provenance")
    bundle = command["source_bundle"]
    bundle_evidence = retained["source_bundle"]
    if (
        not isinstance(bundle_evidence, dict)
        or frozenset(bundle_evidence)
        != {"path", "expected_sha256", "observed_sha256", "incident"}
        or Path(str(bundle_evidence["path"])).resolve()
        != Path(str(bundle.get("path", ""))).resolve()
        or bundle_evidence["expected_sha256"] != bundle.get("sha256")
        or bundle_evidence["observed_sha256"] == bundle_evidence["expected_sha256"]
    ):
        raise ValueError("submitted cancellation source-bundle evidence differs")
    require_file_sha256(
        Path(str(bundle_evidence["path"])).resolve(),
        bundle_evidence["observed_sha256"],
        "observed corrupt source bundle",
    )
    incident = bundle_evidence["incident"]
    if not isinstance(incident, dict) or frozenset(incident) != {"path", "sha256"}:
        raise ValueError("submitted cancellation incident binding is invalid")
    incident_path = require_accounting_evidence_path(
        paths, incident["path"], "source-bundle corruption incident"
    )
    if paths["root"].resolve() == DEFAULT_ROOT.resolve():
        policy = CANONICAL_SOURCE_BUNDLE_RECOVERY
        if (
            incident_path
            != (paths["root"] / str(policy["incident_relative"])).resolve()
            or incident["sha256"] != policy["incident_sha256"]
            or bundle_evidence["expected_sha256"]
            != policy["source_bundle_expected_sha256"]
            or bundle_evidence["observed_sha256"]
            != policy["source_bundle_observed_sha256"]
        ):
            raise ValueError("canonical source-bundle recovery evidence is not promoted")
    require_immutable_recovery_evidence(
        incident_path, incident["sha256"], "source-bundle corruption incident"
    )
    incident_value = json.loads(incident_path.read_text(encoding="utf-8"))
    if (
        not isinstance(incident_value, dict)
        or incident_value.get("schema_version") != 1
        or incident_value.get("record_type")
        != "stage-i-source-bundle-corruption-incident"
        or incident_value.get("execution_epoch") != EXECUTION_EPOCH
        or incident_value.get("expected", {}).get("sha256")
        != bundle_evidence["expected_sha256"]
        or incident_value.get("observed", {}).get("path")
        != str(Path(str(bundle_evidence["path"])).resolve())
        or incident_value.get("observed", {}).get("sha256")
        != bundle_evidence["observed_sha256"]
        or incident_value.get("scheduler_containment", {}).get("job_id") != job_id
        or incident_value.get("scheduler_containment", {}).get("reason")
        != "JobHeldUser"
        or incident_value.get("scheduler_containment", {}).get("released") is not False
    ):
        raise ValueError("source-bundle corruption incident schema or binding differs")
    expected_audit_sha256: object = (
        sha256(evidence["publication_audit"])
        if evidence["publication_audit"].is_file() else None
    )
    if paths["root"].resolve() == DEFAULT_ROOT.resolve():
        expected_audit_sha256 = CANONICAL_SOURCE_BUNDLE_RECOVERY[
            "publication_audit_sha256"
        ]
    validate_recovery_publication_audit(
        evidence, str(expected_sha256), expected_audit_sha256
    )
    return authorization, str(expected_sha256)


def submitted_manifest_before_cancellation(
    cancelled_manifest: dict[str, object],
) -> dict[str, object]:
    """Reconstruct the immutable submitted-state manifest for replay checks."""

    submitted = json.loads(json.dumps(cancelled_manifest))
    if submitted.get("state") != "cancelled":
        raise ValueError("submitted cancellation replay requires cancelled state")
    submitted["state"] = "submitted"
    submitted.pop("cancellation", None)
    submitted.pop("cancellation_notes", None)
    return submitted


def validate_submitted_cancellation_metadata(
    paths: dict[str, Path], manifest_path: Path, manifest: dict[str, object]
) -> None:
    """Authenticate retained no-start scheduler cancellation evidence."""

    cancellation = manifest.get("cancellation")
    required = {
        "cancelled_utc",
        "mode",
        "notes",
        "job_id",
        "scheduler",
        "scheduler_evidence",
        "queue_absence_evidence",
        "recovery_authorization",
    }
    if not isinstance(cancellation, dict) or frozenset(cancellation) != required:
        raise ValueError("submitted cancellation metadata schema is invalid")
    parse_utc_timestamp(cancellation["cancelled_utc"], "submitted cancellation UTC")
    if (
        cancellation["mode"] != "authenticated-submitted-no-start"
        or cancellation["job_id"] != manifest.get("job_id")
        or not isinstance(cancellation["notes"], str)
    ):
        raise ValueError("submitted cancellation metadata identity is invalid")
    scheduler = cancellation["scheduler"]
    if (
        not isinstance(scheduler, dict)
        or scheduler.get("job_id") != manifest.get("job_id")
        or scheduler.get("job_name") != expected_job_name(manifest)
        or scheduler.get("state") != "CANCELLED"
        or scheduler.get("exit_code") != "0:0"
        or scheduler.get("elapsed_seconds") != "0"
        or scheduler.get("nodes") != "0"
    ):
        raise ValueError("submitted cancellation scheduler evidence is invalid")
    job_id = require_numeric_job_id(str(manifest["job_id"]))
    fixed = submitted_cancellation_evidence_paths(paths, job_id)
    bindings = {
        "scheduler_evidence": (
            "post_cancel_sacct", "submitted cancellation scheduler evidence"
        ),
        "queue_absence_evidence": (
            "post_cancel_queue", "submitted cancellation queue-absence evidence"
        ),
        "recovery_authorization": (
            "authorization", "submitted cancellation recovery authorization"
        ),
    }
    retained_paths = {}
    for key, (fixed_key, label) in bindings.items():
        binding = cancellation[key]
        if not isinstance(binding, dict) or frozenset(binding) != {"path", "sha256"}:
            raise ValueError(f"{label} binding is invalid")
        evidence_path = require_accounting_evidence_path(paths, binding["path"], label)
        if evidence_path != fixed[fixed_key]:
            raise ValueError(f"{label} path differs from the fixed recovery namespace")
        require_immutable_recovery_evidence(evidence_path, binding["sha256"], label)
        retained_paths[key] = evidence_path
    source_bundle_cancellation_authorization(
        paths,
        manifest_path,
        submitted_manifest_before_cancellation(manifest),
    )
    parsed = parse_sacct(
        retained_paths["scheduler_evidence"].read_text(encoding="utf-8"),
        str(manifest["job_id"]),
    )
    if parsed != scheduler:
        raise ValueError("submitted cancellation scheduler evidence changed")
    if retained_paths["queue_absence_evidence"].read_text(encoding="utf-8").strip():
        raise ValueError("submitted cancellation queue-absence evidence is nonempty")


@locked_manifest_action
def cancel_submitted(args: argparse.Namespace) -> int:
    """Close a terminally cancelled submitted reservation that never started."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "submitted segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = initialize(root) if offline_local_root else layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    if manifest.get("state") != "submitted":
        raise ValueError("only a submitted segment may use submitted cancellation")
    require_numeric_job_id(args.job_id)
    if manifest.get("job_id") != args.job_id:
        raise ValueError("job ID does not match the submitted segment")
    reservations = read_reservations(paths)
    require_active_reservation_policy(reservations)
    require_reservation_budget(
        read_ledger(paths), reservations, "submitted cancellation baseline"
    )
    require_retained_r17_policy(paths, reservations)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if (
        reservation.get("state") != "submitted"
        or reservation.get("job_id") != args.job_id
    ):
        raise ValueError("reservation does not match the submitted segment")
    require_reservation_matches_manifest(reservation, manifest)
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    authorization, authorization_sha256 = source_bundle_cancellation_authorization(
        paths, manifest_path, manifest
    )
    evidence = submitted_cancellation_evidence_paths(paths, args.job_id)
    for key, label in (
        ("post_cancel_sacct", "submitted cancellation scheduler evidence"),
        ("post_cancel_queue", "submitted cancellation queue-absence evidence"),
    ):
        require_immutable_recovery_evidence(
            evidence[key],
            sha256(evidence[key]) if evidence[key].is_file() else None,
            label,
        )
    scheduler_output = cancellation_sacct_output(args, offline_local_root)
    if scheduler_output != evidence["post_cancel_sacct"].read_text(encoding="utf-8"):
        raise ValueError("live scheduler cancellation evidence differs from retained bytes")
    scheduler = parse_sacct(scheduler_output, args.job_id)
    if (
        scheduler["job_name"] != expected_job_name(manifest)
        or scheduler["state"] != "CANCELLED"
        or scheduler["exit_code"] != "0:0"
        or scheduler["elapsed_seconds"] != "0"
        or scheduler["nodes"] != "0"
    ):
        raise ValueError(
            "submitted cancellation requires a terminal CANCELLED no-start job"
        )
    queue_output = cancelled_job_queue_output(args, offline_local_root)
    if queue_output.strip():
        raise ValueError("cancelled submitted job remains present in squeue")
    if queue_output != evidence["post_cancel_queue"].read_text(encoding="utf-8"):
        raise ValueError("live queue-absence evidence differs from retained bytes")
    reservation["state"] = "cancelled"
    reservation["notes"] = args.notes
    manifest["state"] = "cancelled"
    manifest["cancellation_notes"] = args.notes
    manifest["cancellation"] = {
        "cancelled_utc": utc_now(),
        "mode": "authenticated-submitted-no-start",
        "notes": args.notes,
        "job_id": args.job_id,
        "scheduler": scheduler,
        "scheduler_evidence": {
            "path": str(evidence["post_cancel_sacct"]),
            "sha256": sha256(evidence["post_cancel_sacct"]),
        },
        "queue_absence_evidence": {
            "path": str(evidence["post_cancel_queue"]),
            "sha256": sha256(evidence["post_cancel_queue"]),
        },
        "recovery_authorization": {
            "path": str(authorization),
            "sha256": authorization_sha256,
        },
    }
    validate_submitted_cancellation_metadata(paths, manifest_path, manifest)
    durable_transition(
        paths, "cancelled_submitted", manifest_path, manifest, reservations
    )
    print(f"Cancelled never-started submitted Stage I reservation: {manifest_path}")
    return 0


@locked_manifest_action
def cancel(args: argparse.Namespace) -> int:
    """Release a segment that was prepared but never submitted."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    if offline_local_root:
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
    if manifest.get("state") != "prepared":
        raise ValueError("only an unsubmitted segment may be cancelled")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation has progressed beyond preparation")
    break_glass = getattr(args, "break_glass_cancel_evidence", "")
    if getattr(args, "confirm_break_glass_cancel", False):
        if not str(break_glass).strip():
            raise ValueError("--break-glass-cancel-evidence must be nonempty")
        cancellation_mode = "break-glass"
    else:
        if break_glass:
            raise ValueError(
                "--confirm-break-glass-cancel is required with break-glass evidence"
            )
        require_reconciled_store_consistency(paths)
        require_reserved_execution_intent(
            reservation, manifest, allow_legacy_local=offline_local_root
        )
        authenticate_prepared_execution(
            manifest, manifest_path, allow_legacy_local=offline_local_root
        )
        cancellation_mode = "authenticated"
    reservation["state"] = "cancelled"
    reservation["notes"] = args.notes
    manifest["state"] = "cancelled"
    manifest["cancellation_notes"] = args.notes
    manifest["cancellation"] = {
        "cancelled_utc": utc_now(),
        "mode": cancellation_mode,
        "notes": args.notes,
        "break_glass_evidence": break_glass or None,
    }
    durable_transition(paths, "cancelled", manifest_path, manifest, reservations)
    print(f"Cancelled Stage I reservation: {manifest_path}")
    return 0


def reconcile_report(root: Path) -> dict[str, object]:
    """Read retained stores and report consistency without modifying them."""

    paths = layout(root)
    issues: list[str] = []
    offline_local_root = root.resolve() != DEFAULT_ROOT.expanduser().resolve()
    qualification = qualification_approval_status(paths)
    if (
        qualification["state"] == "invalid"
        or (
            not offline_local_root
            and qualification["state"] != "approved"
        )
    ):
        issues.append(
            "E03 corrected-build Frontier qualification is "
            f"{qualification['state']}: {qualification.get('reason', qualification['path'])}"
        )
    if not paths["transactions"].is_dir():
        issues.append(f"transaction store is missing: {paths['transactions']}")
    transactions = pending_transaction_paths(paths)
    for transaction_path in transactions:
        try:
            transaction = read_transaction(paths, transaction_path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            issues.append(f"cannot read transaction {transaction_path}: {error}")
        else:
            issues.append(
                f"pending {transaction.get('kind')} transaction requires recovery: "
                f"{transaction_path}"
            )
    if paths["ledger"].is_file():
        try:
            ledger = read_ledger(paths)
        except (OSError, ValueError) as error:
            ledger = []
            issues.append(f"cannot read ledger: {error}")
    else:
        ledger = []
        issues.append(f"ledger is missing: {paths['ledger']}")
    if paths["reservations"].is_file():
        try:
            reservations = read_reservations(paths)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            reservations = []
            issues.append(f"cannot read reservations: {error}")
    else:
        reservations = []
        issues.append(f"reservation store is missing: {paths['reservations']}")

    manifests: dict[Path, dict[str, object]] = {}
    if paths["runs"].is_dir():
        for manifest_path in sorted(
            paths["runs"].glob("*/*/manifest/prepared_run.json")
        ):
            resolved = manifest_path.resolve()
            try:
                manifests[resolved] = read_manifest(resolved)
            except (OSError, ValueError, json.JSONDecodeError) as error:
                issues.append(f"cannot read manifest {resolved}: {error}")
    else:
        issues.append(f"run store is missing: {paths['runs']}")
    for run_dir in orphaned_segment_run_directories(paths):
        issues.append(f"orphaned segment run directory: {run_dir}")

    active = [
        reservation for reservation in reservations
        if isinstance(reservation, dict)
        and reservation.get("state") in {"prepared", "submitted"}
    ]
    try:
        require_active_reservation_policy(reservations)
    except ValueError as error:
        issues.append(f"active reservation policy is invalid: {error}")
    try:
        require_retained_r17_policy(paths, reservations)
    except ValueError as error:
        issues.append(f"retained R17 policy is invalid: {error}")
    try:
        require_reservation_budget(ledger, reservations, "active Stage I reservation")
    except ValueError as error:
        issues.append(str(error))

    reservations_by_manifest: dict[Path, list[dict[str, object]]] = {}
    for reservation in reservations:
        if not isinstance(reservation, dict):
            issues.append(f"reservation record is invalid: {reservation!r}")
            continue
        try:
            validate_reservation_record(paths, reservation)
        except (KeyError, TypeError, ValueError) as error:
            issues.append(f"reservation record is invalid: {error}")
            continue
        manifest_path = Path(str(reservation.get("manifest", ""))).resolve()
        reservations_by_manifest.setdefault(manifest_path, []).append(reservation)
        if reservation.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"reservation has wrong execution epoch: {manifest_path}")
        manifest = manifests.get(manifest_path)
        if manifest is None:
            issues.append(f"reservation lacks retained manifest: {manifest_path}")
            continue
        run = manifest.get("run")
        allocation = manifest.get("allocation")
        if not isinstance(run, dict) or not isinstance(allocation, dict):
            issues.append(f"manifest lacks run or allocation metadata: {manifest_path}")
            continue
        for key in ("case_id", "case_name", "segment"):
            if reservation.get(key) != run.get(key):
                issues.append(f"reservation {key} differs from manifest: {manifest_path}")
        for key in ("nodes", "requested_walltime"):
            if reservation.get(key) != allocation.get(key):
                issues.append(
                    f"reservation {key} differs from allocation: {manifest_path}"
                )
        try:
            reserved_difference = abs(
                float(reservation["reserved_node_hours"])
                - float(allocation["reserved_node_hours"])
            )
        except (KeyError, TypeError, ValueError):
            issues.append(f"reservation node-hours are invalid: {manifest_path}")
        else:
            if (
                not math.isfinite(reserved_difference)
                or reserved_difference > 5.0e-12
            ):
                issues.append(
                    f"reservation node-hours differ from allocation: {manifest_path}"
                )
        if reservation.get("state") != manifest.get("state"):
            issues.append(f"reservation state differs from manifest: {manifest_path}")
        try:
            require_reserved_execution_intent(
                reservation, manifest, allow_legacy_local=offline_local_root
            )
        except (KeyError, TypeError, ValueError) as error:
            issues.append(
                f"reservation execution intent differs from manifest "
                f"{manifest_path}: {error}"
            )
        if (
            reservation.get("state") in {"submitted", "recorded", "cancelled"}
            and reservation.get("job_id") != manifest.get("job_id")
        ):
            issues.append(f"reservation job ID differs from manifest: {manifest_path}")

    ledger_by_job: dict[str, list[dict[str, str]]] = {}
    for row in ledger:
        job_id = row.get("job_id", "")
        ledger_by_job.setdefault(job_id, []).append(row)
        if row.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"ledger row has wrong execution epoch: {job_id}")
        try:
            require_numeric_job_id(job_id)
        except ValueError:
            issues.append(f"ledger row has invalid job ID: {job_id!r}")
    for job_id, rows in ledger_by_job.items():
        if len(rows) != 1:
            issues.append(f"ledger has {len(rows)} rows for job {job_id}")

    manifests_by_job: dict[str, list[Path]] = {}
    recorded_manifests_by_job: dict[str, list[Path]] = {}
    for manifest_path, manifest in manifests.items():
        if manifest.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"manifest has wrong execution epoch: {manifest_path}")
        matches = reservations_by_manifest.get(manifest_path, [])
        if len(matches) != 1:
            issues.append(
                f"manifest has {len(matches)} reservation records: {manifest_path}"
            )
        state = manifest.get("state")
        if state not in {"prepared", "submitted", "recorded", "cancelled"}:
            issues.append(f"manifest has invalid state {state!r}: {manifest_path}")
        cancellation = manifest.get("cancellation")
        submitted_cancellation_fields = {
            "scheduler",
            "scheduler_evidence",
            "queue_absence_evidence",
            "recovery_authorization",
        }
        canonical_recovery_manifest = (
            paths["root"].resolve() == DEFAULT_ROOT.resolve()
            and manifest_path
            == (
                paths["root"]
                / str(CANONICAL_SOURCE_BUNDLE_RECOVERY["manifest_relative"])
            ).resolve()
        )
        submitted_cancellation_candidate = (
            state == "cancelled"
            and (
                canonical_recovery_manifest
                or manifest.get("job_id") is not None
                or any(match.get("job_id") is not None for match in matches)
                or (
                    isinstance(cancellation, dict)
                    and bool(submitted_cancellation_fields.intersection(cancellation))
                )
            )
        )
        if submitted_cancellation_candidate:
            try:
                validate_submitted_cancellation_metadata(paths, manifest_path, manifest)
            except (OSError, ValueError, KeyError, TypeError) as error:
                issues.append(
                    f"submitted cancellation evidence drift for {manifest_path}: {error}"
                )
        elif state in {"prepared", "submitted", "recorded"}:
            try:
                command = manifest.get("command", {})
                if (
                    not offline_local_root
                    or (
                        isinstance(command, dict)
                        and isinstance(command.get("overrides"), list)
                    )
                ):
                    prepared_time_tlim_target(manifest)
                authenticate_prepared_execution(
                    manifest, manifest_path,
                    allow_legacy_local=offline_local_root,
                )
            except (OSError, ValueError, KeyError, TypeError) as error:
                issues.append(f"prepared artifact drift for {manifest_path}: {error}")
        if state == "prepared":
            if manifest.get("job_id") is not None:
                issues.append(f"prepared manifest unexpectedly has a job ID: {manifest_path}")
            continue
        if state not in {"submitted", "recorded"}:
            continue
        job_id = str(manifest.get("job_id", ""))
        try:
            require_numeric_job_id(job_id)
        except ValueError:
            issues.append(f"manifest has invalid job ID: {manifest_path}")
        manifests_by_job.setdefault(job_id, []).append(manifest_path)
        if state != "recorded":
            continue
        recorded_manifests_by_job.setdefault(job_id, []).append(manifest_path)
        accounting = manifest.get("accounting")
        if not isinstance(accounting, dict) or accounting.get("job_id") != job_id:
            issues.append(f"recorded manifest lacks matching accounting: {manifest_path}")
        rows = ledger_by_job.get(job_id, [])
        if len(rows) != 1:
            issues.append(f"recorded manifest lacks one ledger row: {manifest_path}")
            continue
        row = rows[0]
        if accounting != row:
            issues.append(f"manifest accounting differs from ledger: {manifest_path}")
        try:
            validate_transaction_ledger_row(row, manifest)
        except (KeyError, TypeError, ValueError) as error:
            issues.append(
                f"recorded ledger provenance differs for {manifest_path}: {error}"
            )
        if not offline_local_root:
            try:
                require_recorded_scheduler_evidence(paths, row, manifest)
            except (OSError, KeyError, TypeError, ValueError) as error:
                issues.append(
                    f"recorded scheduler evidence differs for {manifest_path}: {error}"
                )
        run = manifest.get("run")
        if not isinstance(run, dict):
            issues.append(f"recorded manifest lacks run metadata: {manifest_path}")
        else:
            for key in ("case_id", "case_name", "segment"):
                if row.get(key) != run.get(key):
                    issues.append(f"ledger {key} differs from manifest: {manifest_path}")
        if len(matches) == 1:
            reservation = matches[0]
            if reservation.get("result") != row.get("result"):
                issues.append(f"reservation result differs from ledger: {manifest_path}")
            try:
                difference = abs(
                    float(reservation["actual_node_hours"])
                    - float(row["actual_node_hours"])
                )
            except (KeyError, TypeError, ValueError):
                issues.append(
                    f"reservation actual node-hours are invalid: {manifest_path}"
                )
            else:
                if not math.isfinite(difference) or difference > 5.0e-7:
                    issues.append(
                        f"reservation actual node-hours differ from ledger: "
                        f"{manifest_path}"
                    )
        if isinstance(accounting, dict) and accounting.get("result") in {
            "accepted", "clean_partial"
        }:
            inspection = manifest.get("scientific_inspection")
            if not isinstance(inspection, dict):
                issues.append(f"recorded manifest lacks inspection: {manifest_path}")
            else:
                try:
                    inspected_target = float(inspection["required_time"])
                    if (
                        not math.isfinite(inspected_target)
                        or abs(inspected_target - prepared_time_tlim_target(manifest))
                        > 1.0e-12
                    ):
                        raise ValueError("inspection target differs from preparation")
                    revalidate_inspection_files(inspection, manifest)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    issues.append(f"inspection drift for {manifest_path}: {error}")

    for job_id, manifest_paths in manifests_by_job.items():
        if len(manifest_paths) != 1:
            issues.append(
                f"job {job_id} is attached to {len(manifest_paths)} manifests"
            )
    for job_id in ledger_by_job:
        matches = recorded_manifests_by_job.get(job_id, [])
        if len(matches) != 1:
            issues.append(f"ledger job {job_id} has {len(matches)} recorded manifests")

    return {
        "execution_epoch": EXECUTION_EPOCH,
        "root": str(root),
        "qualification": qualification,
        "consistent": not issues,
        "counts": {
            "transactions": len(transactions),
            "reservations": len(reservations),
            "active_reservations": len(active),
            "ledger_rows": len(ledger),
            "manifests": len(manifests),
        },
        "issues": issues,
    }


@locked_root_action
def reconcile(args: argparse.Namespace) -> int:
    """Print a read-only reservations, ledger, and manifest consistency report."""

    root = require_root(Path(args.root), args.allow_local_root)
    report = reconcile_report(root)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["consistent"] else 1


def parser() -> argparse.ArgumentParser:
    """Build command-line parsing."""

    command = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "E03 production qualification remains pending until the canonical "
            "corrected-build approval token exists. Create it only after review "
            "with approve-qualification."
        ),
    )
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument(
        "--allow-local-root", action="store_true",
        help="Permit a non-project root only for offline validation.",
    )
    actions = command.add_subparsers(dest="action", required=True)
    validate = actions.add_parser("validate-matrix")
    validate.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    validate.add_argument("--source-dir", default=str(ROOT_DIR))
    actions.add_parser("init")
    approved = actions.add_parser(
        "approve-qualification",
        help=(
            "Atomically create the E03 corrected-build qualification token "
            "after Frontier review."
        ),
    )
    approved.add_argument("--executable", required=True)
    approved.add_argument("--build-manifest", required=True)
    approved.add_argument("--approved-by", required=True)
    approved.add_argument("--review-notes", required=True)
    approved.add_argument(
        "--confirm-corrected-build-frontier-qualified", action="store_true"
    )
    approved.add_argument("--replace-existing-approval", action="store_true")
    prepare_parser = actions.add_parser(
        "prepare",
        help=(
            "Prepare one E03 segment; canonical production rejects until the "
            "qualification token exists."
        ),
    )
    prepare_parser.add_argument("--case-id", required=True)
    prepare_parser.add_argument("--segment", type=require_safe_segment, required=True)
    prepare_parser.add_argument("--acceptance-criterion", required=True)
    prepare_parser.add_argument("--executable", required=True)
    prepare_parser.add_argument("--build-manifest", required=True)
    prepare_parser.add_argument("--source-dir", default=str(ROOT_DIR))
    prepare_parser.add_argument(
        "--source-bundle",
        help=(
            "Retained Git bundle beneath the CGL root containing the input "
            "and executable revisions."
        ),
    )
    prepare_parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    prepare_parser.add_argument("--restart-file")
    prepare_parser.add_argument("--nodes", type=int, required=True)
    prepare_parser.add_argument("--walltime", required=True)
    prepare_parser.add_argument("--athena-walltime", required=True)
    prepare_parser.add_argument("--ranks-per-node", type=int, default=8)
    prepare_parser.add_argument("--cpus-per-task", type=int, default=7)
    prepare_parser.add_argument("--override", action="append", default=[])
    prepare_parser.add_argument(
        "--allow-missing-time-target", action="store_true",
        help=(
            "Permit offline local-root preparation without a time/tlim "
            "override. This is never allowed for real production."
        ),
    )
    prepare_parser.add_argument(
        "--allow-missing-restart-time-marker", action="store_true",
        help=(
            "Permit an offline legacy continuation fixture without explicit "
            "time/restart_time markers. This is never allowed for production."
        ),
    )
    checked = actions.add_parser("check-submit")
    submitted_atomically = actions.add_parser("submit")
    for submit_parser in (checked, submitted_atomically):
        submit_parser.add_argument("--manifest", required=True)
        submit_parser.add_argument("--squeue-file")
        submit_parser.add_argument("--skip-slurm-test", action="store_true")
        submit_parser.add_argument(
            "--allow-shared-root-campaign", action="append", default=[],
            help=(
                "Acknowledge one reviewed top-level CGL-root campaign record. "
                "Queued user jobs still fail closed."
            ),
        )
    clearance_refresh = actions.add_parser(
        "render-shared-root-clearance-refresh",
        help=(
            "Print a read-only exact-binding shared-root clearance candidate "
            "for independent review and immutable publication."
        ),
    )
    clearance_refresh.add_argument("--checkpoint", type=int, required=True)
    submitted_atomically.add_argument(
        "--sbatch-output-file",
        help="Use retained sbatch --parsable output only for offline validation.",
    )
    submitted = actions.add_parser("mark-submitted")
    submitted.add_argument("--manifest", required=True)
    submitted.add_argument("--job-id", type=require_numeric_job_id, required=True)
    recovered_submit = actions.add_parser("recover-submit")
    recovered_submit.add_argument("--manifest", required=True)
    recovered_submit.add_argument(
        "--job-id", type=require_numeric_job_id, required=True
    )
    cleared_submit = actions.add_parser("clear-submit-pending")
    cleared_submit.add_argument("--manifest", required=True)
    cleared_submit.add_argument("--notes", required=True)
    cleared_submit.add_argument("--confirm-no-job-submitted", action="store_true")
    cleared_submit.add_argument(
        "--scheduler-absence-evidence-file",
        help="Use machine-readable absence evidence only for offline local fixtures.",
    )
    cleared_submit.add_argument("--break-glass-clear-evidence")
    cleared_submit.add_argument("--confirm-break-glass-clear", action="store_true")
    inspected = actions.add_parser("inspect-segment")
    inspected.add_argument("--manifest", required=True)
    inspected.add_argument("--required-time", type=float, required=True)
    inspected.add_argument(
        "--allow-missing-restart-time-marker", action="store_true",
        help=(
            "Permit offline inspection of legacy restart fixtures without "
            "explicit time/restart_time markers."
        ),
    )
    recorded = actions.add_parser("record")
    recorded.add_argument("--manifest", required=True)
    recorded.add_argument("--job-id", type=require_numeric_job_id, required=True)
    recorded.add_argument(
        "--result",
        choices=("accepted", "clean_partial", "rejected", "failed", "aborted"),
        required=True,
    )
    recorded.add_argument("--notes", default="")
    recorded.add_argument("--sacct-file")
    bundled = actions.add_parser("bundle-case")
    bundled.add_argument("--case-id", required=True)
    bundled.add_argument(
        "--required-final-time", type=positive_finite_float_arg,
        default=REQUIRED_CASE_FINAL_TIME
    )
    bundled.add_argument("--source-dir", default=str(ROOT_DIR))
    bundled.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    bundled.add_argument("--output-dir")
    bundled.add_argument("--replace", action="store_true")
    campaign = actions.add_parser("bundle-campaign")
    campaign.add_argument(
        "--required-final-time", type=positive_finite_float_arg,
        default=REQUIRED_CASE_FINAL_TIME
    )
    campaign.add_argument("--source-dir", default=str(ROOT_DIR))
    campaign.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    campaign.add_argument("--output-dir")
    campaign.add_argument("--replace", action="store_true")
    cancelled = actions.add_parser("cancel")
    cancelled.add_argument("--manifest", required=True)
    cancelled.add_argument("--notes", required=True)
    cancelled.add_argument("--break-glass-cancel-evidence")
    cancelled.add_argument("--confirm-break-glass-cancel", action="store_true")
    cancelled_submitted = actions.add_parser("cancel-submitted")
    cancelled_submitted.add_argument("--manifest", required=True)
    cancelled_submitted.add_argument(
        "--job-id", type=require_numeric_job_id, required=True
    )
    cancelled_submitted.add_argument("--notes", required=True)
    cancelled_submitted.add_argument("--sacct-file")
    cancelled_submitted.add_argument("--squeue-file")
    actions.add_parser("summary")
    actions.add_parser("reconcile")
    actions.add_parser("recover-transactions")
    return command


def main() -> int:
    """Command-line entry point."""

    try:
        authenticate_controller_runtime(sys.argv[1:])
        args = parser().parse_args()
        if args.action == "validate-matrix":
            matrix = validate_matrix(
                Path(args.matrix).expanduser().resolve(),
                Path(args.source_dir).expanduser().resolve(),
            )
            print(f"Validated {len(matrix['cases'])} mapped Stage I cases.")
            return 0
        root = require_root(Path(args.root), args.allow_local_root)
        if args.action == "init":
            initialize(root)
            print(f"Initialized Stage I production accounting beneath {root}.")
            return 0
        if args.action == "approve-qualification":
            return approve_qualification(args)
        if args.action == "prepare":
            prepare(args)
            return 0
        if args.action == "check-submit":
            return check_submit(args)
        if args.action == "submit":
            return submit(args)
        if args.action == "render-shared-root-clearance-refresh":
            return render_managed_shared_root_clearance_refresh(args)
        if args.action == "mark-submitted":
            return mark_submitted(args)
        if args.action == "recover-submit":
            return recover_submit(args)
        if args.action == "clear-submit-pending":
            return clear_submit_pending(args)
        if args.action == "inspect-segment":
            return inspect_segment(args)
        if args.action == "record":
            return record(args)
        if args.action == "bundle-case":
            return bundle_case(args)
        if args.action == "bundle-campaign":
            return bundle_campaign(args)
        if args.action == "cancel":
            return cancel(args)
        if args.action == "cancel-submitted":
            return cancel_submitted(args)
        if args.action == "summary":
            with canonical_root_lock(root):
                if is_offline_local_root(root, args.allow_local_root):
                    paths = initialize(root)
                else:
                    paths = layout(root)
                    require_existing_layout(paths)
                    require_no_pending_transactions(paths)
                    require_no_orphaned_segment_runs(paths)
                    require_reconciled_store_consistency(
                        paths, allow_absent_qualification=True
                    )
                refresh_summary(paths)
            print(f"Wrote {paths['summary']}")
            return 0
        if args.action == "reconcile":
            return reconcile(args)
        if args.action == "recover-transactions":
            return recover_transactions(args)
        raise ValueError(f"unsupported action: {args.action}")
    except (KeyError, OSError, TypeError, ValueError,
            subprocess.CalledProcessError) as error:
        print(f"Stage I production utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
