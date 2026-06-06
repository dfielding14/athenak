#!/usr/bin/env python3
"""Generate deterministic, staged-only Stage I recost recommendation evidence.

The generator consumes one checksum-bound independently reviewed request.  Its
output is always non-authorizing evidence: an exact independent review and
publication audit are required before any controller may consume a
recommendation.  It never submits, promotes, or mutates campaign accounting.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from decimal import Decimal, InvalidOperation, ROUND_CEILING, ROUND_HALF_EVEN
import errno
import fcntl
import hashlib
from importlib.machinery import SourceFileLoader
import importlib.util
import io
import json
import math
import os
from pathlib import Path, PurePosixPath
import pwd
import re
import stat
import subprocess
import sys
import tempfile
from typing import Callable, Iterator


DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
RECOST_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i_recost.py")
STAGE_I_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i.py")
MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
F113_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F113_controller_transition_evidence.json"
)
F113_PUBLICATION_AUDIT_RELATIVE = Path(f"{F113_RELATIVE}.publication_audit.json")
F115_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F115_source_bundle_recovery_supersession_evidence.json"
)
F115_PUBLICATION_AUDIT_RELATIVE = Path(f"{F115_RELATIVE}.publication_audit.json")
F115_PROVENANCE_REVIEW_RELATIVE = Path(f"{F115_RELATIVE}.provenance_security_review.json")
F115_PLASMA_REVIEW_RELATIVE = Path(f"{F115_RELATIVE}.plasma_scientific_review.json")
F115_CANONICAL_SHA256 = {
    "evidence": "cb50beb064678a9446ac33801a0023547d06c432d8c59b6bd3fbe34b11cf0391",
    "publication_audit": "5923e3872b1d4a84a147bcd1d81fddc20ee79b72bfc410683d083039781f5b1f",
    "provenance_review": "6fcd19f9267f36332742f1f103968098216bd8e6b42fa9821964ab8df704bacd",
    "plasma_review": "a78357ed90e593809b1a82a641d2b40d651b16569ec943fc1781940e569b8440",
}
F116_RELATIVE = Path(
    "accounting/"
    "mks24_stage_i_E03_forcing_policy_F116_current_source_authority_supersession_evidence.json"
)
F116_PUBLICATION_AUDIT_RELATIVE = Path(f"{F116_RELATIVE}.publication_audit.json")
F116_PROVENANCE_REVIEW_RELATIVE = Path(f"{F116_RELATIVE}.provenance_security_review.json")
F116_PLASMA_REVIEW_RELATIVE = Path(f"{F116_RELATIVE}.plasma_scientific_review.json")
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
F116_REQUIRED_COMMITTED_TOOLS = {
    "scripts/frontier/cgl_lf_stage_i.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_checkpoint.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_qualification.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_recost.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_source_authority.py": "0755",
    "scripts/frontier/cgl_lf_stage_i_validate_segment.py": "0644",
    "scripts/frontier/cgl_lf_stage_i_wave_plan.py": "0644",
}
QUALIFICATION_APPROVAL_RELATIVE = Path(
    "accounting/mks24_stage_i_E03_forcing_policy_qualification_approval.json"
)
R17_READINESS_RELATIVE = Path(
    "accounting/mks24_stage_i_E03_forcing_policy_R17_readiness_evidence.json"
)
F113_CANONICAL_SHA256 = (
    "a54c52b7c96e974d634cd1c1894559442b37163a6f51ebab8d0d3a155b6a0dea"
)
LEGACY_F114_RECOST_NAME = (
    "mks24_stage_i_E03_forcing_policy_R03_s00_clean_partial_recost_evidence.json"
)
LEGACY_F114_RECOST_SHA256 = (
    "27dd154b8594693afe0308f4de63a7faaed6daf81a396ab31166bb70a8c43b86"
)
LEGACY_F114_PUBLICATION_AUDIT_SHA256 = (
    "3bf50ef1359f7f1798da945185d42b61488cb98779d5a6942839c03e3eb6bc73"
)
AUTHORIZED_CASE_IDS = frozenset(f"R{number:02d}" for number in range(2, 18))
R17_CASE_ID = "R17"
R17_PREDECESSOR_CASE_IDS = tuple(f"R{number:02d}" for number in range(2, 17))
R12_HISTORICAL_JOB_ID = "4766856"
R12_HISTORICAL_SEGMENT = "s00_rankio_t0_t0p25"
R12_HISTORICAL_INVENTORY_IDENTITY = (
    R12_HISTORICAL_JOB_ID,
    "R12",
    R12_HISTORICAL_SEGMENT,
    "clean_partial",
)
R12_FRESH_RERUN_SEGMENT = "s01_rankio_t0_t0p12"
R12_FRESH_RERUN_NODES = 4
R12_FRESH_RERUN_RANKS_PER_NODE = 8
R12_FRESH_RERUN_TOTAL_RANKS = 32
R12_FRESH_RERUN_WALLTIME = "02:00:00"
R12_FRESH_RERUN_ATHENA_WALLTIME = "01:50:00"
R12_FRESH_RERUN_TARGET = 0.12
REQUIRED_CASE_FINAL_TIME = 10.0
MAX_WAVE_NODES = 10
MAX_WALLTIME_SECONDS = 2 * 60 * 60
MINIMUM_SHUTDOWN_MARGIN_SECONDS = 10 * 60
EXPECTED_RANKS_PER_NODE = 8
EXPECTED_CPUS_PER_TASK = 7
ACCOUNT = "AST207"
PARTITION = "batch"
R17_MINIMUM_RETAINED_BYTES = 958_271_710_272
R17_MAX_NORMALIZED_CT_DIVB_TEXT = "1e-12"
R17_PHYSICS_MEASUREMENT_KEYS = {
    "finite_rank_outputs",
    "mass_relative_drift_max",
    "mhd_user_mass_mismatch_max",
    "lf_bad_counts_total",
    "normalized_ct_divb_max",
    "normalized_ct_divb_threshold",
    "normalized_ct_divb_below_threshold",
}
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
R17_SCIENTIFIC_EVIDENCE_KEYS = {
    "schema_version",
    "record_type",
    "case_id",
    "nodes",
    "target_time",
    "scientific_policy",
    "executable",
    "execution_intent_sha256",
    "execution_contract_sha256",
    "forcing_closure_normalized_residual",
    "mhd_history",
    "user_history",
    "snapshots",
    "snapshot_times",
    "restarts",
    "restart_times",
    "restart_load_smoke",
    "terminal_rank_local_outputs",
    "terminal_rank_local_output_inventory_sha256",
    "r17_decomposition",
    "terminal_rank_local_restarts",
    "terminal_rank_local_restart_inventory_sha256",
    "physics_measurements",
    "scientific_agreement_signature",
    "checks",
    "accepted_for_profile_selection",
    "accepted_for_operational_qualification",
}
MINIMUM_PROFILE_STORAGE_BYTES = 1024
STORAGE_EVIDENCE_MAX_AGE = timedelta(hours=24)
AUTHORIZATION_MAX_LIFETIME = timedelta(hours=24)
SCHEDULER_TIME_TOLERANCE = timedelta(minutes=5)
MAX_NORMALIZED_CT_DIVB = Decimal("1e-12")
STORAGE_PROJECTION_METHOD = "observed-stage-i-output-byte-rate-v1"
NODE_HOUR_PROJECTION_METHOD = "observed-stage-i-scoped-node-hour-rate-v2"
F113_REVIEW_AUTHORITY = "campaign-authorized-independent-review"
F113_AUTHORIZED_REVIEWERS = frozenset(
    {
        "Codex execution agent with independent Turing, Hubble, Darwin, and Linnaeus audits",
    }
)
SCIENTIFIC_RESULTS = frozenset({"accepted", "clean_partial"})
RECORDED_RESULTS = frozenset(
    {"accepted", "clean_partial", "rejected", "failed", "aborted"}
)
TERMINAL_SCHEDULER_STATES = frozenset(
    {
        "BOOT_FAIL",
        "CANCELLED",
        "COMPLETED",
        "DEADLINE",
        "FAILED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "REVOKED",
        "SPECIAL_EXIT",
        "TIMEOUT",
    }
)
CONTINUATION_RUNTIME_SAFETY_FACTOR = Decimal("0.90")
DEFAULT_CONTINUATION_TARGET_ALIGNMENT = Decimal("0.25")
CONTINUATION_TARGET_ALIGNMENTS = {"R12": Decimal("0.02")}
NODE_HOUR_QUANTUM = Decimal("0.000001")
CONTROLLER_CLEAN_PARTIAL_SCHEMA_VERSION = 4
CONTROLLER_STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
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
CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS = {
    "schema_version",
    "execution_epoch",
    "inspected_utc",
    "manifest",
    "job_id",
    "case_id",
    "segment",
    "required_time",
    "final_time",
    "maximum_strict_failure_counts",
    "checks",
    "accepted",
    "clean_for_continuation",
    "mhd_history",
    "user_history",
    "plasma_continuation_policy",
    "plasma_continuation_evidence",
    "snapshots",
    "snapshot_times",
    "restarts",
    "restart_times",
    "restart_time_marker_modes",
    "terminal_restart",
    "terminal_restart_time",
    "restart_time_marker_bypass",
    "final_hardwall_projection_count",
}
FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS = (
    CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS
    - {"plasma_continuation_policy", "plasma_continuation_evidence"}
)
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
R17_CT_LIMITATION = FROZEN_E03_CT_DIVERGENCE_REASON
REQUEST_REVIEW_IDENTITY_ASSURANCE = "declared-process-independence-non-cryptographic"
REQUEST_REVIEW_IDENTITY_LIMITATION = (
    "Reviewer identity and process independence are declared evidence, "
    "not cryptographically proven."
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
PLAN_INITIAL_TARGETS = {
    **{f"R{number:02d}": 0.25 for number in range(3, 6)},
    "R06": 0.50,
    **{f"R{number:02d}": 0.25 for number in range(7, 12)},
    "R12": R12_FRESH_RERUN_TARGET,
    **{f"R{number:02d}": 0.25 for number in range(13, 16)},
    "R16": 1.50,
    "R17": 0.25,
}
PLAN_MAX_INCREMENTS = {
    "R03": 0.50,
    **{f"R{number:02d}": 1.0 for number in range(4, 16)},
    "R16": 3.0,
    "R17": 0.25,
}
ACCEPTANCE_POLICIES = {
    **{case_id: "U+A+H" for case_id in (
        "R02", "R03", "R04", "R05", "R10", "R11", "R12", "R13", "R16", "R17"
    )},
    **{case_id: "U+P+H" for case_id in ("R06", "R07", "R08", "R09")},
    **{case_id: "U+A+F" for case_id in ("R14", "R15")},
}
ACCEPTANCE_POLICY_TEXT = {
    "U": (
        "Require the exact endpoint; finite synchronized MHD and user histories; "
        "relative mass drift and MHD/user mass mismatch <= 1e-12; zero lf_dfloor, "
        "lf_pfloor, lf_nonfin, lf_nonpos, lf_hardbd, and hard_vol at every retained "
        "row; positive interval lf_nstage and lf_qface with bounded cap increments; "
        "finite LF heat and pressure-work ledgers; complete ranked products; and one "
        "loadable exact terminal snapshot and restart group. Frozen E03 does not emit "
        "an independently auditable CT-divergence diagnostic, so this recost makes no "
        "CT-divergence bound claim."
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
EXPECTED_PROMOTED_PROFILES = {
    "R03": {1},
    **{f"R{number:02d}": {1, 2, 4} for number in range(4, 16)},
    "R16": {1, 2},
    "R17": {8},
}
SOLE_COMPATIBLE_PROFILE_KEYS = (
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
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
SAFE_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}")
JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
OUTPUT_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}\.json\.staged")
RECOST_ARTIFACT_PATTERN = re.compile(
    r"mks24_stage_i_E03_forcing_policy_F([0-9]+)_recost_evidence\.json"
)
SEGMENT_PATTERN = re.compile(
    r"s(?P<index>[0-9]+)_rankio_t(?P<start>[0-9]+(?:p[0-9]+)?)"
    r"_t(?P<target>[0-9]+(?:p[0-9]+)?)"
)
CANONICAL_NODE_HOURS_PATTERN = re.compile(r"(?:0|[1-9][0-9]*)\.[0-9]{6}")
SELF_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_RECOST_DESCRIPTOR"
PYTHON_DESCRIPTOR_ENV = "_CGL_LF_STAGE_I_RECOST_PYTHON_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_STAGE_I_RECOST_SOURCE"
REPOSITORY_ROOT_ENV = "_CGL_LF_STAGE_I_RECOST_REPOSITORY_ROOT"
CGL_JOB_NAME_PREFIX = "cgl_"
SQUEUE = Path("/usr/bin/squeue")
SACCT = Path("/usr/bin/sacct")
GIT = Path("/usr/lib/git/git")
GIT_EXEC_PATH = Path("/usr/lib/git")
TRUSTED_SYSTEM_PATH = "/usr/bin:/bin"
ACCOUNT_SACCT_FIELDS = (
    "JobIDRaw",
    "JobName",
    "State",
    "ExitCode",
    "NNodes",
    "ElapsedRaw",
    "Submit",
    "Start",
    "End",
    "Partition",
    "Account",
    "User",
)
ACCOUNT_SCHEDULER_HEADER = "|".join(ACCOUNT_SACCT_FIELDS)
ACCOUNT_MISSING_TIMESTAMPS = frozenset({"", "Unknown", "N/A", "None"})
LIVE_JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*(?:_[0-9]+|\+[0-9]+)?")
COUNT_KEYS = (
    "transactions",
    "reservations",
    "active_reservations",
    "ledger_rows",
    "manifests",
)
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


def sha256_bytes(payload: bytes) -> str:
    """Return one SHA-256 digest."""

    return hashlib.sha256(payload).hexdigest()


def sha256_descriptor(descriptor: int) -> str:
    """Hash one descriptor without changing its caller-visible offset."""

    offset = os.lseek(descriptor, 0, os.SEEK_CUR)
    digest = hashlib.sha256()
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                return digest.hexdigest()
            digest.update(block)
    finally:
        os.lseek(descriptor, offset, os.SEEK_SET)


def require_sha256(value: object, label: str) -> str:
    """Require a lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def sha256_arg(value: str) -> str:
    """Parse one SHA-256 CLI argument."""

    try:
        return require_sha256(value, "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def require_revision(value: object, label: str) -> str:
    """Require one full Git commit identifier."""

    if not isinstance(value, str) or REVISION_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a full lowercase Git commit identifier")
    return value


def require_exact_keys(value: object, keys: set[str], label: str) -> dict[str, object]:
    """Require one object with an exact schema."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    retained = set(value)
    if retained != keys:
        raise ValueError(
            f"{label} schema differs; expected {sorted(keys)}, found {sorted(retained)}"
        )
    return value


def require_nonempty_string(value: object, label: str) -> str:
    """Require one nonempty, trimmed string."""

    if not isinstance(value, str) or not value or value.strip() != value:
        raise ValueError(f"{label} must be a nonempty trimmed string")
    return value


def require_safe_id(value: object, label: str) -> str:
    """Require one filesystem-safe identifier."""

    retained = require_nonempty_string(value, label)
    if SAFE_ID_PATTERN.fullmatch(retained) is None:
        raise ValueError(f"{label} is not a safe identifier")
    return retained


def require_job_id(value: object, label: str) -> str:
    """Require one top-level numeric Slurm job identifier."""

    retained = require_nonempty_string(value, label)
    if JOB_ID_PATTERN.fullmatch(retained) is None:
        raise ValueError(f"{label} must be a numeric Slurm job ID")
    return retained


def require_integer(value: object, label: str, *, minimum: int = 0) -> int:
    """Require one non-boolean integer."""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{label} must be an integer >= {minimum}")
    return value


def require_decimal(value: object, label: str, *, positive: bool = False) -> Decimal:
    """Require one finite decimal encoded as a string."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be encoded as a decimal string")
    try:
        retained = Decimal(value)
    except InvalidOperation as error:
        raise ValueError(f"{label} must be a decimal string") from error
    if not retained.is_finite() or (positive and retained <= 0):
        raise ValueError(f"{label} must be finite" + (" and positive" if positive else ""))
    return retained


def decimal_string(value: Decimal) -> str:
    """Return a non-exponent decimal representation."""

    return format(value, "f")


def quantized_node_hours(value: Decimal, label: str) -> Decimal:
    """Return one finite nonnegative node-hour value quantized to the ledger unit."""

    if not value.is_finite() or value < 0:
        raise ValueError(f"{label} must be finite and nonnegative")
    try:
        return value.quantize(NODE_HOUR_QUANTUM, rounding=ROUND_HALF_EVEN)
    except InvalidOperation as error:
        raise ValueError(f"{label} cannot be represented as canonical node-hours") from error


def canonical_node_hours(value: Decimal, label: str) -> str:
    """Return the exact six-decimal serialization used by Stage I accounting."""

    return format(quantized_node_hours(value, label), ".6f")


def require_canonical_node_hours(
    value: object, label: str, *, positive: bool = False
) -> Decimal:
    """Require one exact nonnegative six-decimal Stage I ledger value."""

    if not isinstance(value, str) or CANONICAL_NODE_HOURS_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a canonical six-decimal node-hour string")
    retained = Decimal(value)
    if positive and retained <= 0:
        raise ValueError(f"{label} must be positive")
    return retained


def decimal_ceiling(value: Decimal) -> int:
    """Return one finite nonnegative decimal rounded upward to an integer."""

    if not value.is_finite() or value < 0:
        raise ValueError("cannot round an invalid projected value")
    return int(value.to_integral_value(rounding=ROUND_CEILING))


def require_finite_float(value: object, label: str) -> float:
    """Require one finite non-boolean number."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{label} must be a finite number")
    retained = float(value)
    if not math.isfinite(retained):
        raise ValueError(f"{label} must be a finite number")
    return retained


def parse_utc_timestamp(value: object, label: str) -> datetime:
    """Parse and require one UTC timestamp."""

    retained = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(retained.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None or parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")
    return parsed


def parse_scheduler_timestamp(value: object, label: str) -> datetime:
    """Parse one scheduler timestamp, treating retained naive sacct time as UTC."""

    retained = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(retained.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    if parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")
    return parsed


def acceptance_criterion(case_id: str, target: float) -> str:
    """Return the exact case-aware Stage I acceptance criterion."""

    policy = ACCEPTANCE_POLICIES[case_id]
    codes = policy.split("+")
    target_text = format(Decimal(str(target)).normalize(), "f")
    clauses = " ".join(f"{code}: {ACCEPTANCE_POLICY_TEXT[code]}" for code in codes)
    return (
        f"Accept {case_id} only at exact t={target_text} under policies {policy}. "
        f"{clauses} Any endpoint, provenance, product-inventory, scheduler-accounting, "
        "or policy failure blocks acceptance and successor packet planning."
    )


def parse_segment(value: str, label: str) -> tuple[int, float, float]:
    """Parse one exact rank-local Stage I segment identifier."""

    match = SEGMENT_PATTERN.fullmatch(value)
    if match is None:
        raise ValueError(f"{label} must be a rank-local Stage I segment")
    try:
        start = float(Decimal(match.group("start").replace("p", ".")))
        target = float(Decimal(match.group("target").replace("p", ".")))
    except InvalidOperation as error:
        raise ValueError(f"{label} times are invalid") from error
    if not math.isfinite(start) or not math.isfinite(target) or start < 0 or target <= start:
        raise ValueError(f"{label} interval is invalid")
    return int(match.group("index")), start, target


def walltime_seconds(value: object, label: str) -> int:
    """Parse one strict HH:MM:SS walltime."""

    retained = require_nonempty_string(value, label)
    match = re.fullmatch(r"([0-9]{2}):([0-5][0-9]):([0-5][0-9])", retained)
    if match is None:
        raise ValueError(f"{label} must use HH:MM:SS")
    hours, minutes, seconds = (int(item) for item in match.groups())
    total = hours * 3600 + minutes * 60 + seconds
    if total <= 0 or total > MAX_WALLTIME_SECONDS:
        raise ValueError(f"{label} must be within the two-hour Stage I limit")
    return total


def require_continuation_target_alignment(case_id: str, target: float) -> None:
    """Require the case-specific exact retained-output continuation cadence."""

    alignment = CONTINUATION_TARGET_ALIGNMENTS.get(
        case_id, DEFAULT_CONTINUATION_TARGET_ALIGNMENT
    )
    if Decimal(str(target)) % alignment:
        if case_id == "R12":
            raise ValueError(
                "next profile R12 target is not aligned to the exact 0.02 "
                "fresh-R12 output cadence"
            )
        raise ValueError(
            f"next profile {case_id} target is not aligned to the exact "
            "quarter-time output cadence"
        )


def expected_self_sha256(argv: list[str]) -> str:
    """Read the externally supplied self digest before normal argument parsing."""

    option = "--expected-generator-sha256"
    matches = [
        argv[index + 1]
        for index, item in enumerate(argv[:-1])
        if item == option
    ]
    if len(matches) != 1:
        raise ValueError(f"{option} must be supplied exactly once")
    return require_sha256(matches[0], option)


def initial_source_path() -> Path:
    """Return the named source path for the initial descriptor re-execution."""

    source = Path(__file__).absolute()
    if source.is_symlink():
        raise ValueError("generator source must not be a symlink")
    return source.resolve()


def initial_repository_root(source: Path) -> Path:
    """Return the named repository root for the initial descriptor re-execution."""

    try:
        return source.parents[2].resolve()
    except IndexError as error:
        raise ValueError(f"generator source path is invalid: {source}") from error


def require_reexec_source_relationship(source: Path, repository: Path) -> None:
    """Bind private reexec metadata to the tool's exact repository-relative path."""

    if source != repository / RECOST_RELATIVE:
        raise ValueError("authenticated generator source/repository relationship differs")


def inherited_reexec_path(name: str, label: str) -> Path:
    """Read one private reexec path only after the self descriptor is authenticated."""

    retained = os.environ.get(name)
    if retained is None:
        raise ValueError(f"authenticated generator reexecution lacks {label}")
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != Path(os.path.normpath(retained))
        or ".." in path.parts
        or len(path.parts) < 2
    ):
        raise ValueError(f"authenticated generator {label} is not normalized and absolute")
    require_no_symlink_components(path, f"authenticated generator {label}")
    return path


def reexec_environment(
    descriptor: int, python_descriptor: int, source: Path, repository: Path
) -> dict[str, str]:
    """Return a reexec environment without caller-controlled interpreter/loader state."""

    return {
        SELF_DESCRIPTOR_ENV: str(descriptor),
        PYTHON_DESCRIPTOR_ENV: str(python_descriptor),
        SELF_SOURCE_ENV: str(source),
        REPOSITORY_ROOT_ENV: str(repository),
        "HOME": "/nonexistent",
        "LC_ALL": "C",
        "PATH": TRUSTED_SYSTEM_PATH,
        "PYTHONDONTWRITEBYTECODE": "1",
        "XDG_CONFIG_HOME": "/nonexistent",
    }


def require_regular_profile(profile: os.stat_result, label: str, *,
                            expected_mode: int | None = None) -> None:
    """Require one owned, single-link regular retained file."""

    if not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"{label} must be a regular file")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"{label} must be owned by the effective user")
    if profile.st_nlink != 1:
        raise ValueError(f"{label} must have exactly one link")
    mode = stat.S_IMODE(profile.st_mode)
    if expected_mode is not None and mode != expected_mode:
        raise ValueError(f"{label} mode is {mode:04o}, expected {expected_mode:04o}")
    if mode & 0o022:
        raise ValueError(f"{label} must not be group- or world-writable")


def require_system_executable_profile(profile: os.stat_result, label: str) -> None:
    """Require one root-owned, single-link, non-writable executable."""

    mode = stat.S_IMODE(profile.st_mode)
    if (
        not stat.S_ISREG(profile.st_mode)
        or profile.st_uid != 0
        or profile.st_nlink != 1
        or mode & 0o022
        or mode & 0o111 == 0
    ):
        raise ValueError(f"{label} does not have the trusted system executable profile")


def same_inode(first: os.stat_result, second: os.stat_result) -> bool:
    """Return whether two profiles identify the same retained filesystem object."""

    return (first.st_dev, first.st_ino) == (second.st_dev, second.st_ino)


def inherited_descriptor(name: str, label: str) -> int:
    """Parse one authenticated inherited descriptor marker."""

    retained = os.environ.get(name)
    if retained is None or re.fullmatch(r"[0-9]+", retained) is None:
        raise ValueError(f"{label} descriptor marker is invalid")
    try:
        return int(retained)
    except ValueError as error:
        raise ValueError(f"{label} descriptor marker is invalid") from error


def require_isolated_python_reexec(descriptor: int) -> None:
    """Bind inherited execution to one authenticated isolated Python interpreter."""

    profile = os.fstat(descriptor)
    require_system_executable_profile(profile, "authenticated Python interpreter")
    running = os.stat("/proc/self/exe", follow_symlinks=True)
    if not same_inode(profile, running):
        raise ValueError("authenticated Python descriptor is not the running interpreter")
    if (
        sys.flags.isolated != 1
        or sys.flags.no_site != 1
        or sys.flags.dont_write_bytecode != 1
    ):
        raise ValueError("authenticated Python reexecution is not isolated")


def authenticate_self(argv: list[str]) -> tuple[Path, Path, str]:
    """Re-execute immutable generator bytes and return retained source metadata."""

    expected = expected_self_sha256(argv)
    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        descriptor = inherited_descriptor(SELF_DESCRIPTOR_ENV, "generator")
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError("generator descriptor marker is not attached to this execution")
        require_regular_profile(
            os.fstat(descriptor), "authenticated generator descriptor",
            expected_mode=0o755,
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("authenticated generator descriptor checksum changed")
        python_descriptor = inherited_descriptor(
            PYTHON_DESCRIPTOR_ENV, "Python interpreter"
        )
        require_isolated_python_reexec(python_descriptor)
        source = inherited_reexec_path(SELF_SOURCE_ENV, "source path")
        repository = inherited_reexec_path(REPOSITORY_ROOT_ENV, "repository root")
        require_reexec_source_relationship(source, repository)
        require_directory(repository, "authenticated generator repository root")
        with absolute_descriptor(
            source, "authenticated generator source pathname", flags=os.O_RDONLY
        ) as named_source:
            require_regular_profile(
                os.fstat(named_source),
                "authenticated generator source pathname",
                expected_mode=0o755,
            )
            if not same_inode(os.fstat(named_source), os.fstat(descriptor)):
                raise ValueError(
                    "authenticated generator source pathname is not the inherited descriptor"
                )
            if sha256_descriptor(named_source) != expected:
                raise ValueError("authenticated generator source pathname checksum changed")
        return source, repository, expected
    orphaned = [
        name
        for name in (PYTHON_DESCRIPTOR_ENV, SELF_SOURCE_ENV, REPOSITORY_ROOT_ENV)
        if name in os.environ
    ]
    if orphaned:
        raise ValueError(
            "private generator reexecution path is forbidden without an authenticated "
            "descriptor"
        )
    source = initial_source_path()
    repository = initial_repository_root(source)
    require_reexec_source_relationship(source, repository)
    with absolute_descriptor(source, "retained generator", flags=os.O_RDONLY) as descriptor:
        require_regular_profile(
            os.fstat(descriptor), "retained generator", expected_mode=0o755
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("retained generator checksum differs")
        python = Path("/proc/self/exe").resolve(strict=True)
        with absolute_descriptor(
            python, "authenticated Python interpreter", flags=os.O_RDONLY
        ) as python_descriptor:
            require_system_executable_profile(
                os.fstat(python_descriptor), "authenticated Python interpreter"
            )
            if not same_inode(
                os.fstat(python_descriptor), os.stat("/proc/self/exe", follow_symlinks=True)
            ):
                raise ValueError("authenticated Python pathname differs from running interpreter")
            os.set_inheritable(descriptor, True)
            os.set_inheritable(python_descriptor, True)
            os.execve(
                f"/proc/self/fd/{python_descriptor}",
                [
                    str(python),
                    "-I",
                    "-S",
                    "-B",
                    f"/proc/self/fd/{descriptor}",
                    *argv,
                ],
                reexec_environment(descriptor, python_descriptor, source, repository),
            )
    raise AssertionError("descriptor re-execution unexpectedly returned")


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
    }


def hardened_child_environment(*forbidden_prefixes: str) -> dict[str, str]:
    """Strip private reexec, loader, interpreter, Git, and caller PATH controls."""

    private = {
        SELF_DESCRIPTOR_ENV,
        PYTHON_DESCRIPTOR_ENV,
        SELF_SOURCE_ENV,
        REPOSITORY_ROOT_ENV,
    }
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in private
        and not key.startswith("LD_")
        and not key.startswith("PYTHON")
        and not key.startswith("GIT_")
        and not any(key.startswith(prefix) for prefix in forbidden_prefixes)
    }
    environment.update({"LC_ALL": "C", "PATH": TRUSTED_SYSTEM_PATH})
    return environment


def hardened_git_arguments(repository: Path, arguments: list[str]) -> list[str]:
    """Pin repository mode/worktree and disable local-config execution surfaces."""

    require_directory(repository, "Git repository")
    metadata = repository / ".git"
    retained = list(arguments)
    initializing = bool(retained) and retained[0] == "init"
    if os.path.lexists(metadata):
        require_directory(metadata, "Git metadata directory")
        repository_location = [
            f"--git-dir={metadata}",
            f"--work-tree={repository}",
        ]
    else:
        repository_location = [f"--git-dir={repository}"]
    if initializing:
        repository_location = []
    if retained and retained[0] == "diff":
        retained[1:1] = ["--no-ext-diff", "--no-textconv"]
    return [
        *repository_location,
        "-c",
        "core.fsmonitor=false",
        "-c",
        f"core.hooksPath={os.devnull}",
        "-c",
        f"core.attributesFile={os.devnull}",
        "-c",
        f"core.excludesFile={os.devnull}",
        "-c",
        f"init.templateDir={os.devnull}",
        "-c",
        "init.defaultObjectFormat=sha1",
        *retained,
    ]


def git_run(repository: Path, arguments: list[str], *,
            capture_output: bool = False,
            pass_fds: tuple[int, ...] = ()) -> subprocess.CompletedProcess:
    """Run one exact authenticated Git descriptor with isolated configuration."""

    with absolute_descriptor(GIT, "authenticated absolute Git binary", flags=os.O_RDONLY) as git:
        require_system_executable_profile(
            os.fstat(git), "authenticated absolute Git binary"
        )
        try:
            return subprocess.run(
                [
                    str(GIT),
                    "--no-replace-objects",
                    "-C",
                    str(repository),
                    *hardened_git_arguments(repository, arguments),
                ],
                executable=f"/proc/self/fd/{git}",
                check=False,
                stdin=subprocess.DEVNULL,
                capture_output=capture_output,
                env=hardened_git_environment(),
                pass_fds=(git, *pass_fds),
                timeout=120,
            )
        except (OSError, subprocess.TimeoutExpired) as error:
            raise ValueError("authenticated absolute Git query failed to execute") from error


def require_committed_file(repository: Path, path: Path, expected_sha256: str,
                           label: str, *, revision: str | None = None,
                           expected_mode: int | None = None) -> str:
    """Require a live file to be tracked, committed, and revision-bound."""

    try:
        relative = path.relative_to(repository).as_posix()
    except ValueError as error:
        raise ValueError(f"{label} is outside the source repository") from error
    if git_run(
        repository,
        ["ls-files", "--error-unmatch", "--", relative],
        capture_output=True,
    ).returncode:
        raise ValueError(f"{label} must be tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(repository, arguments).returncode:
            raise ValueError(f"{label} must be committed before recost generation")
    retained = git_run(repository, ["show", f"HEAD:{relative}"], capture_output=True)
    if retained.returncode or sha256_bytes(retained.stdout) != expected_sha256:
        raise ValueError(f"{label} committed HEAD bytes differ from live bytes")
    if revision is not None:
        revision = require_revision(revision, f"{label} revision")
        retained = git_run(
            repository, ["show", f"{revision}:{relative}"], capture_output=True
        )
        if retained.returncode or sha256_bytes(retained.stdout) != expected_sha256:
            raise ValueError(f"{label} requested revision bytes differ")
    with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
        require_regular_profile(
            os.fstat(descriptor), label, expected_mode=expected_mode
        )
        if sha256_descriptor(descriptor) != expected_sha256:
            raise ValueError(f"{label} checksum differs")
    head = git_run(repository, ["rev-parse", "HEAD"], capture_output=True)
    if head.returncode:
        raise ValueError("cannot resolve repository HEAD")
    return head.stdout.decode().strip()


def require_revision_file(repository: Path, relative: PurePosixPath, revision: str,
                          expected_sha256: str, label: str) -> None:
    """Require exact retained bytes at one authenticated historical revision."""

    revision = require_revision(revision, f"{label} revision")
    expected = require_sha256(expected_sha256, f"{label} SHA-256")
    retained = git_run(
        repository, ["show", f"{revision}:{relative.as_posix()}"], capture_output=True
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError(f"{label} historical revision bytes differ")


def authenticate_f116_committed_tools(
    value: object, repository: Path, expected_head: str
) -> list[dict[str, object]]:
    """Authenticate the exact seven committed tools published by F116."""

    if not isinstance(value, list) or len(value) != len(F116_REQUIRED_COMMITTED_TOOLS):
        raise ValueError("F116 committed_tools must contain exactly seven tools")
    retained = []
    for index, (relative, expected_mode) in enumerate(
        sorted(F116_REQUIRED_COMMITTED_TOOLS.items())
    ):
        record = require_exact_keys(
            value[index],
            {"path", "revision", "sha256", "mode"},
            f"F116 committed tool {index}",
        )
        digest = require_sha256(record["sha256"], f"F116 committed tool {relative} SHA-256")
        if record != {
            "path": relative,
            "revision": expected_head,
            "sha256": digest,
            "mode": expected_mode,
        }:
            raise ValueError("F116 committed_tools identity, order, revision, or mode differs")
        require_revision_file(
            repository,
            PurePosixPath(relative),
            expected_head,
            digest,
            f"F116 committed tool {relative}",
        )
        tree = git_run(
            repository,
            ["ls-tree", expected_head, "--", relative],
            capture_output=True,
        )
        git_mode = "100755" if expected_mode == "0755" else "100644"
        try:
            tree_line = tree.stdout.decode("ascii").strip()
        except UnicodeDecodeError as error:
            raise ValueError(f"F116 committed tool {relative} mode is not ASCII") from error
        if (
            tree.returncode
            or not tree_line.startswith(f"{git_mode} blob ")
            or not tree_line.endswith(f"\t{relative}")
        ):
            raise ValueError(f"F116 committed tool {relative} Git mode differs")
        retained.append(dict(record))
    return retained


def require_relative_path(value: object, label: str) -> PurePosixPath:
    """Require one normalized relative POSIX path."""

    retained = require_nonempty_string(value, label)
    path = PurePosixPath(retained)
    if path.is_absolute() or path.as_posix() != retained or ".." in path.parts:
        raise ValueError(f"{label} must be a normalized relative path")
    return path


def require_no_symlink_components(path: Path, label: str, *,
                                  include_leaf: bool = True) -> None:
    """Reject symlinks in every existing component of one absolute path."""

    if not path.is_absolute():
        raise ValueError(f"{label} path must be absolute")
    parts = path.parts
    current = Path(parts[0])
    limit = len(parts) if include_leaf else len(parts) - 1
    for part in parts[1:limit]:
        current = current / part
        try:
            profile = os.lstat(current)
        except FileNotFoundError:
            continue
        if stat.S_ISLNK(profile.st_mode):
            raise ValueError(f"{label} path contains a symlink: {current}")


@contextmanager
def absolute_descriptor(path: Path, label: str, *, flags: int) -> Iterator[int]:
    """Open one absolute path without accepting symlinks in any component."""

    path = path.absolute()
    if (
        not path.is_absolute()
        or len(path.parts) < 2
        or path != Path(os.path.normpath(str(path)))
        or ".." in path.parts
    ):
        raise ValueError(f"{label} path is invalid: {path}")
    parent = os.open("/", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [parent]
    try:
        for component in path.parts[1:-1]:
            parent = os.open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            opened.append(parent)
        descriptor = os.open(path.parts[-1], flags | os.O_NOFOLLOW, dir_fd=parent)
        opened.append(descriptor)
        yield descriptor
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)


def root_path(root: Path, value: object, label: str) -> Path:
    """Return one normalized root-relative path without symlink traversal."""

    relative = require_relative_path(value, label)
    path = root.joinpath(*relative.parts)
    require_no_symlink_components(path, label)
    return path


def repository_path(repository: Path, value: object, label: str) -> Path:
    """Return one normalized repository-relative path."""

    relative = require_relative_path(value, label)
    path = repository.joinpath(*relative.parts)
    require_no_symlink_components(path, label)
    return path


@dataclass(frozen=True)
class TrackedInput:
    """One immutable input to reauthenticate before output creation."""

    path: Path
    sha256: str
    label: str
    expected_mode: int | None


@dataclass(frozen=True)
class DirectoryMeasurement:
    """One exact directory-size observation retained for final revalidation."""

    path: Path
    size_bytes: int
    label: str


class InputTracker:
    """Authenticate inputs and recheck every retained byte before publication."""

    def __init__(self) -> None:
        self._inputs: dict[Path, TrackedInput] = {}
        self._directory_measurements: dict[Path, DirectoryMeasurement] = {}

    def read(self, path: Path, expected_sha256: object, label: str, *,
             expected_mode: int | None = None) -> bytes:
        expected = require_sha256(expected_sha256, f"{label} SHA-256")
        require_no_symlink_components(path, label)
        with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
            require_regular_profile(
                os.fstat(descriptor), label, expected_mode=expected_mode
            )
            chunks = []
            digest = hashlib.sha256()
            while True:
                block = os.read(descriptor, 1024 * 1024)
                if not block:
                    break
                digest.update(block)
                chunks.append(block)
            if digest.hexdigest() != expected:
                raise ValueError(f"{label} checksum differs: {path}")
        retained = TrackedInput(path, expected, label, expected_mode)
        existing = self._inputs.get(path)
        if existing is not None and (
            existing.sha256 != retained.sha256
            or existing.expected_mode != retained.expected_mode
        ):
            raise ValueError(f"{label} has conflicting retained bindings")
        if existing is None:
            self._inputs[path] = retained
        return b"".join(chunks)

    def authenticate(self, path: Path, expected_sha256: object, label: str, *,
                     expected_mode: int | None = None) -> None:
        expected = require_sha256(expected_sha256, f"{label} SHA-256")
        require_no_symlink_components(path, label)
        with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
            require_regular_profile(
                os.fstat(descriptor), label, expected_mode=expected_mode
            )
            if sha256_descriptor(descriptor) != expected:
                raise ValueError(f"{label} checksum differs: {path}")
        retained = TrackedInput(path, expected, label, expected_mode)
        existing = self._inputs.get(path)
        if existing is not None and (
            existing.sha256 != retained.sha256
            or existing.expected_mode != retained.expected_mode
        ):
            raise ValueError(f"{label} has conflicting retained bindings")
        if existing is None:
            self._inputs[path] = retained

    def discover(self, path: Path, label: str, *, expected_mode: int | None = None
                 ) -> tuple[bytes, str]:
        """Read and bind one discovered retained file by its observed digest."""

        require_no_symlink_components(path, label)
        with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
            require_regular_profile(
                os.fstat(descriptor), label, expected_mode=expected_mode
            )
            chunks = []
            digest = hashlib.sha256()
            while True:
                block = os.read(descriptor, 1024 * 1024)
                if not block:
                    break
                chunks.append(block)
                digest.update(block)
        retained_digest = digest.hexdigest()
        retained = TrackedInput(path, retained_digest, label, expected_mode)
        existing = self._inputs.get(path)
        if existing is not None and (
            existing.sha256 != retained.sha256
            or existing.expected_mode != retained.expected_mode
        ):
            raise ValueError(f"{label} has conflicting retained bindings")
        if existing is None:
            self._inputs[path] = retained
        return b"".join(chunks), retained_digest

    def reauthenticate_all(self) -> None:
        """Recheck every bound path immediately before output creation."""

        retained = list(self._inputs.values())
        self._inputs = {}
        for item in retained:
            self.authenticate(
                item.path,
                item.sha256,
                item.label,
                expected_mode=item.expected_mode,
            )
        require_directory_measurement_boundaries(
            tuple(self._directory_measurements.values())
        )

    def retain_directory_measurements(
        self, measurements: tuple[DirectoryMeasurement, ...]
    ) -> None:
        """Retain exact directory-size observations for final revalidation."""

        for measurement in measurements:
            existing = self._directory_measurements.get(measurement.path)
            if existing is not None and existing != measurement:
                raise ValueError(
                    f"{measurement.label} has conflicting retained measurement"
                )
            self._directory_measurements[measurement.path] = measurement


@dataclass(frozen=True)
class BuildResult:
    """Retain deterministic bytes and final live-revalidation bindings."""

    payload: bytes
    tracker: InputTracker
    reconcile: dict[str, object]
    helper_path: Path
    helper_sha256: str
    helper_revision: str
    matrix_path: Path
    matrix_sha256: str
    matrix_revision: str
    storage_available_bytes: int
    storage_retained_stage_i_bytes: int
    storage_required_safety_bytes: int
    projected_storage_bytes: int
    storage_measurements: tuple[DirectoryMeasurement, ...]


def parse_json(payload: bytes, label: str) -> object:
    """Decode one retained JSON input."""

    try:
        return json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error


def input_binding(value: object, label: str) -> tuple[PurePosixPath, str]:
    """Parse one exact root-relative retained-file binding."""

    retained = require_exact_keys(value, {"path", "sha256"}, label)
    return (
        require_relative_path(retained["path"], f"{label} path"),
        require_sha256(retained["sha256"], f"{label} SHA-256"),
    )


def declared_file_binding(value: object, label: str) -> dict[str, object]:
    """Parse one exact published path/hash/mode/link binding."""

    retained = require_exact_keys(value, {"path", "sha256", "mode", "links"}, label)
    return {
        "path": require_nonempty_string(retained["path"], f"{label} path"),
        "sha256": require_sha256(retained["sha256"], f"{label} SHA-256"),
        "mode": require_nonempty_string(retained["mode"], f"{label} mode"),
        "links": require_integer(retained["links"], f"{label} links", minimum=1),
    }


def require_declared_binding(
    value: object,
    path: Path,
    digest: str,
    label: str,
    *,
    mode: str = "0444",
) -> None:
    """Require one declared publication binding to exact retained bytes."""

    expected = {"path": str(path), "sha256": digest, "mode": mode, "links": 1}
    if declared_file_binding(value, label) != expected:
        raise ValueError(f"{label} differs from exact retained publication")


def parse_request_independent_review(
    value: object,
    request_path: Path,
    request_sha256: str,
    request_timestamp: datetime,
    request_author: str,
) -> dict[str, object]:
    """Require independent approval of the exact non-authorizing recost request."""

    retained = require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "reviewed_utc",
            "decision",
            "reviewer",
            "candidate",
            "scope",
        },
        "recost request independent review",
    )
    if (
        retained["schema_version"] != 1
        or retained["record_type"] != "stage-i-recost-request-independent-review"
        or retained["execution_epoch"] != EXECUTION_EPOCH
        or retained["decision"] != "approved-for-evidence-generation"
    ):
        raise ValueError("recost request independent review identity differs")
    reviewed = parse_utc_timestamp(
        retained["reviewed_utc"], "recost request independent review timestamp"
    )
    if reviewed < request_timestamp or reviewed > datetime.now(timezone.utc) + timedelta(minutes=5):
        raise ValueError("recost request independent review chronology differs")
    reviewer = require_exact_keys(
        retained["reviewer"],
        {
            "agent_id",
            "role",
            "declared_process_independence",
            "identity_assurance",
            "identity_assurance_limitation",
        },
        "recost request reviewer",
    )
    reviewer_agent = require_nonempty_string(
        reviewer["agent_id"], "recost request reviewer agent ID"
    )
    if (
        reviewer["role"] != "independent recost request reviewer"
        or reviewer["declared_process_independence"] is not True
        or reviewer["identity_assurance"] != REQUEST_REVIEW_IDENTITY_ASSURANCE
        or reviewer["identity_assurance_limitation"] != REQUEST_REVIEW_IDENTITY_LIMITATION
        or reviewer_agent == request_author
    ):
        raise ValueError(
            "recost request reviewer lacks the exact declared non-cryptographic "
            "process-independence assurance"
        )
    if retained["candidate"] != {
        "path": str(request_path),
        "sha256": request_sha256,
    }:
        raise ValueError("recost request independent review candidate binding differs")
    if retained["scope"] != {"non_authorizing": True}:
        raise ValueError("recost request independent review over-authorizes")
    return retained


def parse_qualification_approval(
    value: object,
    root: Path,
    executable: Path,
    executable_revision: str,
    executable_sha256: str,
    build_manifest: Path,
) -> dict[str, object]:
    """Authenticate the promoted corrected-build qualification approval."""

    retained = require_exact_keys(
        value,
        {
            "schema_version",
            "execution_epoch",
            "approval_scope",
            "approved_by",
            "approved_utc",
            "approved_executable",
            "approved_executable_revision",
            "approved_executable_sha256",
            "build_manifest",
            "review_notes",
        },
        "qualification approval",
    )
    if (
        retained["schema_version"] != 1
        or retained["execution_epoch"] != EXECUTION_EPOCH
        or retained["approved_executable"] != str(executable)
        or retained["approved_executable_revision"] != executable_revision
        or retained["approved_executable_sha256"] != executable_sha256
        or retained["build_manifest"] != str(build_manifest)
    ):
        raise ValueError("qualification approval build binding differs")
    require_nonempty_string(retained["approval_scope"], "qualification approval scope")
    require_nonempty_string(retained["approved_by"], "qualification approver")
    require_nonempty_string(retained["review_notes"], "qualification review notes")
    parse_utc_timestamp(retained["approved_utc"], "qualification approval timestamp")
    if executable != root / executable.relative_to(root):
        raise ValueError("qualification approval executable is outside the campaign root")
    return retained


def parse_historical_f115_source_authority(
    root: Path,
    bindings: object,
    ceiling_path: Path,
    ceiling_sha256: str,
    tracker: InputTracker,
) -> dict[str, object]:
    """Authenticate immutable historical F115 without selecting its retired bundle."""

    retained = require_exact_keys(
        bindings,
        {"evidence", "publication_audit", "provenance_review", "plasma_review"},
        "F115 source authority bindings",
    )
    paths = {
        "evidence": root / F115_RELATIVE,
        "publication_audit": root / F115_PUBLICATION_AUDIT_RELATIVE,
        "provenance_review": root / F115_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / F115_PLASMA_REVIEW_RELATIVE,
    }
    loaded: dict[str, tuple[dict[str, object], str]] = {}
    for key, path in paths.items():
        relative, expected = input_binding(retained[key], f"F115 {key} binding")
        if root_path(root, relative.as_posix(), f"F115 {key}") != path:
            raise ValueError(f"F115 {key} is not the exact retained path")
        if root == DEFAULT_ROOT and expected != F115_CANONICAL_SHA256[key]:
            raise ValueError(f"F115 {key} differs from canonical historical authority")
        value = parse_json(
            tracker.read(path, expected, f"F115 {key}", expected_mode=0o444),
            f"F115 {key}",
        )
        if not isinstance(value, dict):
            raise ValueError(f"F115 {key} must be an object")
        loaded[key] = (value, expected)
    evidence, evidence_sha256 = loaded["evidence"]
    implementation = evidence.get("implementation")
    predecessors = evidence.get("predecessors")
    if (
        evidence.get("schema_version") != 1
        or evidence.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-evidence"
        or evidence.get("checkpoint") != "F-115"
        or evidence.get("execution_epoch") != EXECUTION_EPOCH
        or not isinstance(implementation, dict)
        or not isinstance(predecessors, dict)
    ):
        raise ValueError("F115 source authority identity differs")
    bundle = implementation.get("source_bundle")
    if not isinstance(bundle, dict) or bundle.get("complete_history") is not True:
        raise ValueError("F115 historical source bundle identity differs")
    historical_bundle = root_path(
        root,
        require_relative_path(bundle.get("path"), "F115 historical source bundle path").as_posix(),
        "F115 historical source bundle",
    )
    if historical_bundle.parent != root / "source-archives":
        raise ValueError("F115 historical source bundle is outside source-archives")
    require_sha256(bundle.get("sha256"), "F115 historical source bundle SHA-256")
    f113 = predecessors.get("f113_controller_transition")
    if not isinstance(f113, dict) or {
        "path": f113.get("path"),
        "sha256": f113.get("sha256"),
    } != {
        "path": ceiling_path.relative_to(root).as_posix(),
        "sha256": ceiling_sha256,
    }:
        raise ValueError("F115 source authority F113 predecessor binding differs")

    audit, audit_sha256 = loaded["publication_audit"]
    if (
        audit.get("schema_version") != 1
        or audit.get("record_type")
        != "stage-i-source-bundle-recovery-supersession-publication-audit"
        or audit.get("checkpoint") != "F-115"
        or audit.get("execution_epoch") != EXECUTION_EPOCH
    ):
        raise ValueError("F115 publication audit identity differs")
    require_declared_binding(
        audit.get("artifact"), paths["evidence"], evidence_sha256, "F115 publication artifact"
    )
    reviews = audit.get("independent_reviews")
    if not isinstance(reviews, dict) or reviews.get(
        "reviews_bind_exact_published_f115_sha256"
    ) != evidence_sha256:
        raise ValueError("F115 publication audit independent-review binding differs")
    reviewers: set[str] = set()
    for key, review_kind, decision in (
        ("provenance_review", "provenance-security", "approved-for-publication"),
        ("plasma_review", "plasma-scientific-continuation", "approved"),
    ):
        review, review_sha256 = loaded[key]
        audit_key = "provenance_security" if key == "provenance_review" else (
            "plasma_scientific_continuation"
        )
        require_declared_binding(
            reviews.get(audit_key), paths[key], review_sha256, f"F115 {key}"
        )
        published = review.get("published_f115")
        reviewer = review.get("reviewer")
        if (
            review.get("schema_version") != 1
            or review.get("record_type")
            != "stage-i-source-bundle-recovery-supersession-independent-review"
            or review.get("checkpoint") != "F-115"
            or review.get("execution_epoch") != EXECUTION_EPOCH
            or review.get("review_kind") != review_kind
            or review.get("decision") != decision
            or published
            != {"path": str(paths["evidence"]), "sha256": evidence_sha256}
            or not isinstance(reviewer, dict)
        ):
            raise ValueError(f"F115 {key} identity differs")
        agent_id = require_nonempty_string(reviewer.get("agent_id"), f"F115 {key} reviewer")
        if agent_id in reviewers:
            raise ValueError("F115 independent reviews do not have distinct reviewers")
        reviewers.add(agent_id)
    authority = audit.get("authority_and_enforcement")
    if not isinstance(authority, dict) or authority.get("direct_sbatch_authorized") is not False:
        raise ValueError("F115 publication audit over-authorizes scheduler action")
    return {
        "evidence_sha256": evidence_sha256,
        "publication_audit_sha256": audit_sha256,
        "provenance_review_sha256": loaded["provenance_review"][1],
        "plasma_review_sha256": loaded["plasma_review"][1],
    }


def parse_current_source_authority(
    root: Path,
    repository: Path,
    bindings: object,
    source_bundle_path: Path,
    source_bundle_sha256: str,
    source_bundle_revisions: list[str],
    ceiling_path: Path,
    ceiling_sha256: str,
    tracker: InputTracker,
) -> tuple[dict[str, object], dict[str, object] | None, str]:
    """Authenticate F116 current-source authority and its historical F115 chain."""

    retained = require_exact_keys(
        bindings,
        {
            "checkpoint",
            "evidence",
            "provenance_review",
            "plasma_review",
            "publication_audit",
            "final_source_bundle",
        },
        "F116 current source authority bindings",
    )
    if retained["checkpoint"] != "F-116":
        raise ValueError("current source authority binding is not F-116")
    final_relative, final_sha256, final_revisions = source_bundle_binding(
        retained["final_source_bundle"], "F116 final source bundle binding"
    )
    if (
        root_path(root, final_relative.as_posix(), "F116 final source bundle")
        != source_bundle_path
        or final_sha256 != source_bundle_sha256
        or final_revisions != source_bundle_revisions
    ):
        raise ValueError("F116 six-part binding does not select the current source bundle")
    paths = {
        "evidence": root / F116_RELATIVE,
        "publication_audit": root / F116_PUBLICATION_AUDIT_RELATIVE,
        "provenance_review": root / F116_PROVENANCE_REVIEW_RELATIVE,
        "plasma_review": root / F116_PLASMA_REVIEW_RELATIVE,
    }
    loaded: dict[str, tuple[dict[str, object], str]] = {}
    for key, path in paths.items():
        relative, expected = input_binding(retained[key], f"F116 {key} binding")
        if root_path(root, relative.as_posix(), f"F116 {key}") != path:
            raise ValueError(f"F116 {key} is not the exact retained path")
        value = parse_json(
            tracker.read(path, expected, f"F116 {key}", expected_mode=0o444),
            f"F116 {key}",
        )
        if not isinstance(value, dict):
            raise ValueError(f"F116 {key} must be an object")
        loaded[key] = (value, expected)

    evidence, evidence_sha256 = loaded["evidence"]
    retained_evidence = require_exact_keys(
        evidence,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "generated_utc",
            "scope",
            "predecessor_authorities",
            "implementation",
            "source_archive_catalog",
            "authorization",
            "validation",
            "publication_requirements",
        },
        "F116 evidence",
    )
    generated = parse_utc_timestamp(retained_evidence["generated_utc"], "F116 generation")
    if (
        retained_evidence["schema_version"] != 1
        or retained_evidence["record_type"]
        != "stage-i-current-source-authority-supersession-evidence"
        or retained_evidence["checkpoint"] != "F-116"
        or retained_evidence["execution_epoch"] != EXECUTION_EPOCH
        or retained_evidence["authorization"] != F116_AUTHORIZATION
    ):
        raise ValueError("F116 current source authority identity differs")
    scope = require_exact_keys(
        retained_evidence["scope"],
        {"relationship", "summary", "preserves", "does_not_authorize"},
        "F116 scope",
    )
    if (
        scope["relationship"] != "current-source-selection-only-supersession"
        or not require_nonempty_string(scope["summary"], "F116 scope summary")
        or not isinstance(scope["preserves"], list)
        or not isinstance(scope["does_not_authorize"], list)
        or any(
            item not in scope["does_not_authorize"]
            for item in ("prepare", "submit", "direct sbatch", "scheduler mutation")
        )
    ):
        raise ValueError("F116 scope differs or broadens authority")
    predecessors = require_exact_keys(
        retained_evidence["predecessor_authorities"],
        {"historical_f115"},
        "F116 predecessor authorities",
    )
    historical_f115 = parse_historical_f115_source_authority(
        root,
        predecessors["historical_f115"],
        ceiling_path,
        ceiling_sha256,
        tracker,
    )
    implementation = require_exact_keys(
        retained_evidence["implementation"],
        {
            "publisher",
            "committed_tools",
            "intermediate_36140_bundle",
            "current_source_bundle",
        },
        "F116 implementation",
    )
    current_bundle = require_exact_keys(
        implementation["current_source_bundle"],
        {
            "candidate_path",
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "subject",
        },
        "F116 current source bundle",
    )
    current_head = require_revision(current_bundle["head"], "F116 current source bundle head")
    if (
        current_bundle["path"] != source_bundle_path.relative_to(root).as_posix()
        or current_bundle["sha256"] != source_bundle_sha256
        or current_bundle["complete_history"] is not True
        or current_bundle["selected_as_current"] is not True
        or current_bundle["verified_revisions"] != source_bundle_revisions
        or current_head not in source_bundle_revisions
    ):
        raise ValueError("F116 current source authority bundle binding differs")
    authenticate_f116_committed_tools(
        implementation["committed_tools"], repository, current_head
    )
    require_nonempty_string(current_bundle["subject"], "F116 current source bundle subject")
    require_exact_keys(
        current_bundle["advertised_tip"], {"revision", "name"}, "F116 advertised tip"
    )
    bridge = require_exact_keys(
        implementation["intermediate_36140_bundle"],
        {
            "path",
            "sha256",
            "complete_history",
            "head",
            "advertised_tip",
            "verified_revisions",
            "selected_as_current",
            "role",
        },
        "F116 bridge bundle",
    )
    if (
        bridge["complete_history"] is not True
        or bridge["selected_as_current"] is not False
        or bridge["role"] != "retained-non-current-bridge"
    ):
        raise ValueError("F116 bridge bundle authority differs")

    verified = {
        "authorization_broadening": False,
        "bridge_selected_as_current": False,
        "corrupt_c7_excluded": True,
        "current_source_selection_only": True,
        "final_bundle_sha256": source_bundle_sha256,
        "final_head": current_head,
        "historical_f115_preserved": True,
    }
    reviewers: set[str] = set()
    reviewed_times: list[datetime] = []
    for key, review_kind, decision in (
        ("provenance_review", "provenance-security", "approved-for-publication"),
        ("plasma_review", "plasma-scientific-continuation", "approved"),
    ):
        review, _ = loaded[key]
        retained_review = require_exact_keys(
            review,
            {
                "schema_version",
                "record_type",
                "checkpoint",
                "execution_epoch",
                "review_kind",
                "decision",
                "reviewed_candidate",
                "published_f116",
                "reviewer",
                "reviewed_utc",
                "findings",
                "limitations",
                "verified",
            },
            f"F116 {key}",
        )
        candidate = require_exact_keys(
            retained_review["reviewed_candidate"],
            {"path", "sha256"},
            f"F116 {key} reviewed candidate",
        )
        require_nonempty_string(candidate["path"], f"F116 {key} reviewed candidate path")
        reviewer = require_exact_keys(
            retained_review["reviewer"], {"agent_id", "identity"}, f"F116 {key} reviewer"
        )
        if (
            retained_review["schema_version"] != 1
            or retained_review["record_type"]
            != "stage-i-current-source-authority-supersession-independent-review"
            or retained_review["checkpoint"] != "F-116"
            or retained_review["execution_epoch"] != EXECUTION_EPOCH
            or retained_review["review_kind"] != review_kind
            or retained_review["decision"] != decision
            or candidate["sha256"] != evidence_sha256
            or retained_review["published_f116"]
            != {"path": str(paths["evidence"]), "sha256": evidence_sha256}
            or retained_review["verified"] != verified
        ):
            raise ValueError(f"F116 {key} identity differs")
        agent_id = require_nonempty_string(reviewer["agent_id"], f"F116 {key} reviewer")
        require_nonempty_string(reviewer["identity"], f"F116 {key} reviewer identity")
        if agent_id in reviewers:
            raise ValueError("F116 independent reviews do not have distinct reviewers")
        reviewers.add(agent_id)
        for field in ("findings", "limitations"):
            if (
                not isinstance(retained_review[field], list)
                or not retained_review[field]
                or any(not isinstance(item, str) or not item for item in retained_review[field])
            ):
                raise ValueError(f"F116 {key} {field} differ")
        reviewed = parse_utc_timestamp(retained_review["reviewed_utc"], f"F116 {key} review")
        if reviewed < generated:
            raise ValueError(f"F116 {key} review predates evidence")
        reviewed_times.append(reviewed)

    audit, audit_sha256 = loaded["publication_audit"]
    retained_audit = require_exact_keys(
        audit,
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "execution_epoch",
            "published_utc",
            "artifact",
            "independent_reviews",
            "historical_f115_authority",
            "source_archive_catalog",
            "authority_and_enforcement",
            "publication",
        },
        "F116 publication audit",
    )
    published = parse_utc_timestamp(retained_audit["published_utc"], "F116 publication")
    if (
        retained_audit["schema_version"] != 1
        or retained_audit["record_type"]
        != "stage-i-current-source-authority-supersession-publication-audit"
        or retained_audit["checkpoint"] != "F-116"
        or retained_audit["execution_epoch"] != EXECUTION_EPOCH
        or retained_audit["authority_and_enforcement"] != F116_AUTHORIZATION
        or retained_audit["historical_f115_authority"] != historical_f115
        or retained_audit["publication"]
        != "recoverable-forward-transaction-with-publication-audit-commit-marker-under-stage-i-lock"
        or published < generated
        or any(published < reviewed for reviewed in reviewed_times)
    ):
        raise ValueError("F116 publication audit identity or authority differs")
    require_declared_binding(
        retained_audit["artifact"], paths["evidence"], evidence_sha256, "F116 publication artifact"
    )
    reviews = require_exact_keys(
        retained_audit["independent_reviews"],
        {
            "reviews_bind_exact_published_f116_sha256",
            "provenance_security",
            "plasma_scientific_continuation",
        },
        "F116 publication reviews",
    )
    if reviews["reviews_bind_exact_published_f116_sha256"] != evidence_sha256:
        raise ValueError("F116 publication audit independent-review binding differs")
    require_declared_binding(
        reviews["provenance_security"],
        paths["provenance_review"],
        loaded["provenance_review"][1],
        "F116 provenance review",
    )
    require_declared_binding(
        reviews["plasma_scientific_continuation"],
        paths["plasma_review"],
        loaded["plasma_review"][1],
        "F116 plasma review",
    )
    catalog = require_exact_keys(
        retained_audit["source_archive_catalog"],
        {
            "readme",
            "sha256sums",
            "bridge_bundle",
            "current_source_bundle",
            "corrupt_c7_absent_from_active_checksum_ledger",
            "sole_current_source_bundle",
        },
        "F116 publication source-archive catalog",
    )
    for key, relative in (
        ("readme", Path("source-archives/README.md")),
        ("sha256sums", Path("source-archives/SHA256SUMS")),
    ):
        binding = declared_file_binding(catalog[key], f"F116 catalog {key}")
        path = root / relative
        if binding["path"] != str(path) or binding["mode"] != "0644" or binding["links"] != 1:
            raise ValueError(f"F116 catalog {key} binding differs")
        tracker.authenticate(path, binding["sha256"], f"F116 catalog {key}", expected_mode=0o644)
    published_current = require_exact_keys(
        catalog["current_source_bundle"],
        {"path", "sha256", "mode", "links", "head", "selected_as_current"},
        "F116 published current source bundle",
    )
    if published_current != {
        "path": str(source_bundle_path),
        "sha256": source_bundle_sha256,
        "mode": "0644",
        "links": 1,
        "head": current_head,
        "selected_as_current": True,
    } or catalog["sole_current_source_bundle"] != str(source_bundle_path):
        raise ValueError("F116 published current source bundle binding differs")
    published_bridge = require_exact_keys(
        catalog["bridge_bundle"],
        {"path", "sha256", "mode", "links", "head", "role", "selected_as_current"},
        "F116 published bridge bundle",
    )
    if (
        published_bridge["selected_as_current"] is not False
        or published_bridge["role"] != "retained-non-current-bridge"
        or catalog["corrupt_c7_absent_from_active_checksum_ledger"] is not True
    ):
        raise ValueError("F116 source-archive catalog authority differs")

    historical_audit_path = root / F113_PUBLICATION_AUDIT_RELATIVE
    adoption = None
    if not os.path.lexists(historical_audit_path):
        adoption = {
            "status": "controlled-exact-transitive-F115-F116-supersession",
            "historical_publication_audit_path": str(historical_audit_path),
            "historical_publication_audit_absent": True,
            "artifact": {
                "path": str(ceiling_path),
                "sha256": ceiling_sha256,
                "mode": "0644",
                "links": 1,
            },
            "f115_evidence_sha256": historical_f115["evidence_sha256"],
            "f116_publication_audit_sha256": audit_sha256,
        }

    binding = {
        "checkpoint": "F-116",
        "evidence": {
            "path": F116_RELATIVE.as_posix(),
            "sha256": evidence_sha256,
        },
        "provenance_review": {
            "path": F116_PROVENANCE_REVIEW_RELATIVE.as_posix(),
            "sha256": loaded["provenance_review"][1],
        },
        "plasma_review": {
            "path": F116_PLASMA_REVIEW_RELATIVE.as_posix(),
            "sha256": loaded["plasma_review"][1],
        },
        "publication_audit": {
            "path": F116_PUBLICATION_AUDIT_RELATIVE.as_posix(),
            "sha256": audit_sha256,
        },
        "final_source_bundle": {
            "path": source_bundle_path.relative_to(root).as_posix(),
            "sha256": source_bundle_sha256,
            "verified_revisions": source_bundle_revisions,
        },
    }
    if binding != retained:
        raise ValueError("F116 current source authority binding is not canonical")
    return binding, adoption, published.isoformat()


def repository_binding(value: object, label: str) -> tuple[PurePosixPath, str, str]:
    """Parse one exact repository-relative retained-file binding."""

    retained = require_exact_keys(value, {"path", "revision", "sha256"}, label)
    return (
        require_relative_path(retained["path"], f"{label} path"),
        require_revision(retained["revision"], f"{label} revision"),
        require_sha256(retained["sha256"], f"{label} SHA-256"),
    )


def source_bundle_binding(value: object, label: str) -> tuple[PurePosixPath, str, list[str]]:
    """Parse one exact source-bundle binding."""

    retained = require_exact_keys(
        value, {"path", "sha256", "verified_revisions"}, label
    )
    revisions = retained["verified_revisions"]
    if not isinstance(revisions, list) or not revisions:
        raise ValueError(f"{label} verified revisions must be a nonempty list")
    parsed = [require_revision(item, f"{label} verified revision") for item in revisions]
    if len(set(parsed)) != len(parsed):
        raise ValueError(f"{label} verified revisions are duplicated")
    return (
        require_relative_path(retained["path"], f"{label} path"),
        require_sha256(retained["sha256"], f"{label} SHA-256"),
        parsed,
    )


def require_directory(path: Path, label: str) -> None:
    """Require one owned directory with a trusted replacement boundary."""

    require_no_symlink_components(path, label)
    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        require_directory_profile(os.fstat(descriptor), label)


def require_directory_profile(profile: os.stat_result, label: str) -> None:
    """Require one owned directory profile no broader than 0755 plus setgid."""

    if not stat.S_ISDIR(profile.st_mode):
        raise ValueError(f"{label} must be a directory")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"{label} must be owned by the effective user")
    mode = stat.S_IMODE(profile.st_mode)
    allowed = 0o755 | stat.S_ISGID
    if mode & ~allowed:
        raise ValueError(f"{label} mode is {mode:04o}, exceeds trusted profile 0755")


def require_lock_profile(profile: os.stat_result, path: Path) -> None:
    """Require the retained Stage I lock profile."""

    if (
        not stat.S_ISREG(profile.st_mode)
        or stat.S_IMODE(profile.st_mode) != 0o644
        or profile.st_uid != os.geteuid()
        or profile.st_nlink != 1
        or profile.st_size != 0
    ):
        raise ValueError(f"Stage I lock profile differs: {path}")


@dataclass(frozen=True)
class MutationLock:
    """Retain the exact locked pathname and inode across every mutation."""

    path: Path
    descriptor: int
    device: int
    inode: int

    def authenticate(self) -> None:
        """Require the named lock and held descriptor to remain the retained inode."""

        opened = os.fstat(self.descriptor)
        require_lock_profile(opened, self.path)
        if (opened.st_dev, opened.st_ino) != (self.device, self.inode):
            raise ValueError("held Stage I lock descriptor identity changed")
        with absolute_descriptor(
            self.path, "Stage I lock", flags=os.O_RDWR
        ) as named_descriptor:
            named = os.fstat(named_descriptor)
        require_lock_profile(named, self.path)
        if not same_inode(named, opened):
            raise ValueError("Stage I lock pathname changed while mutation lock is held")


def authenticate_mutation_lock(lock: MutationLock | None) -> None:
    """Reauthenticate a production mutation lock when one is supplied."""

    if lock is not None:
        lock.authenticate()


@contextmanager
def stage_i_lock(root: Path) -> Iterator[MutationLock]:
    """Acquire the existing Stage I mutation lock without creating it."""

    path = root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock"
    require_no_symlink_components(path, "Stage I lock")
    with absolute_descriptor(path, "Stage I lock", flags=os.O_RDWR) as descriptor:
        acquired = False
        require_lock_profile(os.fstat(descriptor), path)
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            acquired = True
        except OSError as error:
            if error.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            raise ValueError(f"another Stage I mutation holds {path}") from error
        opened = os.fstat(descriptor)
        lock = MutationLock(path, descriptor, opened.st_dev, opened.st_ino)
        lock.authenticate()
        try:
            yield lock
        finally:
            try:
                lock.authenticate()
            finally:
                if acquired:
                    fcntl.flock(descriptor, fcntl.LOCK_UN)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse staged generation and its authenticated prerequisite actions."""

    parser = argparse.ArgumentParser(prog="cgl_lf_stage_i_recost.py")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--allow-local-root", action="store_true")
    parser.add_argument(
        "action",
        nargs="?",
        choices=(
            "generate",
            "install-f117-draft-packet",
            "draft-request",
            "install-request-review",
            "retain-generator",
        ),
        default="generate",
    )
    parser.add_argument("--request", type=Path)
    parser.add_argument("--expected-request-sha256", type=sha256_arg)
    parser.add_argument("--review", type=Path)
    parser.add_argument("--expected-review-sha256", type=sha256_arg)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--packet", type=Path)
    parser.add_argument("--expected-packet-sha256", type=sha256_arg)
    parser.add_argument("--expected-generator-sha256", type=sha256_arg, required=True)
    parser.add_argument(
        "--squeue-file",
        type=Path,
        help="local-fixture-only machine-readable queue evidence",
    )
    args = parser.parse_args(argv)
    if args.action == "generate":
        if (
            args.request is None
            or args.expected_request_sha256 is None
            or args.output is None
        ):
            parser.error(
                "generate requires --request, --expected-request-sha256, and --output"
            )
        if any(
            value is not None
            for value in (
                args.packet,
                args.expected_packet_sha256,
                args.review,
                args.expected_review_sha256,
            )
        ):
            parser.error("generate forbids draft-packet and review-install arguments")
    elif args.action in {"install-f117-draft-packet", "draft-request"}:
        if args.packet is None or args.expected_packet_sha256 is None:
            parser.error(
                f"{args.action} requires --packet and --expected-packet-sha256"
            )
        if (
            args.request is not None
            or args.expected_request_sha256 is not None
            or args.review is not None
            or args.expected_review_sha256 is not None
            or args.output is not None
        ):
            parser.error(f"{args.action} forbids request/review/output arguments")
    elif args.action == "install-request-review":
        if any(
            value is None
            for value in (
                args.request,
                args.expected_request_sha256,
                args.review,
                args.expected_review_sha256,
            )
        ):
            parser.error(
                "install-request-review requires --request, --expected-request-sha256, "
                "--review, and --expected-review-sha256"
            )
        if (
            args.output is not None
            or args.packet is not None
            or args.expected_packet_sha256 is not None
        ):
            parser.error("install-request-review forbids output and packet arguments")
    elif any(
        value is not None
        for value in (
            args.request,
            args.expected_request_sha256,
            args.review,
            args.expected_review_sha256,
            args.output,
            args.packet,
            args.expected_packet_sha256,
        )
    ):
        parser.error("retain-generator forbids request, output, and packet arguments")
    return args


def validate_root(args: argparse.Namespace) -> Path:
    """Validate the common campaign-root trust boundary."""

    root = args.root.absolute()
    require_no_symlink_components(root, "campaign root")
    require_directory(root, "campaign root")
    canonical = root == DEFAULT_ROOT
    if not canonical and not args.allow_local_root:
        raise ValueError("--allow-local-root is required outside the canonical root")
    if canonical and args.allow_local_root:
        raise ValueError("--allow-local-root is forbidden for canonical execution")
    if canonical and args.squeue_file is not None:
        raise ValueError("--squeue-file is forbidden for canonical execution")
    require_directory(root / "accounting", "accounting directory")
    return root


def validate_root_and_output(args: argparse.Namespace) -> tuple[Path, Path]:
    """Validate the root trust boundary and explicit staged-only output path."""

    root = validate_root(args)
    accounting = root / "accounting"
    assert args.output is not None
    output = args.output.absolute()
    if output.parent != accounting:
        raise ValueError("staged recost output must be a direct child of accounting")
    if OUTPUT_NAME_PATTERN.fullmatch(output.name) is None:
        raise ValueError("staged recost output must end in a safe .json.staged name")
    require_no_symlink_components(output, "staged recost output", include_leaf=False)
    return root, output


def require_output_namespace_empty(output: Path) -> None:
    """Require an empty artifact namespace or one exact final commit marker."""

    accounting = output.parent
    canonical_name = output.name.removesuffix(".staged")
    with absolute_descriptor(
        accounting, "accounting directory", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        retained = sorted(
            entry.name
            for entry in os.scandir(descriptor)
            if (
                entry.name == canonical_name
                or entry.name.startswith(f"{canonical_name}.")
                or entry.name == f".{canonical_name}"
                or entry.name.startswith(f".{canonical_name}.")
            )
        )
        if retained == [output.name]:
            recovery = os.stat(output.name, dir_fd=descriptor, follow_symlinks=False)
            try:
                require_regular_profile(
                    recovery,
                    "staged recost output final commit marker",
                    expected_mode=0o444,
                )
            except ValueError:
                pass
            else:
                return
    if retained:
        raise ValueError(
            "staged recost output namespace is not empty: " + ", ".join(retained)
        )


def queue_account_name() -> str:
    """Return the scheduler account name for the effective user."""

    try:
        return pwd.getpwuid(os.geteuid()).pw_name
    except KeyError as error:
        raise ValueError("effective UID has no scheduler account name") from error


def parse_queue_rows(payload: str) -> list[str]:
    """Parse strict queue rows and require the Stage I CGL queue to be drained."""

    rows = []
    for line in payload.splitlines():
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != 3 or any(not field or field.strip() != field for field in fields):
            raise ValueError(f"squeue output row is malformed: {line}")
        if fields[1].startswith(CGL_JOB_NAME_PREFIX):
            rows.append(line)
    if rows:
        raise ValueError("active CGL scheduler queue is not drained: " + "; ".join(rows))
    return rows


def require_drained_queue(args: argparse.Namespace, root: Path) -> None:
    """Require no Stage I CGL workflow job to remain queued or running."""

    if args.squeue_file is not None:
        queue = args.squeue_file.absolute()
        if root == DEFAULT_ROOT:
            raise ValueError("queue fixture is forbidden for canonical execution")
        require_no_symlink_components(queue, "queue fixture")
        with absolute_descriptor(queue, "queue fixture", flags=os.O_RDONLY) as descriptor:
            require_regular_profile(os.fstat(descriptor), "queue fixture")
            payload = os.read(descriptor, 1024 * 1024)
            if os.read(descriptor, 1):
                raise ValueError("queue fixture is unexpectedly large")
        try:
            parse_queue_rows(payload.decode())
        except UnicodeDecodeError as error:
            raise ValueError("queue fixture is not UTF-8") from error
        return
    environment = hardened_child_environment("SLURM_")
    try:
        completed = subprocess.run(
            [str(SQUEUE), "-h", "-u", queue_account_name(), "-o", "%i|%j|%T"],
            check=True,
            capture_output=True,
            text=True,
            env=environment,
        )
    except (OSError, subprocess.CalledProcessError) as error:
        raise ValueError("squeue is unavailable; refusing recost generation") from error
    parse_queue_rows(completed.stdout)


def parse_request(payload: bytes) -> dict[str, object]:
    """Parse the strict reviewed recost request schema."""

    request = require_exact_keys(
        parse_json(payload, "recost request"),
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "requested_by",
            "scope",
            "barrier",
            "inputs",
            "recommendations",
        },
        "recost request",
    )
    if request["schema_version"] != 2:
        raise ValueError("recost request schema version must be 2")
    if request["record_type"] != "stage-i-recost-recommendation-request":
        raise ValueError("recost request record type differs")
    checkpoint = require_safe_id(request["checkpoint"], "recost checkpoint")
    artifact_name = require_nonempty_string(
        request["artifact_name"], "recost artifact name"
    )
    match = RECOST_ARTIFACT_PATTERN.fullmatch(artifact_name)
    if match is None or checkpoint != f"F-{match.group(1)}":
        raise ValueError("recost artifact name and checkpoint ID differ")
    if request["execution_epoch"] != EXECUTION_EPOCH:
        raise ValueError("recost request execution epoch differs")
    parse_utc_timestamp(request["generated_utc"], "recost request generation timestamp")
    parse_utc_timestamp(request["expires_utc"], "recost request expiry timestamp")
    require_nonempty_string(request["requested_by"], "recost request author")
    require_nonempty_string(request["scope"], "recost request scope")
    require_exact_keys(request["barrier"], {"recorded_segments"}, "barrier")
    require_exact_keys(
        request["inputs"],
        {
            "reconciliation",
            "ledger",
            "reservations",
            "manifests",
            "scheduler_evidence",
            "storage_evidence",
            "source_bundle",
            "matrix",
            "stage_i_helper",
            "ceiling_evidence",
            "ceiling_publication_audit",
            "source_authority",
            "qualification_approval",
            "predecessor_recost",
            "predecessor_recost_independent_review",
            "predecessor_recost_publication_audit",
            "r17_readiness_evidence",
            "r17_readiness_independent_review",
            "r17_readiness_publication_audit",
        },
        "recost inputs",
    )
    require_exact_keys(
        request["recommendations"],
        {"mode", "max_wave_nodes", "profiles"},
        "recost recommendations",
    )
    return request


def run_authenticated_reconcile(
    helper_path: Path,
    helper_sha256: str,
    root: Path,
    *,
    stage_i_lock_held: bool = False,
) -> dict[str, object]:
    """Run the authenticated committed helper and return its live reconcile report."""

    with absolute_descriptor(
        helper_path, "Stage I helper", flags=os.O_RDONLY
    ) as descriptor:
        require_regular_profile(
            os.fstat(descriptor), "Stage I helper", expected_mode=0o644
        )
        if sha256_descriptor(descriptor) != helper_sha256:
            raise ValueError("Stage I helper checksum changed before reconciliation")
        os.lseek(descriptor, 0, os.SEEK_SET)
        if stage_i_lock_held:
            module_name = (
                f"_cgl_lf_stage_i_locked_reconcile_{os.getpid()}_{descriptor}_{id(root)}"
            )
            loader = SourceFileLoader(module_name, f"/proc/self/fd/{descriptor}")
            spec = importlib.util.spec_from_loader(module_name, loader)
            if spec is None or spec.loader is None:
                raise ValueError("Stage I helper module loader is unavailable")
            module = importlib.util.module_from_spec(spec)
            sys.modules[module_name] = module
            try:
                spec.loader.exec_module(module)
                reconcile_report = getattr(module, "reconcile_report", None)
                if not callable(reconcile_report):
                    raise ValueError("Stage I helper lacks read-only reconcile_report")
                report = reconcile_report(root)
            finally:
                sys.modules.pop(module_name, None)
            if sha256_descriptor(descriptor) != helper_sha256:
                raise ValueError("Stage I helper checksum changed during reconciliation")
            if not isinstance(report, dict):
                raise ValueError("live Stage I reconciliation is not an object")
            return report
        environment = hardened_child_environment("SLURM_")
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        completed = subprocess.run(
            [
                sys.executable,
                "-B",
                f"/proc/self/fd/{descriptor}",
                "--root",
                str(root),
                "reconcile",
            ],
            pass_fds=(descriptor,),
            check=False,
            capture_output=True,
            env=environment,
        )
        if sha256_descriptor(descriptor) != helper_sha256:
            raise ValueError("Stage I helper checksum changed during reconciliation")
    if completed.returncode:
        raise ValueError("authenticated Stage I helper reconciliation failed")
    return parse_json(completed.stdout, "live Stage I reconciliation")


def parse_reconciliation(value: object, root: Path) -> dict[str, object]:
    """Validate exact clean canonical reconciliation evidence."""

    if not isinstance(value, dict):
        raise ValueError("reconciliation evidence must be an object")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("reconciliation execution epoch differs")
    if value.get("root") != str(root):
        raise ValueError("reconciliation root differs")
    if value.get("consistent") is not True or value.get("issues") != []:
        raise ValueError("reconciliation is not clean")
    counts = value.get("counts")
    if not isinstance(counts, dict):
        raise ValueError("reconciliation counts must be an object")
    selected = {key: counts.get(key) for key in COUNT_KEYS}
    if any(isinstance(item, bool) or not isinstance(item, int) or item < 0
           for item in selected.values()):
        raise ValueError("reconciliation counts are invalid")
    if selected["transactions"] != 0:
        raise ValueError("reconciliation retains active transactions")
    if selected["active_reservations"] != 0:
        raise ValueError("reconciliation retains active reservations")
    return value


def parse_ledger(payload: bytes, request_timestamp: datetime) -> list[dict[str, str]]:
    """Parse and validate the retained Stage I ledger."""

    try:
        with io.StringIO(payload.decode(), newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != list(LEDGER_COLUMNS):
                raise ValueError("ledger columns differ from Stage I accounting")
            rows = list(reader)
    except UnicodeDecodeError as error:
        raise ValueError("ledger is not UTF-8") from error
    if not rows:
        raise ValueError("ledger must retain at least one recorded segment")
    jobs: set[str] = set()
    cumulative = Decimal("0.000000")
    for index, row in enumerate(rows):
        label = f"ledger row {index}"
        if set(row) != set(LEDGER_COLUMNS) or any(value is None for value in row.values()):
            raise ValueError(f"{label} columns differ from Stage I accounting")
        if row.get("execution_epoch") != EXECUTION_EPOCH:
            raise ValueError(f"{label} execution epoch differs")
        job_id = require_job_id(row.get("job_id"), f"{label} job ID")
        if job_id in jobs:
            raise ValueError(f"ledger duplicates job {job_id}")
        jobs.add(job_id)
        case_id = require_safe_id(row.get("case_id"), f"{label} case ID")
        if case_id not in AUTHORIZED_CASE_IDS:
            raise ValueError(f"{label} has unknown case {case_id}")
        require_safe_id(row.get("segment"), f"{label} segment")
        result = require_nonempty_string(row.get("result"), f"{label} result")
        if result not in RECORDED_RESULTS:
            raise ValueError(f"{label} result is not a supported terminal outcome")
        state = require_nonempty_string(row.get("state"), f"{label} scheduler state")
        require_nonempty_string(row.get("exit_code"), f"{label} exit code")
        if state not in TERMINAL_SCHEDULER_STATES:
            raise ValueError(f"{label} scheduler state is not terminal")
        if result in SCIENTIFIC_RESULTS and (
            state != "COMPLETED" or row.get("exit_code") != "0:0"
        ):
            raise ValueError(f"{label} scientific result is not a successful completed job")
        nodes = require_integer(int(row["nodes"]), f"{label} nodes", minimum=1)
        elapsed = require_integer(int(row["elapsed_seconds"]), f"{label} elapsed seconds")
        submitted = parse_scheduler_timestamp(row.get("submitted_utc"), f"{label} submit time")
        completed = parse_scheduler_timestamp(
            row.get("completed_utc"), f"{label} completion time"
        )
        if completed <= submitted:
            raise ValueError(f"{label} chronology is invalid")
        if completed > request_timestamp + SCHEDULER_TIME_TOLERANCE:
            raise ValueError(f"{label} completes after the request")
        if elapsed > int((completed - submitted).total_seconds()) + int(
            SCHEDULER_TIME_TOLERANCE.total_seconds()
        ):
            raise ValueError(f"{label} elapsed time exceeds its chronology")
        walltime = walltime_seconds(row.get("requested_walltime"), f"{label} walltime")
        reserved_text = row.get("reserved_node_hours")
        require_canonical_node_hours(reserved_text, f"{label} reserved use", positive=True)
        expected_reserved = canonical_node_hours(
            Decimal(nodes * walltime) / Decimal(3600),
            f"{label} expected reserved use",
        )
        if reserved_text != expected_reserved:
            raise ValueError(f"{label} reserved use differs from allocation arithmetic")
        actual_text = row.get("actual_node_hours")
        require_canonical_node_hours(actual_text, f"{label} actual use")
        expected_actual_text = canonical_node_hours(
            Decimal(nodes * elapsed) / Decimal(3600),
            f"{label} expected actual use",
        )
        if actual_text != expected_actual_text:
            raise ValueError(f"{label} actual use differs from scheduler arithmetic")
        expected_actual = Decimal(expected_actual_text)
        cumulative += expected_actual
        retained_cumulative_text = row.get("cumulative_stage_i_node_hours")
        require_canonical_node_hours(
            retained_cumulative_text, f"{label} cumulative use"
        )
        expected_cumulative = canonical_node_hours(
            cumulative, f"{label} expected cumulative use"
        )
        if retained_cumulative_text != expected_cumulative:
            raise ValueError(f"{label} cumulative use differs from ledger arithmetic")
    return rows


def parse_reservations(value: object) -> list[dict[str, object]]:
    """Validate the retained drained reservation store."""

    if not isinstance(value, list) or not value:
        raise ValueError("reservations must be a nonempty list")
    identities: set[str] = set()
    for index, item in enumerate(value):
        if not isinstance(item, dict):
            raise ValueError(f"reservation {index} must be an object")
        manifest = require_nonempty_string(
            item.get("manifest"), f"reservation {index} manifest"
        )
        raw_job_id = item.get("job_id")
        job_id = (
            None if raw_job_id is None else require_job_id(raw_job_id, f"reservation {index} job ID")
        )
        identity = f"job:{job_id}" if job_id is not None else f"manifest:{manifest}"
        if identity in identities:
            raise ValueError(f"reservations duplicate identity {identity}")
        identities.add(identity)
        label = f"reservation {job_id}" if job_id is not None else f"reservation {index}"
        if item.get("execution_epoch") != EXECUTION_EPOCH:
            raise ValueError(f"{label} execution epoch differs")
        case_id = require_safe_id(item.get("case_id"), f"{label} case ID")
        if case_id not in AUTHORIZED_CASE_IDS:
            raise ValueError(f"{label} has unknown case {case_id}")
        require_safe_id(item.get("segment"), f"{label} segment")
        require_nonempty_string(item.get("case_name"), f"{label} case name")
        require_integer(item.get("nodes"), f"{label} nodes", minimum=1)
        walltime_seconds(
            item.get("requested_walltime"), f"{label} walltime"
        )
        if item.get("state") in {"prepared", "submitted"}:
            raise ValueError(f"{label} remains active")
        if item.get("state") == "recorded":
            if job_id is None:
                raise ValueError(f"{label} recorded state lacks a job ID")
            result = require_nonempty_string(
                item.get("result"), f"{label} result"
            )
            if result not in RECORDED_RESULTS:
                raise ValueError(f"{label} result is unsupported")
            try:
                actual = Decimal(str(item.get("actual_node_hours")))
            except InvalidOperation as error:
                raise ValueError(f"{label} actual use is invalid") from error
            if not actual.is_finite() or actual < 0:
                raise ValueError(f"{label} actual use is invalid")
        elif item.get("state") == "cancelled":
            require_nonempty_string(item.get("notes"), f"{label} cancellation notes")
            if "result" in item or "actual_node_hours" in item:
                raise ValueError(f"{label} fabricates recorded accounting")
        else:
            raise ValueError(f"{label} is not terminal at the barrier")
    return value


def manifest_inventory(root: Path) -> set[str]:
    """Return the exact non-symlink canonical manifest inventory."""

    runs = root / "runs/mks24-stage-i" / EXECUTION_EPOCH
    require_directory(runs, "Stage I run store")
    retained: set[str] = set()
    for path in runs.glob("*/*/manifest/prepared_run.json"):
        require_no_symlink_components(path, "Stage I manifest")
        relative = path.relative_to(root).as_posix()
        retained.add(relative)
    return retained


def parse_manifests(root: Path, bindings: object, tracker: InputTracker
                    ) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Authenticate and parse every retained canonical manifest."""

    if not isinstance(bindings, list) or not bindings:
        raise ValueError("manifest bindings must be a nonempty list")
    parsed_bindings: list[dict[str, object]] = []
    manifests: list[dict[str, object]] = []
    paths: set[str] = set()
    for index, binding in enumerate(bindings):
        relative, expected = input_binding(binding, f"manifest binding {index}")
        relative_string = relative.as_posix()
        if relative_string in paths:
            raise ValueError("manifest bindings contain a duplicate path")
        paths.add(relative_string)
        path = root_path(root, relative_string, f"manifest binding {index}")
        if not relative_string.endswith("/manifest/prepared_run.json"):
            raise ValueError("manifest binding does not name prepared_run.json")
        payload = tracker.read(path, expected, f"manifest {relative_string}", expected_mode=0o644)
        manifest = parse_json(payload, f"manifest {relative_string}")
        if not isinstance(manifest, dict):
            raise ValueError(f"manifest {relative_string} must be an object")
        if manifest.get("execution_epoch") != EXECUTION_EPOCH:
            raise ValueError(f"manifest {relative_string} execution epoch differs")
        if manifest.get("state") not in {"recorded", "cancelled"}:
            raise ValueError(f"manifest {relative_string} is not terminal at the barrier")
        manifest["_manifest_path"] = str(path)
        manifest["_manifest_sha256"] = expected
        parsed_bindings.append({"path": relative_string, "sha256": expected})
        manifests.append(manifest)
    inventory = manifest_inventory(root)
    if paths != inventory:
        missing = sorted(inventory - paths)
        extra = sorted(paths - inventory)
        raise ValueError(
            f"manifest bindings differ from canonical inventory; missing={missing}, extra={extra}"
        )
    return parsed_bindings, manifests


def parse_controller_history(payload: bytes, label: str) -> dict[str, list[float]]:
    """Reproduce the production controller's labeled Athena history parser."""

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
                candidate = [column for _, column in found]
                if (
                    indices != list(range(1, len(candidate) + 1))
                    or len(candidate) != len(set(candidate))
                    or (labels and candidate != labels)
                ):
                    raise ValueError(f"{label} labels are ambiguous or invalid")
                labels = candidate
            continue
        if line.strip():
            try:
                rows.append([float(value) for value in line.split()])
            except ValueError as error:
                raise ValueError(f"{label} contains a nonnumeric value") from error
    if not labels or not rows:
        raise ValueError(f"{label} is missing labels or data")
    if any(len(row) != len(labels) for row in rows):
        raise ValueError(f"{label} row width does not match labels")
    if any(not math.isfinite(value) for row in rows for value in row):
        raise ValueError(f"{label} contains a non-finite value")
    return {
        column: [row[index] for row in rows]
        for index, column in enumerate(labels)
    }


def reproduce_continuation_plasma_evidence(
    case_id: str,
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
    *,
    allow_frozen_e03_no_max_ndiv: bool = False,
) -> dict[str, object]:
    """Reproduce the current controller's complete clean-partial plasma policy."""

    mhd_required = {
        "time", "mass", "tot-E", "lf_nstage", "lf_qface", "lf_qprcap",
        "lf_qpr10", "lf_qpecap", "lf_qpe10", "lf_qprwrk", "lf_qpewrk",
        "lf_cpwrk", "lf_cawrk", "lf_hwproj", *CONTROLLER_STRICT_LF_FAILURE_COLUMNS,
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
        name: max(mhd[name]) for name in CONTROLLER_STRICT_LF_FAILURE_COLUMNS
    }
    if any(
        value != 0.0
        for name in CONTROLLER_STRICT_LF_FAILURE_COLUMNS
        for value in mhd[name]
    ):
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
        if any(
            value != 0.0
            for name in ("lf_cpwrk", "lf_cawrk")
            for value in mhd[name]
        ):
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
    """Read one exact owner-controlled controller migration artifact."""

    path = root_path(DEFAULT_ROOT, binding.get("path"), f"{label} path")
    mode = require_integer(binding.get("mode"), f"{label} mode")
    size_bytes = require_integer(binding.get("size_bytes"), f"{label} size", minimum=1)
    expected_sha256 = require_sha256(binding.get("sha256"), f"{label} SHA-256")
    with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
        profile = os.fstat(descriptor)
        require_regular_profile(profile, label, expected_mode=mode)
        chunks = []
        digest = hashlib.sha256()
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
            digest.update(block)
        retained = b"".join(chunks)
        if profile.st_size != size_bytes or digest.hexdigest() != expected_sha256:
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
    """Read one exact JSON controller migration artifact."""

    _, payload, public_binding = read_frozen_e03_migration_binding(binding, label)
    value = parse_json(payload, label)
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value, public_binding


def validate_frozen_e03_independent_validation(
    contract: dict[str, object],
    validation: dict[str, object],
    review: dict[str, object] | None,
    bindings: dict[str, dict[str, object] | None],
) -> None:
    """Require the controller's exact independent migration validation."""

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
        name: 0 for name in CONTROLLER_STRICT_LF_FAILURE_COLUMNS
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
    inspection: dict[str, object],
    manifest: dict[str, object],
    job_id: str,
    case_id: str,
    segment: str,
    expected_ranks: int,
    mhd: dict[str, list[float]],
    user: dict[str, list[float]],
) -> dict[str, object]:
    """Reproduce the current controller's exact retained R03/R12 migration."""

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
        manifest.get("project_root") != str(DEFAULT_ROOT)
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
        or allocation.get("nodes") * EXPECTED_RANKS_PER_NODE != expected_ranks
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
    evidence = reproduce_continuation_plasma_evidence(
        case_id, mhd, user, allow_frozen_e03_no_max_ndiv=True
    )
    evidence["policy"] = FROZEN_E03_NO_MAX_NDIV_MIGRATION_POLICY
    evidence["continuation_authorized"] = False
    evidence["continuation_eligible"] = False
    evidence["eligibility_only"] = True
    evidence["normalized_ct_divb_evidence"] = normalized_ct_divb_evidence
    evidence["ct_divergence_claimed"] = False
    evidence["ct_divergence_reason"] = FROZEN_E03_CT_DIVERGENCE_REASON
    evidence["authorization_basis"] = (
        "none; historical frozen-E03 evidence is inventory-only and cannot "
        "authorize continuation"
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


def authenticate_controller_history(
    value: object,
    output_dir: Path,
    suffix: str,
    label: str,
) -> dict[str, list[float]]:
    """Authenticate and parse one exact retained controller history binding."""

    record = require_exact_keys(value, {"path", "size_bytes", "sha256"}, label)
    path = absolute_root_member(output_dir, record["path"], f"{label} path")
    if path.parent != output_dir or not path.name.endswith(suffix):
        raise ValueError(f"{label} path differs from the controller output inventory")
    expected_size = require_integer(record["size_bytes"], f"{label} size", minimum=1)
    expected_sha256 = require_sha256(record["sha256"], f"{label} SHA-256")
    with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
        profile = os.fstat(descriptor)
        require_regular_profile(profile, label, expected_mode=0o644)
        if profile.st_size != expected_size:
            raise ValueError(f"{label} size differs")
        payload = b""
        digest = hashlib.sha256()
        chunks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
            digest.update(block)
        payload = b"".join(chunks)
        if digest.hexdigest() != expected_sha256:
            raise ValueError(f"{label} checksum differs")
    return parse_controller_history(payload, label)


def require_controller_product_schema(
    value: object,
    output_dir: Path,
    expected_ranks: int,
    label: str,
) -> list[Path]:
    """Require one exact retained shared or per-rank controller product schema."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    storage = value.get("storage")
    keys = {"path", "size_bytes", "sha256", "storage"}
    if storage == "per_rank":
        keys.add("rank_files")
    record = require_exact_keys(value, keys, label)

    def file_binding(item: object, item_label: str) -> Path:
        binding = require_exact_keys(item, {"path", "size_bytes", "sha256"}, item_label)
        path = absolute_root_member(output_dir, binding["path"], f"{item_label} path")
        require_integer(binding["size_bytes"], f"{item_label} size", minimum=1)
        require_sha256(binding["sha256"], f"{item_label} SHA-256")
        return path

    if storage == "shared_mpiio":
        return [file_binding(
            {key: record[key] for key in ("path", "size_bytes", "sha256")}, label
        )]
    if storage != "per_rank":
        raise ValueError(f"{label} storage differs")
    rank_files = record["rank_files"]
    if not isinstance(rank_files, list) or len(rank_files) != expected_ranks:
        raise ValueError(f"{label} rank inventory differs")
    paths = [
        file_binding(item, f"{label} rank {index}")
        for index, item in enumerate(rank_files)
    ]
    if (
        [path.parent.name for path in paths]
        != [f"rank_{rank:08d}" for rank in range(expected_ranks)]
        or len(set(paths)) != expected_ranks
        or len({path.parent.parent for path in paths}) != 1
        or len({path.name for path in paths}) != 1
    ):
        raise ValueError(f"{label} rank inventory differs")
    first = rank_files[0]
    assert isinstance(first, dict)
    if any(record[key] != first[key] for key in ("path", "size_bytes", "sha256")):
        raise ValueError(f"{label} primary binding differs")
    return paths


def require_clean_partial_controller_evidence(
    inspection: dict[str, object],
    manifest: dict[str, object],
    job_id: str,
    final_time_value: float,
    expected_ranks: int,
) -> None:
    """Reproduce and authenticate the controller's complete schema-4 plasma proof."""

    if (
        require_integer(
            inspection.get("schema_version"),
            f"manifest {job_id} inspection schema version",
            minimum=1,
        )
        != CONTROLLER_CLEAN_PARTIAL_SCHEMA_VERSION
    ):
        raise ValueError(
            f"manifest {job_id} clean-partial continuation lacks controller schema-4 evidence"
        )
    inspection_keys = set(inspection)
    frozen_e03_migration = (
        inspection_keys == FROZEN_E03_NO_MAX_NDIV_INSPECTION_KEYS
    )
    if (
        inspection_keys != CONTROLLER_CLEAN_PARTIAL_INSPECTION_KEYS
        and not frozen_e03_migration
    ):
        raise ValueError(f"manifest {job_id} clean-partial inspection schema differs")
    expected_checks = (
        FROZEN_E03_NO_MAX_NDIV_CHECKS
        if frozen_e03_migration
        else CONTROLLER_CLEAN_PARTIAL_CHECKS
    )
    checks = require_exact_keys(
        inspection["checks"],
        set(expected_checks),
        f"manifest {job_id} clean-partial checks",
    )
    failures = require_exact_keys(
        inspection["maximum_strict_failure_counts"],
        set(CONTROLLER_STRICT_LF_FAILURE_COLUMNS),
        f"manifest {job_id} strict LF failure counts",
    )
    case_id = require_safe_id(inspection["case_id"], f"manifest {job_id} inspection case")
    segment = require_safe_id(inspection["segment"], f"manifest {job_id} inspection segment")
    _, _, segment_target = parse_segment(segment, f"manifest {job_id} inspection segment")
    paths = manifest.get("paths")
    if not isinstance(paths, dict):
        raise ValueError(f"manifest {job_id} lacks output provenance")
    project_root = absolute_root_member(
        Path("/"), manifest.get("project_root"), f"manifest {job_id} project root"
    )
    output_dir = absolute_root_member(
        project_root, paths.get("output_dir"), f"manifest {job_id} output directory"
    )
    manifest_path = absolute_root_member(
        project_root, manifest.get("_manifest_path"), f"manifest {job_id} path"
    )
    mhd = authenticate_controller_history(
        inspection["mhd_history"], output_dir, ".mhd.hst", f"manifest {job_id} MHD history"
    )
    user = authenticate_controller_history(
        inspection["user_history"], output_dir, ".user.hst", f"manifest {job_id} user history"
    )
    if frozen_e03_migration:
        current_plasma = require_frozen_e03_no_max_ndiv_migration(
            inspection,
            manifest,
            job_id,
            case_id,
            segment,
            expected_ranks,
            mhd,
            user,
        )
        plasma_binding_valid = True
    else:
        current_plasma = reproduce_continuation_plasma_evidence(case_id, mhd, user)
        plasma_binding_valid = (
            inspection["plasma_continuation_policy"] == CONTINUATION_PLASMA_POLICY
            and inspection["plasma_continuation_evidence"] == current_plasma
        )
    marker_modes = inspection["restart_time_marker_modes"]
    snapshots = inspection["snapshots"]
    snapshot_times = inspection["snapshot_times"]
    restarts = inspection["restarts"]
    restart_times = inspection["restart_times"]
    if not isinstance(snapshots, list) or not snapshots:
        raise ValueError(f"manifest {job_id} clean-partial snapshots differ")
    if not isinstance(restarts, list) or not restarts:
        raise ValueError(f"manifest {job_id} clean-partial restarts differ")
    snapshot_paths = [
        require_controller_product_schema(
            record, output_dir, expected_ranks, f"manifest {job_id} snapshot {index}"
        )
        for index, record in enumerate(snapshots)
    ]
    restart_paths = [
        require_controller_product_schema(
            record, output_dir, expected_ranks, f"manifest {job_id} restart {index}"
        )
        for index, record in enumerate(restarts)
    ]
    if (
        not isinstance(snapshot_times, list)
        or len(snapshot_times) != len(snapshot_paths)
        or not isinstance(restart_times, list)
        or len(restart_times) != len(restart_paths)
        or max(
            require_finite_float(value, f"manifest {job_id} snapshot time")
            for value in snapshot_times
        ) < final_time_value - 1.0e-10
        or not isinstance(marker_modes, list)
        or len(marker_modes) != len(restart_paths)
        or any(
            not isinstance(group, list)
            or len(group) != expected_ranks
            or any(
                mode not in {"full_precision", "legacy_default_precision"}
                for mode in group
            )
            for group in marker_modes
        )
    ):
        raise ValueError(
            f"manifest {job_id} clean-partial retained-product evidence differs"
        )
    parsed_restart_times = [
        require_finite_float(value, f"manifest {job_id} restart time")
        for value in restart_times
    ]
    terminal_restart_time = require_finite_float(
        inspection["terminal_restart_time"], f"manifest {job_id} terminal restart time"
    )
    terminal_paths = require_controller_product_schema(
        inspection["terminal_restart"],
        output_dir,
        expected_ranks,
        f"manifest {job_id} terminal restart",
    )
    matching_restarts = [
        index
        for index, value in enumerate(parsed_restart_times)
        if abs(value - final_time_value) <= 1.0e-10
    ]
    history_failures = {
        key: max(mhd[key]) for key in CONTROLLER_STRICT_LF_FAILURE_COLUMNS
    }
    if (
        inspection["execution_epoch"] != EXECUTION_EPOCH
        or inspection["job_id"] != job_id
        or manifest.get("run", {}).get("case_id") != case_id
        or manifest.get("run", {}).get("segment") != segment
        or inspection["manifest"] != str(manifest_path)
        or parse_utc_timestamp(
            inspection["inspected_utc"], f"manifest {job_id} inspection timestamp"
        )
        > datetime.now(timezone.utc) + SCHEDULER_TIME_TOLERANCE
        or abs(
            require_finite_float(inspection["required_time"], f"manifest {job_id} required time")
            - segment_target
        )
        > 1.0e-12
        or abs(mhd["time"][-1] - final_time_value) > 1.0e-12
        or checks != expected_checks
        or failures != history_failures
        or failures != {key: 0.0 for key in CONTROLLER_STRICT_LF_FAILURE_COLUMNS}
        or inspection["accepted"] is not False
        or inspection["clean_for_continuation"] is not True
        or not plasma_binding_valid
        or (
            frozen_e03_migration
            and (
                current_plasma.get("continuation_authorized") is not False
                or current_plasma.get("continuation_eligible") is not False
                or current_plasma.get("eligibility_only") is not True
            )
        )
        or (
            not frozen_e03_migration
            and current_plasma["continuation_authorized"] is not True
        )
        or inspection["restart_time_marker_bypass"] is not False
        or matching_restarts != [len(restart_paths) - 1]
        or restart_paths[matching_restarts[0]] != terminal_paths
        or abs(terminal_restart_time - final_time_value) > 1.0e-12
        or require_finite_float(
            inspection["final_hardwall_projection_count"],
            f"manifest {job_id} final hardwall projection count",
        )
        != mhd["lf_hwproj"][-1]
    ):
        raise ValueError(
            f"manifest {job_id} clean-partial continuation physics evidence differs"
        )


def manifest_identity(manifest: dict[str, object]) -> tuple[str | None, str, str, str]:
    """Return and validate one terminal manifest identity and normalized outcome."""

    run = manifest.get("run")
    allocation = manifest.get("allocation")
    accounting = manifest.get("accounting")
    inspection = manifest.get("scientific_inspection")
    if not all(isinstance(item, dict) for item in (run, allocation)):
        raise ValueError("terminal manifest lacks run/allocation")
    assert isinstance(run, dict)
    assert isinstance(allocation, dict)
    case_id = require_safe_id(run.get("case_id"), "manifest case ID")
    segment = require_safe_id(run.get("segment"), "manifest segment")
    _, segment_start, segment_target = parse_segment(segment, "manifest segment")
    if case_id not in AUTHORIZED_CASE_IDS:
        raise ValueError(f"manifest has unknown case {case_id}")
    require_integer(allocation.get("nodes"), f"manifest {case_id}/{segment} nodes", minimum=1)
    if (
        allocation.get("ranks_per_node") != EXPECTED_RANKS_PER_NODE
        or allocation.get("cpus_per_task") != EXPECTED_CPUS_PER_TASK
    ):
        raise ValueError(f"manifest {case_id}/{segment} rank/task allocation differs")
    walltime_seconds(
        allocation.get("requested_walltime"), f"manifest {case_id}/{segment} walltime"
    )
    if manifest.get("state") == "cancelled":
        if accounting is not None or inspection is not None:
            raise ValueError(f"cancelled manifest {case_id}/{segment} fabricates accounting")
        cancellation = manifest.get("cancellation")
        if not isinstance(cancellation, dict):
            raise ValueError(f"cancelled manifest {case_id}/{segment} lacks cancellation evidence")
        require_nonempty_string(
            cancellation.get("notes"), f"cancelled manifest {case_id}/{segment} notes"
        )
        parse_utc_timestamp(
            cancellation.get("cancelled_utc"),
            f"cancelled manifest {case_id}/{segment} timestamp",
        )
        raw_job_id = manifest.get("job_id")
        job_id = None if raw_job_id is None else require_job_id(raw_job_id, "manifest job ID")
        if cancellation.get("job_id") not in {None, job_id}:
            raise ValueError(f"cancelled manifest {case_id}/{segment} job ID differs")
        return job_id, case_id, segment, "cancelled"
    if manifest.get("state") != "recorded" or not isinstance(accounting, dict):
        raise ValueError(f"manifest {case_id}/{segment} is not a recorded terminal outcome")
    job_id = require_job_id(manifest.get("job_id"), "manifest job ID")
    result = require_nonempty_string(accounting.get("result"), "manifest result")
    if accounting.get("job_id") != job_id:
        raise ValueError(f"manifest {job_id} accounting job ID differs")
    if accounting.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(f"manifest {job_id} accounting execution epoch differs")
    if accounting.get("case_id") != case_id or accounting.get("segment") != segment:
        raise ValueError(f"manifest {job_id} accounting identity differs")
    if result not in RECORDED_RESULTS:
        raise ValueError(f"manifest {job_id} result is unsupported")
    state = require_nonempty_string(accounting.get("state"), f"manifest {job_id} scheduler state")
    require_nonempty_string(accounting.get("exit_code"), f"manifest {job_id} exit code")
    if state not in TERMINAL_SCHEDULER_STATES:
        raise ValueError(f"manifest {job_id} scheduler state is not terminal")
    if result in SCIENTIFIC_RESULTS and (
        state != "COMPLETED" or accounting.get("exit_code") != "0:0"
    ):
        raise ValueError(f"manifest {job_id} scientific result is not successful")
    if result not in SCIENTIFIC_RESULTS:
        if inspection is not None:
            raise ValueError(f"manifest {job_id} non-scientific result fabricates inspection")
        return job_id, case_id, segment, result
    if not isinstance(inspection, dict):
        raise ValueError(f"manifest {job_id} scientific result lacks inspection")
    assert isinstance(inspection, dict)
    if inspection.get("case_id") != case_id or inspection.get("segment") != segment:
        raise ValueError(f"manifest {job_id} inspection identity differs")
    if (
        inspection.get("job_id") != job_id
        or inspection.get("execution_epoch") != EXECUTION_EPOCH
    ):
        raise ValueError(f"manifest {job_id} inspection provenance differs")
    observed_final_time = require_finite_float(
        inspection.get("final_time"), f"manifest {job_id} final time"
    )
    if (
        observed_final_time <= segment_start
        or observed_final_time > REQUIRED_CASE_FINAL_TIME
        or observed_final_time > segment_target + 1.0e-10
    ):
        raise ValueError(f"manifest {job_id} final time is outside its segment interval")
    if result == "accepted" and abs(observed_final_time - segment_target) > 1.0e-10:
        raise ValueError(f"manifest {job_id} accepted final time differs from its endpoint")
    if result == "accepted" and inspection.get("accepted") is not True:
        raise ValueError(f"manifest {job_id} accepted result differs from inspection")
    if (
        result == "clean_partial"
        and (
            inspection.get("accepted") is not False
            or inspection.get("clean_for_continuation") is not True
        )
    ):
        raise ValueError(f"manifest {job_id} clean-partial result differs from inspection")
    if result == "clean_partial":
        require_clean_partial_controller_evidence(
            inspection,
            manifest,
            job_id,
            observed_final_time,
            require_integer(allocation["nodes"], f"manifest {job_id} nodes", minimum=1)
            * EXPECTED_RANKS_PER_NODE,
        )
    return job_id, case_id, segment, result


def cross_validate_accounting(
    rows: list[dict[str, str]],
    reservations: list[dict[str, object]],
    manifests: list[dict[str, object]],
    counts: dict[str, object],
) -> dict[str, dict[str, object]]:
    """Cross-bind every recorded or cancelled terminal accounting identity."""

    expected = {
        "transactions": 0,
        "reservations": len(reservations),
        "active_reservations": 0,
        "ledger_rows": len(rows),
        "manifests": len(manifests),
    }
    selected = {key: counts.get(key) for key in COUNT_KEYS}
    if selected != expected:
        raise ValueError(f"reconciliation counts differ from live stores: {selected} != {expected}")
    ledger_by_job = {row["job_id"]: row for row in rows}
    reservations_by_manifest = {str(item["manifest"]): item for item in reservations}
    if len(reservations_by_manifest) != len(reservations):
        raise ValueError("reservations duplicate a manifest path")
    manifests_by_job: dict[str, dict[str, object]] = {}
    recorded_jobs: set[str] = set()
    manifest_paths: set[str] = set()
    for manifest in manifests:
        job_id, case_id, segment, result = manifest_identity(manifest)
        manifest_path = str(manifest.get("_manifest_path"))
        if manifest_path in manifest_paths:
            raise ValueError("manifests duplicate a retained path")
        manifest_paths.add(manifest_path)
        reservation = reservations_by_manifest.get(manifest_path)
        if reservation is None:
            raise ValueError(f"manifest {case_id}/{segment} lacks reservation evidence")
        common_identity = {"case_id": case_id, "segment": segment}
        if any(reservation.get(key) != value for key, value in common_identity.items()):
            raise ValueError(f"reservation identity differs for manifest {case_id}/{segment}")
        if reservation.get("job_id") != job_id:
            raise ValueError(f"reservation job ID differs for manifest {case_id}/{segment}")
        if manifest.get("state") == "cancelled":
            if reservation.get("state") != "cancelled":
                raise ValueError(f"cancelled manifest {case_id}/{segment} reservation differs")
            if job_id is not None and job_id in ledger_by_job:
                raise ValueError(f"cancelled manifest {case_id}/{segment} has a ledger row")
            continue
        if job_id is None:
            raise ValueError(f"recorded manifest {case_id}/{segment} lacks a job ID")
        if job_id in manifests_by_job:
            raise ValueError(f"manifests duplicate job {job_id}")
        manifests_by_job[job_id] = manifest
        recorded_jobs.add(job_id)
        row = ledger_by_job.get(job_id)
        if row is None:
            raise ValueError(f"manifest {job_id} lacks ledger evidence")
        expected_identity = {"case_id": case_id, "segment": segment, "result": result}
        if any(row.get(key) != value for key, value in expected_identity.items()):
            raise ValueError(f"ledger identity differs for job {job_id}")
        if any(reservation.get(key) != value for key, value in expected_identity.items()):
            raise ValueError(f"reservation identity differs for job {job_id}")
        if reservation.get("state") != "recorded":
            raise ValueError(f"recorded reservation state differs for job {job_id}")
        if int(row["nodes"]) != reservation.get("nodes"):
            raise ValueError(f"reservation nodes differ for job {job_id}")
        if reservation.get("requested_walltime") != row.get("requested_walltime"):
            raise ValueError(f"reservation walltime differs for job {job_id}")
        try:
            reservation_actual = Decimal(str(reservation["actual_node_hours"]))
        except (KeyError, InvalidOperation) as error:
            raise ValueError(f"reservation actual use is invalid for job {job_id}") from error
        if canonical_node_hours(
            reservation_actual, f"reservation actual use for job {job_id}"
        ) != row["actual_node_hours"]:
            raise ValueError(f"reservation actual use differs for job {job_id}")
        accounting = manifest["accounting"]
        allocation = manifest["allocation"]
        command = manifest.get("command")
        paths = manifest.get("paths")
        run = manifest.get("run")
        assert isinstance(accounting, dict)
        assert isinstance(allocation, dict)
        if not all(isinstance(item, dict) for item in (command, paths, run)):
            raise ValueError(f"manifest provenance is incomplete for job {job_id}")
        assert isinstance(command, dict)
        assert isinstance(paths, dict)
        assert isinstance(run, dict)
        case_name = require_nonempty_string(run.get("case_name"), f"manifest {job_id} case name")
        if (
            row.get("case_name") != case_name
            or reservation.get("case_name") != case_name
            or accounting.get("case_name") != case_name
        ):
            raise ValueError(f"case name differs for job {job_id}")
        if (
            allocation.get("nodes") != int(row["nodes"])
            or allocation.get("requested_walltime") != row["requested_walltime"]
        ):
            raise ValueError(f"manifest allocation differs for job {job_id}")
        for key in (
            "state",
            "exit_code",
            "nodes",
            "elapsed_seconds",
            "actual_node_hours",
            "submitted_utc",
            "completed_utc",
            "result",
        ):
            if str(accounting.get(key)) != row.get(key):
                raise ValueError(f"manifest accounting {key} differs for job {job_id}")
        for key in (
            "executable_revision",
            "executable_sha256",
            "input_revision",
            "input_file",
        ):
            if (
                str(accounting.get(key)) != row.get(key)
                or str(command.get(key)) != row.get(key)
            ):
                raise ValueError(f"manifest provenance {key} differs for job {job_id}")
        if str(accounting.get("output_dir")) != row.get("output_dir") or str(
            paths.get("output_dir")
        ) != row.get("output_dir"):
            raise ValueError(f"manifest provenance output_dir differs for job {job_id}")
    if recorded_jobs != set(ledger_by_job):
        raise ValueError("ledger and recorded-manifest job sets differ")
    if manifest_paths != set(reservations_by_manifest):
        raise ValueError("reservation and manifest path sets differ")
    return manifests_by_job


def require_empty_transaction_stores(root: Path) -> None:
    """Require no Stage I or recost transaction journal."""

    accounting = root / "accounting"
    for name, label in (
        (f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions", "Stage I transaction store"),
        (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions",
            "recost transaction store",
        ),
    ):
        path = accounting / name
        require_directory(path, label)
        with absolute_descriptor(
            path, label, flags=os.O_RDONLY | os.O_DIRECTORY
        ) as descriptor:
            entries = sorted(item.name for item in os.scandir(descriptor))
        if entries:
            raise ValueError(f"{label} is not empty: {', '.join(entries)}")


def parse_barrier_segments(value: object, manifests_by_job: dict[str, dict[str, object]]
                           ) -> list[dict[str, str]]:
    """Validate the exact just-recorded segment or drained-wave membership."""

    if not isinstance(value, list) or not value:
        raise ValueError("barrier recorded segments must be a nonempty list")
    retained = []
    jobs: set[str] = set()
    for index, item in enumerate(value):
        segment = require_exact_keys(
            item, {"case_id", "segment", "job_id", "result"},
            f"barrier segment {index}",
        )
        parsed = {
            "case_id": require_safe_id(segment["case_id"], f"barrier segment {index} case"),
            "segment": require_safe_id(segment["segment"], f"barrier segment {index} segment"),
            "job_id": require_job_id(segment["job_id"], f"barrier segment {index} job"),
            "result": require_nonempty_string(segment["result"], f"barrier segment {index} result"),
        }
        if parsed["job_id"] in jobs:
            raise ValueError("barrier recorded segments duplicate a job")
        jobs.add(parsed["job_id"])
        manifest = manifests_by_job.get(parsed["job_id"])
        if manifest is None:
            raise ValueError(f"barrier job {parsed['job_id']} lacks a retained manifest")
        job_id, case_id, name, result = manifest_identity(manifest)
        if parsed != {
            "case_id": case_id,
            "segment": name,
            "job_id": job_id,
            "result": result,
        }:
            raise ValueError(f"barrier job {job_id} identity differs from its manifest")
        retained.append(parsed)
    return retained


def parse_scheduler_evidence(
    root: Path,
    bindings: object,
    barrier: list[dict[str, str]],
    ledger_suffix: list[dict[str, str]],
    manifests_by_job: dict[str, dict[str, object]],
    tracker: InputTracker,
    request_timestamp: datetime,
) -> list[dict[str, object]]:
    """Authenticate one strict sacct row for every barrier job."""

    if not isinstance(bindings, list) or not bindings:
        raise ValueError("scheduler evidence bindings must be a nonempty list")
    barrier_jobs = [item["job_id"] for item in barrier]
    suffix_jobs = [row["job_id"] for row in ledger_suffix]
    if barrier_jobs != suffix_jobs:
        raise ValueError("barrier order differs from the exact recorded ledger suffix")
    ledger_by_job = {row["job_id"]: row for row in ledger_suffix}
    retained = []
    jobs: set[str] = set()
    for index, binding in enumerate(bindings):
        relative, expected = input_binding(binding, f"scheduler binding {index}")
        path = root_path(root, relative.as_posix(), f"scheduler binding {index}")
        if path.parent != root / "accounting":
            raise ValueError("scheduler evidence must be a direct child of accounting")
        payload = tracker.read(
            path,
            expected,
            f"scheduler evidence {relative}",
            expected_mode=0o644,
        )
        try:
            lines = payload.decode().splitlines()
        except UnicodeDecodeError as error:
            raise ValueError("scheduler evidence is not UTF-8") from error
        if len(lines) != 1:
            raise ValueError("scheduler evidence must contain exactly one row")
        fields = lines[0].split("|")
        if len(fields) != 8 or any(not field or field.strip() != field for field in fields):
            raise ValueError("scheduler evidence row is malformed")
        job_id, job_name, state, exit_code, nodes, elapsed, submitted, completed = fields
        require_job_id(job_id, f"scheduler evidence {index} job ID")
        if index >= len(barrier_jobs) or job_id != barrier_jobs[index]:
            raise ValueError("scheduler evidence order differs from the ledger suffix")
        if path.name != f"{job_id}.stage_i.sacct.txt":
            raise ValueError("scheduler evidence filename differs from its job ID")
        if job_id in jobs:
            raise ValueError(f"scheduler evidence duplicates job {job_id}")
        jobs.add(job_id)
        manifest = manifests_by_job.get(job_id)
        if manifest is None or job_id not in set(barrier_jobs):
            raise ValueError(f"scheduler evidence job {job_id} is not in the barrier")
        _, case_id, segment, _ = manifest_identity(manifest)
        accounting = manifest["accounting"]
        allocation = manifest["allocation"]
        run = manifest["run"]
        assert isinstance(accounting, dict)
        assert isinstance(allocation, dict)
        assert isinstance(run, dict)
        expected_name = f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_{case_id}_{segment}"
        expected_fields = [
            job_id,
            expected_name,
            str(accounting["state"]),
            str(accounting["exit_code"]),
            str(allocation["nodes"]),
            str(accounting["elapsed_seconds"]),
            str(accounting["submitted_utc"]),
            str(accounting["completed_utc"]),
        ]
        if fields != expected_fields:
            raise ValueError(f"scheduler evidence differs from manifest for job {job_id}")
        ledger = ledger_by_job[job_id]
        ledger_fields = [
            ledger["job_id"],
            expected_name,
            ledger["state"],
            ledger["exit_code"],
            ledger["nodes"],
            ledger["elapsed_seconds"],
            ledger["submitted_utc"],
            ledger["completed_utc"],
        ]
        if fields != ledger_fields:
            raise ValueError(f"scheduler evidence differs from ledger suffix for job {job_id}")
        if state not in TERMINAL_SCHEDULER_STATES:
            raise ValueError(f"scheduler evidence job {job_id} is not terminal")
        submitted_time = parse_scheduler_timestamp(
            submitted, f"scheduler evidence {job_id} submit time"
        )
        completed_time = parse_scheduler_timestamp(
            completed, f"scheduler evidence {job_id} completion time"
        )
        if completed_time <= submitted_time:
            raise ValueError(f"scheduler evidence job {job_id} chronology is invalid")
        elapsed_seconds = require_integer(
            int(elapsed), f"scheduler evidence {job_id} elapsed seconds"
        )
        observed_seconds = int((completed_time - submitted_time).total_seconds())
        if elapsed_seconds > observed_seconds + int(SCHEDULER_TIME_TOLERANCE.total_seconds()):
            raise ValueError(f"scheduler evidence job {job_id} elapsed exceeds chronology")
        if completed_time > request_timestamp + SCHEDULER_TIME_TOLERANCE:
            raise ValueError(f"scheduler evidence job {job_id} completes after the request")
        if root == DEFAULT_ROOT:
            environment = hardened_child_environment("SLURM_")
            try:
                live = subprocess.run(
                    [
                        str(SACCT),
                        "-X",
                        "-j",
                        job_id,
                        "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                        "ElapsedRaw,Submit,End",
                        "-n",
                        "-P",
                    ],
                    check=True,
                    capture_output=True,
                    text=True,
                    env=environment,
                ).stdout
            except (OSError, subprocess.CalledProcessError) as error:
                raise ValueError(
                    f"live scheduler evidence is unavailable for job {job_id}"
                ) from error
            live_rows = [
                line.rstrip("|").split("|")
                for line in live.splitlines()
                if line.rstrip("|").split("|")[0] == job_id
            ]
            if live_rows != [fields]:
                raise ValueError(f"live scheduler evidence differs for job {job_id}")
        retained.append(
            {
                "path": relative.as_posix(),
                "sha256": expected,
                "job_id": job_id,
                "job_name": job_name,
                "state": state,
                "exit_code": exit_code,
                "nodes": int(nodes),
                "elapsed_seconds": elapsed_seconds,
                "submitted_utc": submitted_time.isoformat(),
                "completed_utc": completed_time.isoformat(),
            }
        )
    if jobs != set(barrier_jobs):
        raise ValueError("scheduler evidence job set differs from the barrier")
    retained_by_job = {str(item["job_id"]): item for item in retained}
    return [retained_by_job[job_id] for job_id in barrier_jobs]


def resolution_cell_count(value: object, label: str) -> int:
    """Parse one exact three-dimensional matrix resolution."""

    retained = require_nonempty_string(value, label)
    match = re.fullmatch(r"([1-9][0-9]*)x([1-9][0-9]*)x([1-9][0-9]*)", retained)
    if match is None:
        raise ValueError(f"{label} must use NxNxN positive-integer dimensions")
    dimensions = [int(item) for item in match.groups()]
    cells = math.prod(dimensions)
    if cells <= 0:
        raise ValueError(f"{label} cell count is invalid")
    return cells


def parse_matrix(value: object) -> tuple[dict[str, dict[str, object]], Decimal]:
    """Validate the frozen Stage I mapped matrix."""

    if not isinstance(value, dict):
        raise ValueError("Stage I matrix must be an object")
    authorization = value.get("authorization")
    cases = value.get("cases")
    if not isinstance(authorization, dict) or not isinstance(cases, list):
        raise ValueError("Stage I matrix lacks authorization or cases")
    try:
        project_ceiling = Decimal(str(authorization["project_budget_node_hours"]))
    except (KeyError, InvalidOperation) as error:
        raise ValueError("Stage I matrix project ceiling is invalid") from error
    if not project_ceiling.is_finite() or project_ceiling <= 0:
        raise ValueError("Stage I matrix project ceiling is invalid")
    rule = authorization.get("execution_rule")
    if not isinstance(rule, str) or "R17 last" not in rule:
        raise ValueError("Stage I matrix does not retain the R17-last rule")
    retained: dict[str, dict[str, object]] = {}
    for item in cases:
        if not isinstance(item, dict):
            raise ValueError("Stage I matrix case must be an object")
        case_id = require_safe_id(item.get("id"), "Stage I matrix case ID")
        if case_id in retained:
            raise ValueError(f"Stage I matrix duplicates {case_id}")
        require_nonempty_string(item.get("name"), f"Stage I matrix {case_id} name")
        require_relative_path(item.get("input"), f"Stage I matrix {case_id} input")
        cell_count = resolution_cell_count(
            item.get("resolution"), f"Stage I matrix {case_id} resolution"
        )
        try:
            estimated_value = Decimal(str(item["estimated_node_hours"]))
        except (KeyError, InvalidOperation) as error:
            raise ValueError(
                f"Stage I matrix {case_id} estimated node-hours is invalid"
            ) from error
        estimated = require_decimal(
            decimal_string(estimated_value),
            f"Stage I matrix {case_id} estimated node-hours",
            positive=True,
        )
        item = dict(item)
        item["_estimated_node_hours"] = estimated
        item["_cell_count"] = cell_count
        retained[case_id] = item
    if set(retained) != AUTHORIZED_CASE_IDS:
        raise ValueError("Stage I matrix case set differs from R02-R17")
    return retained, project_ceiling


def expand_promoted_profiles(value: object) -> dict[str, set[int]]:
    """Expand F113-style exact/range node-profile keys."""

    if not isinstance(value, dict) or not value:
        raise ValueError("promoted node profiles must be a nonempty object")
    retained: dict[str, set[int]] = {}
    for key, nodes in value.items():
        if not isinstance(key, str) or not isinstance(nodes, list) or not nodes:
            raise ValueError("promoted node profile entry is invalid")
        parsed_nodes = {
            require_integer(item, f"promoted node profile {key}", minimum=1)
            for item in nodes
        }
        match = re.fullmatch(r"(R[0-9]{2})(?:-(R[0-9]{2}))?", key)
        if match is None:
            raise ValueError(f"promoted node profile key is invalid: {key}")
        first = int(match.group(1)[1:])
        last = int(match.group(2)[1:]) if match.group(2) else first
        if first > last:
            raise ValueError(f"promoted node profile range is reversed: {key}")
        for number in range(first, last + 1):
            case_id = f"R{number:02d}"
            if case_id not in AUTHORIZED_CASE_IDS or case_id in retained:
                raise ValueError(f"promoted node profile case is invalid or duplicated: {case_id}")
            retained[case_id] = parsed_nodes
    return retained


def parse_ceiling_evidence(
    value: object,
) -> tuple[Decimal, Decimal, int, dict[str, set[int]], str, str]:
    """Validate F113-style promoted envelope and concurrency controls."""

    if not isinstance(value, dict):
        raise ValueError("ceiling evidence must be an object")
    if value.get("record_type") != "stage-i-controller-transition-evidence":
        raise ValueError("ceiling evidence record type differs")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("ceiling evidence execution epoch differs")
    if value.get("checkpoint") != "F-113":
        raise ValueError("ceiling evidence checkpoint differs from F-113")
    parse_utc_timestamp(value.get("generated_utc"), "F113 generation timestamp")
    implementation = value.get("implementation")
    controls = value.get("promoted_controls")
    if not isinstance(implementation, dict) or not isinstance(controls, dict):
        raise ValueError("ceiling evidence lacks implementation or promoted controls")
    promoted_revision = require_revision(
        implementation.get("promoted_commit"), "F113 promoted controller revision"
    )
    helper = implementation.get("stage_i_helper")
    if not isinstance(helper, dict):
        raise ValueError("ceiling evidence lacks its historical Stage I helper binding")
    promoted_helper_sha256 = require_sha256(
        helper.get("sha256"), "F113 promoted Stage I helper SHA-256"
    )
    if helper.get("path") != STAGE_I_RELATIVE.as_posix():
        raise ValueError("ceiling evidence historical Stage I helper path differs")
    try:
        envelope = Decimal(str(controls["campaign_budget_node_hours"]))
        project = Decimal(str(controls["project_budget_node_hours"]))
    except (KeyError, InvalidOperation) as error:
        raise ValueError("ceiling evidence budgets are invalid") from error
    if not envelope.is_finite() or not project.is_finite() or envelope <= 0 or project <= 0:
        raise ValueError("ceiling evidence budgets must be finite and positive")
    if envelope > project:
        raise ValueError("promoted Stage I envelope exceeds project ceiling")
    lane_limit = require_integer(
        controls.get("standard_active_lane_limit"), "promoted lane limit", minimum=1
    )
    if lane_limit != 4:
        raise ValueError("promoted lane limit must remain exactly four")
    if require_integer(
        controls.get("prepared_packet_limit"), "promoted prepared-packet limit", minimum=1
    ) != 1:
        raise ValueError("promoted prepared-packet limit must remain one")
    r17_policy = controls.get("r17_policy")
    if (
        not isinstance(r17_policy, str)
        or "exclusive" not in r17_policy.lower()
        or "last" not in r17_policy.lower()
    ):
        raise ValueError("ceiling evidence does not preserve R17 exclusive/last")
    profiles = expand_promoted_profiles(controls.get("standard_multi_node_profiles"))
    if profiles != EXPECTED_PROMOTED_PROFILES:
        raise ValueError("promoted node profiles differ from the reviewed F113 policy")
    return (
        envelope,
        project,
        lane_limit,
        profiles,
        promoted_revision,
        promoted_helper_sha256,
    )


def parse_transition_publication_audit(
    value: object,
    root: Path,
    transition_path: Path,
    transition_sha256: str,
    request_timestamp: datetime,
) -> dict[str, object]:
    """Authenticate the exact reviewed publication of the promoted F113 controls."""

    retained = require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "published_utc",
            "artifact",
            "review",
        },
        "F113 publication audit",
    )
    if retained["schema_version"] != 1:
        raise ValueError("F113 publication audit schema version must be 1")
    if retained["record_type"] != "stage-i-controller-transition-publication-audit":
        raise ValueError("F113 publication audit record type differs")
    if retained["execution_epoch"] != EXECUTION_EPOCH:
        raise ValueError("F113 publication audit execution epoch differs")
    published = parse_utc_timestamp(
        retained["published_utc"], "F113 publication timestamp"
    )
    if published > request_timestamp:
        raise ValueError("recost request predates F113 publication")
    artifact = require_exact_keys(
        retained["artifact"], {"path", "sha256", "mode"}, "F113 publication artifact"
    )
    if artifact != {
        "path": str(transition_path),
        "sha256": transition_sha256,
        "mode": "0644",
    }:
        raise ValueError("F113 publication audit artifact binding differs")
    review = require_exact_keys(
        retained["review"],
        {"status", "reviewed_by", "authority"},
        "F113 publication review",
    )
    if review["status"] != "approved":
        raise ValueError("F113 publication audit is not approved")
    reviewer = require_nonempty_string(review["reviewed_by"], "F113 publication reviewer")
    if (
        review["authority"] != F113_REVIEW_AUTHORITY
        or reviewer not in F113_AUTHORIZED_REVIEWERS
    ):
        raise ValueError("F113 publication reviewer lacks committed campaign authority")
    if transition_path != root / F113_RELATIVE:
        raise ValueError("F113 transition path differs from the exact promoted path")
    return retained


def parse_storage_evidence(
    value: object,
    root: Path,
    request_timestamp: datetime,
    projected_growth_bytes: int,
    profile_projections_sha256: str,
) -> dict[str, object]:
    """Validate reviewed storage headroom for the exact authorized wave."""

    retained = require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "root",
            "measured_utc",
            "available_bytes",
            "retained_stage_i_bytes",
            "required_safety_bytes",
            "projected_authorized_wave_growth_bytes",
            "projection_method",
            "profile_projections_sha256",
        },
        "storage evidence",
    )
    if retained["schema_version"] != 1:
        raise ValueError("storage evidence schema version must be 1")
    if retained["record_type"] != "stage-i-storage-evidence":
        raise ValueError("storage evidence record type differs")
    if retained["execution_epoch"] != EXECUTION_EPOCH or retained["root"] != str(root):
        raise ValueError("storage evidence campaign binding differs")
    measured = parse_utc_timestamp(retained["measured_utc"], "storage evidence timestamp")
    now = datetime.now(timezone.utc)
    if measured > now + timedelta(minutes=5):
        raise ValueError("storage evidence timestamp is in the future")
    if now - measured > STORAGE_EVIDENCE_MAX_AGE:
        raise ValueError("storage evidence is stale")
    if request_timestamp < measured:
        raise ValueError("recost request predates storage evidence")
    available = require_integer(retained["available_bytes"], "available storage bytes")
    with absolute_descriptor(
        root, "campaign root", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as root_descriptor:
        live_storage = os.fstatvfs(root_descriptor)
    live_available = live_storage.f_bavail * live_storage.f_frsize
    if available > live_available:
        raise ValueError("storage evidence exceeds live available bytes")
    retained_bytes = require_integer(
        retained["retained_stage_i_bytes"], "retained Stage I bytes"
    )
    run_store = root / "runs/mks24-stage-i" / EXECUTION_EPOCH
    require_directory(run_store, "Stage I run store")
    live_retained = directory_tree_regular_bytes(run_store, "retained Stage I storage")
    if retained_bytes != live_retained:
        raise ValueError("storage evidence retained-byte count differs from live state")
    safety = require_integer(
        retained["required_safety_bytes"], "storage safety bytes", minimum=1
    )
    reviewed_growth = require_integer(
        retained["projected_authorized_wave_growth_bytes"],
        "reviewed authorized-wave growth bytes",
    )
    if reviewed_growth != projected_growth_bytes:
        raise ValueError("storage evidence growth differs from authorized profiles")
    if (
        retained["projection_method"] != STORAGE_PROJECTION_METHOD
        or retained["profile_projections_sha256"] != profile_projections_sha256
    ):
        raise ValueError("storage evidence does not bind computed observed-rate projections")
    if available < safety + reviewed_growth:
        raise ValueError("storage headroom is exhausted for the authorized wave")
    return retained


def require_live_storage_boundary(
    root: Path,
    available_bytes: int,
    retained_stage_i_bytes: int,
    required_safety_bytes: int,
    projected_growth_bytes: int,
) -> None:
    """Revalidate exact retained bytes and headroom at the final staging barrier."""

    run_store = root / "runs/mks24-stage-i" / EXECUTION_EPOCH
    if directory_tree_regular_bytes(run_store, "retained Stage I storage") != (
        retained_stage_i_bytes
    ):
        raise ValueError("retained Stage I storage changed before staged output creation")
    with absolute_descriptor(
        root, "campaign root", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as root_descriptor:
        live_storage = os.fstatvfs(root_descriptor)
    live_available = live_storage.f_bavail * live_storage.f_frsize
    if live_available < available_bytes:
        raise ValueError("live available storage fell below the reviewed evidence")
    if live_available < required_safety_bytes + projected_growth_bytes:
        raise ValueError("live storage headroom is exhausted at the staging barrier")


def require_directory_measurement_boundaries(
    measurements: tuple[DirectoryMeasurement, ...],
) -> None:
    """Revalidate every directory-size observation used by a projection."""

    for measurement in measurements:
        observed = directory_tree_regular_bytes(measurement.path, measurement.label)
        if observed != measurement.size_bytes:
            raise ValueError(
                f"{measurement.label} changed before staged output creation"
            )


def directory_tree_regular_bytes(path: Path, label: str) -> int:
    """Return regular bytes while validating excluded R02 bundle symlinks."""

    allowed_symlink_root = path / "bundles/R02/cases"

    def validate_excluded_symlink(
        descriptor: int, current: Path, entry: os.DirEntry[str], profile: os.stat_result
    ) -> None:
        symlink_path = current / entry.name
        if (
            not symlink_path.is_relative_to(allowed_symlink_root)
            or profile.st_uid != os.geteuid()
        ):
            raise ValueError(f"{label} contains an unauthorized symlink")
        target_text = os.readlink(entry.name, dir_fd=descriptor)
        target = Path(target_text)
        if (
            not target.is_absolute()
            or target != Path(os.path.normpath(str(target)))
            or not target.is_relative_to(path)
        ):
            raise ValueError(f"{label} contains an unsafe R02 bundle symlink")
        with absolute_descriptor(
            target, f"{label} R02 bundle symlink target", flags=os.O_RDONLY
        ) as target_descriptor:
            require_regular_profile(
                os.fstat(target_descriptor), f"{label} R02 bundle symlink target"
            )
        if os.readlink(entry.name, dir_fd=descriptor) != target_text:
            raise ValueError(f"{label} R02 bundle symlink changed during validation")

    def walk(descriptor: int, current: Path) -> int:
        total = 0
        for entry in os.scandir(descriptor):
            profile = entry.stat(follow_symlinks=False)
            if stat.S_ISREG(profile.st_mode):
                total += profile.st_size
            elif stat.S_ISDIR(profile.st_mode):
                child = os.open(
                    entry.name,
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                    dir_fd=descriptor,
                )
                try:
                    require_directory_profile(os.fstat(child), label)
                    total += walk(child, current / entry.name)
                finally:
                    os.close(child)
            elif stat.S_ISLNK(profile.st_mode):
                validate_excluded_symlink(descriptor, current, entry, profile)
            else:
                raise ValueError(f"{label} contains a special file")
        return total

    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        return walk(descriptor, path)


def absolute_root_member(root: Path, value: object, label: str) -> Path:
    """Require one absolute path lexically beneath the campaign root."""

    retained = require_nonempty_string(value, label)
    path = Path(retained)
    if (
        not path.is_absolute()
        or path != path.absolute()
        or path != Path(os.path.normpath(str(path)))
        or ".." in path.parts
    ):
        raise ValueError(f"{label} must be an absolute normalized path")
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath the campaign root") from error
    require_no_symlink_components(path, label)
    return path


def final_time(manifest: dict[str, object]) -> float:
    """Return one manifest's finite final inspection time."""

    inspection = manifest["scientific_inspection"]
    assert isinstance(inspection, dict)
    return require_finite_float(inspection.get("final_time"), "manifest final time")


def terminal_restart_binding(
    root: Path, manifest: dict[str, object], label: str
) -> dict[str, object]:
    """Validate and return one manifest's complete terminal rank-local restart group."""

    inspection = manifest.get("scientific_inspection")
    allocation = manifest.get("allocation")
    paths = manifest.get("paths")
    if not all(isinstance(item, dict) for item in (inspection, allocation, paths)):
        raise ValueError(f"{label} lacks inspection, allocation, or output provenance")
    assert isinstance(inspection, dict)
    assert isinstance(allocation, dict)
    assert isinstance(paths, dict)
    terminal = inspection.get("terminal_restart")
    if not isinstance(terminal, dict):
        raise ValueError(f"{label} lacks a terminal restart group")
    terminal_time = require_finite_float(
        inspection.get("terminal_restart_time"), f"{label} terminal restart time"
    )
    if abs(terminal_time - final_time(manifest)) > 1.0e-10:
        raise ValueError(f"{label} terminal restart time differs from final time")
    if terminal.get("storage") != "per_rank":
        raise ValueError(f"{label} terminal restart storage is not rank-local")
    rank_files = terminal.get("rank_files")
    expected_ranks = (
        require_integer(allocation.get("nodes"), f"{label} nodes", minimum=1)
        * EXPECTED_RANKS_PER_NODE
    )
    if not isinstance(rank_files, list) or len(rank_files) != expected_ranks:
        raise ValueError(f"{label} terminal restart rank inventory differs")
    output_dir = absolute_root_member(root, paths.get("output_dir"), f"{label} output")
    retained_paths: set[Path] = set()
    ordered_paths: list[Path] = []
    for index, rank_file in enumerate(rank_files):
        if not isinstance(rank_file, dict):
            raise ValueError(f"{label} terminal restart rank {index} is invalid")
        path = absolute_root_member(
            root, rank_file.get("path"), f"{label} terminal restart rank {index}"
        )
        try:
            path.relative_to(output_dir)
        except ValueError as error:
            raise ValueError(f"{label} terminal restart rank is outside parent output") from error
        if path in retained_paths:
            raise ValueError(f"{label} terminal restart rank inventory is duplicated")
        retained_paths.add(path)
        ordered_paths.append(path)
        require_sha256(
            rank_file.get("sha256"), f"{label} terminal restart rank {index} SHA-256"
        )
        require_integer(
            rank_file.get("size_bytes"),
            f"{label} terminal restart rank {index} size",
            minimum=1,
        )
    expected_names = [f"rank_{rank:08d}" for rank in range(expected_ranks)]
    if [path.parent.name for path in ordered_paths] != expected_names:
        raise ValueError(f"{label} terminal restart rank names differ")
    if (
        len({path.parent.parent for path in ordered_paths}) != 1
        or len({path.name for path in ordered_paths}) != 1
        or any(path.suffix != ".rst" for path in ordered_paths)
    ):
        raise ValueError(f"{label} terminal restart rank group layout differs")
    first = rank_files[0]
    assert isinstance(first, dict)
    if (
        terminal.get("path") != first.get("path")
        or terminal.get("sha256") != first.get("sha256")
        or terminal.get("size_bytes") != first.get("size_bytes")
    ):
        raise ValueError(f"{label} terminal restart primary binding differs")
    return terminal


def authenticate_terminal_restart(
    root: Path,
    manifest: dict[str, object],
    tracker: InputTracker,
    label: str,
) -> dict[str, object]:
    """Authenticate every file in one terminal rank-local restart group."""

    terminal = terminal_restart_binding(root, manifest, label)
    rank_files = terminal["rank_files"]
    assert isinstance(rank_files, list)
    for index, rank_file in enumerate(rank_files):
        assert isinstance(rank_file, dict)
        path = absolute_root_member(
            root, rank_file["path"], f"{label} terminal restart rank {index}"
        )
        tracker.authenticate(
            path,
            rank_file["sha256"],
            f"{label} terminal restart rank {index}",
            expected_mode=0o644,
        )
        with absolute_descriptor(
            path, f"{label} terminal restart rank {index}", flags=os.O_RDONLY
        ) as descriptor:
            if os.fstat(descriptor).st_size != rank_file["size_bytes"]:
                raise ValueError(f"{label} terminal restart rank {index} size differs")
    return terminal


def publication_timestamp(value: dict[str, object], label: str) -> datetime:
    """Return one observed or adopted publication timestamp."""

    record_type = value.get("record_type")
    if record_type == "observed-publication":
        key = "published_utc"
    elif record_type == "legacy-canonical-adoption":
        key = "adopted_utc"
    else:
        raise ValueError(f"{label} record type is not a promoted publication")
    return parse_utc_timestamp(value.get(key), f"{label} timestamp")


def parse_predecessor_recost(
    root: Path,
    recost_path: Path,
    recost_sha256: str,
    recost: object,
    review_path: Path | None,
    review_sha256: str | None,
    review: object,
    audit_path: Path,
    audit_sha256: str,
    audit: object,
    request_timestamp: datetime,
    current_artifact_name: str,
    current_checkpoint: str,
    tracker: InputTracker,
) -> dict[str, object]:
    """Authenticate the latest independently reviewed non-authorizing predecessor."""

    if recost_path.parent != root / "accounting" or recost_path.name.endswith(".staged"):
        raise ValueError("predecessor recost must be promoted under accounting")
    if recost_path.name == current_artifact_name:
        raise ValueError("predecessor recost collides with the current artifact")
    if audit_path != recost_path.with_name(
        f"{recost_path.name}.publication_audit.json"
    ):
        raise ValueError("predecessor recost publication-audit path differs")
    legacy_bootstrap = (
        recost_path == root / "accounting" / LEGACY_F114_RECOST_NAME
        and recost_sha256 == LEGACY_F114_RECOST_SHA256
        and audit_sha256 == LEGACY_F114_PUBLICATION_AUDIT_SHA256
    )
    if legacy_bootstrap:
        if (
            current_checkpoint != "F-117"
            or current_artifact_name
            != "mks24_stage_i_E03_forcing_policy_F117_recost_evidence.json"
        ):
            raise ValueError("legacy F114 bootstrap is reserved exactly for first schema-2 F-117")
        if review_path is not None or review_sha256 is not None or review is not None:
            raise ValueError("legacy F114 bootstrap must not fabricate an independent review")
    elif review_path != recost_path.with_name(
        f"{recost_path.name}.independent_review.json"
    ):
        raise ValueError("predecessor recost independent-review path differs")
    if not isinstance(recost, dict):
        raise ValueError("predecessor recost must be an object")
    if recost.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("predecessor recost execution epoch differs")
    generated = parse_utc_timestamp(
        recost.get("generated_utc"), "predecessor recost generation timestamp"
    )
    if generated > request_timestamp:
        raise ValueError("recost request predates its predecessor recost")
    predecessor_checkpoint = require_safe_id(
        recost.get("checkpoint"), "predecessor recost checkpoint ID"
    )
    predecessor_match = re.fullmatch(r"F-([0-9]+)", predecessor_checkpoint)
    current_match = re.fullmatch(r"F-([0-9]+)", current_checkpoint)
    predecessor_name_match = RECOST_ARTIFACT_PATTERN.fullmatch(recost_path.name)
    if (
        predecessor_name_match is not None
        and predecessor_checkpoint != f"F-{predecessor_name_match.group(1)}"
    ):
        raise ValueError("predecessor recost artifact name and checkpoint ID differ")
    if recost.get("artifact_name") not in {None, recost_path.name}:
        raise ValueError("predecessor recost embedded artifact name differs")
    if (
        predecessor_match is None
        or current_match is None
        or int(predecessor_match.group(1)) >= int(current_match.group(1))
    ):
        raise ValueError("predecessor recost checkpoint does not precede the request")
    if legacy_bootstrap:
        if (
            recost.get("schema_version") != 1
            or recost.get("record_type") != "stage-i-clean-partial-recost-checkpoint"
            or predecessor_checkpoint != "F-114"
        ):
            raise ValueError("legacy F114 bootstrap artifact identity differs")
        if not isinstance(audit, dict) or (
            audit.get("schema_version") != 1
            or audit.get("record_type") != "observed-publication"
            or audit.get("execution_epoch") != EXECUTION_EPOCH
            or audit.get("artifact")
            != {
                "path": str(recost_path),
                "sha256": recost_sha256,
                "mode": "0644",
                "links": 1,
            }
        ):
            raise ValueError("legacy F114 bootstrap publication audit differs")
        published = publication_timestamp(audit, "legacy F114 bootstrap publication")
        if published > request_timestamp or published < generated:
            raise ValueError("legacy F114 bootstrap publication chronology differs")
    elif (
        recost.get("schema_version") != 2
        or recost.get("record_type") != "stage-i-recost-recommendation-evidence"
        or recost.get("artifact_name") != recost_path.name
        or recost.get("authority")
        != {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("predecessor recost is not strict non-authorizing recommendation evidence")
    if legacy_bootstrap:
        reviewed = None
    else:
        assert review_path is not None
        assert review_sha256 is not None
        retained_review = require_exact_keys(
            review,
            {
                "schema_version",
                "record_type",
                "execution_epoch",
                "reviewed_utc",
                "decision",
                "reviewer",
                "candidate",
                "scope",
            },
            "predecessor recost independent review",
        )
        if (
            retained_review["schema_version"] != 1
            or retained_review["record_type"]
            != "stage-i-recost-recommendation-independent-review"
            or retained_review["execution_epoch"] != EXECUTION_EPOCH
            or retained_review["decision"] != "approved-for-publication"
            or retained_review["candidate"]
            != {"path": str(recost_path), "sha256": recost_sha256}
            or retained_review["scope"] != {"non_authorizing": True}
        ):
            raise ValueError("predecessor recost independent review differs")
        reviewer = require_exact_keys(
            retained_review["reviewer"],
            {"agent_id", "independent_from_generator"},
            "predecessor recost reviewer",
        )
        require_nonempty_string(reviewer["agent_id"], "predecessor recost reviewer agent ID")
        if reviewer["independent_from_generator"] is not True:
            raise ValueError("predecessor recost review is not independent")
        reviewed = parse_utc_timestamp(
            retained_review["reviewed_utc"], "predecessor recost review timestamp"
        )
        minimal_audit_keys = {
            "schema_version",
            "record_type",
            "execution_epoch",
            "published_utc",
            "artifact",
            "independent_review",
            "authority",
        }
        checkpoint_audit_keys = {
            *minimal_audit_keys,
            "transaction_id",
            "recost_recommendations",
            "counts",
            "generator",
            "scheduler_evidence",
            "source_bundle",
            "stage_i_helper",
            "utility",
            "forensic_copy",
            "publication",
            "generalized_publication_context",
        }
        if not isinstance(audit, dict) or frozenset(audit) not in {
            frozenset(minimal_audit_keys),
            frozenset(checkpoint_audit_keys),
        }:
            raise ValueError("predecessor recost publication audit schema differs")
        retained_audit = audit
        if (
            retained_audit["schema_version"] != 1
            or retained_audit["record_type"]
            != "stage-i-recost-recommendation-publication-audit"
            or retained_audit["execution_epoch"] != EXECUTION_EPOCH
            or retained_audit["authority"]
            != {
                "action_authority": False,
                "scheduler_mutation_authorized": False,
                "canonical_mutation_authorized": False,
            }
        ):
            raise ValueError("predecessor recost publication audit differs")
        if set(retained_audit) == checkpoint_audit_keys:
            if (
                retained_audit["recost_recommendations"] != recost.get("recommendations")
                or retained_audit["generalized_publication_context"].get("authority")
                != recost.get("authority")
            ):
                raise ValueError("predecessor checkpoint publication context differs")
        require_declared_binding(
            retained_audit["artifact"],
            recost_path,
            recost_sha256,
            "predecessor recost publication artifact",
        )
        require_declared_binding(
            retained_audit["independent_review"],
            review_path,
            review_sha256,
            "predecessor recost publication review",
        )
        published = parse_utc_timestamp(
            retained_audit["published_utc"], "predecessor recost publication timestamp"
        )
        if published > request_timestamp:
            raise ValueError("recost request predates predecessor recost publication")
        if reviewed < generated or published < reviewed:
            raise ValueError("predecessor recost publication predates its artifact")

    accounting = root / "accounting"
    with absolute_descriptor(
        accounting, "accounting directory", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        audit_names = sorted(
            entry.name
            for entry in os.scandir(descriptor)
            if entry.name.endswith("_recost_evidence.json.publication_audit.json")
        )
    publications: list[tuple[datetime, Path]] = []
    for name in audit_names:
        candidate = accounting / name
        with absolute_descriptor(
            candidate, f"recost publication audit inventory {name}", flags=os.O_RDONLY
        ) as descriptor:
            require_regular_profile(
                os.fstat(descriptor), f"recost publication audit inventory {name}"
            )
            payload = b""
            while True:
                block = os.read(descriptor, 1024 * 1024)
                if not block:
                    break
                payload += block
        value = parse_json(payload, f"recost publication audit inventory {name}")
        if (
            not isinstance(value, dict)
            or value.get("record_type")
            != "stage-i-recost-recommendation-publication-audit"
            or value.get("execution_epoch") != EXECUTION_EPOCH
        ):
            continue
        tracker.read(
            candidate,
            sha256_bytes(payload),
            f"strict recost publication audit inventory {name}",
            expected_mode=0o444,
        )
        candidate_published = parse_utc_timestamp(
            value.get("published_utc"), f"recost publication audit inventory {name}"
        )
        candidate_artifact = candidate.with_name(
            name.removesuffix(".publication_audit.json")
        )
        artifact_binding = value.get("artifact")
        if not isinstance(artifact_binding, dict):
            raise ValueError(f"recost publication audit lacks artifact binding: {name}")
        with absolute_descriptor(
            candidate_artifact,
            f"promoted recost inventory {candidate_artifact.name}",
            flags=os.O_RDONLY,
        ) as descriptor:
            require_regular_profile(
                os.fstat(descriptor),
                f"promoted recost inventory {candidate_artifact.name}",
                expected_mode=0o444,
            )
            artifact_payload = b""
            while True:
                block = os.read(descriptor, 1024 * 1024)
                if not block:
                    break
                artifact_payload += block
        artifact_sha256 = sha256_bytes(artifact_payload)
        tracker.read(
            candidate_artifact,
            artifact_sha256,
            f"strict promoted recost inventory {candidate_artifact.name}",
            expected_mode=0o444,
        )
        if not artifact_payload:
            raise ValueError(f"promoted recost inventory is empty: {candidate_artifact.name}")
        if {
            "path": artifact_binding.get("path"),
            "sha256": artifact_binding.get("sha256"),
            "mode": artifact_binding.get("mode"),
            "links": artifact_binding.get("links"),
        } != {
            "path": str(candidate_artifact),
            "sha256": artifact_sha256,
            "mode": "0444",
            "links": 1,
        }:
            raise ValueError(f"recost publication audit artifact binding differs: {name}")
        publications.append((candidate_published, candidate))
    if legacy_bootstrap and publications:
        raise ValueError("legacy F114 bootstrap was already consumed by schema-2 publication")
    if legacy_bootstrap:
        return {
            "path": str(recost_path),
            "artifact_name": recost_path.name,
            "sha256": recost_sha256,
            "publication_audit_path": str(audit_path),
            "publication_audit_sha256": audit_sha256,
            "independent_review_path": None,
            "independent_review_sha256": None,
            "checkpoint": predecessor_checkpoint,
            "generated_utc": generated.isoformat(),
            "published_utc": published.isoformat(),
            "bootstrap": "exact-retained-legacy-F114-once",
        }
    if not publications:
        raise ValueError("recost publication audit inventory is empty")
    latest_timestamp = max(item[0] for item in publications)
    latest = [item[1] for item in publications if item[0] == latest_timestamp]
    if len(latest) != 1:
        raise ValueError("latest promoted recost publication is ambiguous")
    if latest[0] != audit_path:
        raise ValueError("predecessor recost is not the latest promoted recost")
    if latest_timestamp > request_timestamp:
        raise ValueError("recost request predates the latest promoted recost")
    return {
        "path": str(recost_path),
        "artifact_name": recost_path.name,
        "sha256": recost_sha256,
        "publication_audit_path": str(audit_path),
        "publication_audit_sha256": audit_sha256,
        "independent_review_path": str(review_path),
        "independent_review_sha256": review_sha256,
        "checkpoint": predecessor_checkpoint,
        "generated_utc": generated.isoformat(),
        "published_utc": published.isoformat(),
    }


def manifest_parent_path(manifest: dict[str, object]) -> Path | None:
    """Return the authenticated parent-manifest path declared by one segment."""

    command = manifest.get("command")
    if not isinstance(command, dict):
        raise ValueError("recorded manifest lacks command provenance")
    parent = command.get("parent_segment")
    if parent is None:
        return None
    if not isinstance(parent, dict):
        raise ValueError("recorded manifest parent segment is invalid")
    path = Path(require_nonempty_string(parent.get("manifest"), "parent manifest path"))
    if not path.is_absolute() or path != path.absolute():
        raise ValueError("parent manifest path must be absolute and normalized")
    job_id, case_id, segment, result = manifest_identity(manifest)
    if (
        parent.get("case_id") != case_id
        or parent.get("segment") == segment
        or parent.get("result") not in {"accepted", "clean_partial"}
        or require_finite_float(parent.get("final_time"), "parent final time")
        >= final_time(manifest)
    ):
        raise ValueError(f"manifest {job_id} parent lineage binding differs")
    return path


def authenticated_case_lineages(
    root: Path,
    manifests: list[dict[str, object]],
) -> dict[str, list[dict[str, object]]]:
    """Build each case's unique latest authenticated restart lineage."""

    by_case: dict[str, list[dict[str, object]]] = {}
    all_by_case: dict[str, list[dict[str, object]]] = {}
    by_path: dict[Path, dict[str, object]] = {}
    for manifest in manifests:
        _, case_id, _, result = manifest_identity(manifest)
        all_by_case.setdefault(case_id, []).append(manifest)
        if result not in SCIENTIFIC_RESULTS:
            continue
        path = Path(str(manifest.get("_manifest_path", ""))).absolute()
        if not path.is_absolute() or path in by_path:
            raise ValueError("manifest lineage path is absent or duplicated")
        by_path[path] = manifest
        by_case.setdefault(case_id, []).append(manifest)

    retained: dict[str, list[dict[str, object]]] = {}
    for case_id, retained_candidates in by_case.items():
        entries_by_index: dict[int, list[dict[str, object]]] = {}
        for item in all_by_case[case_id]:
            index, _, _ = parse_segment(
                manifest_identity(item)[2], f"{case_id} recorded segment"
            )
            entries_by_index.setdefault(index, []).append(item)

        def indexes_are_non_scientific(start: int, stop: int) -> bool:
            return all(
                len(entries_by_index.get(index, [])) == 1
                and manifest_identity(entries_by_index[index][0])[3]
                not in SCIENTIFIC_RESULTS
                for index in range(start, stop)
            )

        candidates = list(retained_candidates)
        r12_historical_inventory = [
            item
            for item in candidates
            if manifest_identity(item) == R12_HISTORICAL_INVENTORY_IDENTITY
        ]
        if (
            case_id == "R12"
            and len(r12_historical_inventory) == 1
            and any(
                manifest_identity(item)[2] == R12_FRESH_RERUN_SEGMENT
                and manifest_parent_path(item) is None
                for item in candidates
            )
        ):
            candidates = [
                item for item in candidates if item is not r12_historical_inventory[0]
            ]
        maximum = max(final_time(item) for item in candidates)
        terminals = [item for item in candidates if final_time(item) == maximum]
        if len(terminals) != 1:
            raise ValueError(f"{case_id} latest recorded lineage is ambiguous")
        lineage = []
        current = terminals[0]
        visited: set[Path] = set()
        while True:
            current_path = Path(str(current["_manifest_path"])).absolute()
            if current_path in visited:
                raise ValueError(f"{case_id} recorded lineage contains a cycle")
            visited.add(current_path)
            lineage.append(current)
            parent_path = manifest_parent_path(current)
            if parent_path is None:
                break
            parent_path = absolute_root_member(
                root, str(parent_path), f"{case_id} recorded lineage parent"
            )
            parent = by_path.get(parent_path)
            if parent is None or manifest_identity(parent)[1] != case_id:
                raise ValueError(f"{case_id} recorded lineage parent is missing")
            parent_record = current["command"]["parent_segment"]
            assert isinstance(parent_record, dict)
            parent_job, _, parent_segment, parent_result = manifest_identity(parent)
            _, _, current_segment, _ = manifest_identity(current)
            current_index, current_start, _ = parse_segment(
                current_segment, f"{case_id} recorded lineage child segment"
            )
            parent_index, _, _ = parse_segment(
                parent_segment, f"{case_id} recorded lineage parent segment"
            )
            parent_command = parent.get("command")
            if not isinstance(parent_command, dict):
                raise ValueError(f"{case_id} recorded lineage parent lacks command provenance")
            terminal = terminal_restart_binding(
                root, parent, f"{case_id} recorded lineage parent"
            )
            terminal_rank_files = terminal["rank_files"]
            assert isinstance(terminal_rank_files, list)
            if (
                parent_record.get("execution_epoch") != EXECUTION_EPOCH
                or parent_record.get("segment") != parent_segment
                or parent_record.get("result") != parent_result
                or parent_record.get("job_id", parent_job) != parent_job
                or parent_record.get("final_time") != final_time(parent)
                or parent_record.get("restart_time") != final_time(parent)
                or parent_record.get("executable_sha256")
                != parent_command.get("executable_sha256")
                or parent_record.get("input_sha256") != parent_command.get("input_sha256")
                or parent_record.get("restart_sha256") != terminal.get("sha256")
                or parent_record.get("restart_files")
                != [item.get("path") for item in terminal_rank_files]
                or current_index <= parent_index
                or not indexes_are_non_scientific(parent_index + 1, current_index)
                or abs(current_start - final_time(parent)) > 1.0e-6
            ):
                raise ValueError(f"{case_id} recorded lineage parent identity differs")
            current = parent
        lineage.reverse()
        root_index, root_start, _ = parse_segment(
            manifest_identity(lineage[0])[2], f"{case_id} recorded lineage root segment"
        )
        r12_fresh_root = (
            case_id == "R12"
            and len(r12_historical_inventory) == 1
            and manifest_identity(lineage[0])[2] == R12_FRESH_RERUN_SEGMENT
            and root_index == 1
            and abs(root_start) <= 1.0e-12
        )
        retry_root = (
            root_index > 0
            and abs(root_start) <= 1.0e-12
            and indexes_are_non_scientific(0, root_index)
        )
        if (
            root_index != 0 or abs(root_start) > 1.0e-12
        ) and not r12_fresh_root and not retry_root:
            raise ValueError(f"{case_id} recorded lineage root is not s00 from t=0")
        if {Path(str(item["_manifest_path"])).absolute() for item in lineage} != {
            Path(str(item["_manifest_path"])).absolute() for item in candidates
        }:
            raise ValueError(f"{case_id} recorded manifests do not form one lineage")
        retained[case_id] = lineage
    return retained


def directory_inventory(
    path: Path, tracker: InputTracker, label: str
) -> list[dict[str, str]]:
    """Authenticate and return one complete flat retained-directory inventory."""

    require_directory(path, label)
    entries = []
    with absolute_descriptor(
        path, label, flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        names = sorted(entry.name for entry in os.scandir(descriptor))
    if not names:
        raise ValueError(f"{label} must not be empty")
    for name in names:
        item = path / name
        require_no_symlink_components(item, label)
        with absolute_descriptor(item, f"{label} file {name}", flags=os.O_RDONLY) as descriptor:
            profile = os.fstat(descriptor)
            if not stat.S_ISREG(profile.st_mode):
                raise ValueError(f"{label} must contain only direct regular files")
            digest = sha256_descriptor(descriptor)
        tracker.authenticate(item, digest, f"{label} file {name}", expected_mode=0o644)
        entries.append(
            {"name": name, "mode": f"{stat.S_IMODE(profile.st_mode):04o}", "sha256": digest}
        )
    return entries


def directory_inventory_sha256(
    path: Path, tracker: InputTracker, label: str
) -> str:
    """Authenticate one flat retained directory and return its inventory digest."""

    entries = directory_inventory(path, tracker, label)
    return sha256_bytes((json.dumps(entries, sort_keys=True) + "\n").encode())


def lineage_summary_sha256(lineages: dict[str, list[dict[str, object]]]) -> str:
    """Return a stable digest for exact retained case lineages."""

    value = {
        case_id: [
            {
                "manifest": item["_manifest_path"],
                "sha256": item["_manifest_sha256"],
                "job_id": manifest_identity(item)[0],
                "segment": manifest_identity(item)[2],
                "result": manifest_identity(item)[3],
                "final_time": final_time(item),
            }
            for item in lineage
        ]
        for case_id, lineage in sorted(lineages.items())
    }
    return sha256_bytes((json.dumps(value, sort_keys=True) + "\n").encode())


def observed_storage_projections(
    root: Path,
    manifests: list[dict[str, object]],
    matrix: dict[str, dict[str, object]],
    profiles: list[dict[str, object]],
) -> tuple[int, list[dict[str, object]], tuple[DirectoryMeasurement, ...]]:
    """Project authorized storage from authenticated observed Stage I output rates."""

    bases: list[dict[str, object]] = []
    measurements: dict[Path, DirectoryMeasurement] = {}
    for manifest in manifests:
        job_id, case_id, segment, result = manifest_identity(manifest)
        if result not in SCIENTIFIC_RESULTS:
            continue
        manifest_path = Path(str(manifest["_manifest_path"])).absolute()
        paths = manifest.get("paths")
        if not isinstance(paths, dict):
            raise ValueError(f"manifest {job_id} lacks output provenance")
        output = absolute_root_member(root, paths.get("output_dir"), f"manifest {job_id} output")
        expected_output = manifest_path.parent.parent / "output"
        if output != expected_output:
            raise ValueError(f"manifest {job_id} output path differs from its run directory")
        _, start, _ = parse_segment(segment, f"manifest {job_id} segment")
        interval = Decimal(str(final_time(manifest))) - Decimal(str(start))
        if interval <= 0:
            raise ValueError(f"manifest {job_id} measured interval is invalid")
        size = directory_tree_regular_bytes(output, f"manifest {job_id} measured output")
        measurements[output] = DirectoryMeasurement(
            output, size, f"manifest {job_id} measured output"
        )
        if size <= 0:
            continue
        cells = require_integer(matrix[case_id]["_cell_count"], f"{case_id} matrix cells", minimum=1)
        bases.append(
            {
                "job_id": job_id,
                "case_id": case_id,
                "segment": segment,
                "output_dir": str(output),
                "observed_bytes": size,
                "observed_cells": cells,
                "observed_simulation_interval": decimal_string(interval),
            }
        )
    if not bases:
        raise ValueError("no nonempty authenticated Stage I output exists for storage projection")

    projections: list[dict[str, object]] = []
    total = 0
    for profile in profiles:
        case_id = str(profile["case_id"])
        segment = str(profile["segment"])
        _, start, target = parse_segment(segment, f"next profile {case_id} segment")
        interval = Decimal(str(target)) - Decimal(str(start))
        cells = require_integer(matrix[case_id]["_cell_count"], f"{case_id} matrix cells", minimum=1)
        candidates: list[tuple[Decimal, dict[str, object]]] = []
        for basis in bases:
            projection = (
                Decimal(int(basis["observed_bytes"]))
                * Decimal(cells)
                * interval
                / Decimal(int(basis["observed_cells"]))
                / Decimal(str(basis["observed_simulation_interval"]))
            )
            candidates.append((projection, basis))
        projected, basis = max(
            candidates, key=lambda item: (item[0], str(item[1]["job_id"]))
        )
        required = max(MINIMUM_PROFILE_STORAGE_BYTES, decimal_ceiling(projected))
        if case_id == R17_CASE_ID:
            required = max(required, R17_MINIMUM_RETAINED_BYTES)
        retained = require_integer(
            profile["estimated_storage_bytes"],
            f"next profile {case_id} estimated storage bytes",
            minimum=1,
        )
        if retained != required:
            raise ValueError(
                f"next profile {case_id} storage estimate differs from observed-rate projection"
            )
        projections.append(
            {
                "case_id": case_id,
                "segment": segment,
                "method": STORAGE_PROJECTION_METHOD,
                "basis_job_id": basis["job_id"],
                "basis_case_id": basis["case_id"],
                "basis_segment": basis["segment"],
                "basis_output_dir": basis["output_dir"],
                "basis_observed_bytes": basis["observed_bytes"],
                "basis_observed_cells": basis["observed_cells"],
                "basis_observed_simulation_interval": basis[
                    "observed_simulation_interval"
                ],
                "projected_cells": cells,
                "projected_simulation_interval": decimal_string(interval),
                "raw_projected_bytes": decimal_string(projected),
                "required_projected_bytes": required,
            }
        )
        total += required
    return (
        total,
        sorted(projections, key=lambda item: (str(item["case_id"]), str(item["segment"]))),
        tuple(sorted(measurements.values(), key=lambda item: str(item.path))),
    )


def compact_json_sha256(value: object) -> str:
    """Return the qualification publisher's canonical compact JSON digest."""

    return sha256_bytes(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    )


def publication_json_sha256(value: object) -> str:
    """Return the controller publisher's canonical indented JSON digest."""

    return sha256_bytes(
        (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode()
    )


def parse_aware_scheduler_timestamp(value: object, label: str) -> datetime:
    """Parse one scheduler timestamp while retaining its explicit offset."""

    retained = require_nonempty_string(value, label)
    try:
        parsed = datetime.fromisoformat(retained.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include an explicit offset")
    return parsed


def parse_optional_account_timestamp(value: str, label: str) -> datetime | None:
    """Parse one optional account-wide scheduler timestamp."""

    if value in ACCOUNT_MISSING_TIMESTAMPS:
        return None
    return parse_aware_scheduler_timestamp(value, label)


def parse_account_scheduler_text(payload: bytes) -> list[dict[str, object]]:
    """Parse retained complete all-users top-level account scheduler evidence."""

    try:
        lines = payload.decode().splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("R17 account scheduler evidence is not UTF-8") from error
    if not lines or lines[0] != ACCOUNT_SCHEDULER_HEADER:
        raise ValueError("R17 account scheduler evidence has an invalid header")
    records: list[dict[str, object]] = []
    seen: set[str] = set()
    for line in lines[1:]:
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != len(ACCOUNT_SACCT_FIELDS):
            raise ValueError("R17 account scheduler evidence has invalid columns")
        (
            job_id,
            job_name,
            state_value,
            exit_code,
            nodes_text,
            elapsed_text,
            submitted,
            started,
            ended,
            partition,
            account,
            owner,
        ) = fields
        state_parts = state_value.split()
        state = state_parts[0].split("+")[0] if state_parts else ""
        if (
            LIVE_JOB_ID_PATTERN.fullmatch(job_id) is None
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
            raise ValueError(
                "R17 account scheduler evidence has invalid numeric values"
            ) from error
        submitted_time = parse_optional_account_timestamp(
            submitted, "R17 account scheduler submit time"
        )
        started_time = parse_optional_account_timestamp(
            started, "R17 account scheduler start time"
        )
        ended_time = parse_optional_account_timestamp(
            ended, "R17 account scheduler end time"
        )
        if (
            nodes < 0
            or elapsed < 0
            or submitted_time is None
            or (
                started_time is None
                and (
                    elapsed != 0
                    or state in {"RUNNING", "COMPLETING", "SUSPENDED", "COMPLETED"}
                    or (ended_time is not None and submitted_time > ended_time)
                )
            )
            or (started_time is not None and nodes <= 0)
            or (started_time is not None and submitted_time > started_time)
            or (
                started_time is not None
                and ended_time is not None
                and (
                    ended_time < started_time
                    or abs((ended_time - started_time).total_seconds() - elapsed) > 1.0
                )
            )
        ):
            raise ValueError("R17 account scheduler evidence chronology differs")
        seen.add(job_id)
        records.append(
            {
                "job_id": job_id,
                "job_name": job_name,
                "state": state,
                "exit_code": exit_code,
                "nodes": nodes,
                "elapsed_seconds": elapsed,
                "submit_utc": submitted,
                "start_utc": started,
                "end_utc": ended,
                "partition": partition,
                "account": account,
                "owner": owner,
            }
        )
    if not records:
        raise ValueError("R17 account scheduler evidence is empty")
    return sorted(records, key=lambda record: str(record["job_id"]))


def expected_account_query_contract(start: datetime, end: datetime) -> dict[str, object]:
    """Return the qualification publisher's exact all-account scheduler query."""

    if start.tzinfo is None or end.tzinfo is None or not start < end:
        raise ValueError("R17 account scheduler query interval differs")
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


def authenticate_r17_decomposition(
    value: object, output_inventory_sha256: str
) -> dict[str, object]:
    """Authenticate schema-2's inline complete 1728-block/64-rank decomposition."""

    evidence = require_exact_keys(
        value,
        {
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
        },
        "R17 inline decomposition evidence",
    )
    inventory = evidence["complete_block_rank_inventory"]
    if not isinstance(inventory, list) or len(inventory) != 64:
        raise ValueError("R17 inline decomposition rank inventory is incomplete")
    expected_locations = {
        (lx1, lx2, lx3, 0)
        for lx1 in range(12)
        for lx2 in range(12)
        for lx3 in range(12)
    }
    observed_locations: list[tuple[int, int, int, int]] = []
    for rank, item in enumerate(inventory):
        record = require_exact_keys(
            item,
            {"rank", "rank_name", "logical_meshblocks"},
            "R17 inline decomposition rank record",
        )
        locations = record["logical_meshblocks"]
        if (
            record["rank"] != rank
            or record["rank_name"] != f"rank_{rank:08d}"
            or not isinstance(locations, list)
            or len(locations) != 27
            or locations != sorted(locations)
            or any(
                not isinstance(location, list)
                or len(location) != 4
                or any(isinstance(index, bool) or not isinstance(index, int) for index in location)
                or tuple(location) not in expected_locations
                for location in locations
            )
            or len({tuple(location) for location in locations}) != 27
        ):
            raise ValueError("R17 inline decomposition does not retain 27 unique blocks per rank")
        observed_locations.extend(tuple(location) for location in locations)
    checks = require_exact_keys(
        evidence["checks"],
        {
            "exact_resolution",
            "exact_rank_count",
            "exact_meshblocks_per_rank",
            "complete_unique_logical_inventory",
        },
        "R17 inline decomposition checks",
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
        or len(set(observed_locations)) != 1728
        or set(observed_locations) != expected_locations
        or checks
        != {
            "exact_resolution": True,
            "exact_rank_count": True,
            "exact_meshblocks_per_rank": True,
            "complete_unique_logical_inventory": True,
        }
    ):
        raise ValueError("R17 inline decomposition evidence differs")
    return evidence


def authenticate_r17_account_exclusivity(
    value: object,
    raw_payload: bytes,
    job_id: str,
    scheduler_fields: list[str],
    completed: datetime,
    request_timestamp: datetime,
) -> dict[str, object]:
    """Reproduce schema-2's retained full-interval account exclusivity proof."""

    evidence = require_exact_keys(
        value,
        {
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
        },
        "R17 account exclusivity evidence",
    )
    records = parse_account_scheduler_text(raw_payload)
    target_records = [record for record in records if record["job_id"] == job_id]
    if len(target_records) != 1:
        raise ValueError("R17 account exclusivity evidence lacks one qualification job")
    target = target_records[0]
    target_without_owner = {
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
    started = parse_aware_scheduler_timestamp(
        target["start_utc"], "R17 account qualification start"
    )
    ended = parse_aware_scheduler_timestamp(
        target["end_utc"], "R17 account qualification end"
    )
    measured = parse_utc_timestamp(
        evidence["measured_utc"], "R17 account exclusivity measurement"
    )
    visibility = require_exact_keys(
        evidence["visibility_contract"],
        {"private_data", "all_users_job_visibility"},
        "R17 account visibility contract",
    )
    private_data = require_nonempty_string(
        visibility["private_data"], "R17 account PrivateData setting"
    )
    overlapping = []
    for record in records:
        record_start = parse_optional_account_timestamp(
            str(record["start_utc"]), "R17 account overlap start"
        )
        record_end = parse_optional_account_timestamp(
            str(record["end_utc"]), "R17 account overlap end"
        )
        if (
            record_start is not None
            and record_start < ended
            and (record_end is None or record_end > started)
        ):
            overlapping.append(str(record["job_id"]))
    if (
        evidence["schema_version"] != 1
        or evidence["record_type"] != "stage-i-r17-account-exclusivity-evidence"
        or evidence["execution_epoch"] != EXECUTION_EPOCH
        or measured < completed
        or measured > request_timestamp
        or evidence["query_contract"] != expected_account_query_contract(started, ended)
        or visibility["all_users_job_visibility"] is not True
        or any(
            setting in {"all", "jobs"}
            for setting in (
                item.strip().casefold() for item in private_data.split(",")
            )
        )
        or evidence["raw_account_scheduler_sha256"] != sha256_bytes(raw_payload)
        or evidence["qualification_job"] != target_without_owner
        or evidence["qualification_job_sha256"] != compact_json_sha256(target_without_owner)
        or evidence["account_jobs"] != records
        or evidence["account_jobs_sha256"] != compact_json_sha256(records)
        or evidence["overlapping_job_ids"] != overlapping
        or overlapping != [job_id]
        or evidence["exclusive_entire_execution_interval"] is not True
        or target["job_name"] != scheduler_fields[1]
        or target["state"] != scheduler_fields[2]
        or target["exit_code"] != scheduler_fields[3]
        or target["nodes"] != int(scheduler_fields[4])
        or target["elapsed_seconds"] != int(scheduler_fields[5])
        or target["submit_utc"] != scheduler_fields[6]
        or target["end_utc"] != scheduler_fields[7]
        or target["partition"] != PARTITION
        or str(target["account"]).casefold() != ACCOUNT.casefold()
    ):
        raise ValueError("R17 account exclusivity evidence differs")
    return evidence


def parse_athinput_payload(payload: bytes, label: str) -> dict[str, str]:
    """Parse one exact Athena input deck for frozen-contract reproduction."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} is not UTF-8") from error
    block = ""
    values: dict[str, str] = {}
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            if not block:
                raise ValueError(f"{label} has an empty parameter block")
            continue
        key, separator, value = line.partition("=")
        if not block or not separator or not key.strip() or not value.strip():
            raise ValueError(f"{label} has an invalid parameter line")
        qualified = f"{block}/{key.strip()}"
        if qualified in values:
            raise ValueError(f"{label} duplicates parameter {qualified}")
        values[qualified] = value.strip()
    if not values:
        raise ValueError(f"{label} is empty")
    return values


def authenticate_r17_absolute_file_record(
    value: object,
    label: str,
    tracker: InputTracker,
    *,
    expected_path: Path | None = None,
    expected_sha256: str | None = None,
    expected_mode: int | None = None,
) -> tuple[Path, dict[str, object]]:
    """Authenticate one producer retained-file record with exact bytes and size."""

    record = require_exact_keys(value, {"path", "sha256", "size_bytes"}, label)
    path = Path(require_nonempty_string(record["path"], f"{label} path"))
    if (
        not path.is_absolute()
        or path != path.absolute()
        or path != Path(os.path.normpath(str(path)))
        or ".." in path.parts
    ):
        raise ValueError(f"{label} path is not absolute and normalized")
    require_no_symlink_components(path, label)
    digest = require_sha256(record["sha256"], f"{label} SHA-256")
    size = require_integer(record["size_bytes"], f"{label} size", minimum=1)
    if expected_path is not None and path != expected_path:
        raise ValueError(f"{label} path differs")
    if expected_sha256 is not None and digest != expected_sha256:
        raise ValueError(f"{label} SHA-256 differs")
    tracker.authenticate(path, digest, label, expected_mode=expected_mode)
    with absolute_descriptor(path, label, flags=os.O_RDONLY) as descriptor:
        if os.fstat(descriptor).st_size != size:
            raise ValueError(f"{label} size differs")
    return path, record


def authenticate_r17_root_binding(
    value: object,
    root: Path,
    label: str,
    tracker: InputTracker,
    *,
    expected_path: Path | None = None,
    expected_mode: int | None = None,
) -> tuple[Path, str, bytes]:
    """Authenticate one exact producer root-relative binding."""

    relative, digest = input_binding(value, f"{label} binding")
    path = root_path(root, relative.as_posix(), label)
    if expected_path is not None and path != expected_path:
        raise ValueError(f"{label} path differs")
    payload = tracker.read(path, digest, label, expected_mode=expected_mode)
    return path, digest, payload


def authenticate_r17_rank_inventory_contract(
    value: object,
    inventory_sha256: object,
    root: Path,
    label: str,
    tracker: InputTracker,
) -> tuple[list[dict[str, str]], str]:
    """Authenticate exactly one ordered terminal file for each R17 rank."""

    if not isinstance(value, list) or len(value) != 64:
        raise ValueError(f"{label} must retain exactly 64 files")
    inventory = []
    ranks: set[str] = set()
    paths: set[Path] = set()
    for index, binding in enumerate(value):
        relative, digest = input_binding(binding, f"{label} record {index}")
        path = root_path(root, relative.as_posix(), f"{label} file {index}")
        rank_parts = [
            part for part in path.relative_to(root).parts
            if re.fullmatch(r"rank_[0-9]{8}", part)
        ]
        if len(rank_parts) != 1 or path in paths:
            raise ValueError(f"{label} rank inventory is incomplete or ambiguous")
        ranks.add(rank_parts[0])
        paths.add(path)
        tracker.authenticate(path, digest, f"{label} file {index}", expected_mode=0o644)
        with absolute_descriptor(path, f"{label} file {index}", flags=os.O_RDONLY) as descriptor:
            if os.fstat(descriptor).st_size <= 0:
                raise ValueError(f"{label} file {index} is empty")
        inventory.append({"path": relative.as_posix(), "sha256": digest})
    digest = sha256_bytes((json.dumps(inventory, sort_keys=True) + "\n").encode())
    if (
        ranks != {f"rank_{rank:08d}" for rank in range(64)}
        or inventory != sorted(inventory, key=lambda item: item["path"])
        or digest != require_sha256(inventory_sha256, f"{label} inventory SHA-256")
    ):
        raise ValueError(f"{label} inventory digest or ordering differs")
    return inventory, digest


def reproduce_r17_frozen_science_build_contract(
    wave: dict[str, object],
    packet: dict[str, object],
    profile: dict[str, object],
    matrix_sha256: str,
    build_manifest_inventory_sha256: str,
    tracker: InputTracker,
) -> tuple[dict[str, object], dict[str, object]]:
    """Reconstruct the producer's exact frozen science/build contract."""

    provenance = require_exact_keys(
        wave["provenance"],
        {
            "source", "source_bundle", "matrix", "executable",
            "build_manifest", "qualification_helper",
        },
        "prepared R17 provenance",
    )
    source = require_exact_keys(
        provenance["source"], {"directory", "revision"}, "prepared R17 source"
    )
    source_revision = require_revision(source["revision"], "prepared R17 source revision")
    bundle = require_exact_keys(
        provenance["source_bundle"],
        {"path", "sha256", "size_bytes", "verified_revisions"},
        "prepared R17 source bundle",
    )
    bundle_record = {
        key: bundle[key] for key in ("path", "sha256", "size_bytes")
    }
    authenticate_r17_absolute_file_record(
        bundle_record,
        "prepared R17 source bundle",
        tracker,
        expected_path=Path(str(profile["source_bundle"])),
        expected_sha256=str(profile["source_bundle_sha256"]),
        expected_mode=0o644,
    )
    matrix = require_exact_keys(
        provenance["matrix"], {"path", "sha256", "size_bytes"}, "prepared R17 matrix"
    )
    authenticate_r17_absolute_file_record(
        matrix,
        "prepared R17 matrix",
        tracker,
        expected_sha256=matrix_sha256,
    )
    executable = require_exact_keys(
        provenance["executable"],
        {"path", "sha256", "size_bytes", "revision"},
        "prepared R17 executable",
    )
    executable_record = {
        key: executable[key] for key in ("path", "sha256", "size_bytes")
    }
    authenticate_r17_absolute_file_record(
        executable_record,
        "prepared R17 executable",
        tracker,
        expected_path=Path(str(profile["executable"])),
        expected_sha256=str(profile["executable_sha256"]),
        expected_mode=0o755,
    )
    helper = require_exact_keys(
        provenance["qualification_helper"],
        {"path", "sha256", "size_bytes", "revision", "committed"},
        "prepared R17 qualification helper",
    )
    helper_record = {
        key: helper[key] for key in ("path", "sha256", "size_bytes")
    }
    authenticate_r17_absolute_file_record(
        helper_record, "prepared R17 qualification helper", tracker
    )
    helper_revision = require_revision(
        helper["revision"], "prepared R17 qualification helper revision"
    )
    revisions = bundle["verified_revisions"]
    if (
        source_revision != profile["executable_revision"]
        or source_revision != profile["input_revision"]
        or executable["revision"] != source_revision
        or helper["committed"] is not True
        or revisions != sorted({source_revision, helper_revision})
    ):
        raise ValueError("prepared R17 source/helper revision binding differs")
    build = require_exact_keys(
        provenance["build_manifest"],
        {"path", "athena_sha256", "environment", "inventory", "inventory_sha256"},
        "prepared R17 build manifest",
    )
    if (
        build["path"] != profile["build_manifest"]
        or build["inventory_sha256"] != build_manifest_inventory_sha256
        or build["inventory"] != profile["_authenticated_build_manifest_inventory"]
    ):
        raise ValueError("prepared R17 build-manifest binding differs")
    intent = require_exact_keys(
        packet["execution_intent"],
        {
            "packet_id", "project_root", "qualification_root", "case_id",
            "case_name", "profile_class", "target_time", "run_basename",
            "job_name", "scientific_policy", "input", "provenance_sha256",
            "allocation", "paths", "execution_contract_sha256",
            "authenticated_commands", "execution_intent_sha256",
        },
        "prepared R17 execution intent",
    )
    intent_sha256 = require_sha256(
        packet["execution_intent_sha256"], "prepared R17 execution-intent SHA-256"
    )
    intent_without_digest = dict(intent)
    embedded_intent_sha256 = intent_without_digest.pop("execution_intent_sha256")
    intent_core = dict(intent_without_digest)
    intent_core.pop("authenticated_commands")
    execution_contract_sha256 = intent_core.pop("execution_contract_sha256")
    allocation = require_exact_keys(
        intent["allocation"],
        {
            "nodes", "walltime", "walltime_seconds", "athena_walltime",
            "athena_walltime_seconds", "ranks_per_node", "cpus_per_task",
        },
        "prepared R17 allocation",
    )
    if (
        packet != {"execution_intent": intent, "execution_intent_sha256": intent_sha256}
        or embedded_intent_sha256 != intent_sha256
        or compact_json_sha256(intent_without_digest) != intent_sha256
        or compact_json_sha256(intent_core) != execution_contract_sha256
        or intent["case_id"] != R17_CASE_ID
        or intent["profile_class"] != "scale_separation_384x384x768"
        or intent["scientific_policy"] != "active_hardwall"
        or intent["target_time"] != 0.25
        or intent["run_basename"] != "qualification_R17_n08"
        or intent["job_name"] != "cglq_R17_n08"
        or intent["provenance_sha256"] != wave["provenance_sha256"]
        or allocation["nodes"] != 8
        or allocation["ranks_per_node"] != 8
        or allocation["cpus_per_task"] != profile["cpus_per_task"]
        or profile["nodes"] != 8
        or profile["ranks_per_node"] != 8
        or profile["time_tlim_target"] != 0.25
    ):
        raise ValueError("prepared R17 execution intent differs")
    input_record = require_exact_keys(
        intent["input"], {"path", "sha256", "size_bytes"}, "prepared R17 input"
    )
    _, authenticated_input = authenticate_r17_absolute_file_record(
        input_record,
        "prepared R17 input",
        tracker,
        expected_sha256=str(profile["input_sha256"]),
    )
    input_path = Path(str(authenticated_input["path"]))
    with absolute_descriptor(input_path, "prepared R17 input", flags=os.O_RDONLY) as descriptor:
        input_payload = b""
        chunks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
        input_payload = b"".join(chunks)
    parameter_contract = parse_athinput_payload(input_payload, "prepared R17 input")
    parameter_contract["job/basename"] = str(intent["run_basename"])
    parameter_contract["time/tlim"] = str(intent["target_time"])
    contract = {
        "case_id": R17_CASE_ID,
        "case_name": intent["case_name"],
        "profile_class": intent["profile_class"],
        "resolution": "384x384x768",
        "mesh_shape": [384, 384, 768],
        "meshblock_shape": [32, 32, 64],
        "target_time": intent["target_time"],
        "scientific_policy": intent["scientific_policy"],
        "run_basename": intent["run_basename"],
        "source_revision": source_revision,
        "source_bundle_sha256": bundle["sha256"],
        "matrix_sha256": matrix["sha256"],
        "input_sha256": input_record["sha256"],
        "provenance_sha256": wave["provenance_sha256"],
        "execution_intent_sha256": intent_sha256,
        "execution_contract_sha256": execution_contract_sha256,
        "parameter_contract": parameter_contract,
        "parameter_contract_sha256": compact_json_sha256(parameter_contract),
        "executable_revision": executable["revision"],
        "executable_sha256": executable["sha256"],
        "build_manifest_inventory_sha256": build_manifest_inventory_sha256,
    }
    return contract, intent


def r17_relative_difference(first: float, second: float) -> float:
    """Reproduce the qualification producer's stable relative difference."""

    return abs(first - second) / max(abs(first), abs(second), 1.0)


def r17_expected_times(target: float, cadence: float) -> list[float]:
    """Reproduce the qualification producer's exact cadence schedule."""

    count = int(math.floor(target / cadence + 1.0e-10))
    retained = [index * cadence for index in range(count + 1)]
    if abs(retained[-1] - target) > 1.0e-10:
        retained.append(target)
    else:
        retained[-1] = target
    return retained


def authenticate_r17_product_groups(
    groups: object,
    times: object,
    root: Path,
    tracker: InputTracker,
    *,
    kind: str,
) -> list[list[dict[str, object]]]:
    """Authenticate complete 64-rank producer snapshot or restart groups."""

    if (
        not isinstance(groups, list)
        or not groups
        or not isinstance(times, list)
        or len(groups) != len(times)
    ):
        raise ValueError(f"R17 {kind} groups or times differ")
    retained_groups = []
    names: set[str] = set()
    parsed_times = [
        require_finite_float(value, f"R17 {kind} time") for value in times
    ]
    if parsed_times != sorted(parsed_times):
        raise ValueError(f"R17 {kind} times are not ordered")
    for group_index, (group_value, group_time) in enumerate(zip(groups, parsed_times)):
        group = require_exact_keys(
            group_value, {"name", "rank_files"}, f"R17 {kind} group {group_index}"
        )
        name = require_nonempty_string(group["name"], f"R17 {kind} group name")
        rank_files = group["rank_files"]
        if name in names or not isinstance(rank_files, list) or len(rank_files) != 64:
            raise ValueError(f"R17 {kind} group rank inventory differs")
        names.add(name)
        retained_files = []
        observed_ranks = set()
        logical_locations = []
        for rank, record_value in enumerate(rank_files):
            record = require_exact_keys(
                record_value,
                {"path", "sha256", "size_bytes", "inspection"},
                f"R17 {kind} group {group_index} rank {rank}",
            )
            path, _ = authenticate_r17_absolute_file_record(
                {key: record[key] for key in ("path", "sha256", "size_bytes")},
                f"R17 {kind} group {group_index} rank {rank}",
                tracker,
                expected_mode=0o644,
            )
            rank_names = [
                part for part in path.relative_to(root).parts
                if re.fullmatch(r"rank_[0-9]{8}", part)
            ]
            if rank_names != [f"rank_{rank:08d}"]:
                raise ValueError(f"R17 {kind} group rank pathname differs")
            observed_ranks.add(rank_names[0])
            inspection = record["inspection"]
            if (
                not isinstance(inspection, dict)
                or abs(
                    require_finite_float(
                        inspection.get("time"), f"R17 {kind} inspection time"
                    )
                    - group_time
                )
                > 1.0e-10
            ):
                raise ValueError(f"R17 {kind} group inspection time differs")
            if kind == "snapshot":
                locations = inspection.get("logical_locations")
                minima = inspection.get("positive_variable_minima")
                if (
                    inspection.get("meshblock_count") != 27
                    or inspection.get("hard_bound_violation_cells") != 0
                    or not isinstance(locations, list)
                    or len(locations) != 27
                    or not isinstance(minima, dict)
                    or any(
                        require_finite_float(
                            minima.get(variable),
                            f"R17 snapshot positive {variable}",
                        )
                        <= 0.0
                        for variable in ("dens", "eint", "p_perp")
                    )
                    or require_finite_float(
                        inspection.get("minimum_mirror_hard_margin"),
                        "R17 snapshot mirror hard margin",
                    )
                    <= 0.0
                    or require_finite_float(
                        inspection.get("minimum_firehose_hard_margin"),
                        "R17 snapshot firehose hard margin",
                    )
                    <= 0.0
                ):
                    raise ValueError("R17 snapshot independently verified physics differs")
                logical_locations.extend(tuple(location) for location in locations)
            retained_files.append(dict(record))
        if observed_ranks != {f"rank_{rank:08d}" for rank in range(64)}:
            raise ValueError(f"R17 {kind} group rank inventory differs")
        if kind == "snapshot" and (
            len(logical_locations) != 1728
            or len(set(logical_locations)) != 1728
        ):
            raise ValueError("R17 snapshot logical meshblock inventory differs")
        retained_groups.append(retained_files)
    return retained_groups


def authenticate_r17_physics_contract(
    scientific: dict[str, object],
    root: Path,
    outputs: list[dict[str, str]],
    restarts: list[dict[str, str]],
    tracker: InputTracker,
) -> dict[str, object]:
    """Independently reproduce the producer's complete R17 physics gate."""

    if require_exact_keys(
        scientific["checks"], set(R17_SCIENTIFIC_CHECKS), "R17 scientific checks"
    ) != R17_SCIENTIFIC_CHECKS:
        raise ValueError("R17 scientific checks differ or contain a failed check")
    histories = {}
    for key, label in (("mhd_history", "R17 MHD history"), ("user_history", "R17 user history")):
        path, record = authenticate_r17_absolute_file_record(
            scientific[key], label, tracker, expected_mode=0o644
        )
        payload = tracker.read(
            path, str(record["sha256"]), label, expected_mode=0o644
        )
        histories[key] = parse_controller_history(payload, label)
    mhd = histories["mhd_history"]
    user = histories["user_history"]
    mhd_required = {
        "time", "mass", "tot-E", "lf_nstage", "lf_qface", "lf_qprcap",
        "lf_qpecap", "lf_qprwrk", "lf_qpewrk", "lf_hwproj", "lf_cpwrk",
        "lf_cawrk", *CONTROLLER_STRICT_LF_FAILURE_COLUMNS,
    }
    user_required = {"time", "mass", "hard_vol", "force_pwr", "force_work", "max_ndiv"}
    if mhd_required - set(mhd) or user_required - set(user):
        raise ValueError("R17 histories lack required independent physics columns")
    expected_times = r17_expected_times(0.25, 0.02)
    if (
        len(mhd["time"]) != len(user["time"])
        or len(mhd["time"]) != len(expected_times)
        or any(
            abs(observed - expected) > 1.0e-10
            for observed, expected in zip(mhd["time"], expected_times)
        )
        or any(
            abs(observed - expected) > 1.0e-10
            for observed, expected in zip(user["time"], expected_times)
        )
    ):
        raise ValueError("R17 histories fail the exact synchronized endpoint schedule")
    mass_relative_drift_max = max(
        [
            r17_relative_difference(mhd["mass"][0], value)
            for value in mhd["mass"]
        ]
        + [
            r17_relative_difference(user["mass"][0], value)
            for value in user["mass"]
        ]
    )
    mass_mismatch_max = max(
        r17_relative_difference(first, second)
        for first, second in zip(mhd["mass"], user["mass"])
    )
    strict_total = math.fsum(
        value
        for name in CONTROLLER_STRICT_LF_FAILURE_COLUMNS
        for value in mhd[name]
    )
    normalized_ct_divb_max = max(user["max_ndiv"])
    if (
        mass_relative_drift_max > 1.0e-12
        or mass_mismatch_max > 1.0e-12
        or strict_total != 0.0
        or any(value != 0.0 for value in user["hard_vol"])
        or any(value < 0.0 for value in user["max_ndiv"])
        or normalized_ct_divb_max >= float(R17_MAX_NORMALIZED_CT_DIVB_TEXT)
    ):
        raise ValueError("R17 independently reproduced mass, LF, hard-bound, or CT gate failed")
    count_columns = ("lf_nstage", "lf_qface", "lf_qprcap", "lf_qpecap", "lf_hwproj")
    if (
        any(
            value < 0.0 or not value.is_integer()
            for name in count_columns
            for value in mhd[name]
        )
        or any(
            later < earlier
            for name in count_columns
            for earlier, later in zip(mhd[name], mhd[name][1:])
        )
        or any(
            later <= earlier
            for name in ("lf_nstage", "lf_qface")
            for earlier, later in zip(mhd[name], mhd[name][1:])
        )
        or any(
            value > qface
            for name in ("lf_qprcap", "lf_qpecap")
            for value, qface in zip(mhd[name], mhd["lf_qface"])
        )
    ):
        raise ValueError("R17 independently reproduced LF interval-count gate failed")
    for name in ("lf_qprcap", "lf_qpecap"):
        if any(
            cap - previous_cap > qface - previous_qface
            for previous_cap, cap, previous_qface, qface in zip(
                mhd[name], mhd[name][1:], mhd["lf_qface"], mhd["lf_qface"][1:]
            )
        ):
            raise ValueError("R17 independently reproduced LF cap-increment gate failed")
    forcing_delta = user["force_work"][-1] - user["force_work"][0]
    energy_delta = mhd["tot-E"][-1] - mhd["tot-E"][0]
    pressure_work_delta = sum(
        abs(mhd[name][-1] - mhd[name][0]) for name in ("lf_cpwrk", "lf_cawrk")
    )
    lf_work_delta = sum(
        abs(mhd[name][-1] - mhd[name][0]) for name in ("lf_qprwrk", "lf_qpewrk")
    )
    if (
        abs(forcing_delta) < 0.25 * 1.0e-4
        or max(abs(value) for value in user["force_pwr"]) < 1.0e-4
        or pressure_work_delta < 0.25 * 1.0e-7
        or lf_work_delta < 0.25 * 1.0e-7
        or abs(user["force_work"][-1] - user["force_work"][-2]) < 0.02 * 1.0e-4
        or sum(
            abs(mhd[name][-1] - mhd[name][-2]) for name in ("lf_cpwrk", "lf_cawrk")
        )
        < pressure_work_delta * 1.0e-4
        or sum(
            abs(mhd[name][-1] - mhd[name][-2]) for name in ("lf_qprwrk", "lf_qpewrk")
        )
        < lf_work_delta * 1.0e-4
        or r17_relative_difference(energy_delta, forcing_delta) > 1.0e-8
    ):
        raise ValueError("R17 independently reproduced activity or closure gate failed")
    snapshot_groups = authenticate_r17_product_groups(
        scientific["snapshots"], scientific["snapshot_times"], root, tracker, kind="snapshot"
    )
    restart_groups = authenticate_r17_product_groups(
        scientific["restarts"], scientific["restart_times"], root, tracker, kind="restart"
    )
    snapshot_times = [float(value) for value in scientific["snapshot_times"]]
    restart_times = [float(value) for value in scientific["restart_times"]]
    if (
        len(snapshot_times) != 2
        or any(
            abs(observed - expected) > 1.0e-10
            for observed, expected in zip(snapshot_times, [0.0, 0.25])
        )
        or sum(abs(value - 0.25) <= 1.0e-10 for value in snapshot_times) != 1
        or sum(abs(value - 0.25) <= 1.0e-10 for value in restart_times) != 1
        or any(value > 0.25 + 1.0e-10 for value in restart_times)
    ):
        raise ValueError("R17 independently reproduced snapshot or restart endpoint gate failed")

    def terminal_inventory(group: list[dict[str, object]]) -> list[dict[str, str]]:
        return [
            {
                "path": Path(str(record["path"])).relative_to(root).as_posix(),
                "sha256": str(record["sha256"]),
            }
            for record in group
        ]

    terminal_snapshot = snapshot_groups[
        next(index for index, value in enumerate(snapshot_times) if abs(value - 0.25) <= 1.0e-10)
    ]
    terminal_restart = restart_groups[
        next(index for index, value in enumerate(restart_times) if abs(value - 0.25) <= 1.0e-10)
    ]
    if terminal_inventory(terminal_snapshot) != outputs or terminal_inventory(terminal_restart) != restarts:
        raise ValueError("R17 terminal product inventories differ from scientific groups")
    measurements = {
        "finite_rank_outputs": len(outputs),
        "mass_relative_drift_max": format(mass_relative_drift_max, ".17g"),
        "mhd_user_mass_mismatch_max": format(mass_mismatch_max, ".17g"),
        "lf_bad_counts_total": 0,
        "normalized_ct_divb_max": format(normalized_ct_divb_max, ".17g"),
        "normalized_ct_divb_threshold": R17_MAX_NORMALIZED_CT_DIVB_TEXT,
        "normalized_ct_divb_below_threshold": True,
    }
    if (
        require_exact_keys(
            scientific["physics_measurements"],
            R17_PHYSICS_MEASUREMENT_KEYS,
            "R17 physics measurements",
        )
        != measurements
        or scientific["forcing_closure_normalized_residual"]
        != r17_relative_difference(energy_delta, forcing_delta)
    ):
        raise ValueError("R17 physics measurements differ from independent reproduction")
    return measurements


def parse_r17_operational_qualification(
    value: object,
    root: Path,
    request_timestamp: datetime,
    profile: dict[str, object],
    matrix_sha256: str,
    qualification_path: Path,
    qualification_sha256: str,
    review_binding: object,
    tracker: InputTracker,
) -> dict[str, object]:
    """Authenticate the full exact producer schema-2 R17 qualification contract."""

    retained = require_exact_keys(
        value,
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
    expected_qualification_path = root / "accounting" / (
        f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_operational_qualification.json"
    )
    job_id = require_job_id(retained["job_id"], "R17 qualification job ID")
    completed = parse_aware_scheduler_timestamp(
        retained["completed_utc"], "R17 qualification completion"
    )
    measured = parse_utc_timestamp(
        retained["measured_utc"], "R17 qualification measurement"
    )
    measured_by = require_nonempty_string(
        retained["measured_by"], "R17 qualification measurement author"
    )
    if (
        qualification_path != expected_qualification_path
        or retained["schema_version"] != 2
        or retained["record_type"] != "stage-i-r17-operational-qualification"
        or retained["execution_epoch"] != EXECUTION_EPOCH
        or retained["state"] != "COMPLETED"
        or retained["exit_code"] != "0:0"
        or retained["nodes"] != 8
        or retained["ranks"] != 64
        or measured < completed.astimezone(timezone.utc)
        or measured > request_timestamp
    ):
        raise ValueError("R17 operational qualification identity differs")

    build_manifest = absolute_root_member(
        root, profile["build_manifest"], "R17 qualification build manifest"
    )
    build_inventory = directory_inventory(
        build_manifest, tracker, "R17 qualification build manifest"
    )
    build_inventory_sha256 = sha256_bytes(
        (json.dumps(build_inventory, sort_keys=True) + "\n").encode()
    )
    if (
        retained["executable_sha256"] != profile["executable_sha256"]
        or retained["build_manifest_inventory"] != build_inventory
        or retained["build_manifest_inventory_sha256"] != build_inventory_sha256
        or build_inventory_sha256 != profile["build_manifest_sha256"]
    ):
        raise ValueError("R17 operational qualification build binding differs")
    profile_with_inventory = dict(profile)
    profile_with_inventory["_authenticated_build_manifest_inventory"] = build_inventory

    wave_path, wave_sha256, wave_payload = authenticate_r17_root_binding(
        retained["prepared_wave"], root, "R17 prepared wave", tracker
    )
    wave = require_exact_keys(
        parse_json(wave_payload, "R17 prepared wave"),
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "policy", "provenance", "provenance_sha256", "waves",
        },
        "R17 prepared wave",
    )
    policy = require_exact_keys(
        wave["policy"],
        {
            "canonical_acceptance_eligible", "submission_action_available",
            "submission_performed", "operational_only", "max_concurrent_nodes",
            "selection_minimum_runtime_savings_seconds",
            "selection_maximum_node_hour_ratio", "selection_runtime_tie_fraction",
        },
        "R17 prepared-wave policy",
    )
    waves = wave["waves"]
    if (
        wave["schema_version"] != 1
        or wave["record_type"] != "cgl_lf_stage_i_qualification_wave"
        or wave["project_root"] != str(root)
        or wave["provenance_sha256"] != compact_json_sha256(wave["provenance"])
        or policy["canonical_acceptance_eligible"] is not False
        or policy["operational_only"] is not True
        or policy["max_concurrent_nodes"] != 8
        or not isinstance(waves, list)
        or len(waves) != 1
    ):
        raise ValueError("R17 prepared wave differs")
    wave_record = require_exact_keys(
        waves[0], {"wave", "total_nodes", "packets"}, "R17 prepared launch wave"
    )
    packets = wave_record["packets"]
    if (
        wave_record["wave"] != 1
        or wave_record["total_nodes"] != 8
        or not isinstance(packets, list)
        or len(packets) != 1
    ):
        raise ValueError("R17 prepared wave is not one exact exclusive packet")
    packet = require_exact_keys(
        packets[0], {"execution_intent", "execution_intent_sha256"}, "R17 prepared packet"
    )
    frozen, intent = reproduce_r17_frozen_science_build_contract(
        wave,
        packet,
        profile_with_inventory,
        matrix_sha256,
        build_inventory_sha256,
        tracker,
    )
    frozen_retained = require_exact_keys(
        retained["frozen_science_build_contract"], set(frozen), "R17 frozen science/build contract"
    )
    if frozen_retained != frozen:
        raise ValueError("R17 frozen science/build binding differs from prepared wave")

    qualification_evidence_path, _, qualification_evidence_payload = (
        authenticate_r17_root_binding(
            retained["qualification_evidence"],
            root,
            "R17 qualification evidence",
            tracker,
        )
    )
    qualification_evidence = require_exact_keys(
        parse_json(qualification_evidence_payload, "R17 qualification evidence"),
        {
            "schema_version", "record_type", "project_root", "qualification_root",
            "prepared_wave", "case_id", "target_time", "results",
        },
        "R17 qualification evidence",
    )
    evidence_wave = require_exact_keys(
        qualification_evidence["prepared_wave"],
        {"path", "sha256", "size_bytes"},
        "R17 qualification evidence prepared wave",
    )
    if (
        qualification_evidence["schema_version"] != 2
        or qualification_evidence["record_type"] != "cgl_lf_stage_i_qualification_evidence"
        or qualification_evidence["project_root"] != str(root)
        or qualification_evidence["case_id"] != R17_CASE_ID
        or qualification_evidence["target_time"] != 0.25
        or evidence_wave["path"] != str(wave_path)
        or evidence_wave["sha256"] != wave_sha256
    ):
        raise ValueError("R17 qualification evidence binding differs")

    outputs, outputs_sha256 = authenticate_r17_rank_inventory_contract(
        retained["rank_local_outputs"],
        retained["rank_local_output_inventory_sha256"],
        root,
        "R17 terminal rank-local outputs",
        tracker,
    )
    restarts, restarts_sha256 = authenticate_r17_rank_inventory_contract(
        retained["rank_local_restarts"],
        retained["rank_local_restart_inventory_sha256"],
        root,
        "R17 terminal rank-local restarts",
        tracker,
    )
    decomposition = authenticate_r17_decomposition(
        retained["decomposition_evidence"], outputs_sha256
    )
    _, _, scientific_payload = authenticate_r17_root_binding(
        retained["scientific_evidence"], root, "R17 scientific evidence", tracker
    )
    scientific = require_exact_keys(
        parse_json(scientific_payload, "R17 scientific evidence"),
        R17_SCIENTIFIC_EVIDENCE_KEYS,
        "R17 scientific evidence",
    )
    if (
        scientific.get("schema_version") != 6
        or scientific.get("record_type")
        != "cgl_lf_stage_i_qualification_scientific_evidence"
        or scientific.get("case_id") != R17_CASE_ID
        or scientific.get("nodes") != 8
        or scientific.get("target_time") != 0.25
        or scientific.get("scientific_policy") != "active_hardwall"
        or scientific.get("executable")
        != {
            "revision": profile["executable_revision"],
            "sha256": profile["executable_sha256"],
        }
        or scientific.get("execution_intent_sha256") != frozen["execution_intent_sha256"]
        or scientific.get("execution_contract_sha256")
        != frozen["execution_contract_sha256"]
        or scientific.get("terminal_rank_local_outputs") != outputs
        or scientific.get("terminal_rank_local_output_inventory_sha256") != outputs_sha256
        or scientific.get("terminal_rank_local_restarts") != restarts
        or scientific.get("terminal_rank_local_restart_inventory_sha256") != restarts_sha256
        or scientific.get("r17_decomposition") != decomposition
        or scientific.get("accepted_for_operational_qualification") is not True
        or scientific.get("accepted_for_profile_selection") is not False
    ):
        raise ValueError("R17 scientific evidence differs from canonical contract")
    physics_measurements = authenticate_r17_physics_contract(
        scientific, root, outputs, restarts, tracker
    )

    evidence_bindings: dict[str, dict[str, object]] = {}
    for key, record_type, expected_measurements in (
        (
            "restart_load_evidence",
            "stage-i-r17-restart-load-evidence",
            {
                "rank_local_restart_inventory_sha256": restarts_sha256,
                "loaded_rank_count": 64,
                "load_state": "COMPLETED",
                "load_exit_code": "0:0",
            },
        ),
        (
            "physics_validation_evidence",
            "stage-i-r17-physics-validation-evidence",
            {
                "rank_local_output_inventory_sha256": outputs_sha256,
                **physics_measurements,
            },
        ),
    ):
        path, digest, payload = authenticate_r17_root_binding(
            retained[key], root, f"R17 qualification {key}", tracker, expected_mode=0o444
        )
        evidence = require_exact_keys(
            parse_json(payload, f"R17 qualification {key}"),
            {
                "schema_version", "record_type", "execution_epoch", "measured_utc",
                "measured_by", "job_id", "executable_sha256",
                "build_manifest_inventory_sha256", "passed", "measurements",
            },
            f"R17 qualification {key}",
        )
        if (
            evidence["schema_version"] != 1
            or evidence["record_type"] != record_type
            or evidence["execution_epoch"] != EXECUTION_EPOCH
            or evidence["measured_utc"] != retained["measured_utc"]
            or evidence["measured_by"] != measured_by
            or evidence["job_id"] != job_id
            or evidence["executable_sha256"] != retained["executable_sha256"]
            or evidence["build_manifest_inventory_sha256"] != build_inventory_sha256
            or evidence["passed"] is not True
            or evidence["measurements"] != expected_measurements
        ):
            raise ValueError(f"R17 qualification {key} differs")
        evidence_bindings[key] = {
            "path": path.relative_to(root).as_posix(),
            "sha256": digest,
            "measured_utc": evidence["measured_utc"],
            "measured_by": evidence["measured_by"],
            "measurements": evidence["measurements"],
        }

    scheduler_path, _, scheduler_payload = authenticate_r17_root_binding(
        retained["scheduler_evidence"],
        root,
        "R17 qualification scheduler evidence",
        tracker,
        expected_path=root / "accounting" / f"{job_id}.r17_qualification.sacct.txt",
        expected_mode=0o444,
    )
    try:
        lines = scheduler_payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("R17 qualification scheduler evidence is not UTF-8") from error
    if len(lines) != 1:
        raise ValueError("R17 qualification scheduler evidence must contain exactly one row")
    fields = lines[0].split("|")
    if len(fields) != 8 or any(not field or field.strip() != field for field in fields):
        raise ValueError("R17 qualification scheduler evidence row is malformed")
    submitted = parse_aware_scheduler_timestamp(
        fields[6], "R17 qualification scheduler submit time"
    )
    scheduler_completed = parse_aware_scheduler_timestamp(
        fields[7], "R17 qualification scheduler completion time"
    )
    if (
        fields[0] != job_id
        or fields[2:5] != ["COMPLETED", "0:0", "8"]
        or require_integer(int(fields[5]), "R17 qualification elapsed seconds", minimum=1) <= 0
        or scheduler_completed != completed
        or scheduler_completed <= submitted
    ):
        raise ValueError("R17 qualification scheduler evidence differs")

    _, _, account_scheduler_payload = authenticate_r17_root_binding(
        retained["account_scheduler_evidence"],
        root,
        "R17 qualification account scheduler evidence",
        tracker,
        expected_path=root / "accounting" / f"{job_id}.r17_qualification.account.sacct.txt",
        expected_mode=0o444,
    )
    account_exclusivity_path, _, account_exclusivity_payload = authenticate_r17_root_binding(
        retained["account_exclusivity_evidence"],
        root,
        "R17 qualification account exclusivity evidence",
        tracker,
        expected_path=root / "accounting" / (
            f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_account_exclusivity_evidence.json"
        ),
        expected_mode=0o444,
    )
    account_exclusivity = authenticate_r17_account_exclusivity(
        parse_json(account_exclusivity_payload, "R17 qualification account exclusivity evidence"),
        account_scheduler_payload,
        job_id,
        fields,
        completed.astimezone(timezone.utc),
        request_timestamp,
    )
    if account_exclusivity["measured_utc"] != retained["measured_utc"]:
        raise ValueError("R17 account exclusivity measurement differs")

    review_contract = require_exact_keys(
        retained["independent_review_contract"],
        {
            "required", "path", "mode", "schema_version", "record_type",
            "execution_epoch", "decision", "candidate_path",
            "candidate_sha256_required", "reviewer_must_differ_from",
            "reviewed_after_utc",
        },
        "R17 independent-review contract",
    )
    expected_review_path = qualification_path.with_name(
        f"{qualification_path.name}.independent_review.json"
    )
    expected_review_contract = {
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
        "reviewed_after_utc": retained["measured_utc"],
    }
    if (
        review_contract != expected_review_contract
        or retained["authority"]
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("R17 independent-review or authority contract differs")
    review_path, _, review_payload = authenticate_r17_root_binding(
        review_binding,
        root,
        "R17 qualification independent review",
        tracker,
        expected_path=expected_review_path,
        expected_mode=0o444,
    )
    review = require_exact_keys(
        parse_json(review_payload, "R17 qualification independent review"),
        {
            "schema_version", "record_type", "execution_epoch", "reviewed_utc",
            "decision", "reviewer", "candidate",
        },
        "R17 qualification independent review",
    )
    reviewer = require_nonempty_string(review["reviewer"], "R17 qualification reviewer")
    reviewed = parse_utc_timestamp(
        review["reviewed_utc"], "R17 qualification review timestamp"
    )
    if (
        review["schema_version"] != review_contract["schema_version"]
        or review["record_type"] != review_contract["record_type"]
        or review["execution_epoch"] != review_contract["execution_epoch"]
        or review["decision"] != review_contract["decision"]
        or review["candidate"] != {"path": str(qualification_path), "sha256": qualification_sha256}
        or reviewer in review_contract["reviewer_must_differ_from"]
        or reviewed < measured
        or reviewed > request_timestamp
    ):
        raise ValueError("R17 qualification independent review differs")

    if root == DEFAULT_ROOT:
        environment = hardened_child_environment("SLURM_")
        try:
            live = subprocess.run(
                [
                    str(SACCT), "-X", "-j", job_id,
                    "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                    "ElapsedRaw,Submit,End", "-n", "-P",
                ],
                check=True,
                capture_output=True,
                text=True,
                env=environment,
            ).stdout
        except (OSError, subprocess.CalledProcessError) as error:
            raise ValueError("live R17 qualification scheduler evidence is unavailable") from error
        live_rows = [
            line.rstrip("|").split("|")
            for line in live.splitlines()
            if line.rstrip("|").split("|")[0] == job_id
        ]
        if live_rows != [fields]:
            raise ValueError("live R17 qualification scheduler evidence differs")
    result = dict(retained)
    result["authenticated_rank_local_outputs"] = outputs
    result["authenticated_rank_local_restarts"] = restarts
    result["authenticated_decomposition_evidence"] = decomposition
    result["authenticated_account_exclusivity_evidence"] = account_exclusivity
    result["authenticated_validation_evidence"] = evidence_bindings
    result["authenticated_frozen_science_build_contract"] = frozen
    result["reviewed_by"] = reviewer
    return result


def parse_r17_readiness(
    value: object,
    root: Path,
    request_timestamp: datetime,
    expires_timestamp: datetime,
    lineage_sha256: str,
    storage_sha256: str,
    projection_sha256: str,
    profile: dict[str, object],
    matrix_sha256: str,
    tracker: InputTracker,
) -> dict[str, object]:
    """Validate exact independently reviewed R17 launch readiness."""

    retained = require_exact_keys(
        value,
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "root",
            "generated_utc",
            "expires_utc",
            "reviewed_by",
            "predecessor_lineages_sha256",
            "storage_evidence_sha256",
            "computed_projection_sha256",
            "executable_sha256",
            "build_manifest_sha256",
            "required_retained_bytes",
            "nodes",
            "ranks",
            "storage_ready",
            "node_hour_ready",
            "rank_64_ready",
            "operational_qualification",
            "operational_qualification_review",
        },
        "R17 readiness evidence",
    )
    if (
        retained["schema_version"] != 1
        or retained["record_type"] != "stage-i-r17-readiness"
        or retained["execution_epoch"] != EXECUTION_EPOCH
        or retained["root"] != str(root)
    ):
        raise ValueError("R17 readiness campaign binding differs")
    generated = parse_utc_timestamp(retained["generated_utc"], "R17 readiness generation")
    expires = parse_utc_timestamp(retained["expires_utc"], "R17 readiness expiry")
    if (
        generated > request_timestamp
        or expires < expires_timestamp
        or expires <= generated
        or expires - generated > AUTHORIZATION_MAX_LIFETIME
    ):
        raise ValueError("R17 readiness lifetime does not cover the recost request")
    readiness_reviewer = require_nonempty_string(
        retained["reviewed_by"], "R17 readiness reviewer"
    )
    for key in (
        "predecessor_lineages_sha256",
        "storage_evidence_sha256",
        "computed_projection_sha256",
        "executable_sha256",
        "build_manifest_sha256",
    ):
        require_sha256(retained[key], f"R17 readiness {key}")
    require_integer(
        retained["required_retained_bytes"], "R17 readiness required retained bytes"
    )
    require_integer(retained["nodes"], "R17 readiness nodes", minimum=1)
    require_integer(retained["ranks"], "R17 readiness ranks", minimum=1)
    if any(
        retained[key] is not True
        for key in ("storage_ready", "node_hour_ready", "rank_64_ready")
    ):
        raise ValueError("R17 readiness flags must be true booleans")
    expected = {
        "predecessor_lineages_sha256": lineage_sha256,
        "storage_evidence_sha256": storage_sha256,
        "computed_projection_sha256": projection_sha256,
        "executable_sha256": profile["executable_sha256"],
        "build_manifest_sha256": profile["build_manifest_sha256"],
        "required_retained_bytes": R17_MINIMUM_RETAINED_BYTES,
        "nodes": 8,
        "ranks": 64,
        "storage_ready": True,
        "node_hour_ready": True,
        "rank_64_ready": True,
    }
    if any(retained.get(key) != expected_value for key, expected_value in expected.items()):
        raise ValueError("R17 readiness reviewed bindings differ")
    qualification_relative, qualification_sha256 = input_binding(
        retained["operational_qualification"], "R17 operational qualification binding"
    )
    qualification_path = root_path(
        root,
        qualification_relative.as_posix(),
        "R17 operational qualification evidence",
    )
    if qualification_path.parent != root / "accounting":
        raise ValueError("R17 operational qualification must be retained under accounting")
    qualification = parse_r17_operational_qualification(
        parse_json(
            tracker.read(
                qualification_path,
                qualification_sha256,
                "R17 operational qualification evidence",
                expected_mode=0o444,
            ),
            "R17 operational qualification evidence",
        ),
        root,
        request_timestamp,
        profile,
        matrix_sha256,
        qualification_path,
        qualification_sha256,
        retained["operational_qualification_review"],
        tracker,
    )
    if qualification["reviewed_by"] == readiness_reviewer:
        raise ValueError("R17 readiness and operational qualification require distinct reviewers")
    result = dict(retained)
    result["operational_qualification_evidence"] = qualification
    return result


def parse_r17_readiness_publication_chain(
    root: Path,
    readiness_path: Path,
    readiness_sha256: str,
    readiness: dict[str, object],
    review_binding: object,
    audit_binding: object,
    request_timestamp: datetime,
    tracker: InputTracker,
) -> dict[str, object]:
    """Authenticate exact independent review and publication of measured R17 readiness."""

    review_relative, review_sha256 = input_binding(
        review_binding, "R17 readiness independent-review binding"
    )
    review_path = root_path(root, review_relative.as_posix(), "R17 readiness independent review")
    if review_path != readiness_path.with_name(f"{readiness_path.name}.independent_review.json"):
        raise ValueError("R17 readiness independent-review path differs")
    review = require_exact_keys(
        parse_json(
            tracker.read(
                review_path,
                review_sha256,
                "R17 readiness independent review",
                expected_mode=0o444,
            ),
            "R17 readiness independent review",
        ),
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "reviewed_utc",
            "decision",
            "reviewer",
            "candidate",
        },
        "R17 readiness independent review",
    )
    if (
        review["schema_version"] != 1
        or review["record_type"] != "stage-i-r17-readiness-independent-review"
        or review["execution_epoch"] != EXECUTION_EPOCH
        or review["decision"] != "approved-for-publication"
        or review["candidate"]
        != {"path": str(readiness_path), "sha256": readiness_sha256}
    ):
        raise ValueError("R17 readiness independent review differs")
    reviewer = require_nonempty_string(review["reviewer"], "R17 readiness reviewer")
    reviewed = parse_utc_timestamp(review["reviewed_utc"], "R17 readiness review timestamp")
    if reviewer != readiness["reviewed_by"]:
        raise ValueError("R17 readiness reviewer differs from independent review")

    audit_relative, audit_sha256 = input_binding(
        audit_binding, "R17 readiness publication-audit binding"
    )
    audit_path = root_path(root, audit_relative.as_posix(), "R17 readiness publication audit")
    if audit_path != readiness_path.with_name(f"{readiness_path.name}.publication_audit.json"):
        raise ValueError("R17 readiness publication-audit path differs")
    audit = require_exact_keys(
        parse_json(
            tracker.read(
                audit_path,
                audit_sha256,
                "R17 readiness publication audit",
                expected_mode=0o444,
            ),
            "R17 readiness publication audit",
        ),
        {
            "schema_version",
            "record_type",
            "execution_epoch",
            "published_utc",
            "artifact",
            "independent_review",
            "authority",
        },
        "R17 readiness publication audit",
    )
    if (
        audit["schema_version"] != 1
        or audit["record_type"] != "stage-i-r17-readiness-publication-audit"
        or audit["execution_epoch"] != EXECUTION_EPOCH
        or audit["authority"]
        != {
            "r17_launch_authorized": False,
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        }
    ):
        raise ValueError("R17 readiness publication audit differs")
    require_declared_binding(
        audit["artifact"], readiness_path, readiness_sha256, "R17 readiness publication artifact"
    )
    require_declared_binding(
        audit["independent_review"],
        review_path,
        review_sha256,
        "R17 readiness publication review",
    )
    published = parse_utc_timestamp(audit["published_utc"], "R17 readiness publication timestamp")
    if published < reviewed or published > request_timestamp:
        raise ValueError("R17 readiness publication chronology differs")
    return {
        "independent_review_sha256": review_sha256,
        "publication_audit_sha256": audit_sha256,
        "reviewed_by": reviewer,
        "published_utc": published.isoformat(),
    }


def require_authoritative_r12_fresh_rerun(
    case_id: str,
    terminal_identity: tuple[str | None, str, str, str] | None,
    segment: str,
    segment_index: int,
    segment_start: float,
    target: float,
    nodes: int,
    ranks_per_node: int,
    walltime: object,
    athena_walltime: object,
    parent_fields: tuple[object, ...],
    expected_segment_index: int,
) -> bool:
    """Require the exact fresh R12 rerun when historical s00 is terminal inventory."""

    if (
        case_id != "R12"
        or terminal_identity != R12_HISTORICAL_INVENTORY_IDENTITY
    ):
        return False
    if (
        segment != R12_FRESH_RERUN_SEGMENT
        or segment_index != expected_segment_index
        or abs(segment_start) > 1.0e-12
        or abs(target - R12_FRESH_RERUN_TARGET) > 1.0e-12
        or abs(target - PLAN_INITIAL_TARGETS["R12"]) > 1.0e-12
        or nodes != R12_FRESH_RERUN_NODES
        or ranks_per_node != R12_FRESH_RERUN_RANKS_PER_NODE
        or nodes * ranks_per_node != R12_FRESH_RERUN_TOTAL_RANKS
        or walltime != R12_FRESH_RERUN_WALLTIME
        or athena_walltime != R12_FRESH_RERUN_ATHENA_WALLTIME
        or any(item is not None for item in parent_fields)
    ):
        raise ValueError(
            "R12 historical clean partial is inventory-only; the authoritative "
            "next profile is fresh R12/s01_rankio_t0_t0p12 on 4 nodes / 32 ranks "
            "with Slurm 02:00:00 and Athena 01:50:00 from t=0 to t=0.12 with "
            "null parent and restart fields"
        )
    return True


def validate_profiles(
    authorization: dict[str, object],
    root: Path,
    repository: Path,
    manifests: list[dict[str, object]],
    matrix: dict[str, dict[str, object]],
    node_profiles: dict[str, set[int]],
    source_bundle_path: Path,
    source_bundle_sha256: str,
    source_bundle_revisions: list[str],
    qualification_approval: dict[str, object],
    qualification_approval_sha256: str,
    lane_limit: int,
    tracker: InputTracker,
) -> tuple[
    str,
    int,
    list[dict[str, object]],
    Decimal,
    int,
    list[dict[str, object]],
    tuple[DirectoryMeasurement, ...],
    dict[str, list[dict[str, object]]],
]:
    """Validate exact sole or bounded-wave next-packet authorization."""

    mode = authorization["mode"]
    if mode not in {"sole-next-profile", "bounded-wave"}:
        raise ValueError("authorization mode must be sole-next-profile or bounded-wave")
    max_wave_nodes = require_integer(
        authorization["max_wave_nodes"], "authorization max wave nodes", minimum=1
    )
    if max_wave_nodes > MAX_WAVE_NODES:
        raise ValueError(f"authorization exceeds the {MAX_WAVE_NODES}-node wave ceiling")
    profiles = authorization["profiles"]
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("authorization profiles must be a nonempty list")
    if mode == "sole-next-profile" and len(profiles) != 1:
        raise ValueError("sole-next-profile mode requires exactly one profile")
    if mode == "bounded-wave" and len(profiles) < 2:
        raise ValueError("bounded-wave mode requires at least two profiles")
    if len(profiles) > lane_limit:
        raise ValueError("authorization profiles exceed the promoted lane limit")

    lineages = authenticated_case_lineages(root, manifests)
    started_r17 = any(manifest_identity(item)[1] == R17_CASE_ID for item in manifests)

    exact_profile_keys = {
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
    retained: list[dict[str, object]] = []
    cases: set[str] = set()
    segments: set[tuple[str, str]] = set()
    total_nodes = 0
    total_reserved = Decimal("0")
    occupied_indexes: dict[str, set[int]] = {}
    for manifest in manifests:
        _, occupied_case, occupied_segment, _ = manifest_identity(manifest)
        occupied_index, _, _ = parse_segment(
            occupied_segment, f"{occupied_case} retained segment"
        )
        occupied_indexes.setdefault(occupied_case, set()).add(occupied_index)
    for index, item in enumerate(profiles):
        profile = require_exact_keys(item, exact_profile_keys, f"next profile {index}")
        case_id = require_safe_id(profile["case_id"], f"next profile {index} case")
        segment = require_safe_id(profile["segment"], f"next profile {index} segment")
        if case_id not in matrix or case_id not in node_profiles:
            raise ValueError(f"next profile {index} has unknown or unpromoted case {case_id}")
        if case_id in cases:
            raise ValueError(f"next profiles duplicate active case {case_id}")
        cases.add(case_id)
        if (case_id, segment) in segments:
            raise ValueError("next profiles duplicate a segment")
        segments.add((case_id, segment))
        if any(
            manifest_identity(manifest)[1:3] == (case_id, segment)
            for manifest in manifests
        ):
            raise ValueError(f"next profile collides with retained segment {case_id}/{segment}")
        nodes = require_integer(profile["nodes"], f"next profile {index} nodes", minimum=1)
        if nodes not in node_profiles[case_id]:
            raise ValueError(f"next profile {case_id} nodes are not promoted")
        if require_integer(
            profile["ranks_per_node"], f"next profile {case_id} ranks per node", minimum=1
        ) != EXPECTED_RANKS_PER_NODE:
            raise ValueError(f"next profile {case_id} ranks per node differ")
        if require_integer(
            profile["cpus_per_task"], f"next profile {case_id} CPUs per task", minimum=1
        ) != EXPECTED_CPUS_PER_TASK:
            raise ValueError(f"next profile {case_id} CPUs per task differ")
        slurm_seconds = walltime_seconds(profile["walltime"], f"next profile {index} walltime")
        if require_integer(
            profile["controller_walltime_max_seconds"],
            f"next profile {case_id} controller walltime maximum",
            minimum=1,
        ) != MAX_WALLTIME_SECONDS:
            raise ValueError(f"next profile {case_id} controller walltime maximum differs")
        athena_seconds = walltime_seconds(
            profile["athena_walltime"], f"next profile {index} Athena walltime"
        )
        if slurm_seconds - athena_seconds < MINIMUM_SHUTDOWN_MARGIN_SECONDS:
            raise ValueError(
                f"next profile {case_id} must retain the ten-minute shutdown margin"
            )
        target = require_finite_float(
            profile["time_tlim_target"], f"next profile {index} target"
        )
        if target <= 0 or target > REQUIRED_CASE_FINAL_TIME:
            raise ValueError(f"next profile {case_id} target is outside Stage I")
        segment_index, segment_start, segment_target = parse_segment(
            segment, f"next profile {case_id} segment"
        )
        if abs(segment_target - target) > 1.0e-10:
            raise ValueError(f"next profile {case_id} segment target differs")
        if profile["source_bundle"] != str(source_bundle_path):
            raise ValueError(f"next profile {case_id} source bundle path differs")
        if profile["source_bundle_sha256"] != source_bundle_sha256:
            raise ValueError(f"next profile {case_id} source bundle SHA-256 differs")
        executable = absolute_root_member(
            root, profile["executable"], f"next profile {case_id} executable"
        )
        executable_sha256 = require_sha256(
            profile["executable_sha256"], f"next profile {case_id} executable SHA-256"
        )
        tracker.authenticate(
            executable,
            executable_sha256,
            f"next profile {case_id} executable",
            expected_mode=0o755,
        )
        executable_revision = require_revision(
            profile["executable_revision"], f"next profile {case_id} executable revision"
        )
        if executable_revision not in source_bundle_revisions:
            raise ValueError(f"next profile {case_id} executable revision is outside the bundle")
        build_manifest = absolute_root_member(
            root, profile["build_manifest"], f"next profile {case_id} build manifest"
        )
        build_manifest_sha256 = directory_inventory_sha256(
            build_manifest, tracker, f"next profile {case_id} build manifest"
        )
        if profile["build_manifest_sha256"] != build_manifest_sha256:
            raise ValueError(f"next profile {case_id} build manifest digest differs")
        parse_qualification_approval(
            qualification_approval,
            root,
            executable,
            executable_revision,
            executable_sha256,
            build_manifest,
        )
        athena_digest_path = build_manifest / "athena.sha256"
        try:
            athena_payload, _ = tracker.discover(
                athena_digest_path,
                f"next profile {case_id} build executable digest",
                expected_mode=0o644,
            )
            athena_digest_line = athena_payload.decode().splitlines()
        except UnicodeDecodeError as error:
            raise ValueError(
                f"next profile {case_id} build executable digest is invalid"
            ) from error
        if athena_digest_line != [f"{executable_sha256}  {executable}"]:
            raise ValueError(f"next profile {case_id} build manifest does not bind executable")
        try:
            environment_payload, _ = tracker.discover(
                build_manifest / "environment.txt",
                f"next profile {case_id} build environment",
                expected_mode=0o644,
            )
            environment_lines = environment_payload.decode().splitlines()
        except UnicodeDecodeError as error:
            raise ValueError(f"next profile {case_id} build environment is invalid") from error
        if [
            line for line in environment_lines if line.startswith("git_revision=")
        ] != [f"git_revision={executable_revision}"]:
            raise ValueError(f"next profile {case_id} build environment revision differs")
        input_relative = require_relative_path(
            profile["input_file"], f"next profile {case_id} input file"
        )
        if input_relative.as_posix() != matrix[case_id]["input"]:
            raise ValueError(f"next profile {case_id} input differs from the matrix")
        input_path = repository_path(
            repository, input_relative.as_posix(), f"next profile {case_id} input"
        )
        input_revision = require_revision(
            profile["input_revision"], f"next profile {case_id} input revision"
        )
        input_sha256 = require_sha256(
            profile["input_sha256"], f"next profile {case_id} input SHA-256"
        )
        require_committed_file(
            repository,
            input_path,
            input_sha256,
            f"next profile {case_id} input",
            revision=input_revision,
        )
        tracker.authenticate(input_path, input_sha256, f"next profile {case_id} input")
        if input_revision not in source_bundle_revisions:
            raise ValueError(f"next profile {case_id} input revision is outside the bundle")
        if profile["output_layout"] != "rank-local":
            raise ValueError(f"next profile {case_id} output layout must be rank-local")
        if profile["acceptance_policy"] != ACCEPTANCE_POLICIES[case_id]:
            raise ValueError(f"next profile {case_id} acceptance policy differs")
        criterion = require_nonempty_string(
            profile["acceptance_criterion"], f"next profile {case_id} acceptance criterion"
        )
        if criterion != acceptance_criterion(case_id, target):
            raise ValueError(f"next profile {case_id} acceptance criterion differs")
        parents = lineages.get(case_id, [])
        expected_segment_index = max(occupied_indexes.get(case_id, {-1})) + 1
        basis = profile["recommendation_basis"]
        parent_fields = (
            profile["parent_job_id"],
            profile["parent_result"],
            profile["parent_segment"],
            profile["restart_file"],
            profile["restart_file_sha256"],
            profile["restart_time"],
        )
        r12_historical_inventory = require_authoritative_r12_fresh_rerun(
            case_id,
            manifest_identity(parents[-1]) if parents else None,
            segment,
            segment_index,
            segment_start,
            target,
            nodes,
            int(profile["ranks_per_node"]),
            profile["walltime"],
            profile["athena_walltime"],
            parent_fields,
            expected_segment_index,
        )
        if not parents or r12_historical_inventory:
            if any(item is not None for item in parent_fields):
                raise ValueError(f"fresh next profile {case_id} has ambiguous parent fields")
            if segment_index != expected_segment_index or abs(segment_start) > 1.0e-12:
                raise ValueError(
                    f"fresh/retry next profile {case_id} segment identity differs"
                )
            if abs(target - PLAN_INITIAL_TARGETS[case_id]) > 1.0e-12:
                raise ValueError(f"fresh next profile {case_id} target differs from the plan")
            if require_exact_keys(
                basis,
                {
                    "kind",
                    "qualification_approval_sha256",
                    "profile_status",
                },
                f"next profile {case_id} recommendation basis",
            ) != {
                "kind": "qualified-initial-calibration",
                "qualification_approval_sha256": qualification_approval_sha256,
                "profile_status": "promoted-by-f113-for-measured-calibration",
            }:
                raise ValueError(
                    f"fresh/retry next profile {case_id} lacks exact qualification basis"
                )
        else:
            if any(item is None for item in parent_fields):
                raise ValueError(f"continuation next profile {case_id} lacks exact parent fields")
            parent_job_id = require_nonempty_string(
                profile["parent_job_id"], f"next profile {case_id} parent job"
            )
            parent_result = require_nonempty_string(
                profile["parent_result"], f"next profile {case_id} parent result"
            )
            parent_segment = require_safe_id(
                profile["parent_segment"], f"next profile {case_id} parent segment"
            )
            parent = parents[-1]
            if manifest_identity(parent) != (
                parent_job_id, case_id, parent_segment, parent_result
            ):
                raise ValueError(f"next profile {case_id} parent is not the lineage terminal")
            maximum = final_time(parent)
            inspection = parent["scientific_inspection"]
            paths = parent.get("paths")
            assert isinstance(inspection, dict)
            if inspection.get("clean_for_continuation") is not True:
                raise ValueError(f"next profile {case_id} parent is not clean for continuation")
            restart_time = require_finite_float(
                profile["restart_time"], f"next profile {case_id} restart time"
            )
            if restart_time != maximum or target <= restart_time:
                raise ValueError(f"next profile {case_id} restart/target times differ")
            if (
                segment_index != expected_segment_index
                or abs(segment_start - restart_time) > 1.0e-6
            ):
                raise ValueError(f"next profile {case_id} segment lineage differs")
            if target - restart_time > PLAN_MAX_INCREMENTS[case_id] + 1.0e-12:
                raise ValueError(f"next profile {case_id} increment exceeds the plan")
            require_continuation_target_alignment(case_id, target)
            if not isinstance(paths, dict):
                raise ValueError(f"next profile {case_id} parent lacks output path")
            output_dir = absolute_root_member(
                root, paths.get("output_dir"), f"next profile {case_id} parent output"
            )
            restart = absolute_root_member(
                root, profile["restart_file"], f"next profile {case_id} restart"
            )
            try:
                restart.relative_to(output_dir)
            except ValueError as error:
                raise ValueError(
                    f"next profile {case_id} restart is outside parent output"
                ) from error
            tracker.authenticate(
                restart,
                profile["restart_file_sha256"],
                f"next profile {case_id} restart",
                expected_mode=0o644,
            )
            terminal = authenticate_terminal_restart(
                root, parent, tracker, f"next profile {case_id} parent"
            )
            if (
                profile["restart_file"] != terminal.get("path")
                or profile["restart_file_sha256"] != terminal.get("sha256")
            ):
                raise ValueError(
                    f"next profile {case_id} restart differs from terminal restart group"
                )
            parent_allocation = parent.get("allocation")
            parent_accounting = parent.get("accounting")
            if not isinstance(parent_allocation, dict) or not isinstance(parent_accounting, dict):
                raise ValueError(f"next profile {case_id} measured parent accounting is absent")
            if parent_allocation.get("nodes") != nodes:
                raise ValueError(
                    f"next profile {case_id} changes nodes without published profile qualification"
                )
            _, parent_start, _ = parse_segment(
                parent_segment, f"next profile {case_id} measured parent segment"
            )
            observed_interval = Decimal(str(maximum)) - Decimal(str(parent_start))
            observed_elapsed = require_integer(
                int(parent_accounting.get("elapsed_seconds")),
                f"next profile {case_id} measured parent elapsed seconds",
                minimum=1,
            )
            recommended_interval = Decimal(str(target)) - Decimal(str(restart_time))
            maximum_interval = (
                observed_interval
                * Decimal(athena_seconds)
                / Decimal(observed_elapsed)
                * CONTINUATION_RUNTIME_SAFETY_FACTOR
            )
            if recommended_interval > maximum_interval + Decimal("0.000000000001"):
                raise ValueError(
                    f"next profile {case_id} interval exceeds measured runtime recommendation"
                )
            expected_basis = {
                "kind": "measured-production-continuation",
                "job_id": parent_job_id,
                "case_id": case_id,
                "nodes": nodes,
                "observed_result": parent_result,
                "observed_simulation_interval": decimal_string(observed_interval),
                "observed_elapsed_seconds": observed_elapsed,
                "athena_walltime_seconds": athena_seconds,
                "runtime_safety_factor": decimal_string(CONTINUATION_RUNTIME_SAFETY_FACTOR),
                "maximum_recommended_interval": decimal_string(maximum_interval),
                "recommended_interval": decimal_string(recommended_interval),
            }
            if require_exact_keys(
                basis, set(expected_basis), f"next profile {case_id} recommendation basis"
            ) != expected_basis:
                raise ValueError(
                    f"next profile {case_id} recommendation basis differs from measured evidence"
                )

        total_nodes += nodes
        total_reserved += Decimal(nodes * slurm_seconds) / Decimal(3600)
        retained.append(dict(profile))

    if total_nodes != max_wave_nodes:
        raise ValueError("authorized profile nodes differ from exact max_wave_nodes")
    if mode == "sole-next-profile" and max_wave_nodes != retained[0]["nodes"]:
        raise ValueError("sole-next-profile max_wave_nodes differs from its profile")
    if R17_CASE_ID in cases:
        if len(retained) != 1 or retained[0]["nodes"] != 8:
            raise ValueError("R17 authorization must be exclusive on eight nodes")
        incomplete = []
        for predecessor in R17_PREDECESSOR_CASE_IDS:
            lineage = lineages.get(predecessor, [])
            if (
                not lineage
                or manifest_identity(lineage[-1])[3] != "accepted"
                or abs(final_time(lineage[-1]) - REQUIRED_CASE_FINAL_TIME) > 1.0e-10
            ):
                incomplete.append(predecessor)
        if incomplete:
            raise ValueError(
                "R17 must remain last; accepted t=10 predecessors are incomplete: "
                + ", ".join(incomplete)
            )
    elif started_r17:
        raise ValueError("R17 has started; lower-resolution authorization is forbidden")
    total_storage, storage_projections, storage_measurements = observed_storage_projections(
        root, manifests, matrix, retained
    )
    return (
        str(mode),
        max_wave_nodes,
        sorted(retained, key=lambda item: (str(item["case_id"]), str(item["segment"]))),
        total_reserved,
        total_storage,
        storage_projections,
        storage_measurements,
        lineages,
    )


def git_descriptor_run(
    repository: Path,
    descriptor: int,
    arguments: list[str],
    *,
    capture_output: bool = False,
) -> subprocess.CompletedProcess:
    """Run Git against an inherited immutable descriptor path."""

    return git_run(
        repository,
        arguments,
        capture_output=capture_output,
        pass_fds=(descriptor,),
    )


def require_valid_git_bundle(
    repository: Path,
    path: Path,
    expected_sha256: str,
    revisions: list[str],
) -> None:
    """Require one descriptor-bound complete bundle covering requested revisions."""

    expected = require_sha256(expected_sha256, "source bundle SHA-256")
    with absolute_descriptor(path, "source bundle", flags=os.O_RDONLY) as descriptor:
        profile = os.fstat(descriptor)
        require_regular_profile(profile, "source bundle", expected_mode=0o644)
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("source bundle checksum differs before verification")
        descriptor_path = f"/proc/self/fd/{descriptor}"
        header = bytearray()
        os.lseek(descriptor, 0, os.SEEK_SET)
        while b"\n\n" not in header:
            block = os.read(descriptor, 4096)
            if not block:
                break
            if len(header) + len(block) > 1024 * 1024:
                raise ValueError("source bundle header is invalid or unexpectedly large")
            header.extend(block)
        try:
            lines = bytes(header).split(b"\n\n", 1)[0].decode().splitlines()
        except UnicodeDecodeError as error:
            raise ValueError("source bundle header is not UTF-8") from error
        if not lines or lines[0] not in {"# v2 git bundle", "# v3 git bundle"}:
            raise ValueError("source bundle header version is invalid")
        if any(line.startswith("-") for line in lines[1:]):
            raise ValueError("source bundle must be self-contained without prerequisites")
        advertised = []
        for line in lines[1:]:
            if line.startswith("@"):
                continue
            match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
            if match is None:
                raise ValueError("source bundle advertised reference is malformed")
            advertised.append(match.group(1))
        if not advertised:
            raise ValueError("source bundle advertises no retained revision")
        if git_descriptor_run(
            repository,
            descriptor,
            ["bundle", "verify", descriptor_path],
            capture_output=True,
        ).returncode:
            raise ValueError("source bundle verification failed")
        for revision in revisions:
            if git_run(repository, ["cat-file", "-e", f"{revision}^{{commit}}"]).returncode:
                raise ValueError(f"source bundle revision is unknown locally: {revision}")
        retained_umask = os.umask(0o077)
        try:
            with tempfile.TemporaryDirectory(prefix="cgl-lf-recost-bundle-") as directory:
                isolated = Path(directory) / "repository.git"
                if git_run(repository, ["init", "--bare", str(isolated)]).returncode:
                    raise ValueError("cannot initialize isolated source-bundle verification")
                if git_descriptor_run(
                    isolated,
                    descriptor,
                    ["bundle", "unbundle", descriptor_path],
                    capture_output=True,
                ).returncode:
                    raise ValueError("source bundle cannot be reconstructed in isolation")
                if git_run(
                    isolated,
                    ["fsck", "--full", "--strict", "--no-reflogs", *advertised],
                    capture_output=True,
                ).returncode:
                    raise ValueError("source bundle advertised history is incomplete")
                for revision in revisions:
                    if git_run(
                        isolated, ["cat-file", "-e", f"{revision}^{{commit}}"]
                    ).returncode:
                        raise ValueError(
                            f"source bundle does not contain requested revision: {revision}"
                        )
                    if not any(
                        git_run(
                            isolated, ["merge-base", "--is-ancestor", revision, head]
                        ).returncode
                        == 0
                        for head in advertised
                    ):
                        raise ValueError(
                            f"source bundle does not cover requested revision: {revision}"
                        )
        finally:
            os.umask(retained_umask)
        with absolute_descriptor(path, "source bundle", flags=os.O_RDONLY) as named:
            named_profile = os.fstat(named)
            if (
                (named_profile.st_dev, named_profile.st_ino)
                != (profile.st_dev, profile.st_ino)
                or sha256_descriptor(named) != expected
            ):
                raise ValueError("source bundle pathname changed during verification")


def calculate_budget(
    rows: list[dict[str, str]],
    lineages: dict[str, list[dict[str, object]]],
    matrix: dict[str, dict[str, object]],
    profiles: list[dict[str, object]],
    next_reserved: Decimal,
    envelope: Decimal,
    project: Decimal,
) -> dict[str, object]:
    """Project remaining use from authenticated observed Stage I node-hour rates."""

    actual = sum(
        (require_decimal(row["actual_node_hours"], "ledger actual node-hours") for row in rows),
        Decimal("0"),
    )
    committed = actual + next_reserved
    if committed > envelope:
        raise ValueError("actual use plus authorized wave exceeds the Stage I envelope")
    if committed > project:
        raise ValueError("actual use plus authorized wave exceeds the project ceiling")

    reserved_by_case: dict[str, Decimal] = {}
    for profile in profiles:
        case_id = str(profile["case_id"])
        seconds = walltime_seconds(profile["walltime"], f"{case_id} projected walltime")
        reserved_by_case[case_id] = (
            reserved_by_case.get(case_id, Decimal("0"))
            + Decimal(int(profile["nodes"]) * seconds) / Decimal(3600)
        )
    if sum(reserved_by_case.values(), Decimal("0")) != next_reserved:
        raise ValueError("authorized-wave budget arithmetic differs")

    rows_by_job = {row["job_id"]: row for row in rows}
    measurements: list[tuple[Decimal, dict[str, object]]] = []
    r12_measurements: list[tuple[Decimal, dict[str, object]]] = []
    historical_r12_inventory_count = 0
    for case_id, lineage in sorted(lineages.items()):
        cells = require_integer(matrix[case_id]["_cell_count"], f"{case_id} matrix cells", minimum=1)
        for manifest in lineage:
            identity = manifest_identity(manifest)
            job_id, _, segment, _ = identity
            if job_id == R12_HISTORICAL_JOB_ID:
                if identity != R12_HISTORICAL_INVENTORY_IDENTITY:
                    raise ValueError("historical R12 inventory identity partially collides")
                historical_r12_inventory_count += 1
                if historical_r12_inventory_count > 1:
                    raise ValueError("historical R12 inventory identity is duplicated")
            row = rows_by_job[job_id]
            actual_job = require_decimal(row["actual_node_hours"], f"ledger job {job_id} actual")
            _, start, _ = parse_segment(segment, f"ledger job {job_id} segment")
            interval = Decimal(str(final_time(manifest))) - Decimal(str(start))
            if interval <= 0:
                raise ValueError(f"ledger job {job_id} measured interval is invalid")
            if actual_job <= 0:
                continue
            normalized_rate = actual_job / Decimal(cells) / interval
            measurement = (
                normalized_rate,
                {
                    "job_id": job_id,
                    "case_id": case_id,
                    "segment": segment,
                    "actual_node_hours": decimal_string(actual_job),
                    "observed_cells": cells,
                    "observed_simulation_interval": decimal_string(interval),
                    "normalized_node_hours_per_cell_per_simulation_time": (
                        decimal_string(normalized_rate)
                    ),
                },
            )
            if case_id == "R12":
                r12_measurements.append(measurement)
            else:
                measurements.append(measurement)
    if not measurements:
        raise ValueError("ledger contains no positive non-R12 observed node-hour measurement")
    normalized_rate, basis = max(
        measurements, key=lambda item: (item[0], str(item[1]["job_id"]))
    )
    r12_rate, r12_basis = (
        max(r12_measurements, key=lambda item: (item[0], str(item[1]["job_id"])))
        if r12_measurements
        else (normalized_rate, basis)
    )

    remaining = Decimal("0")
    breakdown: dict[str, dict[str, object]] = {}
    for case_id, case in sorted(matrix.items()):
        estimated = case["_estimated_node_hours"]
        assert isinstance(estimated, Decimal)
        lineage = lineages.get(case_id, [])
        progress = Decimal("0")
        if (
            lineage
            and manifest_identity(lineage[-1]) != R12_HISTORICAL_INVENTORY_IDENTITY
        ):
            progress = Decimal(str(final_time(lineage[-1]))) / Decimal(
                str(REQUIRED_CASE_FINAL_TIME)
            )
        progress = min(Decimal("1"), max(Decimal("0"), progress))
        cells = require_integer(case["_cell_count"], f"{case_id} matrix cells", minimum=1)
        remaining_time = Decimal(str(REQUIRED_CASE_FINAL_TIME)) * (
            Decimal("1") - progress
        )
        projection_rate = r12_rate if case_id == "R12" else normalized_rate
        projection_basis = r12_basis if case_id == "R12" else basis
        observed_projection = projection_rate * Decimal(cells) * remaining_time
        authorized_reserved = reserved_by_case.get(case_id, Decimal("0"))
        projected_remaining = max(observed_projection, authorized_reserved)
        remaining += projected_remaining
        breakdown[case_id] = {
            "matrix_full_case_node_hours_reference_only": decimal_string(estimated),
            "authenticated_progress_fraction": decimal_string(progress),
            "remaining_simulation_time": decimal_string(remaining_time),
            "projected_cells": str(cells),
            "projection_measurement_basis": projection_basis,
            "observed_rate_projected_remaining_node_hours": decimal_string(
                observed_projection
            ),
            "authorized_profile_reserved_node_hours": decimal_string(authorized_reserved),
            "projected_remaining_node_hours": decimal_string(projected_remaining),
        }
    projected = actual + remaining
    if projected < committed:
        raise ValueError("computed Stage I projection is below actual plus authorized wave")
    r12_lineage = lineages.get("R12", [])
    historical_r12_requires_fresh_calibration = (
        len(r12_lineage) == 1
        and manifest_identity(r12_lineage[0]) == R12_HISTORICAL_INVENTORY_IDENTITY
    )
    fresh_r12_profiles = []
    for profile in profiles:
        try:
            target = Decimal(str(profile.get("time_tlim_target")))
        except InvalidOperation:
            target = None
        if (
            profile.get("case_id") == "R12"
            and profile.get("segment") == R12_FRESH_RERUN_SEGMENT
            and profile.get("nodes") == R12_FRESH_RERUN_NODES
            and profile.get("ranks_per_node") == R12_FRESH_RERUN_RANKS_PER_NODE
            and profile.get("walltime") == R12_FRESH_RERUN_WALLTIME
            and profile.get("athena_walltime") == R12_FRESH_RERUN_ATHENA_WALLTIME
            and target == Decimal(str(R12_FRESH_RERUN_TARGET))
        ):
            fresh_r12_profiles.append(profile)
    provisional_r12_calibration = (
        historical_r12_requires_fresh_calibration and len(fresh_r12_profiles) == 1
    )
    if projected > envelope and not provisional_r12_calibration:
        raise ValueError("computed Stage I projection exceeds the promoted envelope")
    if projected > project:
        raise ValueError("computed Stage I projection exceeds the project ceiling")
    return {
        "method": NODE_HOUR_PROJECTION_METHOD,
        "measurement_basis": basis,
        "actual_stage_i_node_hours": decimal_string(actual),
        "authorized_wave_reserved_node_hours": decimal_string(next_reserved),
        "actual_plus_authorized_wave_node_hours": decimal_string(committed),
        "computed_remaining_stage_i_node_hours": decimal_string(remaining),
        "computed_stage_i_total_node_hours": decimal_string(projected),
        "promoted_stage_i_envelope_node_hours": decimal_string(envelope),
        "project_ceiling_node_hours": decimal_string(project),
        "computed_stage_i_margin_node_hours": decimal_string(envelope - projected),
        "case_breakdown": breakdown,
    }


def build_payload(
    args: argparse.Namespace,
    root: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    *,
    stage_i_lock_held: bool = False,
    require_independent_request_review: bool = True,
) -> BuildResult:
    """Authenticate all request inputs and build deterministic output bytes."""

    tracker = InputTracker()
    assert args.request is not None
    assert args.expected_request_sha256 is not None
    assert args.output is not None
    request_path = args.request.absolute()
    try:
        request_path.relative_to(root)
    except ValueError as error:
        raise ValueError("recost request must be retained beneath the campaign root") from error
    if request_path.parent != root / "accounting":
        raise ValueError("recost request must be a direct child of accounting")
    request_payload = tracker.read(
        request_path,
        args.expected_request_sha256,
        "recost request",
        expected_mode=0o644,
    )
    request = parse_request(request_payload)
    request_review_sha256 = None
    request_review = None
    artifact_name = args.output.name.removesuffix(".staged")
    if request["artifact_name"] != artifact_name:
        raise ValueError("recost request artifact name differs from the output namespace")
    request_timestamp = parse_utc_timestamp(
        request["generated_utc"], "recost request generation timestamp"
    )
    if require_independent_request_review:
        request_review_path = request_path.with_name(
            f"{request_path.name}.independent_review.json"
        )
        request_review_payload, request_review_sha256 = tracker.discover(
            request_review_path,
            "recost request independent review",
            expected_mode=0o444,
        )
        request_review = parse_request_independent_review(
            parse_json(request_review_payload, "recost request independent review"),
            request_path,
            args.expected_request_sha256,
            request_timestamp,
            require_nonempty_string(request["requested_by"], "recost request author"),
        )
    elif os.path.lexists(request_path.with_name(f"{request_path.name}.independent_review.json")):
        raise ValueError("draft request must not retain a pre-existing independent review")
    expires_timestamp = parse_utc_timestamp(
        request["expires_utc"], "recost request expiry timestamp"
    )
    now = datetime.now(timezone.utc)
    if request_timestamp > now + timedelta(minutes=5):
        raise ValueError("recost request generation timestamp is in the future")
    if now - request_timestamp > STORAGE_EVIDENCE_MAX_AGE:
        raise ValueError("recost request is stale")
    if expires_timestamp <= request_timestamp:
        raise ValueError("recost request expiry must follow generation")
    if expires_timestamp - request_timestamp > AUTHORIZATION_MAX_LIFETIME:
        raise ValueError("recost request lifetime exceeds 24 hours")
    if now > expires_timestamp:
        raise ValueError("recost request has expired")
    inputs = request["inputs"]
    assert isinstance(inputs, dict)

    generator_revision = require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    tracker.authenticate(
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )

    helper_binding = require_exact_keys(
        inputs["stage_i_helper"], {"revision", "sha256"}, "Stage I helper binding"
    )
    helper_revision = require_revision(
        helper_binding["revision"], "Stage I helper revision"
    )
    helper_sha256 = require_sha256(
        helper_binding["sha256"], "Stage I helper SHA-256"
    )
    helper_path = repository / STAGE_I_RELATIVE
    require_committed_file(
        repository,
        helper_path,
        helper_sha256,
        "Stage I helper",
        revision=helper_revision,
        expected_mode=0o644,
    )
    tracker.authenticate(helper_path, helper_sha256, "Stage I helper", expected_mode=0o644)

    matrix_relative, matrix_revision, matrix_sha256 = repository_binding(
        inputs["matrix"], "Stage I matrix binding"
    )
    if matrix_relative.as_posix() != MATRIX_RELATIVE.as_posix():
        raise ValueError("Stage I matrix binding is not the canonical matrix path")
    matrix_path = repository_path(repository, matrix_relative.as_posix(), "Stage I matrix")
    require_committed_file(
        repository,
        matrix_path,
        matrix_sha256,
        "Stage I matrix",
        revision=matrix_revision,
    )
    matrix_payload = tracker.read(matrix_path, matrix_sha256, "Stage I matrix")
    matrix, matrix_project_ceiling = parse_matrix(
        parse_json(matrix_payload, "Stage I matrix")
    )

    bundle_relative, bundle_sha256, bundle_revisions = source_bundle_binding(
        inputs["source_bundle"], "source bundle binding"
    )
    bundle_path = root_path(root, bundle_relative.as_posix(), "source bundle")
    source_archives = root / "source-archives"
    require_directory(source_archives, "source archive directory")
    if bundle_path.parent != source_archives:
        raise ValueError("source bundle must be a direct child of source-archives")
    tracker.authenticate(bundle_path, bundle_sha256, "source bundle")
    require_valid_git_bundle(repository, bundle_path, bundle_sha256, bundle_revisions)
    for revision, label in (
        (generator_revision, "Stage I recost generator"),
        (helper_revision, "Stage I helper"),
        (matrix_revision, "Stage I matrix"),
    ):
        if revision not in bundle_revisions:
            raise ValueError(f"source bundle verified revisions omit {label}")

    reconcile_relative, reconcile_sha256 = input_binding(
        inputs["reconciliation"], "reconciliation binding"
    )
    reconcile_path = root_path(root, reconcile_relative.as_posix(), "reconciliation evidence")
    if reconcile_path.parent != root / "accounting":
        raise ValueError("reconciliation evidence must be a direct child of accounting")
    reconcile = parse_reconciliation(
        parse_json(
            tracker.read(
                reconcile_path,
                reconcile_sha256,
                "reconciliation evidence",
                expected_mode=0o644,
            ),
            "reconciliation evidence",
        ),
        root,
    )
    live_reconcile = run_authenticated_reconcile(
        helper_path,
        helper_sha256,
        root,
        stage_i_lock_held=stage_i_lock_held,
    )
    if live_reconcile != reconcile:
        raise ValueError("supplied reconciliation differs from authenticated live report")
    counts = reconcile["counts"]
    assert isinstance(counts, dict)

    ledger_relative, ledger_sha256 = input_binding(inputs["ledger"], "ledger binding")
    ledger_path = root_path(root, ledger_relative.as_posix(), "Stage I ledger")
    if ledger_path != (
        root / "accounting" / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv"
    ):
        raise ValueError("ledger binding is not the canonical Stage I ledger")
    rows = parse_ledger(
        tracker.read(ledger_path, ledger_sha256, "Stage I ledger", expected_mode=0o644),
        request_timestamp,
    )

    reservations_relative, reservations_sha256 = input_binding(
        inputs["reservations"], "reservations binding"
    )
    reservations_path = root_path(
        root, reservations_relative.as_posix(), "Stage I reservations"
    )
    if reservations_path != (
        root / "accounting" / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json"
    ):
        raise ValueError("reservations binding is not the canonical Stage I store")
    reservations = parse_reservations(
        parse_json(
            tracker.read(
                reservations_path,
                reservations_sha256,
                "Stage I reservations",
                expected_mode=0o644,
            ),
            "Stage I reservations",
        )
    )

    manifest_bindings, manifests = parse_manifests(root, inputs["manifests"], tracker)
    manifests_by_job = cross_validate_accounting(rows, reservations, manifests, counts)
    require_empty_transaction_stores(root)

    barrier_value = request["barrier"]
    assert isinstance(barrier_value, dict)
    barrier = parse_barrier_segments(barrier_value["recorded_segments"], manifests_by_job)
    ledger_suffix = rows[-len(barrier):]
    if [item["job_id"] for item in barrier] != [
        row["job_id"] for row in ledger_suffix
    ]:
        raise ValueError("barrier jobs are not the exact recorded ledger suffix")
    scheduler = parse_scheduler_evidence(
        root,
        inputs["scheduler_evidence"],
        barrier,
        ledger_suffix,
        manifests_by_job,
        tracker,
        request_timestamp,
    )

    ceiling_relative, ceiling_sha256 = input_binding(
        inputs["ceiling_evidence"], "ceiling evidence binding"
    )
    ceiling_path = root_path(root, ceiling_relative.as_posix(), "ceiling evidence")
    if ceiling_path != root / F113_RELATIVE:
        raise ValueError("ceiling evidence is not the exact promoted F113 artifact")
    if root == DEFAULT_ROOT and ceiling_sha256 != F113_CANONICAL_SHA256:
        raise ValueError("canonical F113 ceiling digest differs from the reviewed digest")
    ceiling_value = parse_json(
        tracker.read(ceiling_path, ceiling_sha256, "ceiling evidence", expected_mode=0o644),
        "ceiling evidence",
    )
    source_authority, adoption, source_authority_published_utc = parse_current_source_authority(
        root,
        repository,
        inputs["source_authority"],
        bundle_path,
        bundle_sha256,
        bundle_revisions,
        ceiling_path,
        ceiling_sha256,
        tracker,
    )
    ceiling_audit_relative, ceiling_audit_sha256 = input_binding(
        inputs["ceiling_publication_audit"], "F113 publication-audit binding"
    )
    ceiling_audit_path = root_path(
        root, ceiling_audit_relative.as_posix(), "F113 publication audit"
    )
    if adoption is None:
        if ceiling_audit_path != root / F113_PUBLICATION_AUDIT_RELATIVE:
            raise ValueError("F113 publication audit is not the exact retained path")
        ceiling_audit = parse_transition_publication_audit(
            parse_json(
                tracker.read(
                    ceiling_audit_path,
                    ceiling_audit_sha256,
                    "F113 publication audit",
                    expected_mode=0o644,
                ),
                "F113 publication audit",
            ),
            root,
            ceiling_path,
            ceiling_sha256,
            request_timestamp,
        )
    else:
        f116_audit = source_authority["publication_audit"]
        assert isinstance(f116_audit, dict)
        if (
            ceiling_audit_path != root / F116_PUBLICATION_AUDIT_RELATIVE
            or ceiling_audit_sha256 != f116_audit["sha256"]
        ):
            raise ValueError(
                "F113 publication-audit supersession is not the exact F116 publication audit"
            )
        ceiling_audit = {
            "schema_version": 1,
            "record_type": "controlled-F116-supersession",
            "execution_epoch": EXECUTION_EPOCH,
            "published_utc": source_authority_published_utc,
            "artifact": adoption["artifact"],
            "supersession": adoption,
        }
    if parse_utc_timestamp(
        ceiling_audit["published_utc"], "F113 publication timestamp"
    ) < parse_utc_timestamp(ceiling_value.get("generated_utc"), "F113 generation timestamp"):
        raise ValueError("F113 publication predates the promoted transition artifact")
    (
        envelope,
        project,
        lane_limit,
        node_profiles,
        f113_helper_revision,
        f113_helper_sha256,
    ) = parse_ceiling_evidence(ceiling_value)
    if f113_helper_revision not in bundle_revisions:
        raise ValueError("source bundle verified revisions omit the F113 controller")
    require_revision_file(
        repository,
        PurePosixPath(STAGE_I_RELATIVE.as_posix()),
        f113_helper_revision,
        f113_helper_sha256,
        "F113 historical Stage I helper",
    )
    if matrix_project_ceiling != project:
        raise ValueError("Stage I matrix project ceiling differs from promoted evidence")

    qualification_relative, qualification_sha256 = input_binding(
        inputs["qualification_approval"], "qualification approval binding"
    )
    qualification_path = root_path(
        root, qualification_relative.as_posix(), "qualification approval"
    )
    if qualification_path != root / QUALIFICATION_APPROVAL_RELATIVE:
        raise ValueError("qualification approval is not the exact retained path")
    qualification_value = parse_json(
        tracker.read(
            qualification_path,
            qualification_sha256,
            "qualification approval",
            expected_mode=0o644,
        ),
        "qualification approval",
    )
    if not isinstance(qualification_value, dict):
        raise ValueError("qualification approval must be an object")

    recommendations = request["recommendations"]
    assert isinstance(recommendations, dict)
    (
        mode,
        max_wave_nodes,
        profiles,
        next_reserved,
        projected_storage,
        storage_projections,
        storage_measurements,
        lineages,
    ) = validate_profiles(
        recommendations,
        root,
        repository,
        manifests,
        matrix,
        node_profiles,
        bundle_path,
        bundle_sha256,
        bundle_revisions,
        qualification_value,
        qualification_sha256,
        lane_limit,
        tracker,
    )
    tracker.retain_directory_measurements(storage_measurements)
    if mode == "sole-next-profile" and len(scheduler) != 1:
        raise ValueError(
            "sole-next-profile mode requires exactly one barrier scheduler row"
        )
    storage_projections_sha256 = sha256_bytes(
        (json.dumps(storage_projections, sort_keys=True) + "\n").encode()
    )

    storage_relative, storage_sha256 = input_binding(
        inputs["storage_evidence"], "storage evidence binding"
    )
    storage_path = root_path(root, storage_relative.as_posix(), "storage evidence")
    if storage_path.parent != root / "accounting":
        raise ValueError("storage evidence must be a direct child of accounting")
    storage = parse_storage_evidence(
        parse_json(
            tracker.read(storage_path, storage_sha256, "storage evidence", expected_mode=0o644),
            "storage evidence",
        ),
        root,
        request_timestamp,
        projected_storage,
        storage_projections_sha256,
    )

    budget = calculate_budget(
        rows, lineages, matrix, profiles, next_reserved, envelope, project
    )
    budget["storage_projection_evidence"] = {
        "method": STORAGE_PROJECTION_METHOD,
        "profile_projections_sha256": storage_projections_sha256,
        "profile_projections": storage_projections,
    }
    projection_sha256 = sha256_bytes(
        (json.dumps(budget, sort_keys=True) + "\n").encode()
    )
    lineage_sha256 = lineage_summary_sha256(lineages)

    predecessor_relative, predecessor_sha256 = input_binding(
        inputs["predecessor_recost"], "predecessor recost binding"
    )
    predecessor_path = root_path(
        root, predecessor_relative.as_posix(), "predecessor recost"
    )
    predecessor_audit_relative, predecessor_audit_sha256 = input_binding(
        inputs["predecessor_recost_publication_audit"],
        "predecessor recost publication-audit binding",
    )
    predecessor_audit_path = root_path(
        root,
        predecessor_audit_relative.as_posix(),
        "predecessor recost publication audit",
    )
    legacy_f114_bootstrap = (
        predecessor_path == root / "accounting" / LEGACY_F114_RECOST_NAME
        and predecessor_sha256 == LEGACY_F114_RECOST_SHA256
        and predecessor_audit_sha256 == LEGACY_F114_PUBLICATION_AUDIT_SHA256
    )
    if legacy_f114_bootstrap:
        if inputs["predecessor_recost_independent_review"] is not None:
            raise ValueError("legacy F114 bootstrap must not bind a fabricated review")
        predecessor_review_path = None
        predecessor_review_sha256 = None
        predecessor_review = None
        predecessor_mode = 0o644
        predecessor_audit_mode = 0o644
    else:
        predecessor_review_relative, predecessor_review_sha256 = input_binding(
            inputs["predecessor_recost_independent_review"],
            "predecessor recost independent-review binding",
        )
        predecessor_review_path = root_path(
            root,
            predecessor_review_relative.as_posix(),
            "predecessor recost independent review",
        )
        predecessor_review = parse_json(
            tracker.read(
                predecessor_review_path,
                predecessor_review_sha256,
                "predecessor recost independent review",
                expected_mode=0o444,
            ),
            "predecessor recost independent review",
        )
        predecessor_mode = 0o444
        predecessor_audit_mode = 0o444
    predecessor = parse_predecessor_recost(
        root,
        predecessor_path,
        predecessor_sha256,
        parse_json(
            tracker.read(
                predecessor_path,
                predecessor_sha256,
                "predecessor recost",
                expected_mode=predecessor_mode,
            ),
            "predecessor recost",
        ),
        predecessor_review_path,
        predecessor_review_sha256,
        predecessor_review,
        predecessor_audit_path,
        predecessor_audit_sha256,
        parse_json(
            tracker.read(
                predecessor_audit_path,
                predecessor_audit_sha256,
                "predecessor recost publication audit",
                expected_mode=predecessor_audit_mode,
            ),
            "predecessor recost publication audit",
        ),
        request_timestamp,
        artifact_name,
        str(request["checkpoint"]),
        tracker,
    )

    r17_profile = next(
        (profile for profile in profiles if profile["case_id"] == R17_CASE_ID), None
    )
    readiness = None
    readiness_sha256 = None
    if r17_profile is None:
        if any(
            inputs[key] is not None
            for key in (
                "r17_readiness_evidence",
                "r17_readiness_independent_review",
                "r17_readiness_publication_audit",
            )
        ):
            raise ValueError("non-R17 recost must not retain R17 readiness publication chain")
    else:
        readiness_relative, readiness_sha256 = input_binding(
            inputs["r17_readiness_evidence"], "R17 readiness binding"
        )
        readiness_path = root_path(root, readiness_relative.as_posix(), "R17 readiness")
        if readiness_path != root / R17_READINESS_RELATIVE:
            raise ValueError("R17 readiness evidence is not the exact retained path")
        readiness = parse_r17_readiness(
            parse_json(
                tracker.read(
                    readiness_path,
                    readiness_sha256,
                    "R17 readiness evidence",
                    expected_mode=0o444,
                ),
                "R17 readiness evidence",
            ),
            root,
            request_timestamp,
            expires_timestamp,
            lineage_sha256,
            storage_sha256,
            projection_sha256,
            r17_profile,
            matrix_sha256,
            tracker,
        )
        readiness["publication_chain"] = parse_r17_readiness_publication_chain(
            root,
            readiness_path,
            readiness_sha256,
            readiness,
            inputs["r17_readiness_independent_review"],
            inputs["r17_readiness_publication_audit"],
            request_timestamp,
            tracker,
        )

    recommendation_output: dict[str, object] = {
        "mode": mode,
        "authorizing": False,
        "recommended_next_profiles": profiles,
        "bounded_concurrency": {
            "max_active_segments": lane_limit,
            "max_wave_nodes": max_wave_nodes,
            "r17_exclusive_and_last": True,
        },
    }
    if mode == "sole-next-profile":
        recommendation_output["controller_consumption_state"] = (
            "non-authorizing until exact independent review and publication audit"
        )
        recommendation_output["sole_next_segment_recommendation"] = {
            key: profiles[0][key] for key in SOLE_COMPATIBLE_PROFILE_KEYS
        }
    else:
        recommendation_output["controller_consumption_state"] = (
            "non-authorizing bounded-wave recommendation pending exact independent "
            "review, publication audit, and controller-mediated consumption"
        )
    recommendation_output["non_authorizing_reason"] = (
        "A recost request and generated candidate are evidence only; neither may "
        "authorize prepare, submit, scheduler, controller, or canonical mutations."
    )

    provenance: dict[str, object] = {
        "request_sha256": args.expected_request_sha256,
        "request_independent_review_sha256": request_review_sha256,
        "generator_sha256": generator_sha256,
        "generator_revision": generator_revision,
        "stage_i_helper_sha256": helper_sha256,
        "stage_i_helper_revision": helper_revision,
        "matrix_sha256": matrix_sha256,
        "matrix_revision": matrix_revision,
        "source_bundle_sha256": bundle_sha256,
        "source_bundle_verified_revisions": bundle_revisions,
        "source_authority": source_authority,
        "qualification_approval_sha256": qualification_sha256,
        "ceiling_evidence_sha256": ceiling_sha256,
        "ceiling_publication_audit_sha256": ceiling_audit_sha256,
        "f113_historical_helper_revision": f113_helper_revision,
        "f113_historical_helper_sha256": f113_helper_sha256,
        "storage_evidence_sha256": storage_sha256,
        "reconciliation_sha256": reconcile_sha256,
        "ledger_sha256": ledger_sha256,
        "reservations_sha256": reservations_sha256,
        "scheduler_evidence": scheduler,
        "predecessor_recost_sha256": predecessor_sha256,
        "predecessor_recost_independent_review_sha256": predecessor_review_sha256,
        "predecessor_recost_publication_audit_sha256": predecessor_audit_sha256,
        "authenticated_lineages_sha256": lineage_sha256,
        "computed_projection_sha256": projection_sha256,
        "r17_readiness_evidence_sha256": readiness_sha256,
    }
    if mode == "sole-next-profile" and len(scheduler) == 1:
        provenance["scheduler_sha256"] = scheduler[0]["sha256"]

    data = {
        "schema_version": 2,
        "record_type": "stage-i-recost-recommendation-evidence",
        "checkpoint": request["checkpoint"],
        "artifact_name": artifact_name,
        "execution_epoch": EXECUTION_EPOCH,
        "generated_utc": request["generated_utc"],
        "expires_utc": request["expires_utc"],
        "requested_by": request["requested_by"],
        "scope": request["scope"],
        "predecessor_recost": predecessor,
        "authority": {
            "authorizing": False,
            "action_authority": "none-until-independent-review-and-publication",
            "scheduler_mutation_authorized": False,
            "canonical_mutation_authorized": False,
        },
        "publication_requirements": {
            "independent_review_required": True,
            "publication_audit_required": True,
            "published_mode": "0444",
            "published_links": 1,
            "controller_consumption_requires_exact_published_sha256": True,
        },
        "recommendations": recommendation_output,
        "barrier": {
            "job_ids": sorted(item["job_id"] for item in barrier),
            "recorded_segments": barrier,
            "scheduler_evidence": scheduler,
        },
        "budget": budget,
        "storage": {
            "available_bytes": storage["available_bytes"],
            "retained_stage_i_bytes": storage["retained_stage_i_bytes"],
            "required_safety_bytes": storage["required_safety_bytes"],
            "projected_authorized_wave_growth_bytes": projected_storage,
            "headroom_after_authorized_wave_and_safety_bytes": (
                storage["available_bytes"]
                - storage["required_safety_bytes"]
                - projected_storage
            ),
        },
        "ledger": {
            "rows": len(rows),
            "sha256": ledger_sha256,
            "cumulative_stage_i_node_hours": rows[-1]["cumulative_stage_i_node_hours"],
        },
        "reservations": {
            "rows": len(reservations),
            "sha256": reservations_sha256,
            "active": 0,
        },
        "manifests": {
            "rows": len(manifest_bindings),
            "bindings": manifest_bindings,
            "authenticated_lineages_sha256": lineage_sha256,
        },
        "r17_readiness": readiness,
        "promoted_f113": {
            "path": str(ceiling_path),
            "sha256": ceiling_sha256,
            "publication_audit_path": str(ceiling_audit_path),
            "publication_audit_sha256": ceiling_audit_sha256,
            "publication_audit": ceiling_audit,
        },
        "reconcile": reconcile,
        "provenance": provenance,
    }
    return BuildResult(
        payload=(json.dumps(data, indent=2, sort_keys=True) + "\n").encode(),
        tracker=tracker,
        reconcile=reconcile,
        helper_path=helper_path,
        helper_sha256=helper_sha256,
        helper_revision=helper_revision,
        matrix_path=matrix_path,
        matrix_sha256=matrix_sha256,
        matrix_revision=matrix_revision,
        storage_available_bytes=int(storage["available_bytes"]),
        storage_retained_stage_i_bytes=int(storage["retained_stage_i_bytes"]),
        storage_required_safety_bytes=int(storage["required_safety_bytes"]),
        projected_storage_bytes=projected_storage,
        storage_measurements=storage_measurements,
    )


def stable_json_bytes(value: object) -> bytes:
    """Return the repository's canonical indented JSON serialization."""

    return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()


def fsync_directory(path: Path) -> None:
    """Persist one trusted directory's metadata."""

    with absolute_descriptor(
        path, "durability directory", flags=os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
    ) as descriptor:
        os.fsync(descriptor)


def file_security_content_binding(profile: os.stat_result) -> tuple[int, ...]:
    """Return the complete mutation-relevant profile of one retained file."""

    return (
        stat.S_IFMT(profile.st_mode),
        stat.S_IMODE(profile.st_mode),
        profile.st_uid,
        profile.st_gid,
        profile.st_nlink,
        profile.st_size,
        profile.st_mtime_ns,
        profile.st_ctime_ns,
    )


def directory_security_binding(profile: os.stat_result) -> tuple[int, int, int, int]:
    """Return the exact authority profile of one retained directory."""

    return (
        stat.S_IFMT(profile.st_mode),
        stat.S_IMODE(profile.st_mode),
        profile.st_uid,
        profile.st_gid,
    )


def require_parent_path_bound(
    path: Path,
    directory: int,
    expected: os.stat_result,
    label: str,
) -> None:
    """Require an open mutation parent and its public pathname to remain exact."""

    opened = os.fstat(directory)
    require_directory_profile(opened, f"{label} parent")
    if (
        not same_inode(opened, expected)
        or directory_security_binding(opened) != directory_security_binding(expected)
    ):
        raise ValueError(f"{label} parent descriptor authority changed")
    with absolute_descriptor(
        path,
        f"{label} parent",
        flags=os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
    ) as named:
        named_profile = os.fstat(named)
        require_directory_profile(named_profile, f"{label} parent")
        if (
            not same_inode(named_profile, expected)
            or directory_security_binding(named_profile)
            != directory_security_binding(expected)
        ):
            raise ValueError(f"{label} parent pathname changed during publication")


def read_bound_exact_file(
    directory: int,
    name: str,
    payload: bytes,
    mode: int,
    label: str,
    *,
    expected_identity: os.stat_result | None = None,
    expected_links: int = 1,
) -> os.stat_result:
    """Read and bind one exact direct-child file without accepting profile drift."""

    observed = os.stat(name, dir_fd=directory, follow_symlinks=False)
    require_publication_link_profile(
        observed, label, expected_mode=mode, expected_links=expected_links
    )
    if expected_identity is not None and not same_inode(observed, expected_identity):
        raise ValueError(f"{label} inode identity differs")
    descriptor = os.open(
        name, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC, dir_fd=directory
    )
    try:
        before = os.fstat(descriptor)
        require_publication_link_profile(
            before, label, expected_mode=mode, expected_links=expected_links
        )
        if (
            not same_inode(before, observed)
            or file_security_content_binding(before)
            != file_security_content_binding(observed)
        ):
            raise ValueError(f"{label} pathname changed before exact verification")
        if sha256_descriptor(descriptor) != sha256_bytes(payload):
            raise ValueError(f"{label} exists with different bytes; refusing to clobber")
        after = os.fstat(descriptor)
        require_publication_link_profile(
            after, label, expected_mode=mode, expected_links=expected_links
        )
        if (
            not same_inode(after, before)
            or file_security_content_binding(after) != file_security_content_binding(before)
        ):
            raise ValueError(f"{label} file profile changed during exact verification")
    finally:
        os.close(descriptor)
    named = os.stat(name, dir_fd=directory, follow_symlinks=False)
    require_publication_link_profile(
        named, label, expected_mode=mode, expected_links=expected_links
    )
    if (
        not same_inode(named, after)
        or file_security_content_binding(named) != file_security_content_binding(after)
    ):
        raise ValueError(f"{label} pathname changed during exact verification")
    return named


def durably_verify_direct_final(
    directory: int,
    name: str,
    payload: bytes,
    mode: int,
    label: str,
    *,
    expected_identity: os.stat_result | None = None,
    mutation_guard: Callable[[], None] | None = None,
) -> os.stat_result:
    """Bind exact final bytes and profile on both sides of a directory fsync."""

    first = read_bound_exact_file(
        directory,
        name,
        payload,
        mode,
        label,
        expected_identity=expected_identity,
    )
    if mutation_guard is not None:
        mutation_guard()
    os.fsync(directory)
    if mutation_guard is not None:
        mutation_guard()
    second = read_bound_exact_file(
        directory,
        name,
        payload,
        mode,
        label,
        expected_identity=first,
    )
    if file_security_content_binding(second) != file_security_content_binding(first):
        raise ValueError(f"{label} file profile changed across durable verification")
    return second


DIRECT_FINAL_ALLOWED_MODES = frozenset({0o444, 0o644, 0o755})
PUBLICATION_PRIVATE_MODE = 0o600
PUBLICATION_PRIVATE_ATTEMPT_LIMIT = 4096
PUBLICATION_PRIVATE_NAME_PATTERN = re.compile(
    r"\.cgl-lf-recost-publication-([0-9a-f]{64})-([0-9]{8})\.private"
)


@dataclass(frozen=True)
class PublicationTransaction:
    """Bind one no-replace publication to exact intent and producer identity."""

    transaction_id: str
    target: str
    payload_sha256: str
    payload_size: int
    final_mode: int
    label: str
    producer: str
    producer_revision: str

    @property
    def private_prefix(self) -> str:
        """Return the deterministic private-attempt namespace for this transaction."""

        return f".cgl-lf-recost-publication-{self.transaction_id}-"

    def private_name(self, attempt: int) -> str:
        """Return one exact deterministic private-attempt name."""

        if not 0 <= attempt < PUBLICATION_PRIVATE_ATTEMPT_LIMIT:
            raise ValueError("publication private attempt is outside the managed range")
        return f"{self.private_prefix}{attempt:08d}.private"


def publication_producer_revision() -> str:
    """Return the exact content revision of the running recost producer."""

    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None and re.fullmatch(r"[0-9]+", inherited) is not None:
        return sha256_descriptor(int(inherited))
    return sha256_bytes(Path(__file__).read_bytes())


def publication_transaction(
    path: Path,
    payload: bytes,
    mode: int,
    label: str,
) -> PublicationTransaction:
    """Construct one content-addressed publication transaction binding."""

    target = str(path.absolute())
    payload_sha256 = sha256_bytes(payload)
    producer_revision = publication_producer_revision()
    binding = {
        "schema_version": 1,
        "record_type": "stage-i-recost-private-publication-transaction",
        "execution_epoch": EXECUTION_EPOCH,
        "target": target,
        "payload_sha256": payload_sha256,
        "payload_size": len(payload),
        "final_mode": f"{mode:04o}",
        "label": label,
        "producer": RECOST_RELATIVE.as_posix(),
        "producer_revision": producer_revision,
    }
    return PublicationTransaction(
        transaction_id=sha256_bytes(stable_json_bytes(binding)),
        target=target,
        payload_sha256=payload_sha256,
        payload_size=len(payload),
        final_mode=mode,
        label=label,
        producer=RECOST_RELATIVE.as_posix(),
        producer_revision=producer_revision,
    )


def require_publication_link_profile(
    profile: os.stat_result,
    label: str,
    *,
    expected_mode: int,
    expected_links: int,
) -> None:
    """Require one owned, exact-mode regular publication inode and link count."""

    if not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"{label} must be a regular file")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"{label} must be owned by the effective user")
    if profile.st_nlink != expected_links:
        raise ValueError(f"{label} must have exactly {expected_links} links")
    retained_mode = stat.S_IMODE(profile.st_mode)
    if retained_mode != expected_mode:
        raise ValueError(
            f"{label} mode is {retained_mode:04o}, expected {expected_mode:04o}"
        )
    if retained_mode & 0o022:
        raise ValueError(f"{label} must not be group- or world-writable")


def require_bound_private_attempt(
    directory: int,
    name: str,
    descriptor: int,
    expected_identity: os.stat_result,
    label: str,
    *,
    expected_profile: os.stat_result | None = None,
    expected_mode: int = PUBLICATION_PRIVATE_MODE,
    expected_links: int = 1,
) -> os.stat_result:
    """Bind one private publication attempt to its exact descriptor and name."""

    opened = os.fstat(descriptor)
    require_publication_link_profile(
        opened,
        label,
        expected_mode=expected_mode,
        expected_links=expected_links,
    )
    if not same_inode(opened, expected_identity):
        raise ValueError(f"{label} private-attempt descriptor identity changed")
    if (
        expected_profile is not None
        and file_security_content_binding(opened)
        != file_security_content_binding(expected_profile)
    ):
        raise ValueError(f"{label} private-attempt profile changed between mutations")
    named = os.stat(name, dir_fd=directory, follow_symlinks=False)
    require_publication_link_profile(
        named,
        label,
        expected_mode=expected_mode,
        expected_links=expected_links,
    )
    if (
        not same_inode(named, opened)
        or file_security_content_binding(named) != file_security_content_binding(opened)
    ):
        raise ValueError(f"{label} private-attempt pathname changed during publication")
    return opened


def authenticate_private_attempt(
    path: Path,
    directory: int,
    parent_profile: os.stat_result,
    private_name: str,
    descriptor: int,
    expected_identity: os.stat_result,
    label: str,
    mutation_lock: MutationLock | None,
    *,
    expected_profile: os.stat_result | None = None,
    expected_mode: int = PUBLICATION_PRIVATE_MODE,
    expected_links: int = 1,
) -> os.stat_result:
    """Revalidate all authority and one private attempt before mutation."""

    require_parent_path_bound(path.parent, directory, parent_profile, label)
    authenticate_mutation_lock(mutation_lock)
    retained = require_bound_private_attempt(
        directory,
        private_name,
        descriptor,
        expected_identity,
        label,
        expected_profile=expected_profile,
        expected_mode=expected_mode,
        expected_links=expected_links,
    )
    require_parent_path_bound(path.parent, directory, parent_profile, label)
    authenticate_mutation_lock(mutation_lock)
    return retained


def create_private_attempt(
    directory: int,
    transaction: PublicationTransaction,
    label: str,
    mutation_guard: Callable[[], None],
) -> tuple[str, int, os.stat_result]:
    """Create one fresh transaction-bound private attempt without reuse."""

    retained_umask = os.umask(0)
    try:
        for attempt in range(PUBLICATION_PRIVATE_ATTEMPT_LIMIT):
            name = transaction.private_name(attempt)
            try:
                mutation_guard()
                descriptor = os.open(
                    name,
                    os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC,
                    PUBLICATION_PRIVATE_MODE,
                    dir_fd=directory,
                )
            except FileExistsError:
                continue
            break
        else:
            raise ValueError(f"{label} exhausted the managed private-attempt namespace")
    finally:
        os.umask(retained_umask)
    try:
        created = os.fstat(descriptor)
        require_publication_link_profile(
            created,
            label,
            expected_mode=PUBLICATION_PRIVATE_MODE,
            expected_links=1,
        )
        retained = require_bound_private_attempt(
            directory,
            name,
            descriptor,
            created,
            label,
            expected_profile=created,
        )
    except BaseException:
        os.close(descriptor)
        raise
    return name, descriptor, retained


def finalize_private_attempt(
    directory: int,
    path: Path,
    parent_profile: os.stat_result,
    private_name: str,
    descriptor: int,
    private_profile: os.stat_result,
    payload: bytes,
    mode: int,
    label: str,
    mutation_lock: MutationLock | None,
) -> os.stat_result:
    """Write, fsync, and finalize one freshly created private inode."""

    identity = private_profile
    retained = authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
        expected_profile=private_profile,
    )
    offset = 0
    while offset < len(payload):
        authenticate_private_attempt(
            path,
            directory,
            parent_profile,
            private_name,
            descriptor,
            identity,
            label,
            mutation_lock,
            expected_profile=retained,
        )
        written = os.write(descriptor, payload[offset:])
        if written <= 0:
            raise ValueError(f"{label} write made no progress")
        offset += written
        retained = authenticate_private_attempt(
            path,
            directory,
            parent_profile,
            private_name,
            descriptor,
            identity,
            label,
            mutation_lock,
        )
    authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
        expected_profile=retained,
    )
    os.fsync(descriptor)
    retained = authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
    )
    if (
        retained.st_size != len(payload)
        or sha256_descriptor(descriptor) != sha256_bytes(payload)
    ):
        raise ValueError(f"{label} private-attempt bytes differ")
    retained = authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
        expected_profile=retained,
    )
    os.fchmod(descriptor, mode)
    finalized = authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
        expected_mode=mode,
    )
    os.fsync(descriptor)
    finalized = authenticate_private_attempt(
        path,
        directory,
        parent_profile,
        private_name,
        descriptor,
        identity,
        label,
        mutation_lock,
        expected_profile=finalized,
        expected_mode=mode,
    )
    if (
        finalized.st_size != len(payload)
        or sha256_descriptor(descriptor) != sha256_bytes(payload)
    ):
        raise ValueError(f"{label} finalized private-attempt bytes differ")
    return finalized


def transaction_private_names(
    directory: int,
    transaction: PublicationTransaction,
) -> list[str]:
    """Return exact private-attempt names for one transaction."""

    retained = []
    for entry in os.scandir(directory):
        match = PUBLICATION_PRIVATE_NAME_PATTERN.fullmatch(entry.name)
        if (
            match is not None
            and match.group(1) == transaction.transaction_id
            and int(match.group(2)) < PUBLICATION_PRIVATE_ATTEMPT_LIMIT
        ):
            retained.append(entry.name)
    return sorted(retained)


def recover_linked_publication(
    directory: int,
    path: Path,
    parent_profile: os.stat_result,
    transaction: PublicationTransaction,
    payload: bytes,
    mode: int,
    label: str,
    mutation_lock: MutationLock | None,
    *,
    expected_identity: os.stat_result | None = None,
) -> os.stat_result | None:
    """Finish exact same-inode private-link cleanup after publication."""

    try:
        target = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
    except FileNotFoundError:
        return None
    require_publication_link_profile(
        target, label, expected_mode=mode, expected_links=2
    )
    if expected_identity is not None and not same_inode(target, expected_identity):
        raise ValueError(f"{label} public target inode differs from the linked private inode")
    matches = []
    for private_name in transaction_private_names(directory, transaction):
        private = os.stat(private_name, dir_fd=directory, follow_symlinks=False)
        if same_inode(private, target):
            require_publication_link_profile(
                private,
                f"{label} transaction-private link",
                expected_mode=mode,
                expected_links=2,
            )
            matches.append((private_name, private))
    if len(matches) != 1:
        raise ValueError(
            f"{label} linked public target lacks one exact transaction-private inode binding"
        )
    private_name, private = matches[0]
    bound_target = read_bound_exact_file(
        directory,
        path.name,
        payload,
        mode,
        label,
        expected_identity=target,
        expected_links=2,
    )
    bound_private = read_bound_exact_file(
        directory,
        private_name,
        payload,
        mode,
        f"{label} transaction-private link",
        expected_identity=private,
        expected_links=2,
    )
    require_parent_path_bound(path.parent, directory, parent_profile, label)
    authenticate_mutation_lock(mutation_lock)
    target = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
    private = os.stat(private_name, dir_fd=directory, follow_symlinks=False)
    require_publication_link_profile(
        target, label, expected_mode=mode, expected_links=2
    )
    require_publication_link_profile(
        private,
        f"{label} transaction-private link",
        expected_mode=mode,
        expected_links=2,
    )
    if (
        not same_inode(target, bound_target)
        or not same_inode(private, bound_private)
        or file_security_content_binding(target)
        != file_security_content_binding(bound_target)
        or file_security_content_binding(private)
        != file_security_content_binding(bound_private)
        or not same_inode(target, private)
    ):
        raise ValueError(f"{label} transaction-private inode binding changed before unlink")
    try:
        os.unlink(private_name, dir_fd=directory)
    except BaseException as unlink_error:
        try:
            return durably_verify_direct_final(
                directory,
                path.name,
                payload,
                mode,
                label,
                expected_identity=target,
                mutation_guard=lambda: (
                    require_parent_path_bound(
                        path.parent, directory, parent_profile, label
                    ),
                    authenticate_mutation_lock(mutation_lock),
                ),
            )
        except BaseException:
            raise unlink_error
    require_parent_path_bound(path.parent, directory, parent_profile, label)
    authenticate_mutation_lock(mutation_lock)
    return durably_verify_direct_final(
        directory,
        path.name,
        payload,
        mode,
        label,
        expected_identity=target,
        mutation_guard=lambda: (
            require_parent_path_bound(path.parent, directory, parent_profile, label),
            authenticate_mutation_lock(mutation_lock),
        ),
    )


def raise_direct_final_failure(
    directory: int,
    path: Path,
    parent_profile: os.stat_result,
    label: str,
    error: BaseException,
    mutation_lock: MutationLock | None,
) -> None:
    """Durably classify a failed direct-final publication without another mutation."""

    failures: list[BaseException] = []
    try:
        require_parent_path_bound(path.parent, directory, parent_profile, label)
        authenticate_mutation_lock(mutation_lock)
    except BaseException as failure:
        failures.append(failure)
    if not failures:
        try:
            os.fsync(directory)
        except BaseException as failure:
            failures.append(failure)
        try:
            require_parent_path_bound(path.parent, directory, parent_profile, label)
            authenticate_mutation_lock(mutation_lock)
        except BaseException as failure:
            failures.append(failure)
    try:
        retained = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
    except FileNotFoundError:
        state = "absent"
    except BaseException:
        state = "unreadable"
    else:
        state = (
            f"occupied inode {retained.st_dev}:{retained.st_ino} "
            f"mode {stat.S_IMODE(retained.st_mode):04o} links {retained.st_nlink}"
        )
    if failures:
        raise ValueError(
            f"{label} direct-final publication failed; deterministic public target is "
            f"{state}; durability or authority revalidation failed and no rollback was attempted"
        ) from failures[0]
    raise ValueError(
        f"{label} direct-final publication failed; deterministic public target is "
        f"{state}; no rollback was attempted"
    ) from error


def write_exact_or_verify(
    path: Path,
    payload: bytes,
    *,
    mode: int,
    label: str,
    mutation_lock: MutationLock | None = None,
) -> bool:
    """Publish one finalized private inode no-replace, or exact-verify the target."""

    if mode not in DIRECT_FINAL_ALLOWED_MODES:
        raise ValueError(f"{label} final mode {mode:04o} is not a managed publication mode")
    transaction = publication_transaction(path, payload, mode, label)
    authenticate_mutation_lock(mutation_lock)
    require_no_symlink_components(path, label, include_leaf=False)
    with absolute_descriptor(
        path.parent,
        f"{label} parent",
        flags=os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
    ) as directory:
        parent_profile = os.fstat(directory)
        require_directory_profile(parent_profile, f"{label} parent")
        require_parent_path_bound(path.parent, directory, parent_profile, label)

        def mutation_guard() -> None:
            require_parent_path_bound(path.parent, directory, parent_profile, label)
            authenticate_mutation_lock(mutation_lock)

        try:
            observed = os.stat(path.name, dir_fd=directory, follow_symlinks=False)
        except FileNotFoundError:
            observed = None
        if observed is not None:
            if (
                stat.S_IMODE(observed.st_mode) == mode
                and observed.st_nlink == 2
            ):
                try:
                    recover_linked_publication(
                        directory,
                        path,
                        parent_profile,
                        transaction,
                        payload,
                        mode,
                        label,
                        mutation_lock,
                        expected_identity=observed,
                    )
                except BaseException as error:
                    raise_direct_final_failure(
                        directory, path, parent_profile, label, error, mutation_lock
                    )
                return True
            durably_verify_direct_final(
                directory,
                path.name,
                payload,
                mode,
                label,
                mutation_guard=mutation_guard,
            )
            mutation_guard()
            return False

        mutation_guard()
        try:
            private_name, descriptor, private_profile = create_private_attempt(
                directory, transaction, label, mutation_guard
            )
        except BaseException as error:
            raise_direct_final_failure(
                directory, path, parent_profile, label, error, mutation_lock
            )

        try:
            finalized_profile = finalize_private_attempt(
                directory,
                path,
                parent_profile,
                private_name,
                descriptor,
                private_profile,
                payload,
                mode,
                label,
                mutation_lock,
            )
            mutation_guard()
            finalized_profile = authenticate_private_attempt(
                path,
                directory,
                parent_profile,
                private_name,
                descriptor,
                private_profile,
                label,
                mutation_lock,
                expected_profile=finalized_profile,
                expected_mode=mode,
            )
            try:
                os.link(
                    private_name,
                    path.name,
                    src_dir_fd=directory,
                    dst_dir_fd=directory,
                    follow_symlinks=False,
                )
            except BaseException as link_error:
                try:
                    recovered = recover_linked_publication(
                        directory,
                        path,
                        parent_profile,
                        transaction,
                        payload,
                        mode,
                        label,
                        mutation_lock,
                        expected_identity=finalized_profile,
                    )
                except BaseException:
                    raise link_error
                if recovered is None:
                    raise link_error
            else:
                recover_linked_publication(
                    directory,
                    path,
                    parent_profile,
                    transaction,
                    payload,
                    mode,
                    label,
                    mutation_lock,
                    expected_identity=finalized_profile,
                )
        except BaseException as error:
            try:
                os.close(descriptor)
            except BaseException:
                pass
            raise_direct_final_failure(
                directory, path, parent_profile, label, error, mutation_lock
            )
        try:
            os.close(descriptor)
        except BaseException as error:
            raise_direct_final_failure(
                directory, path, parent_profile, label, error, mutation_lock
            )

        mutation_guard()
    return True


def parse_draft_packet(payload: bytes) -> dict[str, object]:
    """Parse one explicit non-authorizing request-draft packet."""

    packet = require_exact_keys(
        parse_json(payload, "recost request draft packet"),
        {
            "schema_version",
            "record_type",
            "checkpoint",
            "artifact_name",
            "execution_epoch",
            "generated_utc",
            "expires_utc",
            "requested_by",
            "scope",
            "barrier",
            "inputs",
            "recommendations",
            "draft_policy",
        },
        "recost request draft packet",
    )
    if (
        packet["schema_version"] != 1
        or packet["record_type"] != "stage-i-recost-request-draft-packet"
        or packet["execution_epoch"] != EXECUTION_EPOCH
    ):
        raise ValueError("recost request draft packet identity differs")
    policy = require_exact_keys(
        packet["draft_policy"],
        {
            "independent_review_created",
            "self_approved",
            "scheduler_mutation_authorized",
            "canonical_mutation_authorized",
            "required_storage_safety_bytes",
        },
        "recost request draft policy",
    )
    if (
        policy["independent_review_created"] is not False
        or policy["self_approved"] is not False
        or policy["scheduler_mutation_authorized"] is not False
        or policy["canonical_mutation_authorized"] is not False
    ):
        raise ValueError("recost request draft packet over-authorizes its workflow")
    require_integer(
        policy["required_storage_safety_bytes"],
        "recost request draft storage safety bytes",
        minimum=1,
    )
    request = {key: value for key, value in packet.items() if key != "draft_policy"}
    request["record_type"] = "stage-i-recost-recommendation-request"
    request["schema_version"] = 2
    inputs = request.get("inputs")
    if not isinstance(inputs, dict):
        raise ValueError("recost request draft inputs must be an object")
    for key in (
        "reconciliation",
        "ledger",
        "reservations",
        "manifests",
        "scheduler_evidence",
        "storage_evidence",
    ):
        if inputs.get(key) is not None:
            raise ValueError(f"recost request draft must leave live binding {key} unset")
    return packet


def draft_output_paths(root: Path, checkpoint: object, artifact_name: object) -> dict[str, Path]:
    """Return canonical direct-child prerequisite paths for one checkpoint."""

    retained_checkpoint = require_safe_id(checkpoint, "recost draft checkpoint")
    retained_artifact = require_nonempty_string(artifact_name, "recost draft artifact name")
    match = RECOST_ARTIFACT_PATTERN.fullmatch(retained_artifact)
    if match is None or retained_checkpoint != f"F-{match.group(1)}":
        raise ValueError("recost draft artifact name and checkpoint ID differ")
    prefix = f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F{match.group(1)}"
    accounting = root / "accounting"
    return {
        "reconciliation": accounting / f"{prefix}_reconciliation_evidence.json",
        "storage": accounting / f"{prefix}_storage_evidence.json",
        "request": accounting / f"{prefix}_recost_request.json",
        "artifact": accounting / retained_artifact,
    }


def discovered_root_binding(
    root: Path,
    path: Path,
    tracker: InputTracker,
    label: str,
    *,
    expected_mode: int,
) -> dict[str, str]:
    """Discover and authenticate one root-relative binding."""

    _, digest = tracker.discover(path, label, expected_mode=expected_mode)
    return {"path": path.relative_to(root).as_posix(), "sha256": digest}


def locked_install_f117_draft_packet(
    args: argparse.Namespace,
    root: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Install or exact-verify the initial F117 packet without out-of-band writes."""

    assert args.packet is not None
    assert args.expected_packet_sha256 is not None
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    tracker = InputTracker()
    packet_source = args.packet.absolute()
    packet_payload = tracker.read(
        packet_source,
        args.expected_packet_sha256,
        "initial F117 recost request draft packet source",
        expected_mode=0o644,
    )
    packet = parse_draft_packet(packet_payload)
    paths = draft_output_paths(root, packet["checkpoint"], packet["artifact_name"])
    expected_artifact = (
        root
        / "accounting"
        / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F117_recost_evidence.json"
    )
    if packet["checkpoint"] != "F-117" or paths["artifact"] != expected_artifact:
        raise ValueError("initial draft-packet installation is restricted to exact F-117")
    target = (
        root
        / "accounting"
        / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_F117_recost_draft_packet.json"
    )
    created = write_exact_or_verify(
        target,
        packet_payload,
        mode=0o644,
        label="managed initial F117 recost request draft packet",
        mutation_lock=mutation_lock,
    )
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    tracker.reauthenticate_all()
    authenticate_mutation_lock(mutation_lock)
    next_action = [
        str(source_path),
        "--root",
        str(root),
        "--packet",
        str(target),
        "--expected-packet-sha256",
        args.expected_packet_sha256,
        "--expected-generator-sha256",
        generator_sha256,
    ]
    if root != DEFAULT_ROOT:
        next_action.append("--allow-local-root")
        if args.squeue_file is not None:
            next_action.extend(["--squeue-file", str(args.squeue_file.absolute())])
    next_action.append("draft-request")
    print(
        json.dumps(
            {
                "action": "install-f117-draft-packet",
                "path": str(target),
                "sha256": args.expected_packet_sha256,
                "mode": "0644",
                "created": created,
                "no_clobber": True,
                "exact_existing_copy_verified": not created,
                "next_action": next_action,
            },
            indent=2,
            sort_keys=True,
        )
    )


def locked_draft_request(
    args: argparse.Namespace,
    root: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Draft exact live prerequisites without creating an independent review."""

    assert args.packet is not None
    assert args.expected_packet_sha256 is not None
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    tracker = InputTracker()
    packet_path = args.packet.absolute()
    if packet_path.parent != root / "accounting":
        raise ValueError("recost request draft packet must be a direct child of accounting")
    packet = parse_draft_packet(
        tracker.read(
            packet_path,
            args.expected_packet_sha256,
            "recost request draft packet",
            expected_mode=0o644,
        )
    )
    paths = draft_output_paths(root, packet["checkpoint"], packet["artifact_name"])
    review_path = paths["request"].with_name(
        f"{paths['request'].name}.independent_review.json"
    )
    if os.path.lexists(review_path):
        raise ValueError("draft request independent review already exists; refusing redraft")
    request_timestamp = parse_utc_timestamp(
        packet["generated_utc"], "recost request draft generation timestamp"
    )
    expires_timestamp = parse_utc_timestamp(
        packet["expires_utc"], "recost request draft expiry timestamp"
    )
    now = datetime.now(timezone.utc)
    if (
        request_timestamp > now + timedelta(minutes=5)
        or now - request_timestamp > STORAGE_EVIDENCE_MAX_AGE
        or expires_timestamp <= request_timestamp
        or expires_timestamp - request_timestamp > AUTHORIZATION_MAX_LIFETIME
        or now > expires_timestamp
    ):
        raise ValueError("recost request draft lifetime is invalid")
    request = {key: value for key, value in packet.items() if key != "draft_policy"}
    request["schema_version"] = 2
    request["record_type"] = "stage-i-recost-recommendation-request"
    inputs = request["inputs"]
    assert isinstance(inputs, dict)

    helper_binding = require_exact_keys(
        inputs["stage_i_helper"], {"revision", "sha256"}, "Stage I helper binding"
    )
    helper_sha256 = require_sha256(helper_binding["sha256"], "Stage I helper SHA-256")
    helper_revision = require_revision(helper_binding["revision"], "Stage I helper revision")
    helper_path = repository / STAGE_I_RELATIVE
    require_committed_file(
        repository,
        helper_path,
        helper_sha256,
        "Stage I helper",
        revision=helper_revision,
        expected_mode=0o644,
    )
    tracker.authenticate(helper_path, helper_sha256, "Stage I helper", expected_mode=0o644)
    reconcile = run_authenticated_reconcile(
        helper_path, helper_sha256, root, stage_i_lock_held=True
    )
    parse_reconciliation(reconcile, root)

    ledger_path = root / "accounting" / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv"
    ledger_binding = discovered_root_binding(
        root, ledger_path, tracker, "Stage I ledger", expected_mode=0o644
    )
    rows = parse_ledger(
        tracker.read(
            ledger_path, ledger_binding["sha256"], "Stage I ledger", expected_mode=0o644
        ),
        request_timestamp,
    )
    reservations_path = (
        root / "accounting" / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json"
    )
    reservations_binding = discovered_root_binding(
        root, reservations_path, tracker, "Stage I reservations", expected_mode=0o644
    )
    reservations = parse_reservations(
        parse_json(
            tracker.read(
                reservations_path,
                reservations_binding["sha256"],
                "Stage I reservations",
                expected_mode=0o644,
            ),
            "Stage I reservations",
        )
    )
    manifest_bindings = [
        discovered_root_binding(
            root,
            root / relative,
            tracker,
            f"Stage I manifest {relative}",
            expected_mode=0o644,
        )
        for relative in sorted(manifest_inventory(root))
    ]
    _, manifests = parse_manifests(root, manifest_bindings, tracker)
    counts = reconcile["counts"]
    assert isinstance(counts, dict)
    manifests_by_job = cross_validate_accounting(rows, reservations, manifests, counts)
    barrier_value = request["barrier"]
    assert isinstance(barrier_value, dict)
    barrier = parse_barrier_segments(barrier_value["recorded_segments"], manifests_by_job)
    ledger_suffix = rows[-len(barrier):]
    scheduler_bindings = [
        discovered_root_binding(
            root,
            root / "accounting" / f"{item['job_id']}.stage_i.sacct.txt",
            tracker,
            f"barrier scheduler evidence {item['job_id']}",
            expected_mode=0o644,
        )
        for item in barrier
    ]
    parse_scheduler_evidence(
        root,
        scheduler_bindings,
        barrier,
        ledger_suffix,
        manifests_by_job,
        tracker,
        request_timestamp,
    )

    matrix_relative, matrix_revision, matrix_sha256 = repository_binding(
        inputs["matrix"], "Stage I matrix binding"
    )
    matrix_path = repository_path(repository, matrix_relative.as_posix(), "Stage I matrix")
    require_committed_file(
        repository,
        matrix_path,
        matrix_sha256,
        "Stage I matrix",
        revision=matrix_revision,
    )
    matrix, _ = parse_matrix(
        parse_json(tracker.read(matrix_path, matrix_sha256, "Stage I matrix"), "Stage I matrix")
    )
    recommendations = request["recommendations"]
    assert isinstance(recommendations, dict)
    profiles = recommendations.get("profiles")
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("recost request draft profiles must be a nonempty list")
    projected_storage, storage_projections, storage_measurements = observed_storage_projections(
        root, manifests, matrix, profiles
    )
    tracker.retain_directory_measurements(storage_measurements)
    storage_projection_sha256 = sha256_bytes(
        (json.dumps(storage_projections, sort_keys=True) + "\n").encode()
    )
    with absolute_descriptor(
        root, "campaign root", flags=os.O_RDONLY | os.O_DIRECTORY
    ) as descriptor:
        storage_profile = os.fstatvfs(descriptor)
    live_available = storage_profile.f_bavail * storage_profile.f_frsize
    required_safety = require_integer(
        packet["draft_policy"]["required_storage_safety_bytes"],
        "recost request draft storage safety bytes",
        minimum=1,
    )
    reviewed_available = required_safety + projected_storage
    if live_available < reviewed_available:
        raise ValueError("live storage headroom is exhausted during request drafting")
    retained_bytes = directory_tree_regular_bytes(
        root / "runs/mks24-stage-i" / EXECUTION_EPOCH,
        "retained Stage I storage",
    )
    storage = {
        "schema_version": 1,
        "record_type": "stage-i-storage-evidence",
        "execution_epoch": EXECUTION_EPOCH,
        "root": str(root),
        "measured_utc": packet["generated_utc"],
        "available_bytes": reviewed_available,
        "retained_stage_i_bytes": retained_bytes,
        "required_safety_bytes": required_safety,
        "projected_authorized_wave_growth_bytes": projected_storage,
        "projection_method": STORAGE_PROJECTION_METHOD,
        "profile_projections_sha256": storage_projection_sha256,
    }
    reconcile_payload = stable_json_bytes(reconcile)
    storage_payload = stable_json_bytes(storage)
    inputs.update(
        {
            "reconciliation": {
                "path": paths["reconciliation"].relative_to(root).as_posix(),
                "sha256": sha256_bytes(reconcile_payload),
            },
            "ledger": ledger_binding,
            "reservations": reservations_binding,
            "manifests": manifest_bindings,
            "scheduler_evidence": scheduler_bindings,
            "storage_evidence": {
                "path": paths["storage"].relative_to(root).as_posix(),
                "sha256": sha256_bytes(storage_payload),
            },
        }
    )
    request_payload = stable_json_bytes(request)
    parse_request(request_payload)

    for path, payload, label in (
        (paths["reconciliation"], reconcile_payload, "draft reconciliation evidence"),
        (paths["storage"], storage_payload, "draft storage evidence"),
        (paths["request"], request_payload, "schema-2 recost request draft"),
    ):
        write_exact_or_verify(
            path,
            payload,
            mode=0o644,
            label=label,
            mutation_lock=mutation_lock,
        )
    build_args = argparse.Namespace(**vars(args))
    build_args.request = paths["request"]
    build_args.expected_request_sha256 = sha256_bytes(request_payload)
    build_args.output = paths["artifact"].with_name(f"{paths['artifact'].name}.staged")
    build = build_payload(
        build_args,
        root,
        source_path,
        repository,
        generator_sha256,
        stage_i_lock_held=True,
        require_independent_request_review=False,
    )
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    tracker.reauthenticate_all()
    build.tracker.reauthenticate_all()
    live_reconcile = run_authenticated_reconcile(
        helper_path, helper_sha256, root, stage_i_lock_held=True
    )
    if live_reconcile != reconcile:
        raise ValueError("live reconciliation changed during request drafting")
    require_live_storage_boundary(
        root,
        build.storage_available_bytes,
        build.storage_retained_stage_i_bytes,
        build.storage_required_safety_bytes,
        build.projected_storage_bytes,
    )
    generation_command = [
        sys.executable,
        str(source_path),
        "--root",
        str(root),
        "--request",
        str(paths["request"]),
        "--expected-request-sha256",
        sha256_bytes(request_payload),
        "--output",
        str(paths["artifact"].with_name(f"{paths['artifact'].name}.staged")),
        "--expected-generator-sha256",
        generator_sha256,
    ]
    if root != DEFAULT_ROOT:
        generation_command.append("--allow-local-root")
        if args.squeue_file is not None:
            generation_command.extend(["--squeue-file", str(args.squeue_file.absolute())])
    request_sha256 = sha256_bytes(request_payload)
    review_install_command = [
        sys.executable,
        str(source_path),
        "--root",
        str(root),
        "--request",
        str(paths["request"]),
        "--expected-request-sha256",
        request_sha256,
        "--review",
        "<external-immutable-review-candidate>",
        "--expected-review-sha256",
        "<external-review-sha256>",
        "--expected-generator-sha256",
        generator_sha256,
    ]
    if root != DEFAULT_ROOT:
        review_install_command.append("--allow-local-root")
        if args.squeue_file is not None:
            review_install_command.extend(
                ["--squeue-file", str(args.squeue_file.absolute())]
            )
    review_install_command.append("install-request-review")
    print(
        json.dumps(
            {
                "action": "draft-request",
                "request": {
                    "path": str(paths["request"]),
                    "sha256": request_sha256,
                    "mode": "0644",
                },
                "reconciliation": {
                    "path": str(paths["reconciliation"]),
                    "sha256": sha256_bytes(reconcile_payload),
                    "mode": "0644",
                },
                "storage": {
                    "path": str(paths["storage"]),
                    "sha256": sha256_bytes(storage_payload),
                    "mode": "0644",
                },
                "independent_review_created": False,
                "self_approved": False,
                "required_independent_review_path": str(review_path),
                "next_required_action": (
                    "externally review the exact request SHA-256, retain the candidate "
                    "immutable, then run the managed review installer"
                ),
                "managed_independent_review_install": {
                    "target": str(review_path),
                    "target_mode": "0444",
                    "command_template": review_install_command,
                },
                "staged_generation_after_external_review": generation_command,
            },
            indent=2,
            sort_keys=True,
        )
    )


def locked_install_request_review(
    args: argparse.Namespace,
    root: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Install exact externally reviewed request approval without canonical hand writes."""

    assert args.request is not None
    assert args.expected_request_sha256 is not None
    assert args.review is not None
    assert args.expected_review_sha256 is not None
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    tracker = InputTracker()
    request_path = args.request.absolute()
    accounting = root / "accounting"
    if request_path.parent != accounting:
        raise ValueError("recost request review candidate must bind a direct child of accounting")
    request_payload = tracker.read(
        request_path,
        args.expected_request_sha256,
        "recost request review candidate",
        expected_mode=0o644,
    )
    request = parse_request(request_payload)
    paths = draft_output_paths(root, request["checkpoint"], request["artifact_name"])
    if request_path != paths["request"]:
        raise ValueError("recost request review candidate path differs from its identity")
    review_target = request_path.with_name(f"{request_path.name}.independent_review.json")
    review_source = args.review.absolute()
    try:
        review_source.relative_to(accounting)
    except ValueError:
        source_is_canonical_accounting = False
    else:
        source_is_canonical_accounting = True
    if review_source == review_target or source_is_canonical_accounting:
        raise ValueError(
            "independent review source must be an external immutable candidate, "
            "not retained beneath canonical accounting"
        )
    review_payload = tracker.read(
        review_source,
        args.expected_review_sha256,
        "external immutable recost request independent review candidate",
        expected_mode=0o444,
    )
    request_timestamp = parse_utc_timestamp(
        request["generated_utc"], "recost request generation timestamp"
    )
    expires_timestamp = parse_utc_timestamp(
        request["expires_utc"], "recost request expiry timestamp"
    )
    now = datetime.now(timezone.utc)
    if (
        request_timestamp > now + timedelta(minutes=5)
        or now - request_timestamp > STORAGE_EVIDENCE_MAX_AGE
        or expires_timestamp <= request_timestamp
        or expires_timestamp - request_timestamp > AUTHORIZATION_MAX_LIFETIME
        or now > expires_timestamp
    ):
        raise ValueError("recost request review candidate lifetime is invalid")
    parse_request_independent_review(
        parse_json(review_payload, "external recost request independent review candidate"),
        request_path,
        args.expected_request_sha256,
        request_timestamp,
        require_nonempty_string(request["requested_by"], "recost request author"),
    )
    created = write_exact_or_verify(
        review_target,
        review_payload,
        mode=0o444,
        label="managed recost request independent review",
        mutation_lock=mutation_lock,
    )
    tracker.authenticate(
        review_target,
        args.expected_review_sha256,
        "installed recost request independent review",
        expected_mode=0o444,
    )
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    tracker.reauthenticate_all()
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    output = paths["artifact"].with_name(f"{paths['artifact'].name}.staged")
    generation_command = [
        sys.executable,
        str(source_path),
        "--root",
        str(root),
        "--request",
        str(request_path),
        "--expected-request-sha256",
        args.expected_request_sha256,
        "--output",
        str(output),
        "--expected-generator-sha256",
        generator_sha256,
    ]
    if root != DEFAULT_ROOT:
        generation_command.append("--allow-local-root")
        if args.squeue_file is not None:
            generation_command.extend(["--squeue-file", str(args.squeue_file.absolute())])
    print(
        json.dumps(
            {
                "action": "install-request-review",
                "request": {
                    "path": str(request_path),
                    "sha256": args.expected_request_sha256,
                    "mode": "0644",
                },
                "external_candidate": {
                    "path": str(review_source),
                    "sha256": args.expected_review_sha256,
                    "mode": "0444",
                },
                "installed_review": {
                    "path": str(review_target),
                    "sha256": args.expected_review_sha256,
                    "mode": "0444",
                },
                "created": created,
                "no_clobber": True,
                "exact_existing_copy_verified": not created,
                "staged_generation": generation_command,
            },
            indent=2,
            sort_keys=True,
        )
    )


def ensure_utilities_directory(
    root: Path, mutation_lock: MutationLock | None = None
) -> Path:
    """Create or authenticate the direct accounting utilities directory."""

    accounting = root / "accounting"
    utilities = accounting / "utilities"
    if not os.path.lexists(utilities):
        authenticate_mutation_lock(mutation_lock)
        with absolute_descriptor(
            accounting,
            "accounting directory",
            flags=os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC,
        ) as descriptor:
            os.mkdir("utilities", mode=0o755, dir_fd=descriptor)
            os.fsync(descriptor)
        authenticate_mutation_lock(mutation_lock)
    require_directory(utilities, "accounting utilities directory")
    return utilities


def locked_retain_generator(
    args: argparse.Namespace,
    root: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Retain exact committed generator bytes without overwriting any existing copy."""

    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    with absolute_descriptor(source_path, "Stage I recost generator", flags=os.O_RDONLY) as descriptor:
        require_regular_profile(
            os.fstat(descriptor), "Stage I recost generator", expected_mode=0o755
        )
        chunks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            chunks.append(block)
    payload = b"".join(chunks)
    if sha256_bytes(payload) != generator_sha256:
        raise ValueError("Stage I recost generator bytes changed before retention")
    target = ensure_utilities_directory(root, mutation_lock) / "cgl_lf_stage_i_recost.py"
    created = write_exact_or_verify(
        target,
        payload,
        mode=0o755,
        label="retained Stage I recost generator",
        mutation_lock=mutation_lock,
    )
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    authenticate_mutation_lock(mutation_lock)
    print(
        json.dumps(
            {
                "action": "retain-generator",
                "path": str(target),
                "sha256": generator_sha256,
                "mode": "0755",
                "checkpoint_generator_relative_path": (
                    "accounting/utilities/cgl_lf_stage_i_recost.py"
                ),
                "created": created,
                "no_clobber": True,
                "exact_existing_copy_verified": not created,
            },
            indent=2,
            sort_keys=True,
        )
    )


def write_staged_output(
    output: Path,
    payload: bytes,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Publish or exact-verify one explicit staged artifact without overwrite."""

    write_exact_or_verify(
        output,
        payload,
        mode=0o444,
        label="staged recost output",
        mutation_lock=mutation_lock,
    )


def locked_main(
    args: argparse.Namespace,
    root: Path,
    output: Path,
    source_path: Path,
    repository: Path,
    generator_sha256: str,
    mutation_lock: MutationLock | None = None,
) -> None:
    """Generate one staged artifact after final barrier revalidation."""

    require_output_namespace_empty(output)
    require_drained_queue(args, root)
    build = build_payload(
        args,
        root,
        source_path,
        repository,
        generator_sha256,
        stage_i_lock_held=True,
    )
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    build.tracker.reauthenticate_all()
    require_committed_file(
        repository,
        source_path,
        generator_sha256,
        "Stage I recost generator",
        expected_mode=0o755,
    )
    require_committed_file(
        repository,
        build.helper_path,
        build.helper_sha256,
        "Stage I helper",
        revision=build.helper_revision,
        expected_mode=0o644,
    )
    require_committed_file(
        repository,
        build.matrix_path,
        build.matrix_sha256,
        "Stage I matrix",
        revision=build.matrix_revision,
    )
    live_reconcile = run_authenticated_reconcile(
        build.helper_path,
        build.helper_sha256,
        root,
        stage_i_lock_held=True,
    )
    if live_reconcile != build.reconcile:
        raise ValueError("live reconciliation changed before staged output creation")
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_output_namespace_empty(output)
    build.tracker.reauthenticate_all()
    require_live_storage_boundary(
        root,
        build.storage_available_bytes,
        build.storage_retained_stage_i_bytes,
        build.storage_required_safety_bytes,
        build.projected_storage_bytes,
    )
    require_directory_measurement_boundaries(build.storage_measurements)
    require_empty_transaction_stores(root)
    require_drained_queue(args, root)
    require_output_namespace_empty(output)
    write_staged_output(output, build.payload, mutation_lock)
    print(output)


def main(argv: list[str] | None = None) -> int:
    """Run one fail-closed recost generation or prerequisite action."""

    retained_argv = sys.argv[1:] if argv is None else argv
    source_path, repository, generator_sha256 = authenticate_self(retained_argv)
    args = parse_args(retained_argv)
    if args.expected_generator_sha256 != generator_sha256:
        raise ValueError("parsed generator SHA-256 differs from authenticated self")
    root = validate_root(args)
    with stage_i_lock(root) as mutation_lock:
        mutation_lock.authenticate()
        if args.action == "generate":
            _, output = validate_root_and_output(args)
            locked_main(
                args,
                root,
                output,
                source_path,
                repository,
                generator_sha256,
                mutation_lock,
            )
        elif args.action == "install-f117-draft-packet":
            locked_install_f117_draft_packet(
                args, root, source_path, repository, generator_sha256, mutation_lock
            )
        elif args.action == "draft-request":
            locked_draft_request(
                args, root, source_path, repository, generator_sha256, mutation_lock
            )
        elif args.action == "install-request-review":
            locked_install_request_review(
                args, root, source_path, repository, generator_sha256, mutation_lock
            )
        else:
            locked_retain_generator(
                args, root, source_path, repository, generator_sha256, mutation_lock
            )
        mutation_lock.authenticate()
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (
        KeyError,
        OSError,
        OverflowError,
        TypeError,
        ValueError,
        subprocess.SubprocessError,
    ) as error:
        print(f"Stage I recost generator failed: {error}", file=sys.stderr)
        raise SystemExit(1)
