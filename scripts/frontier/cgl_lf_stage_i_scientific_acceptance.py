#!/usr/bin/env python3
"""Evaluate non-authorizing scientific acceptance evidence for CGL-LF Stage I.

This utility is deliberately independent of ``analyze_cgl_lf_paper.py``.  It
consumes retained histories and exact replay-verified scientific products,
performs the reviewed statistical and plasma-physics gates, and emits evidence
either to stdout or to an explicitly named candidate file.  It never mutates
canonical campaign state and no result from this utility is publication
authorization.
"""

from __future__ import annotations

import argparse
from array import array
from contextlib import contextmanager
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import random
import re
import stat
import struct
import sys
from typing import Iterable, Iterator


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
UTILITY_RELATIVE_PATH = Path(
    "scripts/frontier/cgl_lf_stage_i_scientific_acceptance.py"
)
DEFAULT_CRITERIA = REPOSITORY_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPOSITORY_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
CANONICAL_CAMPAIGN_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
CANONICAL_EXECUTION_EPOCH = "E03-forcing-policy"
CANONICAL_RUN_ROOT = (
    CANONICAL_CAMPAIGN_ROOT / "runs/mks24-stage-i" / CANONICAL_EXECUTION_EPOCH
)
CANONICAL_CONTROLLER = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-df/scripts/frontier/cgl_lf_stage_i.py"
)
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
CASE_ID_PATTERN = re.compile(r"R(?:0[2-9]|1[0-7])")
RANK_DIRECTORY_PATTERN = re.compile(r"rank_(\d{8})")
HISTORY_LABEL_PATTERN = re.compile(r"\[(\d+)\]=([^\s]+)")
EVIDENCE_DIGEST_METHOD = "sha256-canonical-json-without-evidence_digest-v1"
NON_AUTHORIZING_STATEMENT = (
    "This candidate is scientific assessment evidence only and does not "
    "authorize simulation execution, canonical publication, or Stage I release."
)
IDENTITY_LIMITATION = (
    "Reviewer identity and process independence are declared evidence, not "
    "cryptographically proven."
)
VALID_RESULTS = frozenset(("pass", "fail", "inconclusive", "blocked_out_of_scope"))
MAX_FINITE_FLOAT = sys.float_info.max
MAX_JSON_BYTES = 512 * 1024 * 1024
MAX_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
MAX_RESTART_STATE_BYTES = 32 * 1024 * 1024
RESTART_HEADER_FORMAT = "<ii9d19i19iddi"
RESTART_HEADER_SIZE = struct.calcsize(RESTART_HEADER_FORMAT)
LOGICAL_LOCATION_SIZE = struct.calcsize("<4i")
MESH_BLOCK_COST_SIZE = struct.calcsize("<f")
IO_WRAPPER_SIZE_FORMAT = "<Q"
IO_WRAPPER_SIZE_BYTES = struct.calcsize(IO_WRAPPER_SIZE_FORMAT)


class AcceptanceError(ValueError):
    """Raised when scientific evidence is malformed, unbound, or insufficient."""


def unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys instead of silently selecting one value."""

    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise AcceptanceError(f"duplicate JSON key: {key}")
        value[key] = item
    return value


def reject_json_constant(value: str) -> object:
    """Reject non-finite JSON constants."""

    raise AcceptanceError(f"invalid JSON numeric constant: {value}")


def canonical_json(value: object) -> bytes:
    """Return the compact canonical JSON representation used for evidence digests."""

    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def stable_json(value: object) -> bytes:
    """Return human-readable deterministic JSON."""

    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def sha256_bytes(value: bytes) -> str:
    """Return a lowercase SHA-256 digest."""

    return hashlib.sha256(value).hexdigest()


def require_sha256(value: object, label: str) -> str:
    """Require one lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise AcceptanceError(f"{label} is not a lowercase SHA-256 digest")
    return value


def require_revision(value: object, label: str) -> str:
    """Require one full lowercase Git revision."""

    if not isinstance(value, str) or REVISION_PATTERN.fullmatch(value) is None:
        raise AcceptanceError(f"{label} is not a lowercase full Git revision")
    return value


def require_case_id(value: object, label: str = "case_id") -> str:
    """Require one Stage I case identifier."""

    if not isinstance(value, str) or CASE_ID_PATTERN.fullmatch(value) is None:
        raise AcceptanceError(f"{label} is not a Stage I R02-R17 case identifier")
    return value


def require_finite(value: object, label: str) -> float:
    """Require one finite non-boolean number."""

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise AcceptanceError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise AcceptanceError(f"{label} must be finite")
    return result


def require_finite_number_or_text(value: object, label: str) -> float:
    """Require one finite number represented numerically or as retained text."""

    if isinstance(value, str):
        try:
            result = float(value)
        except ValueError as error:
            raise AcceptanceError(f"{label} must be numeric") from error
        if not math.isfinite(result):
            raise AcceptanceError(f"{label} must be finite")
        return result
    return require_finite(value, label)


def require_positive(value: object, label: str) -> float:
    """Require one finite positive number."""

    result = require_finite(value, label)
    if result <= 0.0:
        raise AcceptanceError(f"{label} must be positive")
    return result


def require_nonnegative(value: object, label: str) -> float:
    """Require one finite nonnegative number."""

    result = require_finite(value, label)
    if result < 0.0:
        raise AcceptanceError(f"{label} must be nonnegative")
    return result


def finite_ratio(numerator: float, denominator: float) -> float:
    """Return a finite nonnegative ratio, saturating exact zero denominators."""

    if numerator == 0.0:
        return 0.0
    if denominator == 0.0:
        return MAX_FINITE_FLOAT
    result = abs(numerator / denominator)
    return result if math.isfinite(result) else MAX_FINITE_FLOAT


def require_int(value: object, label: str, *, minimum: int = 0) -> int:
    """Require one exact integer."""

    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise AcceptanceError(f"{label} must be an integer >= {minimum}")
    return value


def require_dict(value: object, label: str) -> dict[str, object]:
    """Require one JSON object."""

    if not isinstance(value, dict):
        raise AcceptanceError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    """Require one JSON list."""

    if not isinstance(value, list):
        raise AcceptanceError(f"{label} must be a list")
    return value


def stable_profile(profile: os.stat_result) -> tuple[int, int, int, int, int, int]:
    """Return identity and mutation fields used for descriptor-bound reads."""

    return (
        profile.st_dev,
        profile.st_ino,
        profile.st_mode,
        profile.st_nlink,
        profile.st_size,
        profile.st_mtime_ns,
    )


@contextmanager
def open_stable_regular(path: Path, label: str) -> Iterator[int]:
    """Open one unchanged regular file without following its leaf symlink."""

    absolute = path.expanduser().absolute()
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(absolute, flags)
    except OSError as error:
        raise AcceptanceError(f"{label} cannot be opened safely: {absolute}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise AcceptanceError(f"{label} is not a regular file: {absolute}")
        if before.st_nlink != 1:
            raise AcceptanceError(f"{label} must have exactly one hard link: {absolute}")
        yield descriptor
        after = os.fstat(descriptor)
        if stable_profile(after) != stable_profile(before):
            raise AcceptanceError(f"{label} changed while it was read: {absolute}")
    finally:
        os.close(descriptor)


def descriptor_sha256(descriptor: int) -> str:
    """Hash one open descriptor without changing its offset."""

    digest = hashlib.sha256()
    offset = 0
    while True:
        block = os.pread(descriptor, 1024 * 1024, offset)
        if not block:
            break
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def read_stable_bytes(path: Path, label: str, *, maximum: int = MAX_JSON_BYTES) -> bytes:
    """Read one bounded unchanged regular file."""

    with open_stable_regular(path, label) as descriptor:
        profile = os.fstat(descriptor)
        if profile.st_size > maximum:
            raise AcceptanceError(f"{label} exceeds the reviewed size limit")
        payload = bytearray()
        offset = 0
        while len(payload) < profile.st_size:
            block = os.pread(descriptor, min(1024 * 1024, profile.st_size - offset), offset)
            if not block:
                raise AcceptanceError(f"{label} ended before its retained size")
            payload.extend(block)
            offset += len(block)
        if len(payload) != profile.st_size:
            raise AcceptanceError(f"{label} size changed while it was read")
        return bytes(payload)


def load_json(path: Path, label: str) -> tuple[dict[str, object], dict[str, object]]:
    """Load and bind one unambiguous JSON object."""

    payload = read_stable_bytes(path, label)
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_json_object,
            parse_constant=reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AcceptanceError(f"{label} is invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise AcceptanceError(f"{label} must be a JSON object: {path}")
    return value, {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "size_bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def load_json_value(path: Path, label: str) -> tuple[object, dict[str, object]]:
    """Load and bind one unambiguous JSON value."""

    payload = read_stable_bytes(path, label)
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_json_object,
            parse_constant=reject_json_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise AcceptanceError(f"{label} is invalid JSON: {path}") from error
    return value, {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "size_bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def regular_file_binding(path: Path, label: str) -> dict[str, object]:
    """Return exact content binding for one unchanged regular file."""

    with open_stable_regular(path, label) as descriptor:
        profile = os.fstat(descriptor)
        return {
            "path": str(path.expanduser().absolute().resolve(strict=True)),
            "size_bytes": profile.st_size,
            "sha256": descriptor_sha256(descriptor),
        }


def resolve_bound_path(path_value: object) -> Path:
    """Resolve an absolute or repository-relative source binding."""

    if not isinstance(path_value, str) or not path_value:
        raise AcceptanceError("source binding path must be a nonempty string")
    path = Path(path_value)
    return path if path.is_absolute() else REPOSITORY_ROOT / path


def verify_declared_binding(binding: object, label: str) -> dict[str, object]:
    """Verify one declared path/SHA/size binding against unchanged bytes."""

    value = require_dict(binding, label)
    path = resolve_bound_path(value.get("path"))
    expected_sha = require_sha256(value.get("sha256"), f"{label} sha256")
    observed = regular_file_binding(path, label)
    if observed["sha256"] != expected_sha:
        raise AcceptanceError(f"{label} SHA-256 differs from the declared binding")
    expected_size = value.get("size_bytes")
    if expected_size is not None and require_int(
        expected_size, f"{label} size_bytes"
    ) != observed["size_bytes"]:
        raise AcceptanceError(f"{label} size differs from the declared binding")
    return observed


def evidence_digest(value: dict[str, object]) -> str:
    """Digest evidence excluding its self-referential digest member."""

    body = dict(value)
    body.pop("evidence_digest", None)
    return sha256_bytes(canonical_json(body))


def seal_evidence(value: dict[str, object]) -> dict[str, object]:
    """Attach one deterministic self-digest to evidence."""

    if "evidence_digest" in value:
        raise AcceptanceError("evidence is already sealed")
    sealed = dict(value)
    sealed["evidence_digest"] = {
        "method": EVIDENCE_DIGEST_METHOD,
        "sha256": evidence_digest(sealed),
    }
    return sealed


def verify_evidence_digest(value: dict[str, object], label: str) -> None:
    """Verify one evidence self-digest."""

    digest = require_dict(value.get("evidence_digest"), f"{label} evidence_digest")
    if digest.get("method") != EVIDENCE_DIGEST_METHOD:
        raise AcceptanceError(f"{label} evidence digest method is unsupported")
    expected = require_sha256(digest.get("sha256"), f"{label} evidence digest sha256")
    if evidence_digest(value) != expected:
        raise AcceptanceError(f"{label} evidence self-digest differs")


def write_candidate(path: Path | None, value: dict[str, object]) -> None:
    """Write one immutable no-clobber candidate, or emit it to stdout."""

    payload = stable_json(value)
    if path is None:
        sys.stdout.buffer.write(payload)
        return
    absolute = path.expanduser().absolute()
    absolute.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    try:
        descriptor = os.open(absolute, flags, 0o444)
    except FileExistsError as error:
        raise AcceptanceError(f"candidate output already exists: {absolute}") from error
    try:
        os.fchmod(descriptor, 0o444)
        written = 0
        while written < len(payload):
            count = os.write(descriptor, payload[written:])
            if count <= 0:
                raise AcceptanceError(f"candidate output write failed: {absolute}")
            written += count
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    directory = os.open(absolute.parent, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def source_binding_records(criteria: dict[str, object]) -> list[tuple[str, object]]:
    """Return every exact source input binding declared by criteria."""

    sources = require_dict(criteria.get("source_bindings"), "source_bindings")
    records: list[tuple[str, object]] = []
    for key in (
        "acceptance_utility",
        "stage_i_manifest",
        "reference_archive_manifest",
        "verified_reference_archive_manifest",
        "current_source_authority_evidence",
        "current_source_authority_provenance_review",
        "current_source_authority_plasma_review",
        "current_source_authority_publication_audit",
        "source_archive_catalog",
        "qualification_approval",
    ):
        records.append((key, sources.get(key)))
    return records


def validate_criteria_payload(
    criteria: dict[str, object],
    criteria_binding: dict[str, object],
) -> dict[str, object]:
    """Validate preregistered policy and exact source bindings."""

    if criteria.get("schema_version") != 1:
        raise AcceptanceError("criteria schema_version must be 1")
    if criteria.get("record_type") != "stage-i-scientific-acceptance-criteria":
        raise AcceptanceError("criteria record_type is invalid")
    if criteria.get("authority") != "non-authorizing-preregistered-scientific-policy":
        raise AcceptanceError("criteria authority statement is invalid")
    if criteria.get("non_authorizing_statement") != NON_AUTHORIZING_STATEMENT:
        raise AcceptanceError("criteria non-authorizing statement differs")

    verified_sources: dict[str, object] = {}
    for label, record in source_binding_records(criteria):
        verified_sources[label] = verify_declared_binding(record, label)

    manifest_path = Path(
        str(require_dict(criteria["source_bindings"], "source_bindings")[
            "stage_i_manifest"
        ]["path"])
    )
    if not manifest_path.is_absolute():
        manifest_path = REPOSITORY_ROOT / manifest_path
    manifest, _ = load_json(manifest_path, "Stage I manifest")
    cases = require_list(manifest.get("cases"), "Stage I manifest cases")
    expected_cases = [f"R{number:02d}" for number in range(2, 18)]
    observed_cases = [
        require_case_id(require_dict(case, "Stage I manifest case").get("id"))
        for case in cases
    ]
    if observed_cases != expected_cases:
        raise AcceptanceError("criteria Stage I manifest does not enumerate R02-R17")
    required_cases = require_list(criteria.get("required_cases"), "required_cases")
    if required_cases != expected_cases:
        raise AcceptanceError("criteria required_cases must be exactly R02-R17")

    windows = require_dict(criteria.get("analysis_windows"), "analysis_windows")
    expected_windows = {
        "full": [4.0, 10.0],
        "early": [4.0, 7.0],
        "late": [7.0, 10.0],
    }
    if windows != expected_windows:
        raise AcceptanceError("criteria analysis windows differ from the preregistration")

    statistics = require_dict(criteria.get("statistics_policy"), "statistics_policy")
    if statistics.get("time_average") != "endpoint-clipped-trapezoidal":
        raise AcceptanceError("criteria time_average policy is invalid")
    if statistics.get("effective_samples") != "initial-positive-sequence":
        raise AcceptanceError("criteria effective-sample policy is invalid")
    if statistics.get("uncertainty") != "deterministic-time-moving-block-bootstrap":
        raise AcceptanceError("criteria uncertainty policy is invalid")
    if require_int(statistics.get("bootstrap_replicates"), "bootstrap_replicates", minimum=1) != 2000:
        raise AcceptanceError("criteria bootstrap_replicates must be 2000")
    if statistics.get("minimum_independent_time_blocks") != {
        "full": 3.0,
        "half": 1.5,
    }:
        raise AcceptanceError(
            "criteria minimum independent time-block counts differ from the preregistration"
        )
    if statistics.get("minimum_effective_samples") != {
        "full": 3.0,
        "half": 1.5,
    }:
        raise AcceptanceError(
            "criteria minimum effective-sample counts differ from the preregistration"
        )
    gap_policy = require_dict(statistics.get("gap_policy"), "statistics gap_policy")
    if gap_policy != {
        "expected_history_cadence": 0.02,
        "maximum_gap_expected_cadence_multiplier": 2.5,
        "maximum_gap_forcing_tcorr_fraction": 0.25,
        "median_cadence_upper_relative_tolerance": 0.25,
        "rule": (
            "median cadence <= expected cadence * (1 + tolerance) and maximum gap "
            "<= max(expected cadence * multiplier, authenticated forcing_tcorr * fraction)"
        ),
    }:
        raise AcceptanceError("criteria gap policy differs from the preregistration")

    canonical = require_dict(
        criteria.get("canonical_campaign_policy"), "canonical_campaign_policy"
    )
    expected_canonical = {
        "campaign_root": str(CANONICAL_CAMPAIGN_ROOT),
        "controller_path": str(CANONICAL_CONTROLLER),
        "execution_epoch": CANONICAL_EXECUTION_EPOCH,
        "ledger_path": str(
            CANONICAL_CAMPAIGN_ROOT
            / "accounting/mks24_stage_i_E03_forcing_policy_node_hours.csv"
        ),
        "reservation_store_path": str(
            CANONICAL_CAMPAIGN_ROOT
            / "accounting/mks24_stage_i_E03_forcing_policy_reservations.json"
        ),
        "run_root": str(CANONICAL_RUN_ROOT),
        "source_archive_root": str(CANONICAL_CAMPAIGN_ROOT / "source-archives"),
    }
    if canonical != expected_canonical:
        raise AcceptanceError("criteria canonical campaign policy differs")

    products = require_dict(
        criteria.get("scientific_products_policy"), "scientific_products_policy"
    )
    generator = require_dict(
        products.get("reviewed_generator_binding"), "reviewed_generator_binding"
    )
    if generator != {
        "path": None,
        "review_path": None,
        "review_sha256": None,
        "review_status": "changes_required",
        "sha256": None,
        "status": "unavailable_pending_companion_generator",
    }:
        raise AcceptanceError(
            "criteria reviewed scientific-products generator binding differs"
        )
    if products.get("contract") != {
        "contract_schema_version": 2,
        "deterministic_replay_verification": (
            "exact-semantic-replay-from-bound-canonical-case-inputs"
        ),
        "hand_authored_or_unreplayed_products": "inconclusive",
        "output_schema_version": 1,
    }:
        raise AcceptanceError("criteria scientific-products replay contract differs")

    manifest_panels = require_dict(manifest.get("panel_status"), "manifest panel_status")
    configured = require_list(manifest_panels.get("panels"), "manifest panels")
    admitted = {
        str(require_dict(panel, "manifest panel").get("id")): panel
        for panel in configured
        if require_dict(panel, "manifest panel").get("disposition") == "comparison"
    }
    blocked = {
        str(require_dict(panel, "manifest panel").get("id"))
        for panel in configured
        if require_dict(panel, "manifest panel").get("disposition")
        in ("blocked_reference", "external_model")
    }
    criteria_panels = require_list(criteria.get("comparison_panels"), "comparison_panels")
    criteria_panel_ids = [
        str(require_dict(panel, "criteria panel").get("id"))
        for panel in criteria_panels
    ]
    if set(criteria_panel_ids) != set(admitted) or len(criteria_panel_ids) != len(
        set(criteria_panel_ids)
    ):
        raise AcceptanceError("criteria comparison panels differ from admitted panels")
    exclusions = require_dict(criteria.get("allowed_exclusions"), "allowed_exclusions")
    declared_blocked = set(
        str(value)
        for value in require_list(
            exclusions.get("blocked_or_external_panels"),
            "allowed_exclusions blocked_or_external_panels",
        )
    )
    if declared_blocked != blocked:
        raise AcceptanceError("criteria allowed exclusions differ from manifest dispositions")

    utility = require_dict(
        require_dict(criteria.get("source_bindings"), "source_bindings").get(
            "acceptance_utility"
        ),
        "acceptance utility binding",
    )
    if utility.get("path") != str(UTILITY_RELATIVE_PATH):
        raise AcceptanceError("criteria acceptance utility path differs")
    if verified_sources["acceptance_utility"]["sha256"] != utility.get("sha256"):
        raise AcceptanceError("criteria acceptance utility digest differs")

    ct_policy = require_dict(criteria.get("ct_divb_policy"), "ct_divb_policy")
    required_state_times = [
        require_finite(value, "required CT state time")
        for value in require_list(
            ct_policy.get("required_state_times"), "required CT state times"
        )
    ]
    if required_state_times != [9.0, 10.0]:
        raise AcceptanceError("criteria CT required state times must be exactly 9 and 10")

    family = require_dict(criteria.get("family_gates"), "family_gates")
    lf_strength = require_dict(family.get("lf_strength"), "lf_strength")
    if lf_strength.get("cases") != ["R12", "R02", "R06", "R13"]:
        raise AcceptanceError("criteria LF-strength cases differ from the preregistration")

    convergence = require_dict(
        criteria.get("resolution_convergence"), "resolution_convergence"
    )
    if convergence.get("cases") != ["R16", "R02", "R17"]:
        raise AcceptanceError("criteria convergence cases differ from the preregistration")
    if convergence.get("common_k_perp_over_pi") != [4.0, 24.0]:
        raise AcceptanceError("criteria convergence interval differs from the preregistration")
    if convergence.get("alignment_shells") != [4, 6, 8, 12, 16, 24]:
        raise AcceptanceError("criteria alignment shells differ from the preregistration")
    if convergence.get("k_coordinate") != "physical_k_equals_k_perp_over_pi_times_exact_pi":
        raise AcceptanceError("criteria convergence coordinate convention differs")
    if require_positive(
        convergence.get("resolved_disagreement_sigma_gt"),
        "resolved convergence disagreement sigma",
    ) != 2.0:
        raise AcceptanceError("criteria convergence resolution threshold differs")

    return {
        "criteria": criteria,
        "criteria_binding": criteria_binding,
        "manifest": manifest,
        "verified_sources": verified_sources,
    }


def validate_criteria_review(
    review: dict[str, object],
    review_binding: dict[str, object],
    criteria_binding: dict[str, object],
    utility_binding: dict[str, object],
) -> dict[str, object]:
    """Validate a pending or independently approved criteria review."""

    if review.get("schema_version") != 1:
        raise AcceptanceError("criteria review schema_version must be 1")
    if review.get("record_type") != "stage-i-scientific-acceptance-criteria-review":
        raise AcceptanceError("criteria review record_type is invalid")
    if review.get("identity_limitation") != IDENTITY_LIMITATION:
        raise AcceptanceError("criteria review identity limitation differs")
    if review.get("non_authorizing_statement") != NON_AUTHORIZING_STATEMENT:
        raise AcceptanceError("criteria review non-authorizing statement differs")
    declared = require_dict(review.get("criteria"), "criteria review binding")
    if (
        require_sha256(declared.get("sha256"), "criteria review criteria sha256")
        != criteria_binding["sha256"]
    ):
        raise AcceptanceError("criteria review does not bind exact criteria bytes")
    if resolve_bound_path(declared.get("path")).resolve(strict=True) != Path(
        str(criteria_binding["path"])
    ):
        raise AcceptanceError("criteria review binds a different criteria path")
    declared_utility = require_dict(
        review.get("acceptance_utility"), "criteria review acceptance utility"
    )
    if (
        resolve_bound_path(declared_utility.get("path")).resolve(strict=True)
        != Path(str(utility_binding["path"]))
        or require_sha256(
            declared_utility.get("sha256"),
            "criteria review acceptance utility sha256",
        )
        != utility_binding["sha256"]
    ):
        raise AcceptanceError("criteria review does not bind the final acceptance utility")
    status_value = review.get("review_status")
    reviewers = require_list(review.get("reviews"), "criteria review reviews")
    retained_ct = require_dict(
        review.get("retained_ct_observation"), "criteria review retained CT observation"
    )
    if retained_ct != {
        "approval_identity_available": False,
        "case_id": "R02",
        "maximum_normalized_ct_divb": 1.743292792971505e-13,
        "meshblocks_sampled": 216,
        "meshblocks_total": 216,
        "rank_count": 8,
        "scope": "terminal retained state t=10 only; not plasma/statistical criteria approval",
        "state_time": 10.0,
        "status": "retained_observation_not_criteria_approval",
    }:
        raise AcceptanceError("criteria review retained CT observation differs")
    if status_value in ("pending_independent_review", "changes_required"):
        if reviewers:
            raise AcceptanceError("unapproved criteria review must not claim reviewers")
        approved = False
    elif status_value == "approved":
        roles: set[str] = set()
        identities: set[str] = set()
        for item in reviewers:
            record = require_dict(item, "criteria reviewer")
            role = record.get("role")
            identity = record.get("reviewer_id")
            decision = record.get("decision")
            if (
                role not in ("plasma_physics", "statistical_methodology")
                or not isinstance(identity, str)
                or not identity
                or decision != "approved"
                or record.get("independent_of_implementation") is not True
            ):
                raise AcceptanceError("approved criteria reviewer is malformed")
            roles.add(str(role))
            identities.add(identity)
        if roles != {"plasma_physics", "statistical_methodology"} or len(identities) != 2:
            raise AcceptanceError("criteria approval lacks two distinct required reviewers")
        approved = True
    else:
        raise AcceptanceError("criteria review status is invalid")
    return {
        "review": review,
        "review_binding": review_binding,
        "approved": approved,
        "review_status": status_value,
    }


def load_validated_policy(criteria_path: Path, review_path: Path) -> dict[str, object]:
    """Load criteria and its exact review binding."""

    criteria, criteria_binding = load_json(criteria_path, "scientific criteria")
    policy = validate_criteria_payload(criteria, criteria_binding)
    review, review_binding = load_json(review_path, "scientific criteria review")
    policy.update(validate_criteria_review(
        review,
        review_binding,
        criteria_binding,
        require_dict(policy["verified_sources"]["acceptance_utility"], "utility binding"),
    ))
    return policy


def parse_history_bytes(payload: bytes, label: str) -> dict[str, list[float]]:
    """Parse one complete finite Athena history with stable indexed labels."""

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise AcceptanceError(f"{label} is not UTF-8") from error
    labels: list[str] | None = None
    rows: list[list[float]] = []
    for line in text.splitlines():
        if line.startswith("#"):
            found = HISTORY_LABEL_PATTERN.findall(line)
            if found:
                indexed = sorted((int(index), name) for index, name in found)
                indices = [index for index, _ in indexed]
                names = [name for _, name in indexed]
                if (
                    not indices
                    or indices[0] not in (0, 1)
                    or indices != list(range(indices[0], indices[0] + len(indices)))
                    or len(names) != len(set(names))
                ):
                    raise AcceptanceError(f"{label} labels are not unique and contiguous")
                if labels is not None and labels != names:
                    raise AcceptanceError(f"{label} contains conflicting label headers")
                labels = names
            continue
        if not line.strip():
            continue
        try:
            row = [float(value) for value in line.split()]
        except ValueError as error:
            raise AcceptanceError(f"{label} contains a nonnumeric row") from error
        if not all(math.isfinite(value) for value in row):
            raise AcceptanceError(f"{label} contains a non-finite row")
        rows.append(row)
    if labels is None or len(rows) < 2:
        raise AcceptanceError(f"{label} lacks labels or sufficient rows")
    if any(len(row) != len(labels) for row in rows):
        raise AcceptanceError(f"{label} row width differs from its labels")
    result = {
        name: [row[index] for row in rows] for index, name in enumerate(labels)
    }
    times = result.get("time")
    if times is None or any(right <= left for left, right in zip(times, times[1:])):
        raise AcceptanceError(f"{label} time is not strictly increasing")
    return result


def load_history(path: Path, label: str) -> tuple[dict[str, list[float]], dict[str, object]]:
    """Load and bind one retained Athena history."""

    payload = read_stable_bytes(path, label)
    return parse_history_bytes(payload, label), {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "size_bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def interpolate_at(times: list[float], values: list[float], target: float) -> float:
    """Linearly interpolate one finite ordered series."""

    if target < times[0] or target > times[-1]:
        raise AcceptanceError("analysis window is outside retained time coverage")
    for index, time in enumerate(times):
        if time == target:
            return values[index]
        if time > target:
            left_time = times[index - 1]
            fraction = (target - left_time) / (time - left_time)
            return values[index - 1] + fraction * (values[index] - values[index - 1])
    return values[-1]


def clipped_series(
    times: list[float], values: list[float], start: float, end: float
) -> tuple[list[float], list[float]]:
    """Clip one series to exact endpoints using linear interpolation."""

    if len(times) != len(values) or len(times) < 2 or start >= end:
        raise AcceptanceError("invalid series or analysis window")
    if any(right <= left for left, right in zip(times, times[1:])):
        raise AcceptanceError("series time is not strictly increasing")
    clipped_times = [start]
    clipped_values = [interpolate_at(times, values, start)]
    for time, value in zip(times, values):
        if start < time < end:
            clipped_times.append(time)
            clipped_values.append(value)
    clipped_times.append(end)
    clipped_values.append(interpolate_at(times, values, end))
    return clipped_times, clipped_values


def trapezoidal_weights(times: list[float]) -> list[float]:
    """Return endpoint-clipped trapezoidal integration weights."""

    if len(times) < 2:
        raise AcceptanceError("trapezoidal statistics require at least two samples")
    weights = [0.0] * len(times)
    for index, (left, right) in enumerate(zip(times, times[1:])):
        width = right - left
        if not math.isfinite(width) or width <= 0.0:
            raise AcceptanceError("trapezoidal time grid is not strictly increasing")
        weights[index] += 0.5 * width
        weights[index + 1] += 0.5 * width
    return weights


def weighted_mean(values: list[float], weights: list[float]) -> float:
    """Return one finite weighted mean."""

    total = sum(weights)
    if total <= 0.0 or len(values) != len(weights):
        raise AcceptanceError("weighted mean inputs are invalid")
    result = sum(weight * value for weight, value in zip(weights, values)) / total
    if not math.isfinite(result):
        raise AcceptanceError("weighted mean is non-finite")
    return result


def effective_sample_size(values: list[float]) -> tuple[float, float]:
    """Estimate N_eff with Geyer's initial-positive-sequence autocorrelation rule."""

    count = len(values)
    if count < 2:
        return 1.0, 1.0
    mean = sum(values) / count
    centered = [value - mean for value in values]
    variance = sum(value * value for value in centered) / count
    if variance == 0.0:
        return float(count), 1.0

    correlations = [1.0]
    for lag in range(1, count):
        covariance = sum(
            centered[index] * centered[index + lag]
            for index in range(count - lag)
        ) / count
        correlations.append(covariance / variance)
    positive_sum = 0.0
    lag = 1
    while lag < count:
        pair = correlations[lag]
        if lag + 1 < count:
            pair += correlations[lag + 1]
        if pair <= 0.0:
            break
        positive_sum += pair
        lag += 2
    tau = max(1.0, 1.0 + 2.0 * positive_sum)
    return max(1.0, min(float(count), count / tau)), tau


def quantile(values: list[float], probability: float) -> float:
    """Return one linearly interpolated finite quantile."""

    if not values or probability < 0.0 or probability > 1.0:
        raise AcceptanceError("quantile inputs are invalid")
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def median(values: list[float]) -> float:
    """Return the finite median of a nonempty sequence."""

    if not values or any(not math.isfinite(value) for value in values):
        raise AcceptanceError("median inputs are empty or non-finite")
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return 0.5 * (ordered[middle - 1] + ordered[middle])


def bootstrap_means(
    times: list[float],
    values: list[float],
    replicates: int,
    seed_text: str,
    block_duration: float,
) -> list[float]:
    """Return deterministic physical-time moving-block bootstrap means."""

    if len(times) != len(values) or len(times) < 2:
        raise AcceptanceError("bootstrap time/value series is invalid")
    replicates = require_int(replicates, "bootstrap replicates", minimum=1)
    intervals = [
        (right_time - left_time, left_value, right_value)
        for left_time, right_time, left_value, right_value
        in zip(times, times[1:], values, values[1:])
    ]
    if any(
        not math.isfinite(width) or width <= 0.0
        or not math.isfinite(left) or not math.isfinite(right)
        for width, left, right in intervals
    ):
        raise AcceptanceError("bootstrap intervals are non-finite or unordered")
    duration = times[-1] - times[0]
    reviewed_block_duration = min(
        duration, require_positive(block_duration, "bootstrap block duration")
    )
    offsets = [time - times[0] for time in times]

    def circular_integral(start_offset: float, length: float) -> float:
        """Integrate a physical-time block through the circular piecewise-linear path."""

        position = start_offset % duration
        remaining = length
        integral = 0.0
        while remaining > 0.0:
            interval_index = len(intervals) - 1
            for index, right in enumerate(offsets[1:]):
                if position < right:
                    interval_index = index
                    break
            left_offset = offsets[interval_index]
            width, left_value, right_value = intervals[interval_index]
            within = position - left_offset
            retained = min(remaining, width - within)
            if retained <= 0.0:
                position = 0.0
                continue
            start_fraction = within / width
            end_fraction = (within + retained) / width
            retained_left = left_value + start_fraction * (right_value - left_value)
            retained_right = left_value + end_fraction * (right_value - left_value)
            integral += retained * 0.5 * (retained_left + retained_right)
            remaining -= retained
            position += retained
            if position >= duration:
                position = 0.0
        return integral

    seed = int.from_bytes(hashlib.sha256(seed_text.encode("utf-8")).digest()[:16], "big")
    generator = random.Random(seed)
    samples: list[float] = []
    for _ in range(replicates):
        sampled_duration = 0.0
        sampled_integral = 0.0
        while sampled_duration < duration:
            retained = min(reviewed_block_duration, duration - sampled_duration)
            sampled_integral += circular_integral(
                generator.random() * duration, retained
            )
            sampled_duration += retained
        mean = sampled_integral / duration
        if not math.isfinite(mean):
            raise AcceptanceError("bootstrap mean is non-finite")
        samples.append(mean)
    return samples


def window_statistics(
    times: list[float],
    values: list[float],
    start: float,
    end: float,
    *,
    replicates: int,
    seed_text: str,
    minimum_block_duration: float = 0.0,
    expected_cadence: float | None = None,
    maximum_gap_expected_cadence_multiplier: float = 2.5,
    maximum_gap_minimum_block_duration_fraction: float = 0.25,
    median_cadence_upper_relative_tolerance: float = 0.25,
) -> dict[str, object]:
    """Compute time-weighted mean, N_eff, and deterministic bootstrap uncertainty."""

    clipped_times, clipped_values = clipped_series(times, values, start, end)
    weights = trapezoidal_weights(clipped_times)
    mean = weighted_mean(clipped_values, weights)
    variance = weighted_mean([(value - mean) ** 2 for value in clipped_values], weights)
    raw_neff, tau = effective_sample_size(clipped_values)
    gaps = [
        right - left for left, right in zip(clipped_times, clipped_times[1:])
    ]
    cadence = median(gaps)
    maximum_gap = max(gaps)
    estimated_autocorrelation_duration = tau * cadence
    minimum_duration = require_nonnegative(
        minimum_block_duration, "minimum block duration"
    )
    required_block_duration = max(
        cadence,
        estimated_autocorrelation_duration,
        minimum_duration,
    )
    window_duration = clipped_times[-1] - clipped_times[0]
    physical_independence_cap = window_duration / required_block_duration
    neff = min(raw_neff, physical_independence_cap)
    reviewed_expected_cadence = (
        cadence
        if expected_cadence is None
        else require_positive(expected_cadence, "expected cadence")
    )
    multiplier = require_positive(
        maximum_gap_expected_cadence_multiplier,
        "maximum gap expected-cadence multiplier",
    )
    tcorr_fraction = require_nonnegative(
        maximum_gap_minimum_block_duration_fraction,
        "maximum gap minimum-block-duration fraction",
    )
    cadence_tolerance = require_nonnegative(
        median_cadence_upper_relative_tolerance,
        "median cadence upper relative tolerance",
    )
    maximum_allowed_gap = max(
        multiplier * reviewed_expected_cadence,
        tcorr_fraction * minimum_duration,
    )
    cadence_upper_limit = reviewed_expected_cadence * (1.0 + cadence_tolerance)
    gap_adequacy = (
        "pass"
        if cadence <= cadence_upper_limit and maximum_gap <= maximum_allowed_gap
        else "inconclusive"
    )
    bootstrap_block_duration = min(window_duration, required_block_duration)
    bootstraps = bootstrap_means(
        clipped_times,
        clipped_values,
        replicates,
        seed_text,
        bootstrap_block_duration,
    )
    bootstrap_mean = sum(bootstraps) / len(bootstraps)
    bootstrap_variance = sum(
        (value - bootstrap_mean) ** 2 for value in bootstraps
    ) / max(1, len(bootstraps) - 1)
    return {
        "window": {"start": start, "end": end},
        "sample_count": len(clipped_values),
        "mean": mean,
        "standard_deviation": math.sqrt(max(0.0, variance)),
        "effective_sample_count": neff,
        "raw_initial_positive_sequence_effective_sample_count": raw_neff,
        "physical_independence_cap": physical_independence_cap,
        "independent_time_block_count": physical_independence_cap,
        "gap_adequacy": gap_adequacy,
        "median_cadence": cadence,
        "maximum_gap": maximum_gap,
        "maximum_allowed_gap": maximum_allowed_gap,
        "integrated_autocorrelation_time": tau,
        "standard_error": math.sqrt(max(0.0, bootstrap_variance)),
        "confidence_interval_95": [
            quantile(bootstraps, 0.025),
            quantile(bootstraps, 0.975),
        ],
        "method": {
            "time_average": "endpoint-clipped-trapezoidal",
            "effective_samples": "initial-positive-sequence",
            "uncertainty": "deterministic-circular-time-moving-block-bootstrap",
            "bootstrap_replicates": replicates,
            "estimated_autocorrelation_duration": estimated_autocorrelation_duration,
            "minimum_block_duration": minimum_block_duration,
            "required_block_duration": required_block_duration,
            "bootstrap_block_duration": bootstrap_block_duration,
            "expected_cadence": reviewed_expected_cadence,
            "median_cadence_upper_limit": cadence_upper_limit,
            "maximum_gap_rule": (
                "max(expected cadence * multiplier, minimum block duration * fraction)"
            ),
            "maximum_gap_expected_cadence_multiplier": multiplier,
            "maximum_gap_minimum_block_duration_fraction": tcorr_fraction,
        },
    }


def combine_occupancy_series(user: dict[str, list[float]]) -> list[float]:
    """Return mirror-plus-firehose volume occupancy."""

    for column in ("mirror_vol", "fire_vol"):
        if column not in user:
            raise AcceptanceError(f"user history lacks required column: {column}")
    return [
        mirror + fire for mirror, fire in zip(user["mirror_vol"], user["fire_vol"])
    ]


def stationarity_result(
    full: dict[str, object],
    early: dict[str, object],
    late: dict[str, object],
    kind: str,
    policy: dict[str, object],
) -> dict[str, object]:
    """Evaluate preregistered early/late stationarity."""

    difference = abs(float(early["mean"]) - float(late["mean"]))
    combined_se = math.hypot(
        float(early["standard_error"]), float(late["standard_error"])
    )
    z_score = finite_ratio(difference, combined_se)
    relative = finite_ratio(difference, abs(float(full["mean"])))
    z_limit = require_nonnegative(policy.get("z_lte"), "stationarity z_lte")
    if kind == "occupancy":
        absolute_limit = require_nonnegative(
            policy.get("occupancy_absolute_change_lte"),
            "occupancy stationarity absolute limit",
        )
        passed = difference <= absolute_limit and z_score <= z_limit
        limits = {
            "absolute_change_lte": absolute_limit,
            "z_lte": z_limit,
        }
    elif kind == "forcing_power":
        relative_limit = require_nonnegative(
            policy.get("forcing_power_relative_change_lte"),
            "forcing stationarity relative limit",
        )
        passed = relative <= relative_limit
        limits = {"relative_change_lte": relative_limit}
    else:
        relative_limit = require_nonnegative(
            policy.get("scalar_relative_change_lte"),
            "scalar stationarity relative limit",
        )
        passed = relative <= relative_limit and z_score <= z_limit
        limits = {"relative_change_lte": relative_limit, "z_lte": z_limit}
    return {
        "result": "pass" if passed else "fail",
        "absolute_change": difference,
        "relative_change": relative,
        "z_score": z_score,
        "limits": limits,
    }


def metric_statistics(
    history: dict[str, list[float]],
    values: list[float],
    metric: str,
    policy: dict[str, object],
    *,
    kind: str,
    minimum_block_duration: float = 0.0,
) -> dict[str, object]:
    """Compute full/half statistics and one stationarity result."""

    times = history["time"]
    windows = require_dict(policy["criteria"].get("analysis_windows"), "analysis_windows")
    statistics = require_dict(
        policy["criteria"].get("statistics_policy"), "statistics_policy"
    )
    replicates = require_int(statistics.get("bootstrap_replicates"), "bootstrap_replicates", minimum=1)
    gap_policy = require_dict(statistics.get("gap_policy"), "statistics gap policy")
    result: dict[str, object] = {}
    for name in ("full", "early", "late"):
        window = require_list(windows[name], f"{name} window")
        result[name] = window_statistics(
            times,
            values,
            require_finite(window[0], f"{name} start"),
            require_finite(window[1], f"{name} end"),
            replicates=replicates,
            seed_text=f"{metric}:{name}:{policy['criteria_binding']['sha256']}",
            minimum_block_duration=minimum_block_duration,
            expected_cadence=require_positive(
                gap_policy.get("expected_history_cadence"), "expected history cadence"
            ),
            maximum_gap_expected_cadence_multiplier=require_positive(
                gap_policy.get("maximum_gap_expected_cadence_multiplier"),
                "maximum gap expected-cadence multiplier",
            ),
            maximum_gap_minimum_block_duration_fraction=require_nonnegative(
                gap_policy.get("maximum_gap_forcing_tcorr_fraction"),
                "maximum gap forcing-tcorr fraction",
            ),
            median_cadence_upper_relative_tolerance=require_nonnegative(
                gap_policy.get("median_cadence_upper_relative_tolerance"),
                "median cadence upper relative tolerance",
            ),
        )
    minimum = require_dict(statistics.get("minimum_effective_samples"), "minimum_effective_samples")
    minimum_blocks = require_dict(
        statistics.get("minimum_independent_time_blocks"),
        "minimum_independent_time_blocks",
    )
    full_minimum = require_positive(minimum.get("full"), "minimum full effective samples")
    half_minimum = require_positive(minimum.get("half"), "minimum half effective samples")
    full_block_minimum = require_positive(
        minimum_blocks.get("full"), "minimum full independent time blocks"
    )
    half_block_minimum = require_positive(
        minimum_blocks.get("half"), "minimum half independent time blocks"
    )
    adequacy = all(
        float(require_dict(result[name], f"{name} stats")["effective_sample_count"])
        >= sample_required
        and require_dict(result[name], f"{name} stats").get("gap_adequacy") == "pass"
        and float(
            require_dict(result[name], f"{name} stats")["independent_time_block_count"]
        )
        >= block_required
        for name, sample_required, block_required in (
            ("full", full_minimum, full_block_minimum),
            ("early", half_minimum, half_block_minimum),
            ("late", half_minimum, half_block_minimum),
        )
    )
    result["sampling_adequacy"] = "pass" if adequacy else "inconclusive"
    result["stationarity"] = stationarity_result(
        require_dict(result["full"], "full stats"),
        require_dict(result["early"], "early stats"),
        require_dict(result["late"], "late stats"),
        kind,
        require_dict(statistics.get("stationarity"), "stationarity policy"),
    )
    return result


def gate(
    name: str,
    result: str,
    *,
    reason: str,
    observations: object = None,
    limits: object = None,
) -> dict[str, object]:
    """Return one normalized acceptance gate."""

    if result not in VALID_RESULTS:
        raise AcceptanceError(f"invalid gate result for {name}: {result}")
    return {
        "name": name,
        "result": result,
        "reason": reason,
        "observations": observations,
        "limits": limits,
    }


def history_delta(
    history: dict[str, list[float]], column: str, start: float, end: float
) -> float:
    """Return an endpoint-interpolated history increment."""

    if column not in history:
        raise AcceptanceError(f"history lacks required column: {column}")
    return interpolate_at(history["time"], history[column], end) - interpolate_at(
        history["time"], history[column], start
    )


def exact_zero_window(
    history: dict[str, list[float]], column: str, start: float, end: float
) -> bool:
    """Return whether every clipped value is exactly zero."""

    if column not in history:
        raise AcceptanceError(f"history lacks required column: {column}")
    _, values = clipped_series(history["time"], history[column], start, end)
    return all(value == 0.0 for value in values)


def case_name(policy: dict[str, object], case_id: str) -> str:
    """Return one case name from the bound Stage I manifest."""

    for item in require_list(policy["manifest"].get("cases"), "manifest cases"):
        record = require_dict(item, "manifest case")
        if record.get("id") == case_id:
            name = record.get("name")
            if not isinstance(name, str) or not name:
                raise AcceptanceError(f"manifest case {case_id} lacks a name")
            return name
    raise AcceptanceError(f"manifest lacks case {case_id}")


def case_manifest_record(policy: dict[str, object], case_id: str) -> dict[str, object]:
    """Return the exact selected case record from the bound Stage I manifest."""

    for item in require_list(policy["manifest"].get("cases"), "manifest cases"):
        record = require_dict(item, "manifest case")
        if record.get("id") == case_id:
            return record
    raise AcceptanceError(f"manifest lacks case {case_id}")


def binding_identity(binding: object, label: str) -> tuple[Path, str]:
    """Return one resolved path and exact digest from a binding."""

    value = require_dict(binding, label)
    return (
        resolve_bound_path(value.get("path")).resolve(strict=True),
        require_sha256(value.get("sha256"), f"{label} sha256"),
    )


def canonical_bundle_path(case_id: str) -> Path:
    """Return the sole canonical whole-case bundle path for one Stage I case."""

    return CANONICAL_RUN_ROOT / "bundles" / case_id / "manifest.json"


def path_is_exact(path: Path, expected: Path) -> bool:
    """Return whether two existing pathnames resolve to the same exact location."""

    try:
        return (
            path.expanduser().absolute().resolve(strict=True)
            == expected.expanduser().absolute().resolve(strict=True)
        )
    except FileNotFoundError:
        return False


def load_csv_records(
    path: Path, label: str
) -> tuple[list[dict[str, str]], dict[str, object]]:
    """Load and bind one finite unambiguous CSV table."""

    payload = read_stable_bytes(path, label)
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise AcceptanceError(f"{label} is not UTF-8") from error
    reader = csv.DictReader(io.StringIO(text))
    if (
        not reader.fieldnames
        or len(reader.fieldnames) != len(set(reader.fieldnames))
        or any(None in row for row in reader)
    ):
        raise AcceptanceError(f"{label} schema is malformed")
    reader = csv.DictReader(io.StringIO(text))
    records = [dict(row) for row in reader]
    if not records:
        raise AcceptanceError(f"{label} is empty")
    return records, {
        "path": str(path.expanduser().absolute().resolve(strict=True)),
        "size_bytes": len(payload),
        "sha256": sha256_bytes(payload),
    }


def source_archive_catalog(policy: dict[str, object]) -> dict[str, str]:
    """Return the exact current published source-archive checksum ledger."""

    binding = require_dict(
        policy["verified_sources"].get("source_archive_catalog"),
        "source archive catalog binding",
    )
    payload = read_stable_bytes(Path(str(binding["path"])), "source archive catalog")
    records: dict[str, str] = {}
    try:
        lines = payload.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise AcceptanceError("source archive checksum catalog is not UTF-8") from error
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64})  ([^/]+)", line)
        if match is None or match.group(2) in records:
            raise AcceptanceError("source archive checksum catalog is malformed")
        records[match.group(2)] = match.group(1)
    if not records:
        raise AcceptanceError("source archive checksum catalog is empty")
    return records


def authenticate_canonical_accounting(
    policy: dict[str, object],
    case_id: str,
    case_name_value: str,
    segment_paths: list[Path],
    segment_records: list[tuple[dict[str, object], dict[str, object]]],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Authenticate canonical controller/source/accounting lineage for segments."""

    canonical = require_dict(
        policy["criteria"].get("canonical_campaign_policy"),
        "canonical campaign policy",
    )
    ledger_path = Path(str(canonical["ledger_path"]))
    reservation_path = Path(str(canonical["reservation_store_path"]))
    ledger, ledger_binding = load_csv_records(ledger_path, "canonical Stage I ledger")
    reservations_value, reservation_binding = load_json_value(
        reservation_path, "canonical Stage I reservation store"
    )
    reservations = [
        require_dict(value, "canonical Stage I reservation")
        for value in require_list(reservations_value, "canonical Stage I reservations")
    ]
    ledger_by_job = {
        str(row.get("job_id")): row for row in ledger if row.get("job_id")
    }
    if len(ledger_by_job) != len(ledger):
        raise AcceptanceError("canonical Stage I ledger job identities are duplicated")
    reservation_by_job = {
        str(record.get("job_id")): record
        for record in reservations
        if record.get("job_id") is not None
    }
    if len(reservation_by_job) != sum(
        record.get("job_id") is not None for record in reservations
    ):
        raise AcceptanceError("canonical Stage I reservation job identities are duplicated")

    controller_lineage: list[dict[str, object]] = []
    retained_bindings: list[dict[str, object]] = [ledger_binding, reservation_binding]
    source_bundle_bindings: dict[Path, dict[str, object]] = {}
    qualification = require_dict(
        policy["verified_sources"].get("qualification_approval"),
        "qualification approval binding",
    )
    archive_catalog = source_archive_catalog(policy)
    for path, (segment, segment_binding) in zip(segment_paths, segment_records):
        accounting = require_dict(segment.get("accounting"), "canonical segment accounting")
        command = require_dict(segment.get("command"), "canonical segment command")
        segment_name = accounting.get("segment")
        if not isinstance(segment_name, str) or not segment_name:
            raise AcceptanceError("canonical segment lacks a normalized segment name")
        expected_path = (
            CANONICAL_RUN_ROOT
            / case_id
            / segment_name
            / "manifest"
            / "prepared_run.json"
        )
        if not path_is_exact(path, expected_path):
            raise AcceptanceError(
                "accepted segment manifest is outside its exact canonical campaign path"
            )
        job_id = accounting.get("job_id")
        if not isinstance(job_id, str) or not job_id:
            raise AcceptanceError("accepted canonical segment lacks a job identity")
        reservation = reservation_by_job.get(job_id)
        ledger_row = ledger_by_job.get(job_id)
        if reservation is None or ledger_row is None:
            raise AcceptanceError(
                "accepted canonical segment is absent from reservation store or ledger"
            )
        expected_common = {
            "case_id": case_id,
            "case_name": case_name_value,
            "execution_epoch": CANONICAL_EXECUTION_EPOCH,
            "job_id": job_id,
            "segment": segment_name,
        }
        if (
            accounting.get("result") != "accepted"
            or accounting.get("state") != "COMPLETED"
            or any(accounting.get(key) != value for key, value in expected_common.items())
            or reservation.get("state") != "recorded"
            or reservation.get("result") != "accepted"
            or any(reservation.get(key) != value for key, value in expected_common.items())
            or resolve_bound_path(reservation.get("manifest")).resolve(strict=True) != path
            or ledger_row.get("state") != "COMPLETED"
            or ledger_row.get("result") != "accepted"
            or any(ledger_row.get(key) != value for key, value in expected_common.items())
            or ledger_row.get("executable_sha256") != accounting.get("executable_sha256")
            or ledger_row.get("input_file") != accounting.get("input_file")
            or ledger_row.get("output_dir") != accounting.get("output_dir")
        ):
            raise AcceptanceError(
                "accepted canonical segment differs from ledger or reservation lineage"
            )
        utility = require_dict(
            command.get("production_utility"), "canonical segment production utility"
        )
        controller_revision = require_revision(
            utility.get("revision"), "canonical segment controller revision"
        )
        controller_sha = require_sha256(
            utility.get("sha256"), "canonical segment controller sha256"
        )
        if (
            utility.get("committed") is not True
            or resolve_bound_path(utility.get("path")).resolve(strict=True)
            != CANONICAL_CONTROLLER.resolve(strict=True)
        ):
            raise AcceptanceError("accepted canonical segment controller identity differs")
        source_bundle = require_dict(
            command.get("source_bundle"), "canonical segment source bundle"
        )
        source_path = resolve_bound_path(source_bundle.get("path")).resolve(strict=True)
        try:
            source_path.relative_to(
                Path(str(canonical["source_archive_root"])).resolve(strict=True)
            )
        except ValueError as error:
            raise AcceptanceError(
                "accepted canonical segment source bundle is outside the source archive"
            ) from error
        verified_revisions = [
            require_revision(value, "source bundle verified revision")
            for value in require_list(
                source_bundle.get("verified_revisions"), "source bundle verified revisions"
            )
        ]
        source_name = re.fullmatch(
            r"athenak-feature-cgl-through-([0-9a-f]{8,40})\.bundle",
            source_path.name,
        )
        if (
            controller_revision not in verified_revisions
            or source_name is None
            or not controller_revision.startswith(source_name.group(1))
        ):
            raise AcceptanceError(
                "accepted canonical segment controller revision is absent from its source bundle"
            )
        if source_path not in source_bundle_bindings:
            source_bundle_bindings[source_path] = verify_declared_binding(
                source_bundle, "canonical segment source bundle"
            )
            if archive_catalog.get(source_path.name) != source_bundle_bindings[source_path][
                "sha256"
            ]:
                raise AcceptanceError(
                    "accepted canonical segment source bundle is absent from the "
                    "current published source-archive catalog"
                )
            retained_bindings.append(source_bundle_bindings[source_path])
        approval = require_dict(
            command.get("qualification_approval"),
            "canonical segment qualification approval",
        )
        observed_approval = verify_declared_binding(
            approval, "canonical segment qualification approval"
        )
        if binding_identity(observed_approval, "observed qualification approval") != (
            Path(str(qualification["path"])),
            str(qualification["sha256"]),
        ):
            raise AcceptanceError("canonical segment qualification approval differs")
        controller_lineage.append({
            "segment_manifest": segment_binding,
            "job_id": job_id,
            "controller_revision": controller_revision,
            "controller_sha256": controller_sha,
            "source_bundle": source_bundle_bindings[source_path],
        })
    current_authority = {
        key: policy["verified_sources"][key]
        for key in (
            "current_source_authority_evidence",
            "current_source_authority_provenance_review",
            "current_source_authority_plasma_review",
            "current_source_authority_publication_audit",
            "source_archive_catalog",
            "qualification_approval",
        )
    }
    return {
        "campaign_root": str(CANONICAL_CAMPAIGN_ROOT),
        "run_root": str(CANONICAL_RUN_ROOT),
        "ledger": ledger_binding,
        "reservation_store": reservation_binding,
        "current_published_source_authority": current_authority,
        "controller_lineage": controller_lineage,
    }, retained_bindings


def authenticate_case_bundle(
    policy: dict[str, object],
    case_id: str,
    bundle_path: Path | None,
    mhd_binding: dict[str, object],
    user_binding: dict[str, object],
) -> tuple[dict[str, object] | None, list[dict[str, object]], dict[str, object]]:
    """Authenticate one accepted case bundle and its complete accepted lineage."""

    if bundle_path is None:
        return None, [], gate(
            "accepted_case_bundle_lineage",
            "inconclusive",
            reason="accepted case bundle manifest is unavailable",
        )
    manifest, bundle_binding = load_json(bundle_path, "accepted case bundle manifest")
    canonical_path = path_is_exact(bundle_path, canonical_bundle_path(case_id))
    expected_case = case_manifest_record(policy, case_id)
    name = str(expected_case["name"])
    if (
        manifest.get("workflow") != "paper-mks24-stage-i-production"
        or manifest.get("status") != "accepted_for_analysis"
        or manifest.get("production_case_id") != case_id
        or require_finite(manifest.get("required_final_time"), "bundle required final time")
        != 10.0
        or require_finite(manifest.get("accepted_final_time"), "bundle accepted final time")
        != 10.0
    ):
        raise AcceptanceError("accepted case bundle identity or final time differs")
    cases = require_list(manifest.get("cases"), "accepted bundle cases")
    if len(cases) != 1:
        raise AcceptanceError("accepted case bundle must contain exactly one case")
    selected = require_dict(cases[0], "accepted bundle case")
    if (
        selected.get("name") != name
        or selected.get("input") != expected_case.get("input")
        or selected.get("status") != "passed"
    ):
        raise AcceptanceError("accepted bundle selected case differs from Stage I manifest")
    outputs = require_dict(selected.get("outputs"), "accepted bundle case outputs")
    bundle_root = bundle_path.expanduser().absolute().resolve(strict=True).parent
    for key, observed in (
        ("mhd_history", mhd_binding),
        ("user_history", user_binding),
    ):
        output = outputs.get(key)
        if not isinstance(output, str) or not output:
            raise AcceptanceError(f"accepted bundle lacks {key}")
        expected_path = (bundle_root / output).resolve(strict=True)
        if expected_path != Path(str(observed["path"])):
            raise AcceptanceError(f"selected {key} is not the accepted bundle product")

    segment_paths = [
        resolve_bound_path(value).resolve(strict=True)
        for value in require_list(
            manifest.get("production_segment_manifests"),
            "production segment manifests",
        )
    ]
    if not segment_paths or len(segment_paths) != len(set(segment_paths)):
        raise AcceptanceError("accepted bundle segment lineage is empty or duplicated")
    segment_records: list[tuple[dict[str, object], dict[str, object]]] = [
        load_json(path, f"accepted segment manifest {index}")
        for index, path in enumerate(segment_paths)
    ]
    executable_sha: str | None = None
    input_sha: str | None = None
    previous_path: Path | None = None
    previous_final = -math.inf
    segment_bindings: list[dict[str, object]] = []
    for path, (segment, binding) in zip(segment_paths, segment_records):
        accounting = require_dict(segment.get("accounting"), "segment accounting")
        command = require_dict(segment.get("command"), "segment command")
        inspection = require_dict(
            segment.get("scientific_inspection"), "segment scientific inspection"
        )
        segment_executable = require_sha256(
            command.get("executable_sha256"), "segment executable sha256"
        )
        segment_input = require_sha256(command.get("input_sha256"), "segment input sha256")
        final_time = require_finite(inspection.get("final_time"), "segment final time")
        if (
            accounting.get("result") != "accepted"
            or accounting.get("case_id") != case_id
            or accounting.get("case_name") != name
            or accounting.get("executable_sha256") != segment_executable
            or inspection.get("accepted") is not True
            or inspection.get("case_id") != case_id
            or resolve_bound_path(inspection.get("manifest")).resolve(strict=True) != path
            or command.get("matrix_sha256")
            != policy["verified_sources"]["stage_i_manifest"]["sha256"]
            or final_time <= previous_final
        ):
            raise AcceptanceError("accepted segment lineage semantics differ")
        parent = command.get("parent_segment")
        if previous_path is None:
            segment_name = accounting.get("segment")
            if (
                parent is not None
                or not isinstance(segment_name, str)
                or not segment_name.startswith("s00")
                or command.get("restart_file") is not None
                or require_list(
                    command.get("restart_files"), "accepted root restart files"
                )
            ):
                raise AcceptanceError("accepted case lineage root is not an s00 t=0 origin")
        else:
            parent_record = require_dict(parent, "accepted segment parent")
            if (
                resolve_bound_path(parent_record.get("manifest")).resolve(strict=True)
                != previous_path
                or parent_record.get("case_id") != case_id
                or parent_record.get("result") != "accepted"
                or parent_record.get("executable_sha256") != segment_executable
                or parent_record.get("input_sha256") != segment_input
                or require_finite(parent_record.get("final_time"), "parent final time")
                != previous_final
            ):
                raise AcceptanceError("accepted segment parent lineage differs")
        if executable_sha is None:
            executable_sha = segment_executable
            input_sha = segment_input
        elif segment_executable != executable_sha or segment_input != input_sha:
            raise AcceptanceError("accepted segment executable or input lineage differs")
        previous_path = path
        previous_final = final_time
        segment_bindings.append(binding)
    if previous_final != 10.0:
        raise AcceptanceError("accepted case lineage does not terminate at t=10")
    terminal = require_dict(
        segment_records[-1][0]["scientific_inspection"], "terminal inspection"
    )
    terminal_checks = require_dict(terminal.get("checks"), "terminal inspection checks")
    if (
        terminal_checks.get("terminal_restart_physical_time_matches_final") is not True
        or require_finite(terminal.get("terminal_restart_time"), "terminal restart time")
        != 10.0
    ):
        raise AcceptanceError("accepted terminal inspection is incomplete")
    model = require_dict(selected.get("model_choices"), "accepted bundle model choices")
    forcing_tcorr = require_finite_number_or_text(
        model.get("forcing_tcorr"), "forcing tcorr"
    )
    if forcing_tcorr <= 0.0:
        raise AcceptanceError("forcing tcorr must be positive")
    authenticated = {
        "bundle_manifest": bundle_binding,
        "segment_manifests": segment_bindings,
        "case_id": case_id,
        "case_name": name,
        "executable_sha256": executable_sha,
        "input_sha256": input_sha,
        "forcing_tcorr": forcing_tcorr,
        "terminal_segment_manifest": segment_bindings[-1],
        "canonical_campaign_authority_eligible": False,
    }
    retained_bindings = [bundle_binding, *segment_bindings]
    if canonical_path:
        canonical_context, canonical_bindings = authenticate_canonical_accounting(
            policy, case_id, name, segment_paths, segment_records
        )
        authenticated["canonical_campaign_authority_eligible"] = True
        authenticated["canonical_context"] = canonical_context
        retained_bindings.extend(canonical_bindings)
        result = "pass"
        reason = (
            "canonical accepted whole-case bundle, controller/source lineage, "
            "ledger, reservations, and complete t=0..10 segment lineage authenticated"
        )
    else:
        result = "inconclusive"
        reason = (
            "bundle lineage is mechanically valid but its arbitrary noncanonical path "
            "is rejected as campaign-authorizing evidence"
        )
    return authenticated, retained_bindings, gate(
        "accepted_case_bundle_lineage",
        result,
        reason=reason,
        observations=authenticated,
    )


def reviewed_scientific_products_available(policy: dict[str, object]) -> bool:
    """Return whether criteria bind an independently reviewed replayable generator."""

    products = require_dict(
        policy["criteria"].get("scientific_products_policy"),
        "scientific products policy",
    )
    generator = require_dict(
        products.get("reviewed_generator_binding"), "reviewed generator binding"
    )
    return (
        generator.get("status") == "reviewed_available"
        and generator.get("review_status") == "approved"
        and isinstance(generator.get("path"), str)
        and isinstance(generator.get("sha256"), str)
        and isinstance(generator.get("review_path"), str)
        and isinstance(generator.get("review_sha256"), str)
    )


def validate_diagnostics_contract(
    policy: dict[str, object],
    case_id: str,
    window_name: str,
    diagnostics: dict[str, object] | None,
    bundle: dict[str, object] | None,
    mhd_binding: dict[str, object],
    user_binding: dict[str, object],
) -> tuple[dict[str, object] | None, dict[str, object]]:
    """Authenticate one exact-window replayed scientific-products artifact."""

    gate_name = f"scientific_products_contract:{window_name}"
    if not reviewed_scientific_products_available(policy):
        return None, gate(
            gate_name,
            "inconclusive",
            reason=(
                "the exact reviewed scientific-products generator and deterministic "
                "replay verifier are unavailable; hand-authored products cannot pass"
            ),
        )
    if diagnostics is None:
        return None, gate(
            gate_name,
            "inconclusive",
            reason=f"{window_name} scientific products are unavailable",
        )
    contract = diagnostics.get("scientific_acceptance_contract")
    if not isinstance(contract, dict):
        return None, gate(
            gate_name,
            "inconclusive",
            reason=f"{window_name} scientific products lack the immutable replay contract",
        )
    if bundle is None:
        return None, gate(
            gate_name,
            "inconclusive",
            reason="scientific products cannot bind an unavailable accepted case bundle",
        )
    try:
        expected_window = require_list(
            require_dict(
                policy["criteria"].get("analysis_windows"), "analysis windows"
            ).get(window_name),
            f"{window_name} analysis window",
        )
        observed_window = require_dict(contract.get("analysis_window"), "analysis window")
        expected_bundle = binding_identity(bundle["bundle_manifest"], "accepted bundle")
        observed_bundle = binding_identity(
            contract.get("accepted_bundle_manifest"), "diagnostics accepted bundle"
        )
        products = require_dict(
            policy["criteria"].get("scientific_products_policy"),
            "scientific products policy",
        )
        expected_generator = binding_identity(
            require_dict(products.get("reviewed_generator_binding"), "reviewed generator"),
            "expected scientific-products generator",
        )
        observed_generator = binding_identity(
            contract.get("generator"), "scientific-products generator"
        )
        observed_mhd = binding_identity(
            contract.get("mhd_history"), "diagnostics MHD history"
        )
        observed_user = binding_identity(
            contract.get("user_history"), "diagnostics user history"
        )
        if (
            contract.get("schema_version") != 2
            or contract.get("record_type")
            != "stage-i-scientific-acceptance-reviewed-products"
            or contract.get("case_id") != case_id
            or contract.get("case_name") != case_name(policy, case_id)
            or observed_window
            != {"time_start": float(expected_window[0]), "time_end": float(expected_window[1])}
            or observed_bundle != expected_bundle
            or observed_generator != expected_generator
            or observed_mhd != binding_identity(mhd_binding, "selected MHD history")
            or observed_user != binding_identity(user_binding, "selected user history")
            or contract.get("stage_i_manifest_sha256")
            != policy["verified_sources"]["stage_i_manifest"]["sha256"]
            or contract.get("deterministic_replay_verification")
            != require_dict(products.get("contract"), "scientific products contract").get(
                "deterministic_replay_verification"
            )
        ):
            raise AcceptanceError("scientific-products acceptance contract differs")
        selected = diagnostics_case_record(diagnostics, case_name(policy, case_id))
        if selected is None or require_dict(
            selected.get("analysis_window"), "diagnostics case analysis window"
        ) != observed_window:
            raise AcceptanceError("scientific-products selected case or exact window differs")
    except AcceptanceError as error:
        return None, gate(gate_name, "fail", reason=str(error))
    return diagnostics, gate(
        gate_name,
        "pass",
        reason=f"{window_name} reviewed scientific-products replay contract authenticated",
        observations={"analysis_window": observed_window},
    )


def diagnostics_case_record(
    diagnostics: dict[str, object], name: str
) -> dict[str, object] | None:
    """Return one analyzer case record when present."""

    cases = diagnostics.get("cases")
    if not isinstance(cases, dict):
        return None
    record = cases.get(name)
    return record if isinstance(record, dict) else None


def reference_comparison_records(
    diagnostics: dict[str, object],
) -> dict[str, dict[str, object]]:
    """Return analyzer reference curve and surface comparison records."""

    collection = diagnostics.get("reference_curve_comparisons")
    if not isinstance(collection, dict) or collection.get("available") is not True:
        return {}
    records: dict[str, dict[str, object]] = {}
    for key in ("comparisons", "surface_comparisons"):
        values = collection.get(key)
        if not isinstance(values, dict):
            continue
        for product_id, record in values.items():
            if isinstance(product_id, str) and isinstance(record, dict):
                if product_id in records:
                    raise AcceptanceError(
                        f"reference product is duplicated across collections: {product_id}"
                    )
                records[product_id] = record
    return records


def finite_vector(value: object, label: str, *, positive: bool = False) -> list[float]:
    """Return one nonempty finite numeric vector."""

    values = require_list(value, label)
    if not values:
        raise AcceptanceError(f"{label} is empty")
    result = [
        require_positive(item, label) if positive else require_finite(item, label)
        for item in values
    ]
    return result


def vectors_match(left: list[float], right: list[float]) -> bool:
    """Return whether two numeric vectors match at retained JSON/CSV precision."""

    return len(left) == len(right) and all(
        math.isclose(a, b, rel_tol=1.0e-12, abs_tol=1.0e-12)
        for a, b in zip(left, right)
    )


def verified_reference_product(
    policy: dict[str, object],
    product_id: str,
    binding: dict[str, object],
) -> dict[str, object]:
    """Load one exact Stage-I-bound reference manifest entry and CSV product."""

    manifest_status = require_dict(policy["manifest"].get("panel_status"), "panel_status")
    reference_manifests = require_dict(
        manifest_status.get("reference_manifests"), "reference_manifests"
    )
    manifest_id = binding.get("reference_manifest")
    if not isinstance(manifest_id, str):
        raise AcceptanceError(f"reference product {product_id} lacks a manifest id")
    manifest_binding = require_dict(
        reference_manifests.get(manifest_id), f"reference manifest {manifest_id}"
    )
    archive_manifest_path = Path(
        str(
            require_dict(
                policy["verified_sources"]["reference_archive_manifest"],
                "reference archive manifest binding",
            )["path"]
        )
    ).resolve(strict=True)
    archive_root = archive_manifest_path.parent
    relative_manifest = manifest_binding.get("path")
    if not isinstance(relative_manifest, str) or not relative_manifest:
        raise AcceptanceError(f"reference manifest {manifest_id} lacks a path")
    manifest_path = (archive_root / relative_manifest).resolve(strict=True)
    try:
        manifest_path.relative_to(archive_root)
    except ValueError as error:
        raise AcceptanceError("reference manifest escapes the bound archive") from error
    observed_manifest = regular_file_binding(
        manifest_path, f"reference manifest {manifest_id}"
    )
    if observed_manifest["sha256"] != require_sha256(
        manifest_binding.get("sha256"), f"reference manifest {manifest_id} sha256"
    ):
        raise AcceptanceError(f"reference manifest {manifest_id} SHA-256 differs")
    manifest, loaded_manifest = load_json(
        manifest_path, f"reference manifest {manifest_id}"
    )
    if loaded_manifest != observed_manifest or manifest.get("schema_version") != 1:
        raise AcceptanceError(f"reference manifest {manifest_id} binding or schema differs")
    kind = binding.get("kind")
    collection_name = "surfaces" if kind == "surface" else "curves" if kind == "curve" else None
    if collection_name is None:
        raise AcceptanceError(f"reference product {product_id} has an invalid kind")
    entries = [
        require_dict(item, f"reference manifest {manifest_id} product")
        for item in require_list(
            manifest.get(collection_name), f"reference manifest {manifest_id} {collection_name}"
        )
        if require_dict(item, f"reference manifest {manifest_id} product").get("id")
        == product_id
    ]
    if len(entries) != 1:
        raise AcceptanceError(f"reference manifest product {product_id} is absent or duplicated")
    entry = entries[0]
    for key in ("case", "product", "data_file", "data_sha256"):
        if entry.get(key) != binding.get(key):
            raise AcceptanceError(f"reference manifest product {product_id} binding differs")
    data_file = entry.get("data_file")
    if not isinstance(data_file, str) or not data_file:
        raise AcceptanceError(f"reference product {product_id} lacks a data file")
    data_path = (manifest_path.parent / data_file).resolve(strict=True)
    try:
        data_path.relative_to(manifest_path.parent)
    except ValueError as error:
        raise AcceptanceError("reference data file escapes its bound manifest directory") from error
    data_binding = regular_file_binding(data_path, f"reference data {product_id}")
    if data_binding["sha256"] != require_sha256(
        entry.get("data_sha256"), f"reference product {product_id} data sha256"
    ):
        raise AcceptanceError(f"reference product {product_id} data SHA-256 differs")
    payload = read_stable_bytes(data_path, f"reference data {product_id}")
    try:
        reader = csv.DictReader(io.StringIO(payload.decode("utf-8")))
        fields = reader.fieldnames
        rows = list(reader)
    except UnicodeDecodeError as error:
        raise AcceptanceError(f"reference product {product_id} CSV is not UTF-8") from error
    required = (
        {"x", "y", "z", "z_uncertainty"}
        if kind == "surface"
        else {"x", "y", "y_uncertainty"}
    )
    if (
        not fields
        or len(fields) != len(set(fields))
        or not required.issubset(fields)
        or not rows
        or any(None in row for row in rows)
    ):
        raise AcceptanceError(f"reference product {product_id} CSV contract differs")

    def column(name: str, *, positive: bool = False) -> list[float]:
        values: list[float] = []
        for row in rows:
            text = row.get(name)
            try:
                value = float(str(text))
            except ValueError as error:
                raise AcceptanceError(
                    f"reference product {product_id} column {name} is nonnumeric"
                ) from error
            values.append(
                require_positive(value, f"reference product {product_id} {name}")
                if positive
                else require_finite(value, f"reference product {product_id} {name}")
            )
        return values

    result: dict[str, object] = {
        "kind": kind,
        "manifest_sha256": observed_manifest["sha256"],
        "data_path": str(data_path),
        "data_sha256": data_binding["sha256"],
        "interpolation": str(entry.get("interpolation", "bilinear" if kind == "surface" else "linear")),
        "x": column("x"),
    }
    if kind == "surface":
        result.update({
            "y": column("y"),
            "reference": column("z"),
            "uncertainty": column("z_uncertainty", positive=True),
        })
    else:
        x = require_list(result["x"], "reference x")
        if any(float(right) <= float(left) for left, right in zip(x, x[1:])):
            raise AcceptanceError(f"reference product {product_id} x is not strictly increasing")
        result.update({
            "reference": column("y"),
            "uncertainty": column("y_uncertainty", positive=True),
        })
    return result


def normalized_reference_metrics(
    policy: dict[str, object],
    product_id: str,
    record: dict[str, object],
    binding: dict[str, object],
    expected_analysis_case: str,
) -> tuple[float, float, list[float]]:
    """Recompute and internally validate one exact reference comparison."""

    source = verified_reference_product(policy, product_id, binding)
    if (
        record.get("available") is not True
        or record.get("kind") != binding.get("kind")
        or record.get("case") != binding.get("case")
        or record.get("analysis_case") != expected_analysis_case
        or record.get("product") != binding.get("product")
        or record.get("reference_data_file") != binding.get("data_file")
        or record.get("data_sha256") != binding.get("data_sha256")
        or record.get("reference_manifest_sha256") != source["manifest_sha256"]
        or resolve_bound_path(record.get("data_file")).resolve(strict=True)
        != Path(str(source["data_path"]))
        or record.get("interpolation") != source["interpolation"]
        or record.get("stage_i_binding_validated") is not True
    ):
        raise AcceptanceError("reference product semantics or provenance differ")
    reference = record.get("reference_y")
    simulated = record.get("simulated_y")
    uncertainty = record.get("reference_y_uncertainty")
    if binding.get("kind") == "surface":
        reference = record.get("reference_z")
        simulated = record.get("simulated_z")
        uncertainty = record.get("reference_z_uncertainty")
    residual = record.get("residual")
    reference_values = finite_vector(reference, "reference values")
    simulated_values = finite_vector(simulated, "simulated values")
    uncertainty_values = finite_vector(
        uncertainty, "reference uncertainty", positive=True
    )
    residual_values = finite_vector(residual, "reference residual")
    if not (
        len(reference_values)
        == len(simulated_values)
        == len(uncertainty_values)
        == len(residual_values)
    ):
        raise AcceptanceError("reference product vector lengths differ")
    if require_int(record.get("sample_count"), "reference sample count", minimum=1) != len(reference_values):
        raise AcceptanceError("reference product sample count differs")
    source_reference = [float(value) for value in require_list(source["reference"], "source reference")]
    source_uncertainty = [float(value) for value in require_list(source["uncertainty"], "source uncertainty")]
    if not vectors_match(reference_values, source_reference) or not vectors_match(
        uncertainty_values, source_uncertainty
    ):
        raise AcceptanceError("reference values or uncertainties differ from bound CSV bytes")
    x_values = finite_vector(record.get("x"), "reference x")
    if not vectors_match(x_values, [float(value) for value in require_list(source["x"], "source x")]):
        raise AcceptanceError("reference x coordinates differ from bound CSV bytes")
    if binding.get("kind") == "surface":
        y_values = finite_vector(record.get("y"), "reference y")
        if not vectors_match(
            y_values, [float(value) for value in require_list(source["y"], "source y")]
        ):
            raise AcceptanceError("reference y coordinates differ from bound CSV bytes")
    recomputed = [
        simulated_value - reference_value
        for simulated_value, reference_value in zip(simulated_values, reference_values)
    ]
    if any(
        not math.isclose(observed, expected, rel_tol=1.0e-12, abs_tol=1.0e-12)
        for observed, expected in zip(residual_values, recomputed)
    ):
        raise AcceptanceError("reference residual vector is internally inconsistent")
    normalized = [
        value / error for value, error in zip(recomputed, uncertainty_values)
    ]
    rms = math.sqrt(sum(value * value for value in normalized) / len(normalized))
    maximum = max(abs(value) for value in normalized)
    raw_rms = math.sqrt(sum(value * value for value in recomputed) / len(recomputed))
    raw_maximum = max(abs(value) for value in recomputed)
    reported_rms = require_nonnegative(
        record.get("rms_normalized_by_reported_uncertainty"),
        "reported normalized RMS",
    )
    if (
        not math.isclose(reported_rms, rms, rel_tol=1.0e-12, abs_tol=1.0e-12)
        or not math.isclose(
            require_nonnegative(record.get("rms_residual"), "reported residual RMS"),
            raw_rms,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
        or not math.isclose(
            require_nonnegative(
                record.get("maximum_absolute_residual"),
                "reported maximum absolute residual",
            ),
            raw_maximum,
            rel_tol=1.0e-12,
            abs_tol=1.0e-12,
        )
    ):
        raise AcceptanceError("reported residual summaries are internally inconsistent")
    return rms, maximum, simulated_values


def panel_product_assessments(
    policy: dict[str, object],
    case_id: str,
    diagnostics: dict[str, dict[str, object] | None],
) -> list[dict[str, object]]:
    """Evaluate every admitted reference product relevant to one case."""

    manifest_status = require_dict(policy["manifest"].get("panel_status"), "panel_status")
    product_bindings = require_dict(
        manifest_status.get("reference_product_bindings"), "reference_product_bindings"
    )
    aliases = require_dict(
        manifest_status.get("analysis_case_aliases"), "analysis_case_aliases"
    )
    name = case_name(policy, case_id)
    criteria = policy["criteria"]
    panels = require_list(criteria.get("comparison_panels"), "comparison_panels")
    products_reviewed = reviewed_scientific_products_available(policy)
    full = reference_comparison_records(
        diagnostics["full"] if isinstance(diagnostics.get("full"), dict) else {}
    )
    early = reference_comparison_records(
        diagnostics["early"] if isinstance(diagnostics.get("early"), dict) else {}
    )
    late = reference_comparison_records(
        diagnostics["late"] if isinstance(diagnostics.get("late"), dict) else {}
    )
    results: list[dict[str, object]] = []
    for panel_value in panels:
        panel = require_dict(panel_value, "criteria panel")
        panel_id = str(panel["id"])
        manifest_panel = next(
            require_dict(item, "manifest panel")
            for item in require_list(manifest_status.get("panels"), "manifest panels")
            if require_dict(item, "manifest panel").get("id") == panel_id
        )
        if case_id not in manifest_panel["required_cases"]:
            continue
        relevant: list[str] = []
        for product_id in manifest_panel["reference_products"]:
            binding = require_dict(product_bindings.get(product_id), f"binding {product_id}")
            bound_case = str(binding.get("case"))
            if bound_case == name or aliases.get(bound_case) == name:
                relevant.append(str(product_id))
        if not relevant:
            results.append({
                "panel_id": panel_id,
                "product_id": None,
                "result": "inconclusive",
                "reason": "no bound reference product maps to this required case",
            })
            continue
        if not products_reviewed:
            results.extend({
                "panel_id": panel_id,
                "product_id": product_id,
                "result": "inconclusive",
                "reason": (
                    "reviewed scientific-products generator and deterministic replay "
                    "verification are unavailable"
                ),
            } for product_id in relevant)
            continue
        rms_limit = require_nonnegative(panel.get("normalized_residual_rms_lte"), "panel RMS limit")
        max_limit = require_nonnegative(panel.get("maximum_absolute_normalized_residual_lte"), "panel maximum limit")
        drift_limit = require_nonnegative(
            panel.get("early_late_vector_drift_rms_lte"),
            "panel early/late drift limit",
        )
        for product_id in relevant:
            record = full.get(product_id)
            if record is None:
                results.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "result": "inconclusive",
                    "reason": "full-window reference comparison is unavailable",
                })
                continue
            binding = require_dict(product_bindings.get(product_id), f"binding {product_id}")
            if record.get("available") is not True:
                results.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "result": "inconclusive",
                    "reason": "full-window reference comparison is unavailable",
                })
                continue
            try:
                rms, maximum, _ = normalized_reference_metrics(
                    policy, product_id, record, binding, name
                )
            except AcceptanceError as error:
                results.append({
                    "panel_id": panel_id,
                    "product_id": product_id,
                    "result": "fail",
                    "reason": str(error),
                })
                continue
            early_record = early.get(product_id)
            late_record = late.get(product_id)
            drift: float | None = None
            if (
                early_record is not None
                and late_record is not None
                and early_record.get("available") is True
                and late_record.get("available") is True
            ):
                try:
                    _, _, early_y = normalized_reference_metrics(
                        policy, product_id, early_record, binding, name
                    )
                    _, _, late_y = normalized_reference_metrics(
                        policy, product_id, late_record, binding, name
                    )
                    uncertainty = record.get("reference_y_uncertainty")
                    if binding.get("kind") == "surface":
                        uncertainty = record.get("reference_z_uncertainty")
                    if (
                        not isinstance(uncertainty, list)
                        or not early_y
                        or not (len(early_y) == len(late_y) == len(uncertainty))
                    ):
                        raise AcceptanceError("early/late reference vector lengths differ")
                    drift = math.sqrt(sum(
                        (
                            (require_finite(left, "early simulated value")
                             - require_finite(right, "late simulated value"))
                            / require_positive(error, "reference uncertainty")
                        ) ** 2
                        for left, right, error in zip(early_y, late_y, uncertainty)
                    ) / len(early_y))
                except AcceptanceError as error:
                    results.append({
                        "panel_id": panel_id,
                        "product_id": product_id,
                        "result": "fail",
                        "reason": str(error),
                    })
                    continue
            if drift is None:
                result = "inconclusive"
                reason = "early/late panel vectors are unavailable"
            elif rms <= rms_limit and maximum <= max_limit and drift <= drift_limit:
                result = "pass"
                reason = "reference residual and panel stationarity gates passed"
            else:
                result = "fail"
                reason = "reference residual or panel stationarity gate failed"
            results.append({
                "panel_id": panel_id,
                "product_id": product_id,
                "result": result,
                "reason": reason,
                "observations": {
                    "normalized_residual_rms": rms,
                    "maximum_absolute_normalized_residual": maximum,
                    "early_late_vector_drift_rms": drift,
                },
                "limits": {
                    "normalized_residual_rms_lte": rms_limit,
                    "maximum_absolute_normalized_residual_lte": max_limit,
                    "early_late_vector_drift_rms_lte": drift_limit,
                },
            })
    return results


def analyzer_metrics(
    diagnostics: dict[str, object] | None,
    name: str,
    policy: dict[str, object],
    minimum_block_duration: float,
) -> dict[str, dict[str, float]]:
    """Independently recompute analyzer scalar metrics from retained raw samples."""

    if diagnostics is None:
        return {}
    record = diagnostics_case_record(diagnostics, name)
    if record is None:
        return {}
    metrics = record.get("scientific_acceptance_metrics")
    if not isinstance(metrics, dict):
        return {}
    windows = require_dict(policy["criteria"].get("analysis_windows"), "analysis windows")
    full = require_list(windows.get("full"), "full analysis window")
    statistics = require_dict(
        policy["criteria"].get("statistics_policy"), "statistics policy"
    )
    gap_policy = require_dict(statistics.get("gap_policy"), "statistics gap policy")
    replicates = require_int(
        statistics.get("bootstrap_replicates"), "bootstrap_replicates", minimum=1
    )
    result: dict[str, dict[str, float]] = {}
    for metric, value in metrics.items():
        label = f"analyzer metric {metric}"
        metric_record = require_dict(value, label)
        if (
            not isinstance(metric_record.get("sample_times"), list)
            or not isinstance(metric_record.get("sample_values"), list)
        ):
            continue
        times = finite_vector(metric_record.get("sample_times"), f"{label} sample times")
        values = finite_vector(metric_record.get("sample_values"), f"{label} sample values")
        if len(times) != len(values) or any(
            right <= left for left, right in zip(times, times[1:])
        ):
            raise AcceptanceError(f"{label} raw sample contract differs")
        recomputed = window_statistics(
            times,
            values,
            require_finite(full[0], "full analysis start"),
            require_finite(full[1], "full analysis end"),
            replicates=replicates,
            seed_text=f"analyzer:{metric}:full:{policy['criteria_binding']['sha256']}",
            minimum_block_duration=minimum_block_duration,
            expected_cadence=require_positive(
                gap_policy.get("expected_history_cadence"), "expected history cadence"
            ),
            maximum_gap_expected_cadence_multiplier=require_positive(
                gap_policy.get("maximum_gap_expected_cadence_multiplier"),
                "maximum gap expected-cadence multiplier",
            ),
            maximum_gap_minimum_block_duration_fraction=require_nonnegative(
                gap_policy.get("maximum_gap_forcing_tcorr_fraction"),
                "maximum gap forcing-tcorr fraction",
            ),
            median_cadence_upper_relative_tolerance=require_nonnegative(
                gap_policy.get("median_cadence_upper_relative_tolerance"),
                "median cadence upper relative tolerance",
            ),
        )
        normalized = {
            "mean": float(recomputed["mean"]),
            "standard_error": float(recomputed["standard_error"]),
            "standard_deviation": float(recomputed["standard_deviation"]),
        }
        for key, observed in normalized.items():
            reported = (
                require_finite(metric_record.get(key), f"{label} {key}")
                if key == "mean"
                else require_nonnegative(metric_record.get(key), f"{label} {key}")
            )
            if not math.isclose(
                reported,
                observed,
                rel_tol=1.0e-12,
                abs_tol=1.0e-12,
            ):
                raise AcceptanceError(f"{label} reported {key} differs from raw samples")
        result[str(metric)] = normalized
    return result


def finite_curve(value: object, label: str) -> dict[str, list[float]]:
    """Normalize one finite ordered x/y curve."""

    record = require_dict(value, label)
    x = require_list(record.get("x"), f"{label} x")
    y = require_list(record.get("y"), f"{label} y")
    if len(x) != len(y) or len(x) < 3:
        raise AcceptanceError(f"{label} requires at least three x/y samples")
    x_values = [require_finite(item, f"{label} x") for item in x]
    y_values = [require_finite(item, f"{label} y") for item in y]
    if any(right <= left for left, right in zip(x_values, x_values[1:])):
        raise AcceptanceError(f"{label} x is not strictly increasing")
    result = {"x": x_values, "y": y_values}
    if "standard_error" in record:
        standard_error = [
            require_nonnegative(item, f"{label} standard_error")
            for item in require_list(record.get("standard_error"), f"{label} standard_error")
        ]
        if len(standard_error) != len(x_values):
            raise AcceptanceError(f"{label} standard_error length differs")
        result["standard_error"] = standard_error
    return result


def convergence_products(
    diagnostics: dict[str, object] | None, name: str
) -> dict[str, dict[str, list[float]]]:
    """Extract preregistered convergence curves from analyzer diagnostics."""

    if diagnostics is None:
        return {}
    record = diagnostics_case_record(diagnostics, name)
    if record is None:
        return {}
    explicit = record.get("scientific_acceptance_convergence")
    products: dict[str, dict[str, list[float]]] = {}
    if isinstance(explicit, dict):
        for product, curve in explicit.items():
            products[str(product)] = finite_curve(curve, f"convergence product {product}")
    ensemble = record.get("snapshot_ensemble")
    if not isinstance(ensemble, dict):
        return products
    spectra = ensemble.get("spectra")
    if isinstance(spectra, dict):
        for source, target in (
            ("velocity", "velocity_spectrum_shape"),
            ("magnetic_fluctuation", "magnetic_fluctuation_spectrum_shape"),
        ):
            spectral = spectra.get(source)
            if isinstance(spectral, dict):
                products[target] = finite_curve(
                    {"x": spectral.get("k"), "y": spectral.get("power_per_dk")},
                    target,
                )
    peak = ensemble.get("alignment_peak")
    if isinstance(peak, dict):
        products["peak_alignment"] = finite_curve(
            {"x": peak.get("k_perp"), "y": peak.get("cos_theta")},
            "peak_alignment",
        )
    return products


def load_optional_diagnostics(
    values: list[str] | None,
) -> tuple[
    dict[str, dict[str, object] | None],
    list[dict[str, object]],
    dict[str, dict[str, object]],
]:
    """Load optional ``window=path`` reviewed scientific-products artifacts."""

    diagnostics: dict[str, dict[str, object] | None] = {
        "full": None,
        "early": None,
        "late": None,
    }
    bindings: list[dict[str, object]] = []
    bindings_by_window: dict[str, dict[str, object]] = {}
    for item in values or []:
        if "=" not in item:
            raise AcceptanceError("diagnostics arguments must use window=path")
        window, path_text = item.split("=", 1)
        if window not in diagnostics or diagnostics[window] is not None:
            raise AcceptanceError(f"diagnostics window is invalid or duplicated: {window}")
        value, binding = load_json(Path(path_text), f"{window} analyzer diagnostics")
        diagnostics[window] = value
        bindings.append(binding)
        bindings_by_window[window] = binding
    return diagnostics, bindings, bindings_by_window


def aggregate_gate_result(gates: Iterable[dict[str, object]]) -> str:
    """Return fail, inconclusive, or pass for a required gate collection."""

    results = [str(gate_value["result"]) for gate_value in gates]
    if "fail" in results:
        return "fail"
    if "inconclusive" in results:
        return "inconclusive"
    if not results or all(result == "blocked_out_of_scope" for result in results):
        return "inconclusive"
    return "pass"


def evaluate_case(
    policy: dict[str, object],
    case_id: str,
    mhd_path: Path,
    user_path: Path,
    diagnostic_values: list[str] | None,
    ct_path: Path | None,
    bundle_path: Path | None = None,
) -> dict[str, object]:
    """Evaluate one case using retained histories and analyzer diagnostics."""

    case_id = require_case_id(case_id)
    mhd, mhd_binding = load_history(mhd_path, "MHD history")
    user, user_binding = load_history(user_path, "user history")
    diagnostics, diagnostic_bindings, diagnostic_bindings_by_window = (
        load_optional_diagnostics(diagnostic_values)
    )
    bundle, bundle_bindings, bundle_gate = authenticate_case_bundle(
        policy, case_id, bundle_path, mhd_binding, user_binding
    )
    trusted_diagnostics: dict[str, dict[str, object] | None] = {}
    diagnostics_gates: list[dict[str, object]] = []
    for window_name in ("full", "early", "late"):
        trusted, contract_gate = validate_diagnostics_contract(
            policy,
            case_id,
            window_name,
            diagnostics[window_name],
            bundle,
            mhd_binding,
            user_binding,
        )
        trusted_diagnostics[window_name] = trusted
        diagnostics_gates.append(contract_gate)
    name = case_name(policy, case_id)
    criteria = policy["criteria"]
    minimum_block_duration = (
        float(bundle["forcing_tcorr"]) if bundle is not None else 0.0
    )
    windows = require_dict(criteria.get("analysis_windows"), "analysis_windows")
    history_window_gate = gate(
        "history_exact_window_coverage",
        "pass",
        reason="selected histories cover every exact preregistered analysis window",
        observations={
            "analysis_windows": windows,
            "mhd_time_range": [mhd["time"][0], mhd["time"][-1]],
            "user_time_range": [user["time"][0], user["time"][-1]],
        },
    )
    for history, label in ((mhd, "MHD"), (user, "user")):
        for window_name in ("full", "early", "late"):
            window = require_list(windows[window_name], f"{window_name} window")
            interpolate_at(history["time"], history["time"], float(window[0]))
            interpolate_at(history["time"], history["time"], float(window[1]))
    metrics_policy = require_dict(criteria.get("case_metrics"), "case_metrics")
    metrics: dict[str, object] = {}
    stationarity_gates: list[dict[str, object]] = []
    for metric, spec_value in metrics_policy.items():
        spec = require_dict(spec_value, f"metric policy {metric}")
        source = spec.get("history")
        column = spec.get("column")
        kind = str(spec.get("stationarity_kind"))
        history = user if source == "user" else mhd if source == "mhd" else None
        if history is None or not isinstance(column, str) or column not in history:
            raise AcceptanceError(f"metric {metric} source or column is unavailable")
        record = metric_statistics(
            history,
            history[column],
            str(metric),
            policy,
            kind=kind,
            minimum_block_duration=minimum_block_duration,
        )
        metrics[str(metric)] = record
        sampling = str(record["sampling_adequacy"])
        stationarity = str(require_dict(record["stationarity"], "stationarity")["result"])
        result = (
            "inconclusive"
            if sampling != "pass"
            else "fail"
            if stationarity == "fail"
            else "pass"
        )
        stationarity_gates.append(gate(
            f"stationarity:{metric}",
            result,
            reason=(
                "early/late stationarity and effective-sample gates passed"
                if result == "pass"
                else "stationarity failed"
                if result == "fail"
                else "effective sample or independent time-block count is insufficient"
            ),
            observations=record,
        ))

    family = require_dict(criteria.get("family_gates"), "family_gates")
    full_start, full_end = (float(value) for value in windows["full"])
    early_start, early_end = (float(value) for value in windows["early"])
    late_start, late_end = (float(value) for value in windows["late"])
    family_gates: list[dict[str, object]] = []

    finite = require_dict(family.get("finite_limiter"), "finite_limiter")
    if case_id in finite["cases"]:
        hw_zero = exact_zero_window(mhd, "lf_hwproj", full_start, full_end)
        nu_late = require_dict(require_dict(metrics["nu_eff"], "nu_eff")["late"], "nu_eff late")
        occupancy = combine_occupancy_series(user)
        occupancy_stats = metric_statistics(
            user,
            occupancy,
            "unstable_occupancy",
            policy,
            kind="occupancy",
            minimum_block_duration=minimum_block_duration,
        )
        both_halves = (
            float(require_dict(occupancy_stats["early"], "occupancy early")["mean"])
            > float(finite["minimum_occupancy"])
            and float(require_dict(occupancy_stats["late"], "occupancy late")["mean"])
            > float(finite["minimum_occupancy"])
        )
        passed = (
            hw_zero
            and float(nu_late["confidence_interval_95"][0]) > 0.0
            and both_halves
        )
        family_gates.append(gate(
            "finite_limiter_semantics",
            "pass" if passed else "fail",
            reason="finite limiter semantic gates passed" if passed else "finite limiter semantic gate failed",
            observations={
                "hardwall_projection_exact_zero": hw_zero,
                "late_nu_eff_lower_95": nu_late["confidence_interval_95"][0],
                "occupancy": occupancy_stats,
                "occupancy_active_both_halves": both_halves,
            },
            limits={"minimum_occupancy": finite["minimum_occupancy"]},
        ))

    hardwall = require_dict(family.get("hardwall"), "hardwall")
    if case_id in hardwall["cases"]:
        early_delta = history_delta(mhd, "lf_hwproj", early_start, early_end)
        late_delta = history_delta(mhd, "lf_hwproj", late_start, late_end)
        occupancy = combine_occupancy_series(user)
        occupancy_stats = metric_statistics(
            user,
            occupancy,
            "unstable_occupancy",
            policy,
            kind="occupancy",
            minimum_block_duration=minimum_block_duration,
        )
        occupancy_mean = float(require_dict(occupancy_stats["full"], "occupancy full")["mean"])
        hard_zero = exact_zero_window(user, "hard_vol", full_start, full_end)
        passed = (
            early_delta > 0.0
            and late_delta > 0.0
            and occupancy_mean > float(hardwall["minimum_occupancy"])
            and hard_zero
        )
        family_gates.append(gate(
            "hardwall_activation",
            "pass" if passed else "fail",
            reason="hardwall activation gates passed" if passed else "hardwall activation gate failed",
            observations={
                "early_projection_increment": early_delta,
                "late_projection_increment": late_delta,
                "occupancy_mean": occupancy_mean,
                "hard_vol_exact_zero": hard_zero,
            },
            limits={"minimum_occupancy": hardwall["minimum_occupancy"]},
        ))

    lf_policy = require_dict(family.get("lf_strength"), "lf_strength")
    if case_id in lf_policy["cases"]:
        qface = history_delta(mhd, "lf_qface", full_start, full_end)
        qpr = history_delta(mhd, "lf_qprwrk", full_start, full_end)
        qpe = history_delta(mhd, "lf_qpewrk", full_start, full_end)
        state_scale = max(
            abs(float(require_dict(require_dict(metrics["kinetic"], "kinetic")["full"], "kinetic full")["mean"])),
            abs(float(require_dict(require_dict(metrics["magnetic"], "magnetic")["full"], "magnetic full")["mean"])),
            1.0,
        )
        absolute = max(abs(qpr), abs(qpe))
        normalized = absolute / state_scale
        passed = (
            qface > 0.0
            and absolute > float(lf_policy["activity_absolute_gt"])
            and normalized > float(lf_policy["activity_state_normalized_gt"])
        )
        family_gates.append(gate(
            "landau_fluid_activity",
            "pass" if passed else "fail",
            reason="LF face and applied-work activity gates passed" if passed else "LF activity gate failed",
            observations={
                "lf_qface_increment": qface,
                "lf_qprwrk_increment": qpr,
                "lf_qpewrk_increment": qpe,
                "maximum_absolute_work": absolute,
                "state_normalized_work": normalized,
            },
            limits={
                "activity_absolute_gt": lf_policy["activity_absolute_gt"],
                "activity_state_normalized_gt": lf_policy["activity_state_normalized_gt"],
            },
        ))

    passive_cases = set(str(value) for value in family["active_passive"]["passive_cases"])
    active_cases = set(str(value) for value in family["active_passive"]["active_cases"])
    if case_id in passive_cases:
        cp_zero = exact_zero_window(mhd, "lf_cpwrk", full_start, full_end)
        ca_zero = exact_zero_window(mhd, "lf_cawrk", full_start, full_end)
        family_gates.append(gate(
            "passive_pressure_work_exact_zero",
            "pass" if cp_zero and ca_zero else "fail",
            reason="passive pressure-work diagnostics are exactly zero" if cp_zero and ca_zero else "passive pressure-work diagnostics advanced",
            observations={"lf_cpwrk_exact_zero": cp_zero, "lf_cawrk_exact_zero": ca_zero},
        ))
    elif case_id in active_cases:
        cp = history_delta(mhd, "lf_cpwrk", full_start, full_end)
        ca = history_delta(mhd, "lf_cawrk", full_start, full_end)
        activity = max(abs(cp), abs(ca))
        limit = float(family["active_passive"]["activity_absolute_gt"])
        family_gates.append(gate(
            "active_pressure_work_activity",
            "pass" if activity > limit else "fail",
            reason="active pressure work exceeded the preregistered floor" if activity > limit else "active pressure work did not exceed the preregistered floor",
            observations={"lf_cpwrk_increment": cp, "lf_cawrk_increment": ca},
            limits={"activity_absolute_gt": limit},
        ))

    forcing = require_dict(family.get("forcing"), "forcing")
    if case_id in forcing["alfvenic_cases"] or case_id in forcing["random_cases"]:
        if "force_prp2" not in user or "force_prl2" not in user:
            raise AcceptanceError("user history lacks force_prp2 or force_prl2")
        fraction = [
            parallel / max(parallel + perpendicular, 1.0e-300)
            for perpendicular, parallel in zip(user["force_prp2"], user["force_prl2"])
        ]
        fraction_stats = metric_statistics(
            user,
            fraction,
            "parallel_forcing_fraction",
            policy,
            kind="scalar",
            minimum_block_duration=minimum_block_duration,
        )
        full_fraction = require_dict(fraction_stats["full"], "forcing fraction full")
        if case_id in forcing["alfvenic_cases"]:
            passed = float(full_fraction["mean"]) <= float(forcing["alfvenic_parallel_fraction_lte"])
            reason = "Alfvenic forcing remained effectively perpendicular"
        else:
            passed = float(full_fraction["confidence_interval_95"][0]) > float(
                forcing["random_parallel_fraction_lower_95_gt"]
            )
            reason = "random forcing retained a resolved parallel component"
        family_gates.append(gate(
            "forcing_geometry",
            "pass" if passed else "fail",
            reason=reason if passed else "forcing geometry gate failed",
            observations=fraction_stats,
            limits={
                "alfvenic_parallel_fraction_lte": forcing["alfvenic_parallel_fraction_lte"],
                "random_parallel_fraction_lower_95_gt": forcing["random_parallel_fraction_lower_95_gt"],
            },
        ))

    products = panel_product_assessments(policy, case_id, trusted_diagnostics)
    panel_gates = [
        gate(
            f"panel_product:{record['panel_id']}:{record.get('product_id')}",
            str(record["result"]),
            reason=str(record["reason"]),
            observations=record.get("observations"),
            limits=record.get("limits"),
        )
        for record in products
    ]

    ct_gate: dict[str, object]
    ct_binding: dict[str, object] | None = None
    if ct_path is None:
        ct_gate = gate(
            "sampled_restart_ct_divb",
            "inconclusive",
            reason="no sampled retained-state CT-divB audit was supplied",
        )
    else:
        ct_value, ct_binding = load_json(ct_path, "CT-divB evidence")
        verify_evidence_policy_binding(policy, ct_value, "CT-divB evidence")
        if (
            ct_value.get("record_type") != "stage-i-scientific-ct-divb-evidence"
            or ct_value.get("case_id") != case_id
        ):
            raise AcceptanceError("CT-divB evidence case or record type differs")
        ct_value = independently_recompute_ct_evidence(
            policy, ct_value, "CT-divB evidence"
        )
        ct_result = str(ct_value.get("result"))
        if ct_result not in ("pass", "fail", "inconclusive"):
            raise AcceptanceError("CT-divB evidence result is invalid")
        ct_gate = gate(
            "sampled_restart_ct_divb",
            ct_result,
            reason=str(ct_value.get("reason")),
            observations={
                "maximum_normalized_ct_divb": ct_value.get("maximum_normalized_ct_divb"),
                "coverage_complete": ct_value.get("coverage_complete"),
                "required_state_times": ct_value.get("required_state_times"),
                "state_count": ct_value.get("state_count"),
                "restart_count": ct_value.get("restart_count"),
            },
        )

    required_gates = [
        bundle_gate,
        history_window_gate,
        *diagnostics_gates,
        *stationarity_gates,
        *family_gates,
        *panel_gates,
        ct_gate,
    ]
    result = aggregate_gate_result(required_gates)
    inputs = [
        policy["criteria_binding"],
        policy["review_binding"],
        mhd_binding,
        user_binding,
        *bundle_bindings,
        *diagnostic_bindings,
    ]
    if ct_binding is not None:
        inputs.append(ct_binding)
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-scientific-case-evidence",
        "authority": "non-authorizing-scientific-assessment",
        "non_authorizing_statement": NON_AUTHORIZING_STATEMENT,
        "case_id": case_id,
        "case_name": name,
        "result": result,
        "criteria_review_status": policy["review_status"],
        "release_authorizing": False,
        "campaign_authority_eligible": bool(
            bundle is not None
            and bundle.get("canonical_campaign_authority_eligible") is True
        ),
        "analysis_windows": windows,
        "evaluation_inputs": {
            "mhd_history": mhd_binding,
            "user_history": user_binding,
            "accepted_bundle_manifest": (
                bundle["bundle_manifest"] if bundle is not None else None
            ),
            "diagnostics": diagnostic_bindings_by_window,
            "ct_evidence": ct_binding,
        },
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "acceptance_utility": regular_file_binding(Path(__file__), "acceptance utility"),
            "inputs": inputs,
        },
        "metrics": metrics,
        "analyzer_metrics": analyzer_metrics(
            trusted_diagnostics["full"]
            if isinstance(trusted_diagnostics.get("full"), dict)
            else None,
            name,
            policy,
            minimum_block_duration,
        ),
        "convergence_products": convergence_products(
            trusted_diagnostics["full"]
            if isinstance(trusted_diagnostics.get("full"), dict)
            else None,
            name,
        ),
        "panel_products": products,
        "gates": required_gates,
    }
    return seal_evidence(evidence)


def normal_two_sided_p(z_score: float) -> float:
    """Return the two-sided standard-normal tail probability."""

    if not math.isfinite(z_score):
        return 0.0 if z_score > 0.0 else 1.0
    return math.erfc(abs(z_score) / math.sqrt(2.0))


def scalar_from_case(
    case: dict[str, object], metric: str, window: str = "full"
) -> dict[str, float] | None:
    """Return one exact-window scalar estimate from case evidence."""

    if metric == "unstable_occupancy":
        return occupancy_from_case(case, window)
    metrics = case.get("metrics")
    if isinstance(metrics, dict) and metric in metrics:
        full = require_dict(
            require_dict(metrics[metric], f"{metric} metrics").get(window),
            f"{metric} {window}",
        )
        return {
            "mean": float(full["mean"]),
            "standard_error": float(full["standard_error"]),
            "standard_deviation": float(full["standard_deviation"]),
        }
    analyzer = case.get("analyzer_metrics")
    if isinstance(analyzer, dict) and metric in analyzer:
        return {
            key: float(value)
            for key, value in require_dict(analyzer[metric], f"analyzer metric {metric}").items()
        }
    return None


def occupancy_from_case(
    case: dict[str, object], window: str = "full"
) -> dict[str, float] | None:
    """Combine mirror and firehose occupancy estimates."""

    mirror = scalar_from_case(case, "mirror_occupancy", window)
    fire = scalar_from_case(case, "firehose_occupancy", window)
    if mirror is None or fire is None:
        return None
    return {
        "mean": mirror["mean"] + fire["mean"],
        "standard_error": math.hypot(mirror["standard_error"], fire["standard_error"]),
        "standard_deviation": math.hypot(
            mirror["standard_deviation"], fire["standard_deviation"]
        ),
    }


def pair_contrast(
    active: dict[str, object],
    passive: dict[str, object],
    policy: dict[str, object],
) -> dict[str, object]:
    """Evaluate one active/passive pair with Holm-corrected scalar contrasts."""

    family = require_dict(
        require_dict(policy["criteria"].get("family_gates"), "family_gates").get(
            "active_passive"
        ),
        "active_passive",
    )
    metrics = [str(value) for value in family["contrast_metrics"]]
    records: list[dict[str, object]] = []
    for metric in metrics:
        left = occupancy_from_case(active) if metric == "unstable_occupancy" else scalar_from_case(active, metric)
        right = occupancy_from_case(passive) if metric == "unstable_occupancy" else scalar_from_case(passive, metric)
        if left is None or right is None:
            records.append({"metric": metric, "available": False})
            continue
        difference = left["mean"] - right["mean"]
        se = math.hypot(left["standard_error"], right["standard_error"])
        z_score = finite_ratio(difference, se)
        pooled = math.sqrt(
            0.5 * (left["standard_deviation"] ** 2 + right["standard_deviation"] ** 2)
        )
        effect = finite_ratio(difference, pooled)
        records.append({
            "metric": metric,
            "available": True,
            "difference": difference,
            "combined_standard_error": se,
            "z_score": z_score,
            "two_sided_p": normal_two_sided_p(z_score),
            "standardized_effect": effect,
        })
    available = [record for record in records if record.get("available")]
    ordered = sorted(available, key=lambda record: float(record["two_sided_p"]))
    alpha = float(family["holm_alpha"])
    continuing = True
    for index, record in enumerate(ordered):
        threshold = alpha / (len(ordered) - index)
        record["holm_threshold"] = threshold
        significant = continuing and float(record["two_sided_p"]) <= threshold
        record["holm_significant"] = significant
        if not significant:
            continuing = False
    passed = any(
        record.get("holm_significant")
        and float(record["standardized_effect"]) >= float(family["minimum_standardized_effect"])
        for record in available
    )
    result = "pass" if passed else ("inconclusive" if len(available) < len(metrics) else "fail")
    return {
        "result": result,
        "reason": (
            "at least one preregistered active/passive contrast is resolved"
            if passed
            else "one or more preregistered contrast metrics are unavailable"
            if result == "inconclusive"
            else "no preregistered active/passive contrast passed"
        ),
        "metrics": records,
    }


def interpolate_curve(curve: dict[str, object], target: float) -> float:
    """Linearly interpolate one normalized convergence curve."""

    x = [float(value) for value in curve["x"]]
    y = [float(value) for value in curve["y"]]
    return interpolate_at(x, y, target)


def curve_sample_coordinates(
    left: dict[str, object],
    right: dict[str, object],
    interval: tuple[float, float],
    required_samples: list[float] | None,
) -> list[float]:
    """Return the exact preregistered physical-k samples for a curve comparison."""

    left_x = [float(value) for value in left["x"]]
    right_x = [float(value) for value in right["x"]]
    lower, upper = interval
    if (
        lower >= upper
        or left_x[0] > lower
        or left_x[-1] < upper
        or right_x[0] > lower
        or right_x[-1] < upper
    ):
        raise AcceptanceError("convergence curves do not cover the full preregistered interval")
    if required_samples is not None:
        if (
            not required_samples
            or required_samples != sorted(set(required_samples))
            or required_samples[0] < lower
            or required_samples[-1] > upper
            or any(
                not any(math.isclose(sample, value, rel_tol=1.0e-12, abs_tol=1.0e-12)
                        for value in left_x)
                or not any(math.isclose(sample, value, rel_tol=1.0e-12, abs_tol=1.0e-12)
                           for value in right_x)
                for sample in required_samples
            )
        ):
            raise AcceptanceError("convergence curves lack the complete preregistered shell list")
        samples = list(required_samples)
    else:
        samples = sorted(
            set([lower, upper, *(
                value for value in left_x + right_x if lower <= value <= upper
            )])
        )
    if len(samples) < 3:
        raise AcceptanceError("convergence curves lack three common resolved samples")
    return samples


def curve_distance(
    left: dict[str, object],
    right: dict[str, object],
    interval: tuple[float, float],
    kind: str,
    required_samples: list[float] | None = None,
) -> float:
    """Return one preregistered curve distance over a common interval."""

    samples = curve_sample_coordinates(left, right, interval, required_samples)
    left_y = [interpolate_curve(left, value) for value in samples]
    right_y = [interpolate_curve(right, value) for value in samples]
    if kind == "alignment":
        return max(abs(a - b) for a, b in zip(left_y, right_y))
    if any(value <= 0.0 for value in left_y + right_y):
        raise AcceptanceError("spectral convergence curves must be positive")
    left_norm = sum(left_y)
    right_norm = sum(right_y)
    return math.sqrt(sum(
        (math.log(a / left_norm) - math.log(b / right_norm)) ** 2
        for a, b in zip(left_y, right_y)
    ) / len(samples))


def curve_distance_resolution(
    left: dict[str, object],
    right: dict[str, object],
    interval: tuple[float, float],
    kind: str,
    required_samples: list[float] | None,
    sigma_gt: float,
) -> dict[str, object]:
    """Assess whether one curve disagreement is resolved by bound uncertainties."""

    distance = curve_distance(left, right, interval, kind, required_samples)
    if not isinstance(left.get("standard_error"), list) or not isinstance(
        right.get("standard_error"), list
    ):
        return {
            "distance": distance,
            "distance_standard_error": None,
            "resolved": False,
            "resolution_applicability": "not_applicable_uncertainty_unavailable",
        }
    samples = curve_sample_coordinates(left, right, interval, required_samples)
    left_error = {
        "x": left["x"],
        "y": left["standard_error"],
    }
    right_error = {
        "x": right["x"],
        "y": right["standard_error"],
    }
    combined = [
        math.hypot(
            interpolate_curve(left_error, sample),
            interpolate_curve(right_error, sample),
        )
        for sample in samples
    ]
    if kind == "alignment":
        differences = [
            abs(interpolate_curve(left, sample) - interpolate_curve(right, sample))
            for sample in samples
        ]
        index = max(range(len(samples)), key=differences.__getitem__)
        distance_se = combined[index]
    else:
        left_y = [interpolate_curve(left, sample) for sample in samples]
        right_y = [interpolate_curve(right, sample) for sample in samples]
        if any(value <= 0.0 for value in left_y + right_y):
            raise AcceptanceError("spectral convergence curves must be positive")
        distance_se = math.sqrt(sum(
            ((error / left_value) ** 2 + (error / right_value) ** 2)
            for error, left_value, right_value in zip(combined, left_y, right_y)
        ) / len(samples))
    resolved = distance > require_positive(sigma_gt, "curve resolution sigma") * distance_se
    return {
        "distance": distance,
        "distance_standard_error": distance_se,
        "resolved": resolved,
        "resolution_applicability": (
            "applicable_resolved" if resolved else "not_applicable_unresolved"
        ),
    }


def path_from_evidence_binding(binding: object, label: str) -> Path:
    """Reverify and return one exact path from case evaluation inputs."""

    observed = verify_declared_binding(binding, label)
    return Path(str(observed["path"]))


def require_same_evidence_semantics(
    expected: dict[str, object],
    observed: dict[str, object],
    label: str,
) -> None:
    """Require exact independently recomputed evidence excluding its self-digest."""

    expected_body = dict(expected)
    expected_body.pop("evidence_digest", None)
    observed_body = dict(observed)
    observed_body.pop("evidence_digest", None)
    if canonical_json(expected_body) != canonical_json(observed_body):
        raise AcceptanceError(f"{label} semantic contents differ from independent re-evaluation")


def independently_recompute_case(
    policy: dict[str, object], value: dict[str, object], label: str
) -> dict[str, object]:
    """Re-evaluate a case from mandatory bound inputs and compare exact semantics."""

    case_id = require_case_id(value.get("case_id"), f"{label} case_id")
    inputs = require_dict(value.get("evaluation_inputs"), f"{label} evaluation_inputs")
    mhd_path = path_from_evidence_binding(inputs.get("mhd_history"), f"{label} MHD history")
    user_path = path_from_evidence_binding(inputs.get("user_history"), f"{label} user history")
    bundle_value = inputs.get("accepted_bundle_manifest")
    bundle_path = (
        path_from_evidence_binding(bundle_value, f"{label} accepted bundle manifest")
        if isinstance(bundle_value, dict)
        else None
    )
    diagnostics = require_dict(inputs.get("diagnostics"), f"{label} diagnostics")
    diagnostic_values: list[str] = []
    for window in ("full", "early", "late"):
        diagnostic = diagnostics.get(window)
        if isinstance(diagnostic, dict):
            diagnostic_values.append(
                f"{window}={path_from_evidence_binding(diagnostic, f'{label} {window} diagnostics')}"
            )
    ct_value = inputs.get("ct_evidence")
    ct_path = (
        path_from_evidence_binding(ct_value, f"{label} CT evidence")
        if isinstance(ct_value, dict)
        else None
    )
    recomputed = evaluate_case(
        policy,
        case_id,
        mhd_path,
        user_path,
        diagnostic_values,
        ct_path,
        bundle_path,
    )
    require_same_evidence_semantics(value, recomputed, label)
    return recomputed


def case_evidence_is_complete(case: dict[str, object]) -> bool:
    """Return whether one recomputed case has every mandatory acceptance input/gate."""

    inputs = case.get("evaluation_inputs")
    if not isinstance(inputs, dict):
        return False
    diagnostics = inputs.get("diagnostics")
    if (
        not isinstance(inputs.get("mhd_history"), dict)
        or not isinstance(inputs.get("user_history"), dict)
        or not isinstance(inputs.get("accepted_bundle_manifest"), dict)
        or not isinstance(inputs.get("ct_evidence"), dict)
        or not isinstance(diagnostics, dict)
        or any(not isinstance(diagnostics.get(window), dict) for window in ("full", "early", "late"))
    ):
        return False
    gates = case.get("gates")
    if not isinstance(gates, list):
        return False
    by_name = {
        str(require_dict(item, "case gate").get("name")): require_dict(item, "case gate")
        for item in gates
    }
    mandatory = {
        "accepted_case_bundle_lineage",
        "history_exact_window_coverage",
        "scientific_products_contract:full",
        "scientific_products_contract:early",
        "scientific_products_contract:late",
        "sampled_restart_ct_divb",
    }
    return (
        case.get("campaign_authority_eligible") is True
        and
        mandatory.issubset(by_name)
        and all(by_name[name].get("result") == "pass" for name in mandatory)
        and all(require_dict(item, "case gate").get("result") == "pass" for item in gates)
        and case.get("result") == "pass"
    )


def evaluate_campaign(
    policy: dict[str, object],
    case_paths: list[Path],
) -> dict[str, object]:
    """Combine independently generated case evidence into campaign gates."""

    cases: dict[str, dict[str, object]] = {}
    bindings: list[dict[str, object]] = []
    for path in case_paths:
        value, binding = load_json(path, "case evidence")
        verify_evidence_policy_binding(policy, value, "case evidence")
        if value.get("record_type") != "stage-i-scientific-case-evidence":
            raise AcceptanceError("campaign input is not case evidence")
        case_id = require_case_id(value.get("case_id"))
        if case_id in cases:
            raise AcceptanceError(f"campaign case evidence is duplicated: {case_id}")
        cases[case_id] = independently_recompute_case(
            policy, value, f"case evidence {case_id}"
        )
        bindings.append(binding)

    required_cases = [str(value) for value in policy["criteria"]["required_cases"]]
    campaign_gates: list[dict[str, object]] = []
    campaign_gates.append(gate(
        "approved_independent_criteria_review",
        "pass" if policy["approved"] else "inconclusive",
        reason=(
            "independent criteria review approved exact criteria and final utility"
            if policy["approved"]
            else "independent criteria review remains changes-required or pending"
        ),
        observations={
            "review_status": policy["review_status"],
            "criteria_sha256": policy["criteria_binding"]["sha256"],
            "acceptance_utility_sha256": policy["verified_sources"]["acceptance_utility"]["sha256"],
        },
    ))
    products_available = reviewed_scientific_products_available(policy)
    campaign_gates.append(gate(
        "reviewed_scientific_products_generator",
        "pass" if products_available else "inconclusive",
        reason=(
            "exact reviewed generator and deterministic replay contract are available"
            if products_available
            else "companion scientific-products generator/review is unavailable or changes-required"
        ),
        observations=policy["criteria"]["scientific_products_policy"],
    ))
    missing = sorted(set(required_cases) - set(cases))
    campaign_gates.append(gate(
        "required_cases",
        "pass" if not missing else "inconclusive",
        reason="all required cases are present" if not missing else "required case evidence is missing",
        observations={"missing_cases": missing},
    ))
    nonpassing = {
        case_id: value["result"]
        for case_id, value in cases.items()
        if value.get("result") != "pass"
    }
    campaign_gates.append(gate(
        "case_scientific_acceptance",
        "pass" if not missing and not nonpassing else (
            "fail" if any(value == "fail" for value in nonpassing.values()) else "inconclusive"
        ),
        reason="every case passed" if not missing and not nonpassing else "one or more cases have not passed",
        observations={"nonpassing_cases": nonpassing},
    ))
    incomplete_cases = sorted(
        case_id for case_id, value in cases.items()
        if not case_evidence_is_complete(value)
    )
    campaign_gates.append(gate(
        "complete_authenticated_case_evidence",
        "pass" if not missing and not incomplete_cases else "inconclusive",
        reason=(
            "all required cases have independently recomputed complete evidence"
            if not missing and not incomplete_cases
            else "one or more cases lack complete authenticated evidence"
        ),
        observations={"incomplete_cases": incomplete_cases, "missing_cases": missing},
    ))
    lineage_failures = sorted(
        case_id for case_id, value in cases.items()
        if not any(
            require_dict(item, "case gate").get("name") == "accepted_case_bundle_lineage"
            and require_dict(item, "case gate").get("result") == "pass"
            for item in require_list(value.get("gates"), "case gates")
        )
    )
    campaign_gates.append(gate(
        "authenticated_accepted_case_bundle_lineages",
        "pass" if not missing and not lineage_failures else "inconclusive",
        reason=(
            "all 16 accepted case bundle lineages authenticated"
            if not missing and not lineage_failures
            else "one or more accepted case bundle lineages are missing or unauthenticated"
        ),
        observations={"lineage_failures": lineage_failures, "missing_cases": missing},
    ))
    noncanonical_cases = sorted(
        case_id for case_id, value in cases.items()
        if value.get("campaign_authority_eligible") is not True
    )
    campaign_gates.append(gate(
        "canonical_campaign_authority",
        "pass" if not missing and not noncanonical_cases else "inconclusive",
        reason=(
            "all 16 cases bind the exact canonical campaign root, whole-case bundles, "
            "controller/source lineage, ledger, reservations, and current source authority"
            if not missing and not noncanonical_cases
            else "one or more cases are local, offline, or lack exact canonical authority"
        ),
        observations={
            "canonical_campaign_root": str(CANONICAL_CAMPAIGN_ROOT),
            "noncanonical_cases": noncanonical_cases,
            "missing_cases": missing,
        },
    ))

    manifest_status = require_dict(policy["manifest"].get("panel_status"), "panel_status")
    manifest_panels = {
        str(require_dict(item, "manifest panel")["id"]): require_dict(item, "manifest panel")
        for item in manifest_status["panels"]
    }
    for panel_value in policy["criteria"]["comparison_panels"]:
        panel = require_dict(panel_value, "criteria panel")
        panel_id = str(panel["id"])
        expected = set(str(value) for value in manifest_panels[panel_id]["reference_products"])
        if not products_available:
            campaign_gates.append(gate(
                f"panel:{panel_id}",
                "inconclusive",
                reason=(
                    "reviewed scientific-products generator and deterministic replay "
                    "verification are unavailable"
                ),
                observations={"missing_products": sorted(expected)},
            ))
            continue
        products: dict[str, dict[str, object]] = {}
        for case in cases.values():
            for product in require_list(case.get("panel_products"), "case panel_products"):
                record = require_dict(product, "case panel product")
                if record.get("panel_id") == panel_id and isinstance(record.get("product_id"), str):
                    if str(record["product_id"]) in products:
                        raise AcceptanceError(
                            f"campaign panel product is duplicated: {record['product_id']}"
                        )
                    products[str(record["product_id"])] = record
        missing_products = sorted(expected - set(products))
        failures = sorted(
            product for product, record in products.items()
            if record.get("result") == "fail"
        )
        inconclusive = sorted(
            product for product, record in products.items()
            if record.get("result") == "inconclusive"
        )
        result = "fail" if failures else (
            "inconclusive" if missing_products or inconclusive else "pass"
        )
        campaign_gates.append(gate(
            f"panel:{panel_id}",
            result,
            reason=(
                "all admitted panel products passed"
                if result == "pass"
                else "one or more panel products failed"
                if result == "fail"
                else "one or more panel products are missing or inconclusive"
            ),
            observations={
                "missing_products": missing_products,
                "failed_products": failures,
                "inconclusive_products": inconclusive,
            },
        ))
    for panel_id in policy["criteria"]["allowed_exclusions"]["blocked_or_external_panels"]:
        campaign_gates.append(gate(
            f"panel:{panel_id}",
            "blocked_out_of_scope",
            reason="panel is explicitly blocked or external in the bound Stage I manifest",
        ))

    active_passive = policy["criteria"]["family_gates"]["active_passive"]
    for active, passive in active_passive["pairs"]:
        if active not in cases or passive not in cases:
            contrast = {"result": "inconclusive", "reason": "pair case evidence is missing"}
        else:
            contrast = pair_contrast(cases[active], cases[passive], policy)
        campaign_gates.append(gate(
            f"active_passive_pair:{active}:{passive}",
            str(contrast["result"]),
            reason=str(contrast["reason"]),
            observations=contrast.get("metrics"),
        ))

    if "R14" in cases and "R15" in cases:
        lower = scalar_from_case(cases["R14"], "nu_eff", "late")
        upper = scalar_from_case(cases["R15"], "nu_eff", "late")
        if lower is None or upper is None:
            ordering_result = "inconclusive"
            ordering_observation = None
        else:
            difference = upper["mean"] - lower["mean"]
            lower_95 = difference - 1.96 * math.hypot(
                upper["standard_error"], lower["standard_error"]
            )
            ordering_result = "pass" if lower_95 > 0.0 else "fail"
            ordering_observation = {
                "R15_minus_R14": difference,
                "lower_95": lower_95,
            }
    else:
        ordering_result = "inconclusive"
        ordering_observation = None
    campaign_gates.append(gate(
        "finite_limiter_ordering:R15_gt_R14",
        ordering_result,
        reason=(
            "R15 effective collisionality exceeds R14 with 95% confidence"
            if ordering_result == "pass"
            else "finite-limiter ordering is unavailable or unresolved"
        ),
        observations=ordering_observation,
    ))

    lf_policy = require_dict(
        require_dict(
            policy["criteria"].get("family_gates"), "family_gates"
        ).get("lf_strength"),
        "lf_strength",
    )
    lf_cases = [str(value) for value in require_list(lf_policy.get("cases"), "LF cases")]
    lf_resolved: list[dict[str, object]] = []
    for case_id in (value for value in lf_cases if value != "R02"):
        left = scalar_from_case(cases[case_id], "peak_alignment") if case_id in cases else None
        right = scalar_from_case(cases["R02"], "peak_alignment") if "R02" in cases else None
        if left is None or right is None:
            lf_resolved.append({
                "pair": [case_id, "R02"],
                "available": False,
                "missing_dependency": "immutable analyzer metric peak_alignment",
            })
            continue
        difference = abs(left["mean"] - right["mean"])
        lower_95 = difference - 1.96 * math.hypot(
            left["standard_error"], right["standard_error"]
        )
        lf_resolved.append({
            "pair": [case_id, "R02"],
            "available": True,
            "absolute_difference": difference,
            "lower_95": lower_95,
            "resolved": lower_95 > 0.0,
        })
    lf_result = (
        "pass"
        if all(record.get("resolved") for record in lf_resolved)
        else "inconclusive"
        if any(not record.get("available") for record in lf_resolved)
        else "fail"
    )
    campaign_gates.append(gate(
        "lf_strength_resolved_response",
        lf_result,
        reason=(
            "R12/R02/R06/R13 peak-alignment responses are resolved"
            if lf_result == "pass"
            else "LF-strength response is unavailable or unresolved"
        ),
        observations=lf_resolved,
    ))

    convergence_policy = require_dict(
        policy["criteria"].get("resolution_convergence"), "resolution_convergence"
    )
    convergence_cases = ("R16", "R02", "R17")
    interval_over_pi = tuple(
        float(value) for value in convergence_policy["common_k_perp_over_pi"]
    )
    interval = tuple(value * math.pi for value in interval_over_pi)
    if not products_available:
        records = []
        convergence_result = "inconclusive"
    elif all(case_id in cases for case_id in convergence_cases):
        convergence_records: list[dict[str, object]] = []
        product_specs = {
            "peak_alignment": ("alignment", float(convergence_policy["alignment_max_abs_lte"])),
            "velocity_spectrum_shape": ("spectrum", float(convergence_policy["spectral_log_rms_lte"])),
            "magnetic_fluctuation_spectrum_shape": ("spectrum", float(convergence_policy["spectral_log_rms_lte"])),
        }
        for product, (kind, limit) in product_specs.items():
            curves = [
                require_dict(cases[case_id].get("convergence_products"), "convergence products").get(product)
                for case_id in convergence_cases
            ]
            if not all(isinstance(curve, dict) for curve in curves):
                convergence_records.append({
                    "product": product,
                    "available": False,
                    "missing_dependency": "immutable analyzer convergence curve",
                })
                continue
            try:
                required_samples = (
                    [
                        float(value) * math.pi
                        for value in convergence_policy["alignment_shells"]
                    ]
                    if kind == "alignment"
                    else None
                )
                low_mid = curve_distance_resolution(
                    curves[0],
                    curves[1],
                    interval,
                    kind,
                    required_samples,
                    float(convergence_policy["resolved_disagreement_sigma_gt"]),
                )
                mid_high = curve_distance(
                    curves[1], curves[2], interval, kind, required_samples
                )
            except AcceptanceError as error:
                convergence_records.append({
                    "product": product,
                    "available": False,
                    "missing_dependency": str(error),
                })
                continue
            improvement_applicable = bool(low_mid["resolved"])
            improved = (
                mid_high
                <= float(convergence_policy["improvement_ratio_lte"])
                * float(low_mid["distance"])
                if improvement_applicable
                else True
            )
            convergence_records.append({
                "product": product,
                "available": True,
                "coordinate": "physical_k",
                "physical_k_interval": list(interval),
                "R16_R02_distance": low_mid["distance"],
                "R16_R02_distance_standard_error": low_mid["distance_standard_error"],
                "resolved_R16_R02_difference": low_mid["resolved"],
                "R02_R17_distance": mid_high,
                "R02_R17_limit": limit,
                "improvement_ratio_lte": convergence_policy["improvement_ratio_lte"],
                "improvement_applicability": low_mid["resolution_applicability"],
                "improved": improved,
                "passed": mid_high <= limit and improved,
            })
        scalar_records: list[dict[str, object]] = []
        for metric in convergence_policy["scalar_metrics"]:
            values = [scalar_from_case(cases[case_id], metric) for case_id in convergence_cases]
            if any(value is None for value in values):
                scalar_records.append({
                    "metric": metric,
                    "available": False,
                    "missing_dependency": "immutable exact-window scalar metric",
                })
                continue
            low, mid, high = values
            mid_high = abs(high["mean"] - mid["mean"])
            relative = finite_ratio(mid_high, abs(mid["mean"]))
            within_uncertainty = mid_high <= 2.0 * math.hypot(
                high["standard_error"], mid["standard_error"]
            )
            low_mid = abs(mid["mean"] - low["mean"])
            resolved_low_mid = low_mid > 2.0 * math.hypot(
                mid["standard_error"], low["standard_error"]
            )
            improved = (
                not resolved_low_mid
                or mid_high <= float(convergence_policy["improvement_ratio_lte"]) * low_mid
            )
            scalar_records.append({
                "metric": metric,
                "available": True,
                "relative_R02_R17_difference": relative,
                "within_two_combined_standard_errors": within_uncertainty,
                "resolved_R16_R02_difference": resolved_low_mid,
                "improved": improved,
                "passed": (
                    (relative <= float(convergence_policy["scalar_relative_difference_lte"]) or within_uncertainty)
                    and improved
                ),
            })
        records = convergence_records + scalar_records
        convergence_result = (
            "pass"
            if records and all(record.get("passed") for record in records)
            else "inconclusive"
            if any(not record.get("available") for record in records)
            else "fail"
        )
    else:
        records = []
        convergence_result = "inconclusive"
    campaign_gates.append(gate(
        "R16_R02_R17_resolution_convergence",
        convergence_result,
        reason="resolution convergence gates passed" if convergence_result == "pass" else "resolution convergence is unavailable or failed",
        observations=records,
        limits=convergence_policy,
    ))

    result = aggregate_gate_result(campaign_gates)
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-scientific-campaign-evidence",
        "authority": "non-authorizing-scientific-assessment",
        "non_authorizing_statement": NON_AUTHORIZING_STATEMENT,
        "result": result,
        "criteria_review_status": policy["review_status"],
        "release_authorizing": False,
        "case_results": {case_id: value["result"] for case_id, value in sorted(cases.items())},
        "gates": campaign_gates,
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "acceptance_utility": regular_file_binding(Path(__file__), "acceptance utility"),
            "case_evidence": bindings,
        },
    }
    return seal_evidence(evidence)


def independently_recompute_campaign(
    policy: dict[str, object], value: dict[str, object], label: str
) -> dict[str, object]:
    """Re-evaluate one campaign from its exact bound case-evidence files."""

    provenance = require_dict(value.get("provenance"), f"{label} provenance")
    case_paths = [
        path_from_evidence_binding(binding, f"{label} case evidence {index}")
        for index, binding in enumerate(
            require_list(provenance.get("case_evidence"), f"{label} case evidence")
        )
    ]
    recomputed = evaluate_campaign(policy, case_paths)
    require_same_evidence_semantics(value, recomputed, label)
    return recomputed


def parse_parameter_dump(payload: bytes) -> tuple[dict[str, str], int]:
    """Parse the bounded text parameter dump at the beginning of a restart."""

    marker = b"<par_end>"
    location = payload.find(marker)
    if location < 0:
        raise AcceptanceError("restart lacks <par_end> within the reviewed prefix")
    end = location + len(marker) + 1
    try:
        text = payload[:end].decode("utf-8")
    except UnicodeDecodeError as error:
        raise AcceptanceError("restart parameter dump is not UTF-8") from error
    block: str | None = None
    values: dict[str, str] = {}
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1]
            continue
        if "=" in line and block is not None:
            key, value = (part.strip() for part in line.split("=", 1))
            qualified = f"{block}/{key}"
            if qualified in values:
                raise AcceptanceError(f"restart parameter is duplicated: {qualified}")
            values[qualified] = value
    return values, end


def parameter_int(parameters: dict[str, str], key: str, *, default: int | None = None) -> int:
    """Return one exact restart integer parameter."""

    value = parameters.get(key)
    if value is None and default is not None:
        return default
    try:
        result = int(str(value))
    except (TypeError, ValueError) as error:
        raise AcceptanceError(f"restart parameter {key} is not an integer") from error
    return result


def parameter_float(parameters: dict[str, str], key: str) -> float:
    """Return one finite restart float parameter."""

    try:
        result = float(parameters[key])
    except (KeyError, ValueError) as error:
        raise AcceptanceError(f"restart parameter {key} is not a float") from error
    if not math.isfinite(result):
        raise AcceptanceError(f"restart parameter {key} is non-finite")
    return result


def pread_exact(descriptor: int, size: int, offset: int, label: str) -> bytes:
    """Read exact bytes at one descriptor offset."""

    payload = os.pread(descriptor, size, offset)
    if len(payload) != size:
        raise AcceptanceError(f"{label} is truncated")
    return payload


def read_real_array(descriptor: int, count: int, offset: int, label: str) -> array:
    """Read one native little-endian binary64 array."""

    payload = pread_exact(descriptor, count * 8, offset, label)
    values = array("d")
    values.frombytes(payload)
    if sys.byteorder != "little":
        values.byteswap()
    if len(values) != count or any(not math.isfinite(value) for value in values):
        raise AcceptanceError(f"{label} is malformed or non-finite")
    return values


def flat_index_3d(k: int, j: int, i: int, nj: int, ni: int) -> int:
    """Return LayoutRight index for one 3-D array."""

    return (k * nj + j) * ni + i


def rank_identity(path: Path) -> int:
    """Return the canonical rank id encoded by one rank-local parent directory."""

    match = RANK_DIRECTORY_PATTERN.fullmatch(path.parent.name)
    if match is None:
        raise AcceptanceError(f"rank-local restart path lacks canonical rank identity: {path}")
    return int(match.group(1))


def validate_region_indices(
    indices: tuple[int, ...],
    label: str,
    *,
    require_three_d: bool,
    validate_coarse: bool,
) -> None:
    """Validate defined AthenaK RegionIndcs active and coarse ranges."""

    if len(indices) != 19:
        raise AcceptanceError(f"{label} RegionIndcs ABI length differs")
    ng, nx1, nx2, nx3, is_, ie, js, je, ks, ke = indices[:10]
    if ng < 0 or min(nx1, nx2, nx3) <= 0:
        raise AcceptanceError(f"{label} dimensions or ghost count are invalid")
    if require_three_d and min(nx1, nx2, nx3) <= 1:
        raise AcceptanceError("CT audit explicitly requires a three-dimensional state")
    if (
        is_ > ie
        or js > je
        or ks > ke
        or ie - is_ + 1 != nx1
        or je - js + 1 != nx2
        or ke - ks + 1 != nx3
        or (is_, js, ks) != (ng, ng, ng)
        or (ie, je, ke) != (ng + nx1 - 1, ng + nx2 - 1, ng + nx3 - 1)
    ):
        raise AcceptanceError(f"{label} active ranges or ghost containment differ")
    if not validate_coarse:
        return
    cnx1, cnx2, cnx3, cis, cie, cjs, cje, cks, cke = indices[10:]
    expected_coarse = (nx1 // 2, nx2 // 2, nx3 // 2)
    if (
        min(cnx1, cnx2, cnx3) <= 0
        or (cnx1, cnx2, cnx3) != expected_coarse
        or (cis, cjs, cks) != (ng, ng, ng)
        or (cie, cje, cke)
        != (ng + cnx1 - 1, ng + cnx2 - 1, ng + cnx3 - 1)
    ):
        raise AcceptanceError(f"{label} coarse ranges or ghost containment differ")


def athenak_rank_assignments(costs: list[float], rank_count: int) -> list[int]:
    """Reproduce AthenaK's contiguous load-balance assignment from retained costs."""

    def f32(value: float) -> float:
        try:
            return struct.unpack("<f", struct.pack("<f", value))[0]
        except OverflowError as error:
            raise AcceptanceError("restart MeshBlock cost accumulation overflowed") from error

    if rank_count <= 0 or len(costs) < rank_count:
        raise AcceptanceError("restart rank count cannot cover retained MeshBlocks")
    if any(not math.isfinite(value) or value <= 0.0 for value in costs):
        raise AcceptanceError("restart MeshBlock costs must be finite and positive")
    total = 0.0
    for value in costs:
        total = f32(total + value)
    target = f32(total / rank_count)
    rank = rank_count - 1
    accumulated = 0.0
    assignments = [0] * len(costs)
    for index in range(len(costs) - 1, -1, -1):
        accumulated = f32(accumulated + costs[index])
        assignments[index] = rank
        if accumulated >= target and rank > 0:
            rank -= 1
            total = f32(total - accumulated)
            accumulated = 0.0
            target = f32(total / (rank + 1))
    if rank != 0 or sorted(set(assignments)) != list(range(rank_count)):
        raise AcceptanceError("restart MeshBlock costs do not assign every canonical rank")
    return assignments


def audit_restart_file(
    path: Path,
    expected_sha256: str,
    threshold: float,
    *,
    expected_rank: int | None = None,
    expected_rank_count: int | None = None,
) -> dict[str, object]:
    """Read one native Stage I rank-local restart and audit normalized CT-divB."""

    expected_sha256 = require_sha256(expected_sha256, "expected restart SHA-256")
    canonical_rank = rank_identity(path)
    if expected_rank is not None and canonical_rank != expected_rank:
        raise AcceptanceError("rank-local restart canonical rank differs from inventory")
    if expected_rank_count is not None and (
        expected_rank_count <= 0 or canonical_rank >= expected_rank_count
    ):
        raise AcceptanceError("rank-local restart rank is outside expected inventory")
    with open_stable_regular(path, "rank-local restart") as descriptor:
        profile = os.fstat(descriptor)
        observed_sha = descriptor_sha256(descriptor)
        if observed_sha != expected_sha256:
            raise AcceptanceError(f"rank-local restart SHA-256 differs: {path}")
        prefix = pread_exact(
            descriptor,
            min(profile.st_size, MAX_PARAMETER_DUMP_BYTES),
            0,
            "restart parameter prefix",
        )
        parameters, parameter_end = parse_parameter_dump(prefix)
        if (
            parameters.get("mhd/eos") != "cgl"
            or parameters.get("output3/single_file_per_rank") != "true"
            or parameter_int(parameters, "mhd/nscalars", default=0) != 0
        ):
            raise AcceptanceError("restart is not a qualified rank-local CGL Stage I state")
        header_bytes = pread_exact(
            descriptor, RESTART_HEADER_SIZE, parameter_end, "restart mesh header"
        )
        unpacked = struct.unpack(RESTART_HEADER_FORMAT, header_bytes)
        nmb_total = unpacked[0]
        root_level = unpacked[1]
        region = unpacked[2:11]
        mesh_indices = unpacked[11:30]
        block_indices = unpacked[30:49]
        time_value = unpacked[49]
        if (
            nmb_total <= 0
            or root_level < 0
            or not math.isfinite(time_value)
            or any(not math.isfinite(value) for value in region)
        ):
            raise AcceptanceError("restart mesh header is invalid")
        validate_region_indices(
            mesh_indices,
            "restart mesh",
            require_three_d=True,
            validate_coarse=False,
        )
        validate_region_indices(
            block_indices,
            "restart MeshBlock",
            require_three_d=True,
            validate_coarse=True,
        )
        if tuple(mesh_indices[:4]) != (
            parameter_int(parameters, "mesh/nghost"),
            parameter_int(parameters, "mesh/nx1"),
            parameter_int(parameters, "mesh/nx2"),
            parameter_int(parameters, "mesh/nx3"),
        ):
            raise AcceptanceError("restart mesh indices differ from its parameter dump")
        ng, nx1, nx2, nx3 = block_indices[:4]
        if (nx1, nx2, nx3) != (
            parameter_int(parameters, "meshblock/nx1"),
            parameter_int(parameters, "meshblock/nx2"),
            parameter_int(parameters, "meshblock/nx3"),
        ):
            raise AcceptanceError("restart MeshBlock indices differ from parameters")
        is_, ie, js, je, ks, ke = block_indices[4:10]
        nout1 = nx1 + 2 * ng
        nout2 = nx2 + 2 * ng
        nout3 = nx3 + 2 * ng
        if min(nout1, nout2, nout3) <= 0:
            raise AcceptanceError("restart output extents are empty")
        cell_count = nout1 * nout2 * nout3
        x1_count = (nout1 + 1) * nout2 * nout3
        x2_count = nout1 * (nout2 + 1) * nout3
        x3_count = nout1 * nout2 * (nout3 + 1)
        expected_data_size = (6 * cell_count + x1_count + x2_count + x3_count) * 8
        if expected_data_size <= 0 or expected_data_size > profile.st_size:
            raise AcceptanceError("restart variable-data size contract is impossible")
        list_start = parameter_end + RESTART_HEADER_SIZE
        list_bytes = nmb_total * (LOGICAL_LOCATION_SIZE + MESH_BLOCK_COST_SIZE)
        list_end = list_start + list_bytes
        if list_end >= profile.st_size:
            raise AcceptanceError("restart logical-location list is truncated")
        locations = pread_exact(
            descriptor,
            nmb_total * LOGICAL_LOCATION_SIZE,
            list_start,
            "restart logical locations",
        )
        global_locations = [
            struct.unpack_from("<4i", locations, index * LOGICAL_LOCATION_SIZE)
            for index in range(nmb_total)
        ]
        if len(set(global_locations)) != len(global_locations):
            raise AcceptanceError("restart global logical locations are duplicated")
        if any(location[3] != root_level for location in global_locations):
            raise AcceptanceError("CT audit currently rejects AMR restart states")
        mesh_nx1, mesh_nx2, mesh_nx3 = mesh_indices[1:4]
        if (
            mesh_nx1 % nx1 != 0
            or mesh_nx2 % nx2 != 0
            or mesh_nx3 % nx3 != 0
        ):
            raise AcceptanceError("restart MeshBlock geometry does not tile the mesh")
        root_counts = (mesh_nx1 // nx1, mesh_nx2 // nx2, mesh_nx3 // nx3)
        expected_root_level = (max(root_counts) - 1).bit_length()
        expected_locations = {
            (lx1, lx2, lx3, root_level)
            for lx3 in range(root_counts[2])
            for lx2 in range(root_counts[1])
            for lx1 in range(root_counts[0])
        }
        if (
            root_level != expected_root_level
            or nmb_total != math.prod(root_counts)
            or set(global_locations) != expected_locations
        ):
            raise AcceptanceError("restart logical locations do not exactly cover the root mesh")
        costs_offset = list_start + nmb_total * LOGICAL_LOCATION_SIZE
        costs_payload = pread_exact(
            descriptor,
            nmb_total * MESH_BLOCK_COST_SIZE,
            costs_offset,
            "restart MeshBlock costs",
        )
        costs = list(struct.unpack(f"<{nmb_total}f", costs_payload))
        if any(not math.isfinite(value) or value <= 0.0 for value in costs):
            raise AcceptanceError("restart MeshBlock costs are malformed")
        search_size = min(MAX_RESTART_STATE_BYTES, profile.st_size - list_end)
        state = pread_exact(descriptor, search_size, list_end, "restart internal state")
        marker = struct.pack(IO_WRAPPER_SIZE_FORMAT, expected_data_size)
        candidates: list[tuple[int, int]] = []
        cursor = state.find(marker)
        while cursor >= 0:
            data_start = list_end + cursor + IO_WRAPPER_SIZE_BYTES
            remaining = profile.st_size - data_start
            if remaining > 0 and remaining % expected_data_size == 0:
                local_blocks = remaining // expected_data_size
                if 0 < local_blocks <= nmb_total:
                    candidates.append((data_start, local_blocks))
            cursor = state.find(marker, cursor + 1)
        if len(candidates) != 1:
            raise AcceptanceError("restart variable-data boundary is ambiguous")
        data_start, local_blocks = candidates[0]
        if data_start + local_blocks * expected_data_size != profile.st_size:
            raise AcceptanceError("restart variable-data size contract differs")
        if expected_rank_count is not None:
            assignments = athenak_rank_assignments(costs, expected_rank_count)
            local_ids = [
                index for index, rank in enumerate(assignments) if rank == canonical_rank
            ]
            expected_local_blocks = len(local_ids)
            if local_blocks != expected_local_blocks:
                raise AcceptanceError("rank-local MeshBlock count differs from canonical partition")
            local_locations = [global_locations[index] for index in local_ids]
        else:
            local_locations = []
        dx1, dx2, dx3 = region[6:9]
        bounds = list(zip(region[:3], region[3:6]))
        expected_dx = [
            (right - left) / cells
            for (left, right), cells in zip(bounds, mesh_indices[1:4])
        ]
        if (
            min(dx1, dx2, dx3) <= 0.0
            or any(right <= left for left, right in bounds)
            or any(
                not math.isclose(observed, expected, rel_tol=1.0e-14, abs_tol=0.0)
                for observed, expected in zip((dx1, dx2, dx3), expected_dx)
            )
        ):
            raise AcceptanceError("restart mesh spacing is invalid")
        bfloor = require_positive(parameter_float(parameters, "mhd/bfloor"), "restart bfloor")
        dx_min = min(dx1, dx2, dx3)
        maximum = 0.0
        cc_bytes = 6 * cell_count * 8
        for block in range(local_blocks):
            block_start = data_start + block * expected_data_size
            x1_offset = block_start + cc_bytes
            x2_offset = x1_offset + x1_count * 8
            x3_offset = x2_offset + x2_count * 8
            x1f = read_real_array(descriptor, x1_count, x1_offset, "x1 face field")
            x2f = read_real_array(descriptor, x2_count, x2_offset, "x2 face field")
            x3f = read_real_array(descriptor, x3_count, x3_offset, "x3 face field")
            for k in range(ks, ke + 1):
                for j in range(js, je + 1):
                    for i in range(is_, ie + 1):
                        x1_left = x1f[flat_index_3d(k, j, i, nout2, nout1 + 1)]
                        x1_right = x1f[flat_index_3d(k, j, i + 1, nout2, nout1 + 1)]
                        x2_left = x2f[flat_index_3d(k, j, i, nout2 + 1, nout1)]
                        x2_right = x2f[flat_index_3d(k, j + 1, i, nout2 + 1, nout1)]
                        x3_left = x3f[flat_index_3d(k, j, i, nout2, nout1)]
                        x3_right = x3f[flat_index_3d(k + 1, j, i, nout2, nout1)]
                        divb = (
                            (x1_right - x1_left) / dx1
                            + (x2_right - x2_left) / dx2
                            + (x3_right - x3_left) / dx3
                        )
                        bx = 0.5 * (x1_left + x1_right)
                        by = 0.5 * (x2_left + x2_right)
                        bz = 0.5 * (x3_left + x3_right)
                        normalized = abs(divb) * dx_min / max(
                            math.sqrt(bx * bx + by * by + bz * bz), bfloor
                        )
                        if not math.isfinite(normalized):
                            raise AcceptanceError("normalized CT-divB is non-finite")
                        maximum = max(maximum, normalized)
        final_sha = descriptor_sha256(descriptor)
        if final_sha != observed_sha or final_sha != expected_sha256:
            raise AcceptanceError("rank-local restart changed during hash/read/hash audit")
        contract = {
            "parameter_sha256": sha256_bytes(canonical_json(parameters)),
            "root_level": root_level,
            "region": list(region),
            # AthenaK does not initialize mesh_indcs coarse fields before writing them.
            "mesh_indices_defined": list(mesh_indices[:10]),
            "block_indices": list(block_indices),
            "dx": [dx1, dx2, dx3],
            "dx_min": dx_min,
            "bfloor": bfloor,
            "data_size_per_meshblock": expected_data_size,
            "meshblock_costs_sha256": sha256_bytes(costs_payload),
        }
        return {
            "path": str(path.expanduser().absolute().resolve(strict=True)),
            "sha256": observed_sha,
            "size_bytes": profile.st_size,
            "rank": canonical_rank,
            "time": time_value,
            "nmb_total": nmb_total,
            "local_meshblocks": local_blocks,
            "data_size_per_meshblock": expected_data_size,
            "global_logical_locations_sha256": sha256_bytes(locations),
            "global_logical_locations": [list(value) for value in global_locations],
            "local_logical_locations": [list(value) for value in local_locations],
            "restart_contract": contract,
            "restart_contract_sha256": sha256_bytes(canonical_json(contract)),
            "maximum_normalized_ct_divb": maximum,
            "passed": maximum < threshold,
        }


def parse_restart_argument(value: str) -> tuple[Path, str]:
    """Parse one ``path=sha256`` restart binding."""

    if "=" not in value:
        raise AcceptanceError("restart arguments must use path=sha256")
    path_text, digest = value.rsplit("=", 1)
    return Path(path_text), require_sha256(digest, "restart argument SHA-256")


def load_bound_json(binding: object, label: str) -> tuple[dict[str, object], dict[str, object]]:
    """Load one JSON object only after authenticating its declared binding."""

    observed = verify_declared_binding(binding, label)
    value, loaded = load_json(Path(str(observed["path"])), label)
    if loaded != observed:
        raise AcceptanceError(f"{label} changed between binding verification and load")
    return value, observed


def validate_ct_inventory(
    policy: dict[str, object],
    case_id: str,
    inventory_path: Path,
) -> tuple[
    dict[str, object],
    dict[str, object],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, object] | None,
]:
    """Authenticate a complete multi-state CT inventory against accepted inspections."""

    inventory, inventory_binding = load_json(inventory_path, "restart CT inventory")
    if (
        inventory.get("schema_version") != 2
        or inventory.get("record_type") != "stage-i-restart-ct-inventory"
        or inventory.get("case_id") != case_id
        or inventory.get("coverage_complete") is not True
    ):
        raise AcceptanceError("restart CT inventory type or case differs")
    ct_policy = require_dict(policy["criteria"].get("ct_divb_policy"), "ct_divb_policy")
    required_times = [
        require_finite(value, "required CT state time")
        for value in require_list(
            ct_policy.get("required_state_times"), "required CT state times"
        )
    ]
    declared_times = [
        require_finite(value, "inventory required state time")
        for value in require_list(
            inventory.get("required_state_times"), "inventory required state times"
        )
    ]
    if declared_times != required_times:
        raise AcceptanceError("restart CT inventory required state times differ")
    executable_sha = require_sha256(
        inventory.get("executable_sha256"), "inventory executable sha256"
    )
    bundle, bundle_binding = load_bound_json(
        inventory.get("accepted_bundle_manifest"), "inventory accepted bundle manifest"
    )
    canonical_bundle = path_is_exact(
        Path(str(bundle_binding["path"])), canonical_bundle_path(case_id)
    )
    if (
        bundle.get("workflow") != "paper-mks24-stage-i-production"
        or bundle.get("status") != "accepted_for_analysis"
        or bundle.get("production_case_id") != case_id
        or require_finite(bundle.get("required_final_time"), "bundle required final time")
        != 10.0
        or require_finite(bundle.get("accepted_final_time"), "bundle accepted final time")
        != 10.0
    ):
        raise AcceptanceError("inventory accepted bundle identity differs")
    selected_cases = require_list(bundle.get("cases"), "inventory accepted bundle cases")
    expected_case = case_manifest_record(policy, case_id)
    if len(selected_cases) != 1:
        raise AcceptanceError("inventory accepted bundle must contain exactly one case")
    selected_case = require_dict(selected_cases[0], "inventory accepted bundle case")
    if (
        selected_case.get("name") != expected_case.get("name")
        or selected_case.get("input") != expected_case.get("input")
        or selected_case.get("status") != "passed"
    ):
        raise AcceptanceError("inventory accepted bundle selected case differs")
    bundle_segments = {
        resolve_bound_path(value).resolve(strict=True)
        for value in require_list(
            bundle.get("production_segment_manifests"), "bundle segment manifests"
        )
    }
    states = require_list(inventory.get("states"), "inventory states")
    if len(states) != len(required_times):
        raise AcceptanceError("restart CT inventory state count differs")
    observed_times: list[float] = []
    state_specs: list[dict[str, object]] = []
    retained_bindings: list[dict[str, object]] = [
        inventory_binding,
        bundle_binding,
    ]
    all_paths: set[Path] = set()
    all_digests: set[str] = set()
    for state_index, state_value in enumerate(states):
        state = require_dict(state_value, f"inventory state {state_index}")
        state_time = require_finite(state.get("time"), "inventory state time")
        observed_times.append(state_time)
        segment, segment_binding = load_bound_json(
            state.get("accepted_segment_manifest"),
            f"inventory state {state_time} accepted segment manifest",
        )
        segment_path = Path(str(segment_binding["path"]))
        if segment_path not in bundle_segments:
            raise AcceptanceError("inventory state segment is absent from accepted bundle lineage")
        accounting = require_dict(segment.get("accounting"), "inventory segment accounting")
        command = require_dict(segment.get("command"), "inventory segment command")
        inspection = require_dict(
            segment.get("scientific_inspection"), "inventory segment inspection"
        )
        if (
            accounting.get("result") != "accepted"
            or accounting.get("case_id") != case_id
            or accounting.get("case_name") != expected_case.get("name")
            or accounting.get("executable_sha256") != executable_sha
            or inspection.get("accepted") is not True
            or inspection.get("case_id") != case_id
            or inspection.get("schema_version") != 3
            or inspection.get("segment") != accounting.get("segment")
            or resolve_bound_path(inspection.get("manifest")).resolve(strict=True)
            != segment_path
            or require_finite(inspection.get("final_time"), "inspection final time")
            != state_time
            or require_finite(inspection.get("required_time"), "inspection required time")
            != state_time
            or require_finite(
                inspection.get("terminal_restart_time"), "inspection terminal restart time"
            )
            != state_time
            or command.get("executable_sha256") != executable_sha
            or command.get("matrix_sha256")
            != policy["verified_sources"]["stage_i_manifest"]["sha256"]
        ):
            raise AcceptanceError("inventory state accepted segment semantics differ")
        terminal = require_dict(
            inspection.get("terminal_restart"), "inspection terminal restart"
        )
        checks = require_dict(inspection.get("checks"), "inspection checks")
        if (
            terminal.get("storage") != "per_rank"
            or checks.get("required_time_reached") is not True
            or checks.get("restart_retained") is not True
            or checks.get("terminal_restart_physical_time_matches_final") is not True
            or state_time not in [
                require_finite(value, "inspection restart time")
                for value in require_list(
                    inspection.get("restart_times"), "inspection restart times"
                )
            ]
            or terminal not in require_list(
                inspection.get("restarts"), "inspection retained restarts"
            )
        ):
            raise AcceptanceError("inventory state terminal inspection contract differs")
        inspected_ranks = require_list(
            terminal.get("rank_files"), "inspection terminal restart rank files"
        )
        ranks = require_list(state.get("rank_files"), "inventory state rank files")
        allocation = require_dict(segment.get("allocation"), "inventory segment allocation")
        expected_rank_count = (
            require_int(allocation.get("nodes"), "segment nodes", minimum=1)
            * require_int(
                allocation.get("ranks_per_node"), "segment ranks per node", minimum=1
            )
        )
        if len(ranks) != expected_rank_count:
            raise AcceptanceError("inventory state rank count differs from accepted allocation")
        normalized: list[dict[str, object]] = []
        for rank_value in ranks:
            rank = require_dict(rank_value, "inventory rank file")
            rank_id = require_int(rank.get("rank"), "inventory rank", minimum=0)
            observed = verify_declared_binding(rank, f"inventory rank {rank_id}")
            path = Path(str(observed["path"]))
            if rank_identity(path) != rank_id:
                raise AcceptanceError("inventory rank identity differs from canonical path")
            if canonical_bundle:
                expected_rank_root = (
                    segment_path.parents[1]
                    / "output"
                    / "rst"
                    / f"rank_{rank_id:08d}"
                ).resolve(strict=True)
                try:
                    path.relative_to(expected_rank_root)
                except ValueError as error:
                    raise AcceptanceError(
                        "canonical inventory rank file is outside its accepted segment output"
                    ) from error
            if path in all_paths or str(observed["sha256"]) in all_digests:
                raise AcceptanceError("inventory contains duplicated path or copied rank bytes")
            all_paths.add(path)
            all_digests.add(str(observed["sha256"]))
            normalized.append({
                **observed,
                "rank": rank_id,
            })
        if sorted(int(record["rank"]) for record in normalized) != list(
            range(expected_rank_count)
        ):
            raise AcceptanceError("inventory state canonical ranks are not exact and contiguous")
        inspected = {
            (
                str(resolve_bound_path(require_dict(item, "inspected rank").get("path")).resolve(strict=True)),
                require_sha256(
                    require_dict(item, "inspected rank").get("sha256"),
                    "inspected rank sha256",
                ),
                require_int(
                    require_dict(item, "inspected rank").get("size_bytes"),
                    "inspected rank size",
                    minimum=1,
                ),
            )
            for item in inspected_ranks
        }
        declared = {
            (str(record["path"]), str(record["sha256"]), int(record["size_bytes"]))
            for record in normalized
        }
        if declared != inspected:
            raise AcceptanceError(
                "inventory state rank files differ from accepted terminal inspection"
            )
        rank_zero = next(record for record in normalized if record["rank"] == 0)
        if (
            resolve_bound_path(terminal.get("path")).resolve(strict=True)
            != Path(str(rank_zero["path"]))
            or require_sha256(terminal.get("sha256"), "terminal restart sha256")
            != rank_zero["sha256"]
            or require_int(
                terminal.get("size_bytes"), "terminal restart size", minimum=1
            )
            != rank_zero["size_bytes"]
        ):
            raise AcceptanceError("terminal restart summary differs from canonical rank zero")
        retained_bindings.append(segment_binding)
        retained_bindings.extend(normalized)
        state_specs.append({
            "time": state_time,
            "accepted_segment_manifest": segment_binding,
            "rank_count": expected_rank_count,
            "rank_files": normalized,
        })
    if observed_times != required_times or len(set(observed_times)) != len(observed_times):
        raise AcceptanceError("restart CT inventory does not enumerate exact required states")
    canonical_context: dict[str, object] | None = None
    if canonical_bundle:
        ordered_bundle_segments = [
            resolve_bound_path(value).resolve(strict=True)
            for value in require_list(
                bundle.get("production_segment_manifests"), "bundle segment manifests"
            )
        ]
        records = [
            load_json(path, f"canonical CT whole-case segment {index}")
            for index, path in enumerate(ordered_bundle_segments)
        ]
        canonical_context, canonical_inputs = authenticate_canonical_accounting(
            policy,
            case_id,
            str(expected_case["name"]),
            ordered_bundle_segments,
            records,
        )
        retained_bindings.extend(binding for _, binding in records)
        retained_bindings.extend(canonical_inputs)
    return (
        inventory,
        inventory_binding,
        state_specs,
        retained_bindings,
        canonical_context,
    )


def validate_ct_state_records(
    records: list[dict[str, object]],
    expected_time: float,
    *,
    require_complete: bool,
) -> dict[str, object]:
    """Validate sibling consistency and logical-MeshBlock coverage for one state."""

    if not records:
        raise AcceptanceError("CT state has no restart records")
    if any(float(record["time"]) != expected_time for record in records):
        raise AcceptanceError("rank-local restart state time differs from inventory")
    ranks = [int(record["rank"]) for record in records]
    if len(ranks) != len(set(ranks)):
        raise AcceptanceError("rank-local restart state has duplicate canonical ranks")
    if require_complete and sorted(ranks) != list(range(len(records))):
        raise AcceptanceError("rank-local restart state canonical ranks are not exact")
    if len({str(record["sha256"]) for record in records}) != len(records):
        raise AcceptanceError("rank-local restart state contains copied duplicate bytes")
    for key, label in (
        ("nmb_total", "global MeshBlock count"),
        ("global_logical_locations_sha256", "global logical locations"),
        ("restart_contract_sha256", "restart ABI/geometry contract"),
    ):
        if len({str(record[key]) for record in records}) != 1:
            raise AcceptanceError(f"rank-local restart siblings differ in {label}")
    nmb_total = int(records[0]["nmb_total"])
    global_locations = [
        tuple(int(value) for value in location)
        for location in require_list(
            records[0].get("global_logical_locations"), "global logical locations"
        )
    ]
    if len(global_locations) != nmb_total or len(set(global_locations)) != nmb_total:
        raise AcceptanceError("global logical MeshBlock identity set is incomplete")
    local_locations = [
        tuple(int(value) for value in location)
        for record in records
        for location in require_list(
            record.get("local_logical_locations"), "local logical locations"
        )
    ]
    complete = (
        len(local_locations) == nmb_total
        and len(set(local_locations)) == nmb_total
        and set(local_locations) == set(global_locations)
    )
    if require_complete and not complete:
        raise AcceptanceError("rank-local restart siblings do not uniquely cover every MeshBlock")
    maximum = max(float(record["maximum_normalized_ct_divb"]) for record in records)
    sampled_meshblocks = (
        len(local_locations)
        if local_locations
        else sum(int(record["local_meshblocks"]) for record in records)
    )
    return {
        "time": expected_time,
        "rank_count": len(records),
        "state_meshblocks": nmb_total,
        "sampled_meshblocks": sampled_meshblocks,
        "meshblock_coverage_complete": complete,
        "global_logical_locations_sha256": records[0]["global_logical_locations_sha256"],
        "restart_contract_sha256": records[0]["restart_contract_sha256"],
        "maximum_normalized_ct_divb": maximum,
        "restart_audits": records,
    }


def audit_ct_divb(
    policy: dict[str, object],
    case_id: str,
    restart_values: list[str],
    inventory_path: Path | None,
) -> dict[str, object]:
    """Audit sampled retained-state CT-divB from native rank-local restarts."""

    case_id = require_case_id(case_id)
    ct_policy = require_dict(policy["criteria"].get("ct_divb_policy"), "ct_divb_policy")
    threshold = require_positive(ct_policy.get("normalized_ct_divb_lt"), "CT-divB threshold")
    supplied = [parse_restart_argument(value) for value in restart_values]
    if len({str(path.expanduser().absolute()) for path, _ in supplied}) != len(supplied):
        raise AcceptanceError("CT-divB restart paths are duplicated")
    inventory_binding: dict[str, object] | None = None
    retained_inputs: list[dict[str, object]] = [
        policy["criteria_binding"],
        policy["review_binding"],
    ]
    state_audits: list[dict[str, object]] = []
    coverage_complete = False
    canonical_context: dict[str, object] | None = None
    if inventory_path is not None:
        (
            _,
            inventory_binding,
            state_specs,
            inventory_inputs,
            canonical_context,
        ) = validate_ct_inventory(policy, case_id, inventory_path)
        declared = {
            (str(rank["path"]), str(rank["sha256"]))
            for state in state_specs
            for rank in require_list(state["rank_files"], "state rank files")
        }
        supplied_set = {
            (str(path.expanduser().absolute().resolve(strict=True)), digest)
            for path, digest in supplied
        }
        if supplied and declared != supplied_set:
            raise AcceptanceError("supplied restarts differ from the bound CT inventory")
        retained_inputs.extend(inventory_inputs)
        for state in state_specs:
            rank_files = require_list(state["rank_files"], "state rank files")
            records = [
                audit_restart_file(
                    Path(str(require_dict(rank, "state rank")["path"])),
                    str(require_dict(rank, "state rank")["sha256"]),
                    threshold,
                    expected_rank=int(require_dict(rank, "state rank")["rank"]),
                    expected_rank_count=int(state["rank_count"]),
                )
                for rank in rank_files
            ]
            state_audits.append(validate_ct_state_records(
                records, float(state["time"]), require_complete=True
            ))
        coverage_complete = True
    else:
        if not supplied:
            raise AcceptanceError("CT-divB audit requires restarts or an explicit inventory")
        records = [
            audit_restart_file(path, digest, threshold) for path, digest in supplied
        ]
        retained_inputs.extend(records)
        by_time: dict[float, list[dict[str, object]]] = {}
        for record in records:
            by_time.setdefault(float(record["time"]), []).append(record)
        state_audits = [
            validate_ct_state_records(group, time, require_complete=False)
            for time, group in sorted(by_time.items())
        ]
    maximum = max(
        float(state["maximum_normalized_ct_divb"]) for state in state_audits
    )
    numeric_pass = maximum < threshold
    campaign_authority_eligible = canonical_context is not None
    result = "pass" if numeric_pass and coverage_complete and campaign_authority_eligible else (
        "fail" if not numeric_pass else "inconclusive"
    )
    restart_count = sum(int(state["rank_count"]) for state in state_audits)
    sampled_meshblocks = sum(int(state["sampled_meshblocks"]) for state in state_audits)
    state_meshblocks = sum(int(state["state_meshblocks"]) for state in state_audits)
    meshblock_coverage_complete = all(
        bool(state["meshblock_coverage_complete"]) for state in state_audits
    )
    evidence = {
        "schema_version": 1,
        "record_type": "stage-i-scientific-ct-divb-evidence",
        "authority": "non-authorizing-read-only-sampled-state-audit",
        "non_authorizing_statement": NON_AUTHORIZING_STATEMENT,
        "case_id": case_id,
        "result": result,
        "reason": (
            "every required inventory-bound accepted retained state is below the CT-divB threshold"
            if result == "pass"
            else "one or more sampled retained states exceed the CT-divB threshold"
            if result == "fail"
            else (
                "numeric CT-divB passed but exact canonical campaign authority or required "
                "multi-state inventory is unbound"
            )
        ),
        "claim_scope": ct_policy["claim_scope"],
        "historical_limitation": ct_policy["historical_limitation"],
        "coverage_complete": coverage_complete,
        "numerical_result": "pass" if numeric_pass else "fail",
        "required_state_times": ct_policy["required_state_times"],
        "state_count": len(state_audits),
        "restart_count": restart_count,
        "sampled_meshblocks": sampled_meshblocks,
        "state_meshblocks": state_meshblocks,
        "sampled_meshblock_fraction": sampled_meshblocks / state_meshblocks,
        "meshblock_coverage_complete": meshblock_coverage_complete,
        "campaign_authority_eligible": campaign_authority_eligible,
        "canonical_campaign_context": canonical_context,
        "maximum_normalized_ct_divb": maximum,
        "normalized_ct_divb_lt": threshold,
        "state_audits": state_audits,
        "release_authorizing": False,
        "evaluation_inputs": {
            "inventory": inventory_binding,
            "restarts": [
                {
                    "path": str(path.expanduser().absolute().resolve(strict=True)),
                    "sha256": digest,
                }
                for path, digest in supplied
            ],
        },
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "acceptance_utility": regular_file_binding(Path(__file__), "acceptance utility"),
            "inputs": retained_inputs,
        },
    }
    return seal_evidence(evidence)


def independently_recompute_ct_evidence(
    policy: dict[str, object], value: dict[str, object], label: str
) -> dict[str, object]:
    """Re-run one CT audit from its exact bound inputs and compare semantics."""

    inputs = require_dict(value.get("evaluation_inputs"), f"{label} evaluation inputs")
    inventory_binding = inputs.get("inventory")
    inventory_path = (
        path_from_evidence_binding(inventory_binding, f"{label} inventory")
        if isinstance(inventory_binding, dict)
        else None
    )
    restart_values = [
        f"{path_from_evidence_binding(item, f'{label} restart {index}')}="
        f"{require_sha256(require_dict(item, f'{label} restart {index}').get('sha256'), f'{label} restart sha256')}"
        for index, item in enumerate(
            require_list(inputs.get("restarts"), f"{label} restarts")
        )
    ]
    recomputed = audit_ct_divb(
        policy,
        require_case_id(value.get("case_id"), f"{label} case_id"),
        restart_values,
        inventory_path,
    )
    require_same_evidence_semantics(value, recomputed, label)
    return recomputed


def verify_input_bindings(value: dict[str, object]) -> int:
    """Reverify every retained input binding present in evidence provenance."""

    provenance = require_dict(value.get("provenance"), "evidence provenance")
    collections: list[object] = []
    for key in ("inputs", "case_evidence"):
        item = provenance.get(key)
        if isinstance(item, list):
            collections.extend(item)
    count = 0
    for index, binding in enumerate(collections):
        verify_declared_binding(binding, f"evidence input binding {index}")
        count += 1
    return count


def verify_evidence_policy_binding(
    policy: dict[str, object],
    value: dict[str, object],
    label: str,
    *,
    reverify_inputs: bool = True,
) -> int:
    """Require evidence to bind the current policy, utility, and retained inputs."""

    verify_evidence_digest(value, label)
    provenance = require_dict(value.get("provenance"), f"{label} provenance")
    if require_dict(provenance.get("criteria"), f"{label} criteria").get(
        "sha256"
    ) != policy["criteria_binding"]["sha256"]:
        raise AcceptanceError(f"{label} binds different criteria")
    if require_dict(
        provenance.get("criteria_review"), f"{label} criteria review"
    ).get("sha256") != policy["review_binding"]["sha256"]:
        raise AcceptanceError(f"{label} binds a different criteria review")
    utility = verify_declared_binding(
        provenance.get("acceptance_utility"), f"{label} acceptance utility"
    )
    if Path(str(utility["path"])) != Path(__file__).resolve(strict=True):
        raise AcceptanceError(f"{label} binds a different acceptance utility")
    expected_utility = require_dict(
        policy["verified_sources"]["acceptance_utility"], "policy acceptance utility"
    )
    if utility["sha256"] != expected_utility["sha256"]:
        raise AcceptanceError(f"{label} acceptance utility differs from criteria binding")
    return verify_input_bindings(value) if reverify_inputs else 0


def verify_evidence(
    policy: dict[str, object], path: Path, expected_sha256: str
) -> dict[str, object]:
    """Verify one exact candidate and all retained input bindings it declares."""

    expected_sha256 = require_sha256(expected_sha256, "expected evidence SHA-256")
    value, binding = load_json(path, "scientific evidence")
    if binding["sha256"] != expected_sha256:
        raise AcceptanceError("scientific evidence file SHA-256 differs")
    inputs = verify_evidence_policy_binding(policy, value, "scientific evidence")
    record_type = value.get("record_type")
    if record_type == "stage-i-scientific-case-evidence":
        independently_recompute_case(policy, value, "scientific evidence")
    elif record_type == "stage-i-scientific-campaign-evidence":
        independently_recompute_campaign(policy, value, "scientific evidence")
    elif record_type == "stage-i-scientific-ct-divb-evidence":
        independently_recompute_ct_evidence(policy, value, "scientific evidence")
    elif record_type == "stage-i-scientific-criteria-validation":
        require_same_evidence_semantics(
            value, validate_criteria_evidence(policy), "scientific evidence"
        )
    else:
        raise AcceptanceError("scientific evidence record type cannot be re-evaluated")
    provenance = require_dict(value.get("provenance"), "scientific evidence provenance")
    utility = require_dict(
        provenance.get("acceptance_utility"), "evidence acceptance utility"
    )
    return seal_evidence({
        "schema_version": 1,
        "record_type": "stage-i-scientific-evidence-verification",
        "authority": "non-authorizing-verification",
        "non_authorizing_statement": NON_AUTHORIZING_STATEMENT,
        "verified": True,
        "release_authorizing": False,
        "evidence": binding,
        "evidence_record_type": value.get("record_type"),
        "evidence_result": value.get("result"),
        "reverified_input_bindings": inputs,
        "criteria_review_status": policy["review_status"],
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "acceptance_utility": utility,
            "inputs": [binding],
        },
    })


def validate_criteria_evidence(policy: dict[str, object]) -> dict[str, object]:
    """Return non-authorizing criteria validation evidence."""

    return seal_evidence({
        "schema_version": 1,
        "record_type": "stage-i-scientific-criteria-validation",
        "authority": "non-authorizing-policy-validation",
        "non_authorizing_statement": NON_AUTHORIZING_STATEMENT,
        "valid": True,
        "release_authorizing": False,
        "criteria_review_status": policy["review_status"],
        "independent_review_complete": bool(policy["approved"]),
        "required_cases": policy["criteria"]["required_cases"],
        "admitted_panel_count": len(policy["criteria"]["comparison_panels"]),
        "provenance": {
            "criteria": policy["criteria_binding"],
            "criteria_review": policy["review_binding"],
            "acceptance_utility": regular_file_binding(Path(__file__), "acceptance utility"),
            "inputs": [
                policy["criteria_binding"],
                policy["review_binding"],
                *policy["verified_sources"].values(),
            ],
        },
    })


def add_common_policy_arguments(parser: argparse.ArgumentParser) -> None:
    """Add common criteria/review/output arguments."""

    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)
    parser.add_argument("--output", type=Path)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser("validate-criteria")
    add_common_policy_arguments(validate)

    case = subparsers.add_parser("evaluate-case")
    add_common_policy_arguments(case)
    case.add_argument("--case-id", required=True)
    case.add_argument("--mhd-history", type=Path, required=True)
    case.add_argument("--user-history", type=Path, required=True)
    case.add_argument(
        "--diagnostics",
        action="append",
        help="optional reviewed scientific products as full=PATH, early=PATH, or late=PATH",
    )
    case.add_argument("--ct-evidence", type=Path)
    case.add_argument("--bundle-manifest", type=Path)

    campaign = subparsers.add_parser("evaluate-campaign")
    add_common_policy_arguments(campaign)
    campaign.add_argument("--case-evidence", type=Path, action="append", required=True)

    ct = subparsers.add_parser("audit-ct-divb")
    add_common_policy_arguments(ct)
    ct.add_argument("--case-id", required=True)
    ct.add_argument(
        "--restart",
        action="append",
        help="rank-local native restart binding as PATH=SHA256",
    )
    ct.add_argument("--inventory", type=Path)

    verify = subparsers.add_parser("verify-evidence")
    add_common_policy_arguments(verify)
    verify.add_argument("--evidence", type=Path, required=True)
    verify.add_argument("--expected-evidence-sha256", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run one non-authorizing scientific acceptance command."""

    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        policy = load_validated_policy(args.criteria, args.criteria_review)
        if args.command == "validate-criteria":
            evidence = validate_criteria_evidence(policy)
        elif args.command == "evaluate-case":
            evidence = evaluate_case(
                policy,
                args.case_id,
                args.mhd_history,
                args.user_history,
                args.diagnostics,
                args.ct_evidence,
                args.bundle_manifest,
            )
        elif args.command == "evaluate-campaign":
            evidence = evaluate_campaign(policy, args.case_evidence)
        elif args.command == "audit-ct-divb":
            evidence = audit_ct_divb(
                policy, args.case_id, args.restart or [], args.inventory
            )
        elif args.command == "verify-evidence":
            evidence = verify_evidence(
                policy, args.evidence, args.expected_evidence_sha256
            )
        else:
            raise AcceptanceError(f"unsupported command: {args.command}")
        write_candidate(args.output, evidence)
    except (AcceptanceError, OSError) as error:
        parser.error(str(error))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
