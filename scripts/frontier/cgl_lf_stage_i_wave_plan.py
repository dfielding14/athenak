#!/usr/bin/env python3
"""Emit an evidence-bound, deterministic, read-only CGL-LF Stage I wave plan.

The planner never launches work and never mutates campaign state.  It accepts
only a fresh review packet retained at the exact canonical campaign root.  The
packet checksum-binds the frozen matrix, promoted controller/helper/source/build
provenance, controller summary, a controller-reconciliation snapshot, canonical
ledger/reservation/manifest/transaction stores, independently published
allocation profiles, any required exact published F-115 R03 continuation
authority chain, live Slurm state for submitted jobs, and evidence-backed R17
readiness artifacts.
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
GIT = Path("/usr/bin/git")
ACCOUNT = "AST207"

CONTROLLER_REVISION = "5834a91e448a69ec0df5d011b7be3fe786666806"
CONTROLLER_SHA256 = "0a6b30a60ba70faf32dae722c4472e96acb38031d8c1a4b52b816fc93e58e53a"
SOURCE_REVISION = "9e07542281e4e6d125582f253df3ad2e3b8b154d"
MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
SOURCE_BUNDLE_SHA256 = "a6aa40f8f3350d65022be6d898d5575a60185b05bba7b72f184c3f2ebd401ec8"
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
STATE_MAX_AGE = timedelta(minutes=20)
R17_READINESS_MAX_AGE = timedelta(hours=24)
FUTURE_SKEW = timedelta(minutes=5)
R17_REQUIRED_RETENTION_BYTES = 958271710272
MAX_RESTART_PARAMETER_DUMP_BYTES = 11 * 4096 + 1
RESTART_MESH_HEADER_SIZE = 252
RESTART_TIME_OFFSET = 232
RESTART_TIME_FORMAT = "<d"
ALLOWED_RESTART_MARKER_MODES = frozenset(("full_precision", "legacy_default_precision"))
ACTIVE_SCHEDULER_STATES = frozenset(
    ("PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "SUSPENDED")
)

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
F115_SOURCE_BUNDLE_REQUIRED_REVISIONS = (
    SOURCE_REVISION,
    R03_F114_CONTROLLER_REVISION,
    "b0d3e8d526d8f3b4e5977333000db24c09d4eab1",
    "38aedd2a65c3c11855721858f5dbcce20bae11e4",
    CONTROLLER_REVISION,
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
ALLOCATION_AUTHORITY_RELATIVE = Path(
    "accounting/mks24_stage_i_E03_forcing_policy_independent_wave_allocation_authority.json"
)

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
JOB_RE = re.compile(r"^[1-9][0-9]*$")
DIGITS_RE = re.compile(r"^[0-9]+$")
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


def summary_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_budget_summary.md")


def ledger_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv")


def reservations_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json")


def qualification_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_qualification_approval.json")


def profile_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reviewed_wave_profile.json")


def reconciliation_path() -> Path:
    return expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reconciliation_snapshot.json")


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


def allocation_authority_path() -> Path:
    """Return the exact independently published allocation-authority path."""

    return expected_path(ALLOCATION_AUTHORITY_RELATIVE.as_posix())


def allocation_authority_audit_path() -> Path:
    """Return the exact allocation-authority publication-audit path."""

    return expected_path(f"{ALLOCATION_AUTHORITY_RELATIVE.as_posix()}.publication_audit.json")


def source_bundle_path() -> Path:
    return expected_path("source-archives/athenak-feature-cgl-through-5834a91e4.bundle")


def controller_repository_path() -> Path:
    """Return the exact repository containing the promoted controller helper."""

    return normalized_path(CONTROLLER_HELPER.parents[2])


def git_read_only(
    repository: Path,
    arguments: list[str],
    *,
    descriptor: int | None = None,
) -> subprocess.CompletedProcess:
    """Run one bounded, environment-sanitized, non-shell Git query."""

    environment = {
        key: value for key, value in os.environ.items() if not key.startswith("GIT_")
    }
    environment.update(
        {
            "GIT_CONFIG_GLOBAL": "/dev/null",
            "GIT_CONFIG_NOSYSTEM": "1",
            "GIT_OPTIONAL_LOCKS": "0",
            "LC_ALL": "C",
        }
    )
    try:
        return subprocess.run(
            [
                str(GIT), "--no-replace-objects", "-C", str(repository), *arguments
            ],
            check=False,
            stdin=subprocess.DEVNULL,
            capture_output=True,
            env=environment,
            pass_fds=(() if descriptor is None else (descriptor,)),
            timeout=60,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ValueError("read-only Git source-bundle validation failed to execute") from error


def validate_source_bundle_coverage(
    expected_evidence: dict[str, object] | None = None,
) -> dict[str, object]:
    """Independently require a usable self-contained bundle covering F-115 revisions."""

    path = source_bundle_path()
    payload, evidence = read_stable_regular_file(path, "promoted source bundle")
    if evidence["sha256"] != SOURCE_BUNDLE_SHA256:
        raise ValueError("promoted source bundle digest changed")
    if expected_evidence is not None and evidence != expected_evidence:
        raise ValueError("promoted source bundle changed before Git validation")
    try:
        header = payload.split(b"\n\n", 1)[0].decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("promoted source bundle header is not UTF-8") from error
    if not header or header[0] not in {"# v2 git bundle", "# v3 git bundle"}:
        raise ValueError("promoted source bundle is not a Git bundle")
    if any(line.startswith("-") for line in header[1:]):
        raise ValueError("promoted source bundle is not self-contained")
    advertised_header = []
    for line in header[1:]:
        if line.startswith("@"):
            continue
        match = re.fullmatch(r"([0-9a-f]{40}) (.+)", line)
        if match is None:
            raise ValueError("promoted source bundle advertised reference is malformed")
        advertised_header.append((match.group(1), match.group(2)))
    if (CONTROLLER_REVISION, "HEAD") not in advertised_header:
        raise ValueError("promoted source bundle does not advertise the exact controller HEAD")

    components_before = component_identity(path, "promoted source bundle Git validation")
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                break
            blocks.append(block)
        if sha256_bytes(b"".join(blocks)) != SOURCE_BUNDLE_SHA256:
            raise ValueError("promoted source bundle descriptor digest changed")
        os.lseek(descriptor, 0, os.SEEK_SET)
        descriptor_path = f"/proc/self/fd/{descriptor}"
        repository = controller_repository_path()
        verify = git_read_only(
            repository, ["bundle", "verify", descriptor_path], descriptor=descriptor
        )
        if verify.returncode:
            raise ValueError("promoted source bundle fails git bundle verify")
        heads = git_read_only(
            repository, ["bundle", "list-heads", descriptor_path], descriptor=descriptor
        )
        if heads.returncode:
            raise ValueError("promoted source bundle heads cannot be listed")
        try:
            advertised = [
                tuple(line.split(" ", 1))
                for line in heads.stdout.decode("utf-8").splitlines()
                if line
            ]
        except UnicodeDecodeError as error:
            raise ValueError("promoted source bundle heads are not UTF-8") from error
        if advertised != advertised_header:
            raise ValueError("promoted source bundle header and Git-advertised heads differ")
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    stable_fields = (
        "st_dev", "st_ino", "st_mode", "st_uid", "st_nlink", "st_size",
        "st_mtime_ns", "st_ctime_ns",
    )
    if any(getattr(before, field) != getattr(after, field) for field in stable_fields):
        raise ValueError("promoted source bundle changed during Git validation")
    if component_identity(path, "promoted source bundle Git validation") != components_before:
        raise ValueError("promoted source bundle path changed during Git validation")
    _, retained = read_stable_regular_file(path, "promoted source bundle")
    if retained != evidence:
        raise ValueError("promoted source bundle changed after Git validation")

    repository = controller_repository_path()
    for revision in F115_SOURCE_BUNDLE_REQUIRED_REVISIONS:
        if git_read_only(
            repository, ["cat-file", "-e", f"{revision}^{{commit}}"]
        ).returncode:
            raise ValueError(f"F-115 required bundle revision is unknown: {revision}")
        if git_read_only(
            repository,
            ["merge-base", "--is-ancestor", revision, CONTROLLER_REVISION],
        ).returncode:
            raise ValueError(f"F-115 required revision is not covered by bundle HEAD: {revision}")
    try:
        helper_relative = CONTROLLER_HELPER.relative_to(repository).as_posix()
    except ValueError as error:
        raise ValueError("promoted controller helper is outside its repository") from error
    committed_helper = git_read_only(
        repository, ["show", f"{CONTROLLER_REVISION}:{helper_relative}"]
    )
    if committed_helper.returncode or sha256_bytes(committed_helper.stdout) != CONTROLLER_SHA256:
        raise ValueError("bundle HEAD does not bind the promoted controller helper bytes")
    return {
        "evidence": evidence,
        "git_bundle_verify": "passed",
        "self_contained": True,
        "advertised_heads": [
            {"revision": revision, "name": name} for revision, name in advertised
        ],
        "required_revisions": list(F115_SOURCE_BUNDLE_REQUIRED_REVISIONS),
        "committed_controller_helper_sha256": CONTROLLER_SHA256,
    }


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
    profile: dict[str, object],
    matrix_evidence: dict[str, object],
    input_evidence: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Authenticate the promoted helper/source/build and qualification records."""

    provenance = profile["provenance"]
    if not isinstance(provenance, dict):
        raise ValueError("profile provenance must be an object")
    expected_provenance = {
        "controller_revision": CONTROLLER_REVISION,
        "controller_sha256": CONTROLLER_SHA256,
        "source_revision": SOURCE_REVISION,
        "matrix_sha256": MATRIX_SHA256,
        "source_bundle_sha256": SOURCE_BUNDLE_SHA256,
        "executable_revision": SOURCE_REVISION,
        "executable_sha256": EXECUTABLE_SHA256,
        "build_manifest": str(build_manifest_path()),
    }
    if provenance != expected_provenance:
        raise ValueError("profile provenance differs from the promoted production identity")

    bindings = profile["evidence"]
    if not isinstance(bindings, dict):
        raise ValueError("profile evidence must be an object")
    require_exact_keys(
        bindings,
        (
            "matrix", "summary", "reconciliation", "ledger", "reservations",
            "qualification", "controller_helper", "source_bundle", "executable",
            "build_manifest", "inputs", "manifests",
        ),
        "profile evidence",
    )
    validate_binding(bindings["matrix"], matrix_evidence, frozen_matrix_path(), "matrix")

    _, helper_evidence = read_stable_regular_file(CONTROLLER_HELPER, "promoted controller helper")
    if helper_evidence["sha256"] != CONTROLLER_SHA256:
        raise ValueError("promoted controller helper digest changed")
    validate_binding(bindings["controller_helper"], helper_evidence, CONTROLLER_HELPER, "controller helper")

    _, bundle_evidence = read_stable_regular_file(source_bundle_path(), "promoted source bundle")
    if bundle_evidence["sha256"] != SOURCE_BUNDLE_SHA256:
        raise ValueError("promoted source bundle digest changed")
    validate_binding(bindings["source_bundle"], bundle_evidence, source_bundle_path(), "source bundle")
    bundle_validation = validate_source_bundle_coverage(bundle_evidence)

    _, executable_evidence = read_stable_regular_file(executable_path(), "qualified executable")
    if executable_evidence["sha256"] != EXECUTABLE_SHA256:
        raise ValueError("qualified executable digest changed")
    validate_binding(bindings["executable"], executable_evidence, executable_path(), "executable")

    build_bindings = bindings["build_manifest"]
    if not isinstance(build_bindings, dict) or set(build_bindings) != set(BUILD_FILE_SHA256):
        raise ValueError("build-manifest evidence is incomplete")
    build_evidence: dict[str, dict[str, object]] = {}
    build_payloads: dict[str, bytes] = {}
    for filename, digest in BUILD_FILE_SHA256.items():
        path = build_manifest_path() / filename
        payload, retained = read_stable_regular_file(path, f"build manifest {filename}")
        if retained["sha256"] != digest:
            raise ValueError(f"build manifest {filename} digest changed")
        validate_binding(build_bindings[filename], retained, path, f"build manifest {filename}")
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
    validate_binding(
        bindings["qualification"], qualification_evidence, qualification_path(), "qualification"
    )
    if (
        qualification.get("schema_version") != 1
        or qualification.get("execution_epoch") != EXECUTION_EPOCH
        or qualification.get("approved_executable") != str(executable_path())
        or qualification.get("approved_executable_revision") != SOURCE_REVISION
        or qualification.get("approved_executable_sha256") != EXECUTABLE_SHA256
        or qualification.get("build_manifest") != str(build_manifest_path())
    ):
        raise ValueError("qualification approval differs from the promoted build")

    input_bindings = bindings["inputs"]
    if not isinstance(input_bindings, dict) or set(input_bindings) != set(ALL_CASES):
        raise ValueError("profile input evidence differs from the exact R02-R17 set")
    for case_id in ALL_CASES:
        validate_binding(
            input_bindings[case_id], input_evidence[case_id], input_path(case_id),
            f"{case_id} frozen input",
        )
    return {
        "controller_helper": {
            "path": str(CONTROLLER_HELPER),
            "revision": CONTROLLER_REVISION,
            "sha256": CONTROLLER_SHA256,
            "evidence": helper_evidence,
        },
        "source_revision": SOURCE_REVISION,
        "source_bundle": {
            "path": str(source_bundle_path()),
            "sha256": SOURCE_BUNDLE_SHA256,
            "verified_revisions": list(F115_SOURCE_BUNDLE_REQUIRED_REVISIONS),
            "evidence": bundle_evidence,
            "independent_validation": bundle_validation,
        },
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
    try:
        completed = subprocess.run(
            [
                str(SQUEUE), "-j", job_id, "-h",
                "-t", ",".join(sorted(ACTIVE_SCHEDULER_STATES)),
                "-o", "%i|%j|%T|%D|%u|%a",
            ],
            check=True,
            capture_output=True,
            text=True,
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
    try:
        completed = subprocess.run(
            [str(SCONTROL), "write", "batch_script", job_id, "-"],
            check=True,
            capture_output=True,
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
        if utility != {
            "committed": True,
            "path": str(CONTROLLER_HELPER),
            "revision": CONTROLLER_REVISION,
            "sha256": CONTROLLER_SHA256,
        }:
            raise ValueError(f"{case_id}/{segment} active helper provenance is not promoted")
        if not isinstance(bundle, dict) or (
            bundle.get("path") != str(source_bundle_path())
            or bundle.get("sha256") != SOURCE_BUNDLE_SHA256
            or bundle.get("verified_revisions") != [SOURCE_REVISION, CONTROLLER_REVISION]
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


def validated_lineage(
    case_id: str, infos: list[dict[str, object]]
) -> tuple[list[dict[str, object]], dict[str, object] | None]:
    """Return one parent-authenticated lineage while preserving cancelled identities."""

    recorded = [info for info in infos if info["state"] == "recorded"]
    active = [info for info in infos if info["state"] in {"prepared", "submitted"}]
    cancelled = [info for info in infos if info["state"] == "cancelled"]
    if len(active) > 1:
        raise ValueError(f"{case_id} has multiple active manifests")
    index_owners: dict[int, dict[str, object]] = {}
    for info in infos:
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
        if result == "clean_partial" and case_id != R03:
            raise ValueError(f"{case_id} clean_partial requires a dedicated reviewed recovery path")
        final_time = decimal_value(inspection.get("final_time"), "inspection final_time")
        required_time = decimal_value(inspection.get("required_time"), "inspection required_time")
        if required_time != info["target"]:
            raise ValueError(f"{case_id}/{info['segment']} inspection target differs from preparation")
        if result == "accepted" and final_time != info["target"]:
            raise ValueError(f"{case_id}/{info['segment']} accepted endpoint is not exact")
        if result == "clean_partial" and not (info["start"] < final_time < info["target"]):
            raise ValueError(f"{case_id}/{info['segment']} clean_partial endpoint is invalid")
        if prior is None:
            if info["index"] != 0 or info["start"] != 0 or manifest["command"].get("parent_segment") is not None:
                raise ValueError(f"{case_id} recorded lineage does not start fresh at index zero")
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
            int(info["index"]) for info in recorded + cancelled
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


def validate_reconciliation_snapshot(
    snapshot: dict[str, object],
    snapshot_evidence: dict[str, object],
    actual_evidence: dict[str, object],
    counts: dict[str, int],
    observed_utc: str,
) -> None:
    """Require a fresh promoted-controller snapshot bound to exact canonical stores."""

    require_exact_keys(
        snapshot,
        ("schema_version", "record_type", "execution_epoch", "generated_utc", "observed_utc", "generator", "evidence", "report"),
        "reconciliation snapshot",
    )
    if (
        snapshot["schema_version"] != 1
        or snapshot["record_type"] != "cgl_lf_stage_i_controller_reconciliation_snapshot"
        or snapshot["execution_epoch"] != EXECUTION_EPOCH
        or snapshot["observed_utc"] != observed_utc
    ):
        raise ValueError("reconciliation snapshot identity is invalid")
    generated = require_fresh(snapshot_evidence, snapshot["generated_utc"], "reconciliation snapshot")
    if generated < require_utc(observed_utc, "observed state"):
        raise ValueError("reconciliation snapshot predates its observed state")
    if snapshot["generator"] != {
        "path": str(CONTROLLER_HELPER),
        "revision": CONTROLLER_REVISION,
        "sha256": CONTROLLER_SHA256,
    }:
        raise ValueError("reconciliation snapshot generator is not the promoted controller")
    evidence = snapshot["evidence"]
    if not isinstance(evidence, dict) or evidence != actual_evidence:
        raise ValueError("reconciliation snapshot does not bind the exact canonical stores")
    report = snapshot["report"]
    if not isinstance(report, dict):
        raise ValueError("reconciliation report must be an object")
    require_exact_keys(
        report, ("execution_epoch", "root", "qualification", "consistent", "counts", "issues"),
        "reconciliation report",
    )
    if (
        report["execution_epoch"] != EXECUTION_EPOCH
        or report["root"] != str(CANONICAL_ROOT)
        or report["consistent"] is not True
        or report["issues"] != []
        or report["counts"] != counts
    ):
        raise ValueError("controller reconciliation report is not clean/current")
    qualification = report["qualification"]
    if not isinstance(qualification, dict) or (
        qualification.get("state") != "approved"
        or qualification.get("path") != str(qualification_path())
        or qualification.get("sha256") != actual_evidence["qualification"]["sha256"]
        or qualification.get("approved_executable_revision") != SOURCE_REVISION
        or qualification.get("approved_executable_sha256") != EXECUTABLE_SHA256
    ):
        raise ValueError("reconciliation qualification does not match the promoted build")


def validate_publication_audit(
    value: dict[str, object],
    artifact_path: Path,
    artifact_sha256: str,
    *,
    record_type: str,
    label: str,
    profile_reviewer: str | None = None,
) -> None:
    """Require one independently reviewed exact-path publication audit."""

    require_exact_keys(
        value,
        ("schema_version", "record_type", "execution_epoch", "published_utc", "artifact", "review"),
        label,
    )
    published = require_utc(value["published_utc"], f"{label} published_utc")
    if (
        value["schema_version"] != 1
        or value["record_type"] != record_type
        or value["execution_epoch"] != EXECUTION_EPOCH
        or published > current_utc() + FUTURE_SKEW
    ):
        raise ValueError(f"{label} identity is invalid")
    artifact = value["artifact"]
    if artifact != {
        "path": str(artifact_path),
        "sha256": artifact_sha256,
        "mode": "0644",
    }:
        raise ValueError(f"{label} artifact binding differs")
    review = value["review"]
    if not isinstance(review, dict):
        raise ValueError(f"{label} review must be an object")
    require_exact_keys(
        review, ("status", "reviewed_by", "independent_from_profile_author"), f"{label} review"
    )
    reviewer = require_nonempty_string(review["reviewed_by"], f"{label} reviewer")
    if (
        review["status"] != "approved"
        or review["independent_from_profile_author"] is not True
        or (profile_reviewer is not None and reviewer == profile_reviewer)
    ):
        raise ValueError(f"{label} lacks independent approval")


def validate_allocation_authority(
    profile: dict[str, object],
    profile_evidence: dict[str, object],
    profiles: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Require profiles to come from a separately published exact authority."""

    binding = profile["allocation_authority"]
    if not isinstance(binding, dict):
        raise ValueError("allocation authority binding must be an object")
    require_exact_keys(binding, ("artifact", "publication_audit"), "allocation authority binding")
    artifact_path = allocation_authority_path()
    artifact, artifact_evidence = read_json_file(artifact_path, "allocation authority")
    validate_binding(binding["artifact"], artifact_evidence, artifact_path, "allocation authority")
    require_exact_keys(
        artifact,
        (
            "schema_version", "record_type", "execution_epoch", "canonical_root",
            "generated_utc", "expires_utc", "matrix_sha256", "controller_revision",
            "controller_sha256", "profiles_sha256", "profiles", "projection",
        ),
        "allocation authority",
    )
    generated = require_fresh(artifact_evidence, artifact["generated_utc"], "allocation authority")
    expires = require_utc(artifact["expires_utc"], "allocation authority expiry")
    ordered_profiles = [profiles[case_id] for case_id in sorted(profiles)]
    serialized_profiles = [
        {**value, "next_increment": decimal_text(value["next_increment"])}
        for value in ordered_profiles
    ]
    if (
        artifact["schema_version"] != 1
        or artifact["record_type"] != "cgl_lf_stage_i_independent_wave_allocation_authority"
        or artifact["execution_epoch"] != EXECUTION_EPOCH
        or artifact["canonical_root"] != str(CANONICAL_ROOT)
        or artifact["matrix_sha256"] != MATRIX_SHA256
        or artifact["controller_revision"] != CONTROLLER_REVISION
        or artifact["controller_sha256"] != CONTROLLER_SHA256
        or expires <= generated
        or expires < current_utc()
        or artifact["profiles"] != serialized_profiles
        or artifact["profiles_sha256"] != sha256_bytes(canonical_json(serialized_profiles))
    ):
        raise ValueError("allocation profiles differ from independent published authority")
    audit_path = allocation_authority_audit_path()
    audit, audit_evidence = read_json_file(audit_path, "allocation authority publication audit")
    validate_binding(
        binding["publication_audit"], audit_evidence, audit_path,
        "allocation authority publication audit",
    )
    validate_publication_audit(
        audit,
        artifact_path,
        artifact_evidence["sha256"],
        record_type="cgl_lf_stage_i_wave_allocation_authority_publication_audit",
        label="allocation authority publication audit",
        profile_reviewer=str(profile["review"]["reviewed_by"]),
    )
    if require_utc(audit["published_utc"], "allocation authority publication") < generated:
        raise ValueError("allocation authority publication predates the authority artifact")
    return {
        "artifact": artifact,
        "artifact_evidence": artifact_evidence,
        "publication_audit": audit,
        "publication_audit_evidence": audit_evidence,
        "profile_evidence": profile_evidence,
    }


def validate_profile(
    profile: dict[str, object],
    profile_evidence: dict[str, object],
    matrix_evidence: dict[str, object],
    input_evidence: dict[str, dict[str, object]],
) -> tuple[dict[str, dict[str, object]], dict[str, object], dict[str, object]]:
    """Validate the fixed-path fresh reviewed profile artifact."""

    if normalized_path(Path(str(profile_evidence["path"]))) != profile_path():
        raise ValueError("allocation profile is not retained at the exact canonical path")
    require_exact_keys(
        profile,
        (
            "schema_version", "record_type", "execution_epoch", "canonical_root",
            "review", "provenance", "evidence", "profiles", "r03_authorization",
            "r17_readiness", "allocation_authority",
        ),
        "allocation profile",
    )
    if (
        profile["schema_version"] != 2
        or profile["record_type"] != "cgl_lf_stage_i_authenticated_wave_profile"
        or profile["execution_epoch"] != EXECUTION_EPOCH
        or profile["canonical_root"] != str(CANONICAL_ROOT)
    ):
        raise ValueError("allocation profile identity/root is invalid")
    review = profile["review"]
    if not isinstance(review, dict):
        raise ValueError("allocation profile review must be an object")
    require_exact_keys(
        review, ("decision", "reviewed_by", "reviewed_utc", "state_observed_utc", "notes"),
        "allocation profile review",
    )
    if review["decision"] != "approved":
        raise ValueError("allocation profile review decision must be approved")
    require_nonempty_string(review["reviewed_by"], "allocation profile reviewer")
    require_nonempty_string(review["notes"], "allocation profile review notes")
    reviewed = require_fresh(profile_evidence, review["reviewed_utc"], "allocation profile")
    observed = require_utc(review["state_observed_utc"], "allocation profile observed state")
    if reviewed < observed:
        raise ValueError("allocation profile review predates its observed state")
    static = validate_static_provenance(profile, matrix_evidence, input_evidence)
    values = profile["profiles"]
    if not isinstance(values, list):
        raise ValueError("allocation profiles must be a list")
    profiles: dict[str, dict[str, object]] = {}
    for index, value in enumerate(values):
        if not isinstance(value, dict):
            raise ValueError(f"allocation profile {index} must be an object")
        require_exact_keys(
            value,
            ("case_id", "nodes", "walltime", "athena_walltime", "next_increment", "estimated_runtime_seconds", "rationale"),
            f"allocation profile {index}",
        )
        case_id = require_nonempty_string(value["case_id"], f"profile {index} case_id")
        if case_id not in PROFILE_CASES or case_id in profiles:
            raise ValueError(f"allocation profile case {case_id} is unsupported or duplicated")
        increment = decimal_value(value["next_increment"], f"{case_id} next_increment")
        if increment <= 0 or increment > MAX_REVIEWED_INCREMENTS[case_id]:
            raise ValueError(f"{case_id} next_increment exceeds the reviewed packet bound")
        estimate = value["estimated_runtime_seconds"]
        if isinstance(estimate, bool) or not isinstance(estimate, int) or not 1 <= estimate <= 6000:
            raise ValueError(f"{case_id} estimated runtime must be in [1, 6000] seconds")
        if isinstance(value["nodes"], bool) or not isinstance(value["nodes"], int):
            raise ValueError(f"{case_id} reviewed profile nodes must be an integer")
        validate_resource_values(
            case_id, value["nodes"], value["walltime"], value["athena_walltime"],
            Decimal(int(value["nodes"]) * walltime_seconds(value["walltime"], "profile walltime")) / Decimal(3600),
            f"{case_id} reviewed profile",
        )
        require_nonempty_string(value["rationale"], f"{case_id} profile rationale")
        profiles[case_id] = {**value, "next_increment": increment}
    if set(profiles) != set(PROFILE_CASES):
        raise ValueError("allocation profile cases differ from the exact R04-R17 set")
    authority = validate_allocation_authority(profile, profile_evidence, profiles)
    return profiles, static, authority


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
    profile: dict[str, object],
    summary: dict[str, object],
    summary_evidence: dict[str, object],
    reconciliation: dict[str, object],
    reconciliation_evidence: dict[str, object],
) -> dict[str, object]:
    """Authenticate exact canonical stores and return validated case lineages."""

    transaction_state = discover_transaction_state()
    bindings = profile["evidence"]
    validate_binding(bindings["summary"], summary_evidence, summary_path(), "summary")
    summary_observed = require_fresh(summary_evidence, summary["updated_utc"], "controller summary")
    if profile["review"]["state_observed_utc"] != summary["updated_utc"]:
        raise ValueError("allocation profile does not review the current controller summary")

    ledger, ledger_evidence = read_ledger(ledger_path())
    validate_binding(bindings["ledger"], ledger_evidence, ledger_path(), "ledger")
    reservations_value, reservations_evidence = read_json_value(reservations_path(), "reservations")
    if not isinstance(reservations_value, list):
        raise ValueError("reservations store must contain a list")
    reservations = reservations_value
    validate_binding(bindings["reservations"], reservations_evidence, reservations_path(), "reservations")
    qualification, qualification_evidence = read_json_file(qualification_path(), "qualification")
    validate_binding(bindings["qualification"], qualification_evidence, qualification_path(), "qualification")

    manifest_bindings = bindings["manifests"]
    if not isinstance(manifest_bindings, list):
        raise ValueError("manifest bindings must be a list")
    discovered = sorted(
        normalized_path(path)
        for path in expected_path(f"runs/mks24-stage-i/{EXECUTION_EPOCH}").glob(
            "*/*/manifest/prepared_run.json"
        )
    )
    bound_paths = sorted(normalized_path(Path(str(item.get("path", "")))) for item in manifest_bindings if isinstance(item, dict))
    if discovered != bound_paths:
        raise ValueError("reviewed manifest bindings differ from the exact canonical manifest store")
    infos: list[dict[str, object]] = []
    manifest_evidence = []
    for binding in manifest_bindings:
        if not isinstance(binding, dict):
            raise ValueError("manifest binding must be an object")
        path = normalized_path(Path(str(binding.get("path", ""))))
        manifest, retained = read_json_file(path, "canonical manifest")
        validate_binding(binding, retained, path, "manifest")
        infos.append(validate_manifest(manifest, retained, cases))
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

    expected_recorded = [
        {
            "job": row["job_id"], "case": row["case_id"], "segment": row["segment"],
            "state": row["state"], "hours": f"{decimal_value(row['actual_node_hours'], 'actual'):.6f}",
            "result": row["result"],
        }
        for row in ledger
    ]
    expected_active = [
        {
            "case": info["case_id"],
            "segment": info["segment"],
            "nodes": str(info["nodes"]),
            "walltime": info["manifest"]["allocation"]["requested_walltime"],
            "hours": f"{info['reserved']:.6f}",
            "state": info["state"],
        }
        for info in active_infos
    ]
    if summary["recorded"] != expected_recorded or summary["active"] != expected_active:
        raise ValueError("controller summary rows differ from canonical ledger/reservations")
    expected_values = {
        "actual": actual,
        "reserved": active_reserved,
        "stage_remaining": max(Decimal("0"), STAGE_I_BUDGET_NODE_HOURS - actual - active_reserved),
        "project_remaining": PROJECT_BUDGET_NODE_HOURS - actual - active_reserved,
    }
    for key, expected in expected_values.items():
        if abs(summary[key] - expected) > Decimal("0.000001"):
            raise ValueError(f"controller summary {key} differs from canonical accounting")

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
        "summary": summary_evidence,
        "ledger": ledger_evidence,
        "reservations": reservations_evidence,
        "qualification": qualification_evidence,
        "manifests": manifest_evidence,
    }
    validate_binding(
        bindings["reconciliation"], reconciliation_evidence, reconciliation_path(),
        "reconciliation snapshot",
    )
    validate_reconciliation_snapshot(
        reconciliation, reconciliation_evidence, actual_evidence, counts, summary["updated_utc"]
    )

    by_case: dict[str, list[dict[str, object]]] = {case_id: [] for case_id in ALL_CASES}
    for info in infos:
        by_case[info["case_id"]].append(info)
    lineages = {}
    for case_id in ALL_CASES:
        lineages[case_id] = validated_lineage(case_id, by_case[case_id])
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
        "next_indexes": next_indexes,
        "summary_observed": summary_observed,
        "evidence": {
            **actual_evidence,
            "transaction_stores": transaction_state,
            "submitted_authentication": submitted_authentication,
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
            "revision": CONTROLLER_REVISION,
            "sha256": CONTROLLER_SHA256,
        }
        or command.get("source_bundle")
        != {
            "path": str(source_bundle_path()),
            "sha256": SOURCE_BUNDLE_SHA256,
            "verified_revisions": [SOURCE_REVISION, CONTROLLER_REVISION],
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
    profile: dict[str, object],
    lineage: list[dict[str, object]],
    active: dict[str, object] | None,
    case_infos: list[dict[str, object]],
) -> tuple[dict[str, object] | None, dict[str, object] | None]:
    """Validate the permanent published F-115 anchor and any exact R03 s02 intent."""

    binding = profile["r03_authorization"]
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
    if binding is None and not any(
        info["segment"] == R03_F115_SEGMENT for info in case_infos
    ):
        raise ValueError("incomplete R03 requires one bound sole-next authorization")
    if binding is not None and not isinstance(binding, dict):
        raise ValueError("R03 F-115 profile binding is malformed")
    path = r03_f115_path()
    if isinstance(binding, dict):
        bound_path = normalized_path(Path(str(binding.get("path", ""))))
        if bound_path != path:
            raise ValueError("R03 authorization is not the exact published F-115 artifact")
    authorization, evidence = read_json_file(path, "R03 sole-next authorization")
    if isinstance(binding, dict):
        validate_binding(binding, evidence, path, "R03 sole-next authorization")
    require_immutable_publication_evidence(evidence, R03_F115_SHA256, "R03 F-115 authority")
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
        "source_bundle": str(source_bundle_path()),
        "source_bundle_sha256": SOURCE_BUNDLE_SHA256,
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
            "to": str(source_bundle_path()),
        },
        "source_bundle_sha256": {
            "from": R03_F114_SOURCE_BUNDLE_SHA256,
            "to": SOURCE_BUNDLE_SHA256,
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
        "head": CONTROLLER_REVISION,
        "links": 1,
        "mode": "0644",
        "path": source_bundle_path().relative_to(CANONICAL_ROOT).as_posix(),
        "sha256": SOURCE_BUNDLE_SHA256,
        "verified_revisions": list(F115_SOURCE_BUNDLE_REQUIRED_REVISIONS),
    }:
        raise ValueError("R03 F-115 source-bundle implementation declaration differs")
    if not isinstance(implementation_helper, dict) or (
        implementation_helper.get("links") != 1
        or implementation_helper.get("mode") != "0644"
        or implementation_helper.get("sha256") != CONTROLLER_SHA256
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
        expected=source_bundle_path(),
        digest=SOURCE_BUNDLE_SHA256,
        mode="0644",
        label="R03 F-115 authoritative source bundle",
    )
    if (
        bundle.get("complete_history") is not True
        or bundle.get("head") != CONTROLLER_REVISION
        or helper
        != {
            "path": "scripts/frontier/cgl_lf_stage_i.py",
            "sha256": CONTROLLER_SHA256,
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


def validate_r17_rank_readiness(value: dict[str, object]) -> dict[str, object]:
    """Require measured 64-rank scheduler, output, restart, and performance evidence."""

    require_exact_keys(
        value,
        (
            "schema_version", "record_type", "execution_epoch", "generated_utc",
            "reviewed_utc", "reviewed_by", "decision", "nodes", "ranks", "job_id",
            "final_time", "scheduler_evidence", "batch_script", "output_inventory",
            "terminal_restart", "performance",
        ),
        "R17 rank readiness",
    )
    job_id = require_nonempty_string(value["job_id"], "R17 readiness job ID")
    if (
        value["record_type"] != "cgl_lf_stage_i_r17_64_rank_readiness"
        or value["nodes"] != 8
        or value["ranks"] != 64
        or JOB_RE.fullmatch(job_id) is None
    ):
        raise ValueError("R17 64-rank readiness identity is incomplete")
    readiness_root = expected_path(
        f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_64_rank_readiness"
    )
    scheduler_path = expected_path(f"accounting/{job_id}.stage_i.sacct.txt")
    scheduler_payload, scheduler_evidence = read_stable_regular_file(
        scheduler_path, "R17 rank-readiness scheduler evidence"
    )
    validate_binding(
        value["scheduler_evidence"], scheduler_evidence, scheduler_path,
        "R17 rank-readiness scheduler evidence",
    )
    scheduler = parse_scheduler_evidence(scheduler_payload, job_id)
    elapsed = scheduler["elapsed_seconds"]
    if (
        scheduler["job_name"] != "cgl_mks24_r17_64_rank_readiness"
        or scheduler["state"] != "COMPLETED"
        or scheduler["exit_code"] != "0:0"
        or scheduler["nodes"] != "8"
        or DIGITS_RE.fullmatch(elapsed) is None
        or int(elapsed) <= 0
    ):
        raise ValueError("R17 64-rank scheduler evidence is not a successful measured run")
    submitted = require_scheduler_timestamp(
        scheduler["submitted_utc"], "R17 rank-readiness submitted time"
    )
    completed = require_scheduler_timestamp(
        scheduler["completed_utc"], "R17 rank-readiness completed time"
    )
    if completed < submitted or completed > current_utc() + FUTURE_SKEW:
        raise ValueError("R17 64-rank scheduler chronology is invalid")
    batch_path = readiness_root / "r17_64_rank_readiness.sbatch"
    batch_payload, batch_evidence = read_stable_regular_file(
        batch_path, "R17 rank-readiness batch script"
    )
    validate_binding(value["batch_script"], batch_evidence, batch_path, "R17 rank-readiness batch script")
    if normalized_batch_script_sha256(batch_payload) != BATCH_SCRIPT_DIGEST_PATTERN.findall(
        batch_payload.decode("utf-8")
    )[0]:
        raise ValueError("R17 rank-readiness batch script self-digest is invalid")
    outputs = value["output_inventory"]
    if not isinstance(outputs, list) or len(outputs) != 64:
        raise ValueError("R17 rank-readiness output inventory is incomplete")
    output_evidence = []
    seen: set[Path] = set()
    for index, binding in enumerate(outputs):
        if not isinstance(binding, dict):
            raise ValueError("R17 rank-readiness output binding is malformed")
        path = normalized_path(Path(str(binding.get("path", ""))))
        if path.parent != readiness_root / "output" / f"rank_{index:08d}":
            raise ValueError("R17 rank-readiness output inventory is not exact per rank")
        if path in seen:
            raise ValueError("R17 rank-readiness output inventory contains duplicates")
        seen.add(path)
        _, retained = read_stable_regular_file(path, f"R17 rank-readiness output {index}")
        validate_binding(binding, retained, path, f"R17 rank-readiness output {index}")
        output_evidence.append(retained)
    final_time = decimal_value(value["final_time"], "R17 rank-readiness final time")
    if final_time <= 0:
        raise ValueError("R17 rank-readiness final time must be positive")
    restart = validate_restart_group(
        value["terminal_restart"],
        expected_root=readiness_root / "restarts",
        expected_count=64,
        expected_time=final_time,
        label="R17 rank readiness",
    )
    performance = value["performance"]
    if not isinstance(performance, dict):
        raise ValueError("R17 rank-readiness performance must be an object")
    require_exact_keys(
        performance,
        ("elapsed_seconds", "node_hours", "simulated_time", "meshblocks_per_rank"),
        "R17 rank-readiness performance",
    )
    expected_node_hours = Decimal(8 * int(elapsed)) / Decimal(3600)
    if (
        performance["elapsed_seconds"] != int(elapsed)
        or abs(decimal_value(performance["node_hours"], "R17 readiness node-hours") - expected_node_hours)
        > Decimal("0.000001")
        or decimal_value(performance["simulated_time"], "R17 readiness simulated time") != final_time
        or isinstance(performance["meshblocks_per_rank"], bool)
        or not isinstance(performance["meshblocks_per_rank"], int)
        or performance["meshblocks_per_rank"] <= 0
    ):
        raise ValueError("R17 rank-readiness performance differs from measured evidence")
    return {
        "scheduler": scheduler_evidence,
        "batch_script": batch_evidence,
        "outputs": output_evidence,
        "restart": restart,
    }


def validate_r17_readiness(
    profile: dict[str, object],
    state: dict[str, object],
    cases: dict[str, dict[str, object]],
    profiles: dict[str, dict[str, object]],
) -> dict[str, object]:
    """Require digest-bound measured and independently reproduced R17 readiness."""

    readiness = profile["r17_readiness"]
    if not isinstance(readiness, dict):
        raise ValueError("R17 readiness must be an object")
    require_exact_keys(
        readiness, ("decision", "reviewed_by", "reviewed_utc", "notes", "evidence"),
        "R17 readiness",
    )
    if readiness["decision"] != "approved":
        raise ValueError("R17 prerequisites are not approved")
    require_nonempty_string(readiness["reviewed_by"], "R17 readiness reviewer")
    require_nonempty_string(readiness["notes"], "R17 readiness notes")
    reviewed = require_utc(readiness["reviewed_utc"], "R17 readiness reviewed_utc")
    now = current_utc()
    if (
        reviewed < state["summary_observed"]
        or reviewed > now + FUTURE_SKEW
        or now - reviewed > R17_READINESS_MAX_AGE
    ):
        raise ValueError("R17 readiness is stale relative to canonical state")
    evidence = readiness["evidence"]
    if not isinstance(evidence, dict):
        raise ValueError("R17 readiness evidence must be an object")
    exact_paths = {
        "storage": expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_storage_readiness.json"),
        "recost": expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_node_hour_recost.json"),
        "rank_readiness": expected_path(f"accounting/mks24_stage_i_{EXECUTION_EPOCH_SLUG}_R17_64_rank_readiness.json"),
    }
    if set(evidence) != set(exact_paths):
        raise ValueError("R17 readiness evidence set is incomplete")
    values = {}
    retained = {}
    for key, path in exact_paths.items():
        value, file_evidence = read_json_file(path, f"R17 {key}")
        validate_binding(evidence[key], file_evidence, path, f"R17 {key}")
        require_fresh(file_evidence, value.get("generated_utc"), f"R17 {key}", R17_READINESS_MAX_AGE)
        artifact_reviewed = require_fresh(
            file_evidence, value.get("reviewed_utc"), f"R17 {key} review",
            R17_READINESS_MAX_AGE,
        )
        if (
            value.get("schema_version") != 1
            or value.get("execution_epoch") != EXECUTION_EPOCH
            or value.get("decision") != "approved"
            or not require_nonempty_string(value.get("reviewed_by"), f"R17 {key} reviewer")
            or artifact_reviewed < state["summary_observed"]
        ):
            raise ValueError(f"R17 {key} evidence is not approved for current state")
        values[key] = value
        retained[key] = file_evidence
    storage = values["storage"]
    if storage.get("record_type") != "cgl_lf_stage_i_r17_storage_readiness":
        raise ValueError("R17 storage evidence record_type is invalid")
    available = storage.get("available_bytes")
    required = storage.get("required_retention_bytes")
    live_available = available_storage_bytes(CANONICAL_ROOT)
    if (
        isinstance(available, bool) or not isinstance(available, int)
        or isinstance(required, bool) or not isinstance(required, int)
        or required < R17_REQUIRED_RETENTION_BYTES
        or available < required
        or live_available < required
    ):
        raise ValueError("R17 storage evidence does not satisfy retained-capacity policy")
    recost = values["recost"]
    if recost.get("record_type") != "cgl_lf_stage_i_r17_node_hour_recost":
        raise ValueError("R17 recost evidence record_type is invalid")
    r17_reserved = Decimal(
        int(profiles[R17]["nodes"])
        * walltime_seconds(profiles[R17]["walltime"], "R17 profile walltime")
    ) / Decimal(3600)
    projection = campaign_projection(cases, state, {R17: r17_reserved})
    if (
        recost.get("projection") != projection
        or recost.get("projection_sha256") != sha256_bytes(canonical_json(projection))
        or decimal_value(
            recost.get("projected_r17_reserved_node_hours"), "R17 projected reserved use"
        )
        != r17_reserved
    ):
        raise ValueError("R17 recost evidence differs from reproduced campaign projection")
    require_projection_within_budget(projection, "R17 recost")
    rank = values["rank_readiness"]
    retained["rank_readiness_support"] = validate_r17_rank_readiness(rank)
    retained["live_available_bytes"] = live_available
    retained["reproduced_recost_projection"] = projection
    return retained


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

    increment = profile["next_increment"]
    if start == 0 and increment != INITIAL_TARGETS[case_id]:
        raise ValueError(f"{case_id} fresh next_increment differs from the reviewed initial target")
    target = min(Decimal("10"), start + increment)
    if target <= start:
        raise ValueError(f"{case_id} has no positive exact next interval")
    codes = policy_codes(case_id)
    segment = f"s{next_index:02d}_rankio_t{time_token(start)}_t{time_token(target)}"
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
            else "planning_only_independent_profile_non_authorizing"
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
            "estimated_runtime_seconds": profile["estimated_runtime_seconds"],
            "review_rationale": profile["rationale"],
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
    profile: dict[str, object],
    profile_evidence: dict[str, object],
    summary: dict[str, object],
    summary_evidence: dict[str, object],
    reconciliation: dict[str, object],
    reconciliation_evidence: dict[str, object],
) -> dict[str, object]:
    """Build one evidence-bound deterministic read-only advisory diagnostic."""

    cases, input_evidence = validate_matrix(matrix, matrix_evidence)
    profiles, static, allocation_authority = validate_profile(
        profile, profile_evidence, matrix_evidence, input_evidence
    )
    state = validate_canonical_state(
        cases, profile, summary, summary_evidence, reconciliation, reconciliation_evidence
    )
    base_projection = campaign_projection(cases, state)
    validate_authority_projection(allocation_authority, base_projection)
    lineages = state["lineages"]
    active_map = {info["case_id"]: info for info in state["active_infos"]}

    r03_lineage, r03_active = lineages[R03]
    r03_sole, r03_authority = validate_r03_authorization(
        profile, r03_lineage, r03_active, state["case_infos"][R03]
    )
    for case_id in PROFILE_CASES:
        lineage, active = lineages[case_id]
        if not lineage and profiles[case_id]["next_increment"] != INITIAL_TARGETS[case_id]:
            raise ValueError(f"{case_id} fresh next_increment differs from the reviewed initial target")
        if lineage and lineage[0]["target"] != INITIAL_TARGETS[case_id]:
            raise ValueError(f"{case_id} first accepted target differs from the reviewed plan")
        validate_active_reviewed_profile(case_id, lineage, active, profiles[case_id])
    lower_state = {
        case_id: {
            "lineage": lineages[case_id][0],
            "active": lineages[case_id][1],
            "endpoint": lineage_endpoint(lineages[case_id][0]),
        }
        for case_id in LOWER_CASES
    }
    predecessor_state = {
        case_id: lineage_endpoint(lineages[case_id][0]) == Decimal("10")
        for case_id in tuple(f"R{index:02d}" for index in range(2, 17))
    }
    r17_lineage, r17_active = lineages[R17]
    r17_endpoint = lineage_endpoint(r17_lineage)
    if r17_lineage or r17_active:
        if not all(predecessor_state.values()):
            raise ValueError("R17 started before every authenticated R02-R16 lineage reached exact t=10")
        if any(info["case_id"] != R17 for info in state["active_infos"]):
            raise ValueError("R17 is not exclusive")

    packets = []
    deferred: list[dict[str, str]] = []
    if r17_endpoint == Decimal("10"):
        if state["active_infos"]:
            raise ValueError("R17 is complete but an active Stage I lane remains")
        wave_status = "campaign_complete"
        r17_readiness = None
    elif r17_active is not None:
        wave_status = "r17_active"
        r17_readiness = validate_r17_readiness(profile, state, cases, profiles)
    elif all(predecessor_state.values()):
        if state["active_infos"]:
            raise ValueError("R17-ready predecessor state still has an active lower lane")
        r17_readiness = validate_r17_readiness(profile, state, cases, profiles)
        continuation = terminal_restart(r17_lineage[-1]) if r17_lineage else None
        packets.append(
            build_packet(
                R17, cases[R17], profiles[R17], int(state["next_indexes"][R17]),
                r17_endpoint,
                static, continuation,
            )
        )
        wave_status = "r17_exclusive"
    else:
        if r17_lineage or r17_active:
            raise ValueError("R17 cannot coexist with incomplete predecessors")
        r17_readiness = None
        candidates: list[str] = []
        if lineage_endpoint(r03_lineage) < Decimal("10") and r03_active is None:
            candidates.append(R03)
        candidates.extend(
            case_id for case_id in PRODUCTION_PRIORITY[1:]
            if lower_state[case_id]["lineage"]
            and lower_state[case_id]["endpoint"] < Decimal("10")
            and lower_state[case_id]["active"] is None
        )
        candidates.extend(
            case_id for case_id in PRODUCTION_PRIORITY[1:]
            if not lower_state[case_id]["lineage"] and lower_state[case_id]["active"] is None
        )
        planned_nodes = state["active_nodes"]
        planned_lanes = len(state["active_infos"])
        planned_reserved = state["active_reserved_node_hours"]
        planned_reserved_by_case: dict[str, Decimal] = {}
        active_nonconcurrent = any(
            info["case_id"] not in CONCURRENT_CASES for info in state["active_infos"]
        )
        for case_id in candidates:
            if active_nonconcurrent:
                deferred.append({"case_id": case_id, "reason": "active_nonconcurrent_case"})
                continue
            if planned_lanes >= MAX_LANES:
                deferred.append({"case_id": case_id, "reason": "four_lane_ceiling"})
                continue
            if case_id == R03:
                if r03_sole is None or r03_authority is None:
                    raise ValueError("incomplete inactive R03 lacks sole-next authorization")
                profile_value = {
                    "nodes": 1,
                    "walltime": r03_sole["walltime"],
                    "athena_walltime": r03_sole["athena_walltime"],
                    "next_increment": decimal_value(r03_sole["time_tlim_target"], "R03 target")
                    - lineage_endpoint(r03_lineage),
                    "estimated_runtime_seconds": 6000,
                    "rationale": (
                        "planning-only reference to the exact published F-115 R03 s02 "
                        "continuation profile"
                    ),
                }
                nodes = 1
            else:
                profile_value = profiles[case_id]
                nodes = int(profile_value["nodes"])
            if planned_nodes + nodes > MAX_NODES:
                deferred.append({"case_id": case_id, "reason": "ten_node_ceiling"})
                continue
            candidate_reserved = Decimal(
                nodes * walltime_seconds(profile_value["walltime"], f"{case_id} profile walltime")
            ) / Decimal(3600)
            candidate_projection = campaign_projection(
                cases,
                state,
                {
                    **planned_reserved_by_case,
                    case_id: planned_reserved_by_case.get(case_id, Decimal("0"))
                    + candidate_reserved,
                },
            )
            if (
                state["actual_node_hours"] + planned_reserved + candidate_reserved
                > STAGE_I_BUDGET_NODE_HOURS
                or state["actual_node_hours"] + planned_reserved + candidate_reserved
                > PROJECT_BUDGET_NODE_HOURS
                or projection_total(candidate_projection) > STAGE_I_BUDGET_NODE_HOURS
                or projection_total(candidate_projection) > PROJECT_BUDGET_NODE_HOURS
            ):
                deferred.append({"case_id": case_id, "reason": "budget_ceiling"})
                continue
            if case_id == R03:
                continuation = r03_authority["continuation"]
                authorization = {
                    **{
                        key: value for key, value in r03_authority.items()
                        if key != "continuation"
                    },
                    "sole_next_segment_profile": r03_sole,
                }
                lineage = r03_lineage
                endpoint = lineage_endpoint(lineage)
            else:
                lineage = lower_state[case_id]["lineage"]
                endpoint = lower_state[case_id]["endpoint"]
                continuation = terminal_restart(lineage[-1]) if lineage else None
                authorization = None
            packets.append(
                build_packet(
                    case_id, cases[case_id], profile_value,
                    int(state["next_indexes"][case_id]), endpoint,
                    static, continuation, authorization,
                )
            )
            planned_nodes += nodes
            planned_lanes += 1
            planned_reserved += candidate_reserved
            planned_reserved_by_case[case_id] = (
                planned_reserved_by_case.get(case_id, Decimal("0")) + candidate_reserved
            )
        wave_status = "lower_resolution"

    total_nodes = state["active_nodes"] + sum(int(packet["allocation"]["nodes"]) for packet in packets)
    total_lanes = len(state["active_infos"]) + len(packets)
    planned_reserved_node_hours = sum(
        (packet_reserved_node_hours(packet) for packet in packets), Decimal("0")
    )
    total_committed_node_hours = (
        state["actual_node_hours"]
        + state["active_reserved_node_hours"]
        + planned_reserved_node_hours
    )
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
    if len(packet_cases) != len(set(packet_cases)) or set(packet_cases) & set(active_map):
        raise ValueError("planned wave contains duplicate case lanes")
    final_projection = campaign_projection(
        cases,
        state,
        {
            packet["case_id"]: packet_reserved_node_hours(packet)
            for packet in packets
        },
    )
    require_projection_within_budget(final_projection, "planned wave")

    core = {
        "schema_version": 3,
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
        ],
        "evidence": {
            "matrix": matrix_evidence,
            "allocation_profile": profile_evidence,
            "controller_summary": summary_evidence,
            "reconciliation": reconciliation_evidence,
            "canonical_state": state["evidence"],
            "allocation_authority": allocation_authority,
            "r03_f115_authority": r03_authority,
            "r17_readiness": r17_readiness,
        },
        "policy": {
            "max_lanes": MAX_LANES,
            "max_nodes": MAX_NODES,
            "max_segment_seconds": MAX_SEGMENT_SECONDS,
            "shutdown_margin_seconds": SHUTDOWN_MARGIN_SECONDS,
            "ranks_per_node": RANKS_PER_NODE,
            "cpus_per_task": CPUS_PER_TASK,
            "initial_lower_wave_priority": list(PRODUCTION_PRIORITY[:4]),
            "rolling_priority": list(ROLLING_PRIORITY),
            "r17_exclusive_last": True,
        },
        "physics_policies": POLICY_TEXT,
        "observed_state": {
            "updated_utc": summary["updated_utc"],
            "reconciliation_counts": state["counts"],
            "actual_node_hours": decimal_text(state["actual_node_hours"]),
            "active_reserved_node_hours": decimal_text(state["active_reserved_node_hours"]),
            "planned_reserved_node_hours": decimal_text(planned_reserved_node_hours),
            "total_committed_node_hours": decimal_text(total_committed_node_hours),
            "credible_remaining_campaign_projection": final_projection,
            "active_lanes": summary["active"],
            "predecessors_accepted_exact_t10": predecessor_state,
        },
        "wave": {
            "status": wave_status,
            "active_lane_count": len(state["active_infos"]),
            "active_nodes": state["active_nodes"],
            "planned_lane_count": len(packets),
            "planned_nodes": sum(int(packet["allocation"]["nodes"]) for packet in packets),
            "total_lane_count": total_lanes,
            "total_nodes": total_nodes,
            "packets": packets,
            "deferred": sorted(deferred, key=lambda item: (item["case_id"], item["reason"])),
        },
    }
    core["evidence"]["global_toctou_boundary"] = revalidate_global_boundary(
        [
            profile,
            matrix_evidence,
            profile_evidence,
            summary_evidence,
            reconciliation_evidence,
            static,
            state,
            allocation_authority,
            r03_authority,
            r17_readiness,
            packets,
        ],
        state["transaction_state"],
        state["submitted_scheduler"],
    )
    return {**core, "plan_sha256": sha256_bytes(canonical_json(core))}


def plan_from_paths(
    matrix_path_value: Path,
    profile_path_value: Path,
    summary_path_value: Path,
    reconciliation_path_value: Path,
) -> dict[str, object]:
    """Read all exact evidence and build one non-authorizing advisory diagnostic."""

    matrix, matrix_evidence = read_json_file(matrix_path_value, "matrix")
    profile, profile_evidence = read_json_file(profile_path_value, "allocation profile")
    summary_payload, summary_evidence = read_stable_regular_file(summary_path_value, "controller summary")
    summary = parse_summary(summary_payload)
    reconciliation, reconciliation_evidence = read_json_file(
        reconciliation_path_value, "reconciliation snapshot"
    )
    return build_wave_plan(
        matrix, matrix_evidence, profile, profile_evidence, summary, summary_evidence,
        reconciliation, reconciliation_evidence,
    )


def parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument("--matrix", required=True, type=Path)
    command.add_argument("--allocation-profile", required=True, type=Path)
    command.add_argument("--summary", required=True, type=Path)
    command.add_argument("--reconciliation", required=True, type=Path)
    return command


def main(argv: list[str] | None = None) -> int:
    """Read exact canonical evidence and print one command-free JSON plan."""

    args = parser().parse_args(argv)
    try:
        if normalized_path(args.matrix) != frozen_matrix_path():
            raise ValueError("CLI matrix path is not the exact frozen matrix")
        if normalized_path(args.allocation_profile) != profile_path():
            raise ValueError("CLI allocation-profile path is not the exact canonical profile")
        if normalized_path(args.summary) != summary_path():
            raise ValueError("CLI summary path is not the exact canonical controller summary")
        if normalized_path(args.reconciliation) != reconciliation_path():
            raise ValueError("CLI reconciliation path is not the exact canonical snapshot")
        plan = plan_from_paths(
            args.matrix, args.allocation_profile, args.summary, args.reconciliation
        )
    except (OSError, ValueError, KeyError, TypeError) as error:
        print(f"Stage I wave planner failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(plan, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
