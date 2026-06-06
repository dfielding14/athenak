#!/usr/bin/env python3
"""Audit accelerated Stage I retained state for available CT evidence.

The direct-fast reporter indexes rank-local AthenaK ``.bin`` snapshots.  These
files retain cell-centered ``bcc1/bcc2/bcc3`` fields, not the staggered
face-centered magnetic state advanced by constrained transport.  This adapter
authenticates an externally SHA-bound assembled inventory and its selected case
provenance, binds and inspects available snapshot bytes, and reports that claim
boundary without converting cell-centered data into false CT-divB evidence.  It
also discovers exact rank-local native restart groups beneath authenticated
selected direct-fast lineage roots and reuses the reviewed native restart
CT-divB parser for exact preregistered t=9 and t=10 states when they exist.
The selected ``accepted_r02_bundle`` lineage is also supported: its exact
externally bound whole-case bundle and accepted production-segment manifests
authenticate the declared t=9 and t=10 rank-local restart groups.

The adapter supports partial campaigns and emits deterministic JSON, CSV, and
Markdown.  Direct-fast numerical CT evidence is non-authorizing and does not
substitute for a canonical accepted-bundle inventory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import stat
import struct
import sys
import tempfile
from typing import Iterable


SCHEMA_VERSION = 2
RECORD_TYPE = "stage-i-direct-fast-ct-audit"
EVIDENCE_DIGEST_METHOD = "sha256-canonical-json-without-evidence-digest"
CASE_ID = re.compile(r"R(?:0[2-9]|1[0-7])")
SHA256 = re.compile(r"[0-9a-f]{64}")
RANK_DIRECTORY = re.compile(r"rank_(\d{8})")
MAX_JSON_BYTES = 256 * 1024 * 1024
MAX_HEADER_BYTES = 16 * 1024 * 1024
CELL_CENTERED_MAGNETIC_FIELDS = ("bcc1", "bcc2", "bcc3")
NATIVE_CT_ACCEPTANCE_PATH = Path(__file__).with_name(
    "cgl_lf_stage_i_scientific_acceptance.py"
)
CT_ABSENCE_REASON = (
    "Authenticated AthenaK rank-local .bin snapshots expose cell-centered "
    "bcc1/bcc2/bcc3, not the staggered face-centered magnetic state advanced "
    "by constrained transport; discrete CT divB is therefore not measurable "
    "from these retained snapshots."
)
DIRECT_FAST_AUTHORITY_LIMITATION = (
    "The selected direct-fast lineage authenticates output roots and rank counts, "
    "but it is not a canonical accepted whole-case bundle; native restart CT-divB "
    "evidence is therefore numerical, sampled-state, and non-authorizing."
)
RESTART_INDEX_LIMITATION = (
    "Direct-fast lineage manifests and snapshots.json do not declare produced "
    "restart groups. Exact groups must be discovered and byte-bound at audit time "
    "beneath authenticated selected direct-fast output roots."
)
ACCEPTED_R02_AUTHORITY_LIMITATION = (
    "The exact selected accepted R02 whole-case bundle and its accepted t=9,t=10 "
    "segment states support sampled numerical CT-divB evidence, but this adapter "
    "does not grant canonical campaign or release authority."
)
ACCEPTED_R02_RESTART_INDEX = (
    "The selected accepted R02 whole-case bundle names the exact accepted production "
    "segment manifests; their scientific inspections declare the exact retained "
    "t=9,t=10 rank-local native restart groups."
)


class FastCtAuditError(RuntimeError):
    """Raised when campaign-level authentication or configuration fails."""


class CaseAuthenticationError(FastCtAuditError):
    """Raised when one assembled case cannot be authenticated."""


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise FastCtAuditError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> object:
    """Reject non-finite JSON constants."""

    raise FastCtAuditError(f"invalid JSON numeric constant: {value}")


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise FastCtAuditError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise FastCtAuditError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise FastCtAuditError(f"{label} must be a nonempty string")
    return value


def require_int(value: object, label: str, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise FastCtAuditError(f"{label} must be an integer >= {minimum}")
    return value


def require_finite(value: object, label: str) -> float:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(float(value))
    ):
        raise FastCtAuditError(f"{label} must be finite")
    return float(value)


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise FastCtAuditError(f"{label} must be a lowercase SHA-256")
    return value


def canonical_json(value: object) -> bytes:
    """Return deterministic compact JSON bytes."""

    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def evidence_digest(value: dict[str, object]) -> str:
    body = dict(value)
    body.pop("evidence_digest", None)
    return hashlib.sha256(canonical_json(body)).hexdigest()


def seal_evidence(value: dict[str, object]) -> dict[str, object]:
    result = dict(value)
    result["evidence_digest"] = {
        "method": EVIDENCE_DIGEST_METHOD,
        "sha256": evidence_digest(result),
    }
    return result


def stable_profile(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def descriptor_sha256(descriptor: int) -> str:
    digest = hashlib.sha256()
    offset = 0
    while True:
        block = os.pread(descriptor, 1024 * 1024, offset)
        if not block:
            break
        digest.update(block)
        offset += len(block)
    return digest.hexdigest()


def regular_file_binding(
    path: Path,
    label: str,
    *,
    expected_sha256: str | None = None,
    expected_size: int | None = None,
) -> dict[str, object]:
    """Bind one unchanged regular file, resolving a declared leaf symlink."""

    declared = path.expanduser().absolute()
    try:
        resolved = declared.resolve(strict=True)
    except OSError as error:
        raise FastCtAuditError(f"{label} does not resolve: {declared}") from error
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(resolved, flags)
    except OSError as error:
        raise FastCtAuditError(f"{label} cannot be opened safely: {resolved}") from error
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise FastCtAuditError(f"{label} is not a regular file: {resolved}")
        digest = descriptor_sha256(descriptor)
        after = os.fstat(descriptor)
        if stable_profile(before) != stable_profile(after):
            raise FastCtAuditError(f"{label} changed while it was hashed")
        if expected_sha256 is not None and digest != expected_sha256:
            raise FastCtAuditError(f"{label} SHA-256 differs from its authority")
        if expected_size is not None and before.st_size != expected_size:
            raise FastCtAuditError(f"{label} size differs from its authority")
        result: dict[str, object] = {
            "path": str(resolved),
            "size_bytes": before.st_size,
            "sha256": digest,
        }
        if declared != resolved:
            result["declared_path"] = str(declared)
        return result
    finally:
        os.close(descriptor)


def read_bound_json(
    path: Path, label: str, *, expected_sha256: str | None = None
) -> tuple[dict[str, object], dict[str, object]]:
    binding = regular_file_binding(path, label, expected_sha256=expected_sha256)
    size = require_int(binding["size_bytes"], f"{label} size")
    if size > MAX_JSON_BYTES:
        raise FastCtAuditError(f"{label} exceeds the JSON size limit")
    resolved = Path(str(binding["path"]))
    payload = resolved.read_bytes()
    if len(payload) != size or hashlib.sha256(payload).hexdigest() != binding["sha256"]:
        raise FastCtAuditError(f"{label} changed between binding and read")
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise FastCtAuditError(f"{label} is invalid JSON") from error
    return require_dict(value, label), binding


def verify_declared_binding(value: object, label: str) -> dict[str, object]:
    record = require_dict(value, label)
    path = Path(require_text(record.get("path"), f"{label} path"))
    expected_sha = require_sha256(record.get("sha256"), f"{label} sha256")
    expected_size = require_int(record.get("size_bytes"), f"{label} size")
    return regular_file_binding(
        path,
        label,
        expected_sha256=expected_sha,
        expected_size=expected_size,
    )


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as stream:
            staged = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def safe_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except ValueError:
        return False


def expand_cases(values: Iterable[str], available: Iterable[str]) -> list[str]:
    requested: list[str] = []
    for value in values:
        for item in value.split(","):
            item = item.strip()
            if not item:
                continue
            match = re.fullmatch(r"R(\d{2})-R(\d{2})", item)
            if match:
                start, end = (int(token) for token in match.groups())
                if start > end:
                    raise FastCtAuditError(f"descending case range is invalid: {item}")
                requested.extend(f"R{number:02d}" for number in range(start, end + 1))
            else:
                requested.append(item)
    if not requested:
        requested = sorted(available)
    invalid = [case for case in requested if CASE_ID.fullmatch(case) is None]
    if invalid:
        raise FastCtAuditError(f"unsupported cases: {', '.join(invalid)}")
    return list(dict.fromkeys(requested))


def load_native_ct_authority() -> tuple[
    object, dict[str, object], list[float], float, dict[str, object]
]:
    """Load the reviewed native restart CT parser and its validated policy."""

    name = "_cgl_lf_stage_i_fast_ct_audit_acceptance"
    module = sys.modules.get(name)
    if module is None:
        spec = importlib.util.spec_from_file_location(name, NATIVE_CT_ACCEPTANCE_PATH)
        if spec is None or spec.loader is None:
            raise FastCtAuditError(
                f"cannot load native restart CT authority: {NATIVE_CT_ACCEPTANCE_PATH}"
            )
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        try:
            spec.loader.exec_module(module)
        except Exception as error:
            sys.modules.pop(name, None)
            raise FastCtAuditError(
                f"native restart CT authority could not be loaded: {error}"
            ) from error
    try:
        policy = module.load_validated_policy(
            module.DEFAULT_CRITERIA, module.DEFAULT_CRITERIA_REVIEW
        )
        criteria = require_dict(policy.get("criteria"), "native CT criteria")
        ct_policy = require_dict(criteria.get("ct_divb_policy"), "native CT policy")
        required_times = [
            require_finite(value, "required native CT state time")
            for value in require_list(
                ct_policy.get("required_state_times"), "required native CT state times"
            )
        ]
        threshold = require_finite(
            ct_policy.get("normalized_ct_divb_lt"), "native CT-divB threshold"
        )
        if threshold <= 0.0 or required_times != [9.0, 10.0]:
            raise FastCtAuditError("reviewed native CT policy differs from Stage I t=9,10")
        sources = {
            "scientific_acceptance": regular_file_binding(
                NATIVE_CT_ACCEPTANCE_PATH, "native CT acceptance utility"
            ),
            "criteria": require_dict(policy.get("criteria_binding"), "criteria binding"),
            "criteria_review": require_dict(
                policy.get("review_binding"), "criteria review binding"
            ),
        }
    except FastCtAuditError:
        raise
    except Exception as error:
        raise FastCtAuditError(
            f"native restart CT authority could not be validated: {error}"
        ) from error
    return module, policy, required_times, threshold, sources


def native_restart_time(acceptance: object, path: Path) -> float:
    """Read one native restart physical time with the reviewed prefix parser."""

    with acceptance.open_stable_regular(path, "direct-fast restart representative") as descriptor:
        profile = os.fstat(descriptor)
        prefix = acceptance.pread_exact(
            descriptor,
            min(profile.st_size, acceptance.MAX_PARAMETER_DUMP_BYTES),
            0,
            "direct-fast restart parameter prefix",
        )
        _, parameter_end = acceptance.parse_parameter_dump(prefix)
        header = acceptance.pread_exact(
            descriptor,
            acceptance.RESTART_HEADER_SIZE,
            parameter_end,
            "direct-fast restart mesh header",
        )
        unpacked = struct.unpack(acceptance.RESTART_HEADER_FORMAT, header)
        return require_finite(unpacked[49], "direct-fast restart physical time")


def authenticate_direct_fast_segment(
    case_id: str,
    case_name: object,
    value: object,
    matrix_sha256: str,
) -> dict[str, object]:
    """Authenticate one selected direct-fast segment as a restart source root."""

    segment = require_dict(value, f"{case_id} direct-fast lineage segment")
    if segment.get("kind") != "fast":
        raise CaseAuthenticationError("lineage segment is not direct-fast")
    order = require_int(segment.get("order"), f"{case_id} direct-fast segment order")
    segment_dir = Path(
        require_text(segment.get("segment_dir"), f"{case_id} direct-fast segment directory")
    ).absolute()
    output = Path(
        require_text(segment.get("output"), f"{case_id} direct-fast output")
    ).absolute()
    if output != segment_dir / "output":
        raise CaseAuthenticationError(
            f"{case_id} direct-fast segment {order} output root differs"
        )
    manifest_record = require_dict(
        segment.get("manifest"), f"{case_id} direct-fast segment {order} manifest"
    )
    manifest_path = Path(
        require_text(
            manifest_record.get("path"), f"{case_id} direct-fast manifest path"
        )
    ).absolute()
    if manifest_path != segment_dir / "manifest" / "fast_run.json":
        raise CaseAuthenticationError(
            f"{case_id} direct-fast segment {order} manifest path differs"
        )
    manifest, manifest_binding = read_bound_json(
        manifest_path,
        f"{case_id} direct-fast segment {order} manifest",
        expected_sha256=require_sha256(
            manifest_record.get("sha256"), f"{case_id} direct-fast manifest SHA-256"
        ),
    )
    if int(manifest_binding["size_bytes"]) != require_int(
        manifest_record.get("size_bytes"), f"{case_id} direct-fast manifest size"
    ):
        raise CaseAuthenticationError(
            f"{case_id} direct-fast segment {order} manifest size differs"
        )
    ranks = require_int(segment.get("ranks"), f"{case_id} direct-fast ranks", 1)
    nodes = require_int(manifest.get("nodes"), f"{case_id} direct-fast nodes", 1)
    ranks_per_node = require_int(
        manifest.get("ranks_per_node"), f"{case_id} direct-fast ranks per node", 1
    )
    expected = {
        "schema_version": 1,
        "case_id": case_id,
        "case_name": case_name,
        "matrix_sha256": matrix_sha256,
        "output_dir": str(output),
        "run_dir": str(segment_dir),
        "ranks": ranks,
    }
    for key, expected_value in expected.items():
        if manifest.get(key) != expected_value:
            raise CaseAuthenticationError(
                f"{case_id} direct-fast segment {order} manifest {key} differs"
            )
    for key in ("input_sha256", "executable_sha256"):
        declared = segment.get(key)
        if declared is not None and manifest.get(key) != declared:
            raise CaseAuthenticationError(
                f"{case_id} direct-fast segment {order} manifest {key} differs"
            )
    if nodes * ranks_per_node != ranks:
        raise CaseAuthenticationError(
            f"{case_id} direct-fast segment {order} allocation rank count differs"
        )
    return {
        "source_kind": "direct_fast",
        "lineage_order": order,
        "segment": segment.get("segment"),
        "output": str(output),
        "rank_count": ranks,
        "manifest": manifest_binding,
    }


def r02_authority_chain_record(
    acceptance: object, policy: dict[str, object]
) -> dict[str, object]:
    """Record the reviewed source-authority result without making it numerical gate."""

    try:
        _, authority, bindings = acceptance.current_source_archive_catalog(policy)
    except Exception as error:
        return {
            "status": "blocked",
            "campaign_authority_eligible": False,
            "release_authorizing": False,
            "blocker": str(error),
        }
    return {
        "status": "authenticated_non_authorizing",
        "campaign_authority_eligible": False,
        "release_authorizing": False,
        "authority": authority,
        "bindings": bindings,
        "blocker": None,
    }


def declared_r02_rank_files(
    case_id: str,
    terminal: dict[str, object],
    segment_dir: Path,
    expected_rank_count: int,
) -> list[dict[str, object]]:
    """Validate one exact accepted-segment restart declaration without reading bytes."""

    if terminal.get("storage") != "per_rank":
        raise CaseAuthenticationError(f"{case_id} terminal restart storage is not per_rank")
    declared = require_list(
        terminal.get("rank_files"), f"{case_id} terminal restart rank files"
    )
    if len(declared) != expected_rank_count:
        raise CaseAuthenticationError(
            f"{case_id} terminal restart rank count differs from accepted allocation"
        )
    records: list[dict[str, object]] = []
    paths: set[str] = set()
    digests: set[str] = set()
    expected_root = segment_dir / "output" / "rst"
    for rank, value in enumerate(declared):
        item = require_dict(value, f"{case_id} terminal restart rank {rank}")
        path = Path(require_text(item.get("path"), "terminal restart rank path")).absolute()
        expected_parent = expected_root / f"rank_{rank:08d}"
        if path.parent != expected_parent or path.suffix != ".rst":
            raise CaseAuthenticationError(
                f"{case_id} terminal restart rank {rank} path is noncanonical"
            )
        digest = require_sha256(
            item.get("sha256"), f"{case_id} terminal restart rank {rank} SHA-256"
        )
        size = require_int(
            item.get("size_bytes"), f"{case_id} terminal restart rank {rank} size", 1
        )
        if str(path) in paths or digest in digests:
            raise CaseAuthenticationError(
                f"{case_id} terminal restart declaration duplicates path or rank bytes"
            )
        paths.add(str(path))
        digests.add(digest)
        records.append({
            "path": str(path),
            "sha256": digest,
            "size_bytes": size,
            "rank": rank,
        })
    representative = {
        "path": require_text(terminal.get("path"), f"{case_id} terminal restart path"),
        "sha256": require_sha256(
            terminal.get("sha256"), f"{case_id} terminal restart SHA-256"
        ),
        "size_bytes": require_int(
            terminal.get("size_bytes"), f"{case_id} terminal restart size", 1
        ),
    }
    if representative != {
        key: records[0][key] for key in ("path", "sha256", "size_bytes")
    }:
        raise CaseAuthenticationError(
            f"{case_id} terminal restart representative differs from rank zero"
        )
    return records


def authenticate_accepted_r02_bundle(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    case_name: object,
    value: object,
    matrix_sha256: str,
    required_times: list[float],
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Authenticate selected accepted-R02 bundle lineage and exact required states."""

    segment = require_dict(value, f"{case_id} accepted R02 lineage segment")
    if case_id != "R02" or segment.get("kind") != "accepted_r02_bundle":
        raise CaseAuthenticationError("lineage segment is not accepted R02 bundle")
    order = require_int(segment.get("order"), f"{case_id} accepted bundle order")
    if segment.get("state") != "accepted_for_analysis":
        raise CaseAuthenticationError(f"{case_id} accepted bundle lineage state differs")
    bundle_root = Path(
        require_text(segment.get("output"), f"{case_id} accepted bundle output")
    ).absolute()
    manifest_record = require_dict(
        segment.get("manifest"), f"{case_id} accepted bundle manifest"
    )
    manifest_path = Path(
        require_text(manifest_record.get("path"), f"{case_id} accepted bundle manifest path")
    ).absolute()
    if manifest_path != bundle_root / "manifest.json":
        raise CaseAuthenticationError(f"{case_id} accepted bundle manifest path differs")
    bundle, bundle_binding = read_bound_json(
        manifest_path,
        f"{case_id} accepted whole-case bundle",
        expected_sha256=require_sha256(
            manifest_record.get("sha256"), f"{case_id} accepted bundle SHA-256"
        ),
    )
    if int(bundle_binding["size_bytes"]) != require_int(
        manifest_record.get("size_bytes"), f"{case_id} accepted bundle size"
    ):
        raise CaseAuthenticationError(f"{case_id} accepted bundle size differs")

    expected_case = acceptance.case_manifest_record(policy, case_id)
    selected_cases = require_list(bundle.get("cases"), f"{case_id} accepted bundle cases")
    if (
        bundle.get("workflow") != "paper-mks24-stage-i-production"
        or bundle.get("status") != "accepted_for_analysis"
        or bundle.get("production_case_id") != case_id
        or require_finite(bundle.get("required_final_time"), "bundle required final time")
        != 10.0
        or require_finite(bundle.get("accepted_final_time"), "bundle accepted final time")
        != 10.0
        or len(selected_cases) != 1
    ):
        raise CaseAuthenticationError(f"{case_id} accepted whole-case bundle identity differs")
    selected_case = require_dict(selected_cases[0], f"{case_id} accepted bundle case")
    if (
        case_name != expected_case.get("name")
        or selected_case.get("name") != expected_case.get("name")
        or selected_case.get("input") != expected_case.get("input")
        or selected_case.get("status") != "passed"
    ):
        raise CaseAuthenticationError(f"{case_id} accepted bundle selected case differs")
    expected_matrix = require_sha256(
        policy["verified_sources"]["stage_i_manifest"]["sha256"],
        "reviewed Stage I manifest SHA-256",
    )
    if matrix_sha256 != expected_matrix:
        raise CaseAuthenticationError(f"{case_id} assembled matrix differs from reviewed policy")
    rank_count = require_int(segment.get("ranks"), f"{case_id} accepted bundle ranks", 1)
    epoch_root = bundle_root.parent.parent
    case_root = epoch_root / case_id
    segment_paths: list[Path] = []
    for item in require_list(
        bundle.get("production_segment_manifests"),
        f"{case_id} accepted production segment manifests",
    ):
        path = Path(require_text(item, f"{case_id} accepted segment path")).absolute()
        try:
            relative = path.relative_to(case_root)
        except ValueError as error:
            raise CaseAuthenticationError(
                f"{case_id} accepted segment path escapes the production case root"
            ) from error
        if len(relative.parts) != 3 or relative.parts[1:] != (
            "manifest", "prepared_run.json"
        ):
            raise CaseAuthenticationError(f"{case_id} accepted segment path is noncanonical")
        segment_paths.append(path)
    if not segment_paths or len(segment_paths) != len(set(segment_paths)):
        raise CaseAuthenticationError(
            f"{case_id} accepted segment inventory is empty or duplicated"
        )

    accepted_segments: list[dict[str, object]] = []
    candidates: list[dict[str, object]] = []
    final_times: list[float] = []
    executables: set[str] = set()
    for segment_path in segment_paths:
        prepared, prepared_binding = read_bound_json(
            segment_path, f"{case_id} accepted production segment"
        )
        accounting = require_dict(
            prepared.get("accounting"), f"{case_id} accepted segment accounting"
        )
        command = require_dict(prepared.get("command"), f"{case_id} accepted segment command")
        inspection = require_dict(
            prepared.get("scientific_inspection"), f"{case_id} segment inspection"
        )
        allocation = require_dict(
            prepared.get("allocation"), f"{case_id} accepted segment allocation"
        )
        final_time = require_finite(
            inspection.get("final_time"), f"{case_id} accepted segment final time"
        )
        executable = require_sha256(
            command.get("executable_sha256"), f"{case_id} accepted executable SHA-256"
        )
        segment_rank_count = require_int(
            allocation.get("nodes"), f"{case_id} accepted segment nodes", 1
        ) * require_int(
            allocation.get("ranks_per_node"), f"{case_id} accepted ranks per node", 1
        )
        if (
            accounting.get("result") != "accepted"
            or accounting.get("case_id") != case_id
            or accounting.get("case_name") != expected_case.get("name")
            or accounting.get("executable_sha256") != executable
            or command.get("matrix_sha256") != matrix_sha256
            or inspection.get("schema_version") != 3
            or inspection.get("accepted") is not True
            or inspection.get("case_id") != case_id
            or Path(require_text(inspection.get("manifest"), "inspection manifest")).absolute()
            != segment_path
            or segment_rank_count != rank_count
        ):
            raise CaseAuthenticationError(
                f"{case_id} whole-case bundle segment is not an accepted CT source"
            )
        segment_name = require_text(
            accounting.get("segment"), f"{case_id} accepted segment name"
        )
        accepted_segments.append({
            "segment": segment_name,
            "time": final_time,
            "rank_count": segment_rank_count,
            "manifest": prepared_binding,
        })
        final_times.append(final_time)
        executables.add(executable)
        if final_time not in required_times:
            continue
        required_time = final_time
        checks = require_dict(inspection.get("checks"), f"{case_id} required state checks")
        terminal = require_dict(
            inspection.get("terminal_restart"), f"{case_id} required terminal restart"
        )
        restart_times = [
            require_finite(item, f"{case_id} retained restart time")
            for item in require_list(
                inspection.get("restart_times"), f"{case_id} retained restart times"
            )
        ]
        if (
            require_finite(inspection.get("required_time"), "required state time")
            != required_time
            or require_finite(
                inspection.get("terminal_restart_time"), "terminal restart time"
            )
            != required_time
            or inspection.get("segment") != segment_name
            or checks.get("required_time_reached") is not True
            or checks.get("restart_retained") is not True
            or checks.get("terminal_restart_physical_time_matches_final") is not True
            or required_time not in restart_times
            or terminal not in require_list(
                inspection.get("restarts"), f"{case_id} retained restarts"
            )
        ):
            raise CaseAuthenticationError(
                f"{case_id} accepted segment lacks exact t={required_time:g} terminal restart"
            )
        candidates.append({
            "source_kind": "accepted_r02_bundle",
            "lineage_order": order,
            "segment": segment_name,
            "name": Path(require_text(terminal.get("path"), "terminal restart path")).name,
            "rank_count": rank_count,
            "time": required_time,
            "declared_rank_files": declared_r02_rank_files(
                case_id, terminal, segment_path.parent.parent, rank_count
            ),
            "accepted_bundle_manifest": bundle_binding,
            "accepted_segment_manifest": prepared_binding,
        })
    if (
        final_times != sorted(final_times)
        or len(final_times) != len(set(final_times))
        or final_times[-1] != 10.0
        or len(executables) != 1
    ):
        raise CaseAuthenticationError(f"{case_id} accepted segment lineage differs")
    for required_time in required_times:
        if sum(float(item["time"]) == required_time for item in candidates) != 1:
            raise CaseAuthenticationError(
                f"{case_id} accepted bundle must declare exactly one t={required_time:g} state"
            )
    return {
        "source_kind": "accepted_r02_bundle",
        "lineage_order": order,
        "segment": segment.get("segment"),
        "output": str(bundle_root),
        "rank_count": rank_count,
        "manifest": bundle_binding,
        "accepted_segment_count": len(accepted_segments),
        "accepted_segments": accepted_segments,
    }, candidates


def compact_native_state(
    state: dict[str, object], source: dict[str, object]
) -> dict[str, object]:
    """Retain deterministic native CT evidence without repeated large location lists."""

    audits: list[dict[str, object]] = []
    for value in require_list(state.get("restart_audits"), "native restart audits"):
        audit = require_dict(value, "native restart audit")
        audits.append({
            key: audit.get(key)
            for key in (
                "path",
                "sha256",
                "size_bytes",
                "rank",
                "time",
                "nmb_total",
                "local_meshblocks",
                "data_size_per_meshblock",
                "global_logical_locations_sha256",
                "restart_contract_sha256",
                "maximum_normalized_ct_divb",
                "passed",
            )
        })
    return {
        key: state.get(key)
        for key in (
            "time",
            "rank_count",
            "state_meshblocks",
            "sampled_meshblocks",
            "meshblock_coverage_complete",
            "global_logical_locations_sha256",
            "restart_contract_sha256",
            "maximum_normalized_ct_divb",
        )
    } | {
        "source": source,
        "restart_audits": audits,
    }


def audit_selected_lineage_restarts(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    case: dict[str, object],
    matrix_sha256: str,
    required_times: list[float],
    threshold: float,
) -> dict[str, object]:
    """Audit exact required native restart states from authenticated selected lineage."""

    criteria = require_dict(policy.get("criteria"), "native CT criteria")
    ct_policy = require_dict(criteria.get("ct_divb_policy"), "native CT policy")
    blockers: list[str] = []
    source_errors: list[str] = []
    fast_segments: list[dict[str, object]] = []
    accepted_bundles: list[dict[str, object]] = []
    accepted_candidates: list[dict[str, object]] = []
    authority_chain: dict[str, object] = {
        "status": "not_applicable",
        "campaign_authority_eligible": False,
        "release_authorizing": False,
        "blocker": None,
    }
    for value in require_list(case.get("lineage"), f"{case_id} selected lineage"):
        segment = require_dict(value, f"{case_id} lineage segment")
        try:
            if segment.get("kind") == "fast":
                fast_segments.append(
                    authenticate_direct_fast_segment(
                        case_id, case.get("case_name"), segment, matrix_sha256
                    )
                )
            elif segment.get("kind") == "accepted_r02_bundle":
                bundle, states = authenticate_accepted_r02_bundle(
                    acceptance,
                    policy,
                    case_id,
                    case.get("case_name"),
                    segment,
                    matrix_sha256,
                    required_times,
                )
                accepted_bundles.append(bundle)
                accepted_candidates.extend(states)
        except (FastCtAuditError, OSError) as error:
            source_errors.append(str(error))
    if accepted_bundles:
        authority_chain = r02_authority_chain_record(acceptance, policy)
    if not fast_segments and not accepted_bundles:
        blockers.append(
            "selected assembled lineage contains no authenticated direct-fast or "
            "accepted-R02 source that can supply native restart state"
        )

    discovered: list[dict[str, object]] = list(accepted_candidates)
    incomplete: list[dict[str, object]] = []
    candidates: dict[float, list[dict[str, object]]] = {
        value: [] for value in required_times
    }
    for record in accepted_candidates:
        candidates[float(record["time"])].append(record)
    for segment in fast_segments:
        rank_count = int(segment["rank_count"])
        root = Path(str(segment["output"])) / "rst"
        if not root.is_dir():
            blockers.append(
                f"direct-fast segment {segment['segment']} has no retained rst directory"
            )
            continue
        actual_dirs = sorted(
            path for path in root.glob("rank_*")
            if path.is_dir() and RANK_DIRECTORY.fullmatch(path.name)
        )
        expected_dirs = [root / f"rank_{rank:08d}" for rank in range(rank_count)]
        if actual_dirs != expected_dirs:
            blockers.append(
                f"direct-fast segment {segment['segment']} restart rank directories "
                "are not exact and contiguous"
            )
            continue
        names = sorted({
            path.name
            for directory in actual_dirs
            for path in directory.glob("*.rst")
            if path.is_file()
        })
        for name in names:
            paths = [directory / name for directory in actual_dirs]
            missing = [str(path) for path in paths if not path.is_file()]
            empty = [
                str(path) for path in paths if path.is_file() and path.stat().st_size == 0
            ]
            identity = {
                "source_kind": "direct_fast",
                "lineage_order": segment["lineage_order"],
                "segment": segment["segment"],
                "name": name,
                "rank_count": rank_count,
            }
            if missing or empty:
                incomplete.append({
                    **identity,
                    "missing_rank_files": missing,
                    "empty_rank_files": empty,
                })
                continue
            try:
                time_value = native_restart_time(acceptance, paths[0])
            except Exception as error:
                incomplete.append({
                    **identity,
                    "missing_rank_files": [],
                    "empty_rank_files": [],
                    "error": str(error),
                })
                continue
            record = {
                **identity,
                "time": time_value,
                "paths": [str(path.absolute()) for path in paths],
            }
            discovered.append(record)
            if time_value in candidates:
                candidates[time_value].append(record)

    state_audits: list[dict[str, object]] = []
    selected_sources: list[dict[str, object]] = []
    audit_errors: list[str] = []
    duplicate_candidates: list[dict[str, object]] = []
    for required_time in required_times:
        matches = candidates[required_time]
        if len(matches) > 1:
            duplicate_candidates.append({
                "time": required_time,
                "candidate_count": len(matches),
                "selected": {
                    key: matches[-1][key]
                    for key in (
                        "source_kind", "lineage_order", "segment", "name", "rank_count"
                    )
                },
            })
            if any(item.get("source_kind") == "accepted_r02_bundle" for item in matches):
                audit_errors.append(
                    f"{case_id} accepted R02 lineage does not supply exactly one "
                    f"t={required_time:g} restart group"
                )
                continue
        if not matches:
            blockers.append(
                f"exact required native restart state t={required_time:g} is not retained "
                "by an authenticated selected-lineage source"
            )
            continue
        selected = matches[-1]
        selected_source = {
            key: selected[key]
            for key in (
                "source_kind", "lineage_order", "segment", "name", "rank_count", "time"
            )
        }
        for key in ("accepted_bundle_manifest", "accepted_segment_manifest"):
            if key in selected:
                selected_source[key] = selected[key]
        try:
            declared = selected.get("declared_rank_files")
            if declared is None:
                bindings = [
                    regular_file_binding(
                        Path(path),
                        f"{case_id} native restart t={required_time:g} rank {rank}",
                    )
                    for rank, path in enumerate(selected["paths"])
                ]
            else:
                bindings = []
                for rank, value in enumerate(
                    require_list(declared, f"{case_id} declared restart ranks")
                ):
                    item = require_dict(value, f"{case_id} declared restart rank {rank}")
                    bindings.append(regular_file_binding(
                        Path(require_text(item.get("path"), "declared restart path")),
                        f"{case_id} accepted native restart t={required_time:g} rank {rank}",
                        expected_sha256=require_sha256(
                            item.get("sha256"), "declared restart SHA-256"
                        ),
                        expected_size=require_int(
                            item.get("size_bytes"), "declared restart size", 1
                        ),
                    ))
            records = [
                acceptance.audit_restart_file(
                    Path(str(binding["path"])),
                    str(binding["sha256"]),
                    threshold,
                    expected_rank=rank,
                    expected_rank_count=int(selected["rank_count"]),
                )
                for rank, binding in enumerate(bindings)
            ]
            state = acceptance.validate_ct_state_records(
                records, required_time, require_complete=True
            )
            state_audits.append(compact_native_state(state, selected_source))
            selected_sources.append(selected_source)
        except Exception as error:
            audit_errors.append(
                f"{case_id} exact native restart state t={required_time:g} audit failed: {error}"
            )
    audited_times = [float(state["time"]) for state in state_audits]
    missing_times = [value for value in required_times if value not in audited_times]
    maximum = (
        max(float(state["maximum_normalized_ct_divb"]) for state in state_audits)
        if state_audits else None
    )
    restart_count = sum(int(state["rank_count"]) for state in state_audits)
    sampled_meshblocks = sum(int(state["sampled_meshblocks"]) for state in state_audits)
    state_meshblocks = sum(int(state["state_meshblocks"]) for state in state_audits)
    numerical = (
        "inconclusive"
        if maximum is None else "pass" if maximum < threshold else "fail"
    )
    authentication_failed = bool(source_errors or audit_errors)
    coverage_complete = audited_times == required_times and not authentication_failed
    result = (
        "fail"
        if numerical == "fail" else "pass"
        if numerical == "pass" and coverage_complete else "inconclusive"
    )
    authority_blockers = (
        [str(authority_chain["blocker"])] if authority_chain.get("blocker") else []
    )
    accepted_source = bool(accepted_bundles)
    return {
        "status": (
            "authentication_failed"
            if authentication_failed else "complete"
            if coverage_complete else "partial"
        ),
        "result": result,
        "numerical_result": numerical,
        "coverage_complete": coverage_complete,
        "campaign_authority_eligible": False,
        "release_authorizing": False,
        "claim_scope": ct_policy.get("claim_scope"),
        "historical_limitation": ct_policy.get("historical_limitation"),
        "required_state_times": required_times,
        "audited_state_times": audited_times,
        "missing_required_state_times": missing_times,
        "restart_count": restart_count,
        "sampled_meshblocks": sampled_meshblocks,
        "state_meshblocks": state_meshblocks,
        "sampled_meshblock_fraction": (
            sampled_meshblocks / state_meshblocks if state_meshblocks else None
        ),
        "meshblock_coverage_complete": bool(state_audits) and all(
            bool(state["meshblock_coverage_complete"]) for state in state_audits
        ),
        "normalized_ct_divb_lt": threshold,
        "maximum_normalized_ct_divb": maximum,
        "ct_evidence_available": bool(state_audits),
        "ct_claim_supported": result in ("pass", "fail"),
        "native_parser": (
            "cgl_lf_stage_i_scientific_acceptance.audit_restart_file plus "
            "validate_ct_state_records(require_complete=True)"
        ),
        "public_audit_ct_divb_without_inventory": {
            "can_consume_rank_file_bindings": True,
            "can_establish_complete_multirank_state_coverage": False,
            "reason": (
                "free-form audit-ct-divb has no expected rank-count authority; this "
                "adapter supplies the authenticated selected-lineage rank count directly "
                "to the same reviewed native parser and coverage validator"
            ),
        },
        "authority_chain": authority_chain,
        "authority_blockers": authority_blockers,
        "authority_limitation": (
            ACCEPTED_R02_AUTHORITY_LIMITATION
            if accepted_source else DIRECT_FAST_AUTHORITY_LIMITATION
        ),
        "restart_inventory_limitation": (
            ACCEPTED_R02_RESTART_INDEX if accepted_source else RESTART_INDEX_LIMITATION
        ),
        "authenticated_fast_segments": fast_segments,
        "authenticated_accepted_bundles": accepted_bundles,
        "discovered_state_times": sorted({
            float(record["time"]) for record in discovered
        }),
        "discovered_complete_group_count": len(discovered),
        "incomplete_groups": incomplete,
        "duplicate_required_state_candidates": duplicate_candidates,
        "selected_state_sources": selected_sources,
        "state_audits": state_audits,
        "blockers": blockers,
        "errors": [*source_errors, *audit_errors],
    }


def header_parameter(header: list[str], block_name: str, key_name: str) -> str | None:
    block = ""
    for line in header:
        if line.startswith("<"):
            block = line.strip()
            continue
        if block != block_name or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == key_name:
            return value.strip()
    return None


def binary_snapshot_metadata(
    path: Path, label: str, expected_size: int
) -> tuple[dict[str, object], dict[str, object]]:
    """Bind and inspect Athena binary metadata without reading field arrays."""

    before = regular_file_binding(path, label, expected_size=expected_size)
    resolved = Path(str(before["path"]))
    try:
        with resolved.open("rb") as stream:
            code_header = stream.readline(4096).split()
            if not code_header or code_header[0] != b"Athena":
                raise FastCtAuditError(f"{label} has an invalid Athena binary header")
            if code_header[-1].split(b"=")[-1] != b"1.1":
                raise FastCtAuditError(f"{label} has an unsupported binary version")
            count_line = stream.readline(4096)
            count = int(count_line.split(b"=")[-1])
            if count < 1 or count > 1024:
                raise FastCtAuditError(f"{label} has an invalid preheader count")
            preheader: dict[str, str] = {}
            for _ in range(count - 1):
                line = stream.readline(4096).decode("utf-8")
                if "=" not in line:
                    raise FastCtAuditError(f"{label} has a malformed preheader")
                key, value = line.split("=", 1)
                preheader[key.strip()] = value.strip()
            nvars = int(stream.readline(4096).split(b"=")[-1])
            variables = [
                value.decode("utf-8") for value in stream.readline(65536).split()[1:]
            ]
            if nvars != len(variables) or len(set(variables)) != len(variables):
                raise FastCtAuditError(f"{label} has an invalid variable inventory")
            header_size = int(stream.readline(4096).split(b"=")[-1])
            if header_size < 0 or header_size > MAX_HEADER_BYTES:
                raise FastCtAuditError(f"{label} has an invalid parameter-header size")
            payload = stream.read(header_size)
            if len(payload) != header_size:
                raise FastCtAuditError(f"{label} has a truncated parameter header")
            header = [
                line.decode("utf-8").split("#", 1)[0].strip()
                for line in payload.split(b"\n")
            ]
            header = [line for line in header if line]
    except (OSError, UnicodeDecodeError, ValueError) as error:
        raise FastCtAuditError(f"{label} metadata parse failed: {error}") from error
    after = regular_file_binding(
        path,
        f"{label} after metadata parse",
        expected_sha256=str(before["sha256"]),
        expected_size=int(before["size_bytes"]),
    )
    if before != after:
        raise FastCtAuditError(f"{label} changed while metadata was parsed")
    metadata: dict[str, object] = {
        "time": require_finite(float(preheader["time"]), f"{label} time"),
        "cycle": int(preheader["cycle"]),
        "variable_count": nvars,
        "variables": variables,
        "magnetic_fields_present": [
            field for field in CELL_CENTERED_MAGNETIC_FIELDS if field in variables
        ],
        "cell_centered_magnetic_triplet_complete": all(
            field in variables for field in CELL_CENTERED_MAGNETIC_FIELDS
        ),
        "face_centered_ct_state_available": False,
    }
    for block, prefix in (("<mesh>", "mesh"), ("<meshblock>", "meshblock")):
        for axis in (1, 2, 3):
            value = header_parameter(header, block, f"nx{axis}")
            if value is not None:
                metadata[f"{prefix}_nx{axis}"] = int(value)
    return metadata, before


def verify_case_lineage(
    case_id: str,
    case: dict[str, object],
    case_dir: Path,
    matrix_sha256: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Authenticate the inventory-selected assembled case provenance."""

    lineage_path = case_dir / "lineage.json"
    lineage, lineage_binding = read_bound_json(lineage_path, f"{case_id} lineage")
    if lineage != case:
        raise CaseAuthenticationError(
            f"{case_id} inventory record differs from bound lineage.json"
        )
    if lineage.get("case_id") != case_id:
        raise CaseAuthenticationError(f"{case_id} lineage case identity differs")
    retained_errors = lineage.get("errors")
    if not isinstance(retained_errors, list):
        raise CaseAuthenticationError(f"{case_id} lineage errors field is malformed")
    if retained_errors:
        raise CaseAuthenticationError(f"{case_id} lineage retains assembly errors")
    verified: list[dict[str, object]] = []
    input_record = lineage.get("input")
    if isinstance(input_record, dict) and input_record.get("sha256") is not None:
        verified.append(verify_declared_binding(input_record, f"{case_id} input"))
    histories = lineage.get("histories")
    if isinstance(histories, dict):
        for kind in ("mhd", "user"):
            history = histories.get(kind)
            if not isinstance(history, dict) or not history.get("available"):
                continue
            binding = history.get("binding")
            if binding is not None:
                verified.append(
                    verify_declared_binding(binding, f"{case_id} merged {kind} history")
                )
    segments = require_list(lineage.get("lineage"), f"{case_id} selected lineage")
    if not segments:
        raise CaseAuthenticationError(f"{case_id} selected lineage is empty")
    for order, value in enumerate(segments):
        segment = require_dict(value, f"{case_id} lineage segment {order}")
        if segment.get("order") != order:
            raise CaseAuthenticationError(f"{case_id} lineage segment order differs")
        manifest = segment.get("manifest")
        if manifest is None:
            raise CaseAuthenticationError(f"{case_id} lineage segment lacks manifest")
        verified.append(
            verify_declared_binding(manifest, f"{case_id} lineage segment {order} manifest")
        )
        declared_matrix = segment.get("matrix_sha256")
        if declared_matrix is not None and declared_matrix != matrix_sha256:
            raise CaseAuthenticationError(
                f"{case_id} lineage segment {order} matrix SHA-256 differs"
            )
    identities = lineage.get("lineage_identities")
    if identities is not None:
        identity_record = require_dict(identities, f"{case_id} lineage identities")
        matrices = require_list(
            identity_record.get("matrix_sha256"), f"{case_id} matrix identities"
        )
        if matrices != [matrix_sha256]:
            raise CaseAuthenticationError(f"{case_id} matrix lineage identity differs")
    return lineage_binding, verified


def validate_snapshot_group(
    case_id: str, value: object, position: int
) -> tuple[dict[str, object], bool]:
    group = require_dict(value, f"{case_id} snapshot {position}")
    complete = group.get("complete")
    if not isinstance(complete, bool):
        raise CaseAuthenticationError(f"{case_id} snapshot {position} complete flag is malformed")
    require_finite(group.get("time"), f"{case_id} snapshot {position} time")
    expected = require_int(
        group.get("expected_ranks"), f"{case_id} snapshot {position} expected ranks", 1
    )
    members = require_list(group.get("rank_files"), f"{case_id} snapshot {position} ranks")
    if len(members) != expected:
        raise CaseAuthenticationError(f"{case_id} snapshot {position} rank count differs")
    paths: list[Path] = []
    ranks: list[int] = []
    for rank, member_value in enumerate(members):
        member = require_dict(member_value, f"{case_id} snapshot {position} rank {rank}")
        path = Path(require_text(member.get("path"), "snapshot rank path")).absolute()
        match = RANK_DIRECTORY.fullmatch(path.parent.name)
        if match is None:
            raise CaseAuthenticationError(
                f"{case_id} snapshot {position} rank path is noncanonical"
            )
        paths.append(path)
        ranks.append(int(match.group(1)))
        size = member.get("size_bytes")
        if size is not None:
            require_int(size, f"{case_id} snapshot {position} rank {rank} size", 1)
    if ranks != list(range(expected)) or len(set(paths)) != expected:
        raise CaseAuthenticationError(
            f"{case_id} snapshot {position} rank inventory is not exact and contiguous"
        )
    representative = Path(
        require_text(group.get("representative"), f"{case_id} snapshot representative")
    ).absolute()
    if representative != paths[0]:
        raise CaseAuthenticationError(
            f"{case_id} snapshot {position} representative is not rank zero"
        )
    return group, complete


def inspect_snapshot_group(
    case_id: str, group: dict[str, object], position: int
) -> dict[str, object]:
    """Authenticate one complete group and classify its available CT evidence."""

    expected_time = require_finite(group.get("time"), f"{case_id} snapshot time")
    metadata_records: list[dict[str, object]] = []
    bindings: list[dict[str, object]] = []
    for rank, value in enumerate(require_list(group["rank_files"], "snapshot ranks")):
        member = require_dict(value, f"{case_id} snapshot rank {rank}")
        size = require_int(member.get("size_bytes"), f"{case_id} snapshot rank size", 1)
        metadata, binding = binary_snapshot_metadata(
            Path(require_text(member.get("path"), "snapshot rank path")),
            f"{case_id} snapshot {position} rank {rank}",
            size,
        )
        if not math.isclose(
            float(metadata["time"]), expected_time, rel_tol=0.0, abs_tol=1.0e-12
        ):
            raise CaseAuthenticationError(
                f"{case_id} snapshot {position} binary time differs from index"
            )
        metadata_records.append(metadata)
        bindings.append(binding)
    reference = metadata_records[0]
    for rank, metadata in enumerate(metadata_records[1:], start=1):
        if metadata != reference:
            raise CaseAuthenticationError(
                f"{case_id} snapshot {position} rank metadata differs at rank {rank}"
            )
    return {
        "time": expected_time,
        "segment": group.get("segment"),
        "lineage_order": group.get("lineage_order"),
        "rank_count": len(bindings),
        "authenticated_bytes": sum(int(binding["size_bytes"]) for binding in bindings),
        "rank_files": bindings,
        "metadata": reference,
        "inspection_scope": "authenticated Athena binary metadata only",
        "ct_evidence_available": False,
        "ct_divb_measurement": None,
        "ct_claim_supported": False,
        "reason": CT_ABSENCE_REASON,
    }


def audit_snapshot_index(
    case_id: str,
    index_path: Path,
    expected_path: Path,
    expected_count: int,
    expected_complete_count: int,
    snapshot_policy: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    if index_path.absolute() != expected_path.absolute():
        raise CaseAuthenticationError(f"{case_id} snapshot index path differs")
    index, index_binding = read_bound_json(index_path, f"{case_id} snapshot index")
    if index.get("schema_version") != 1:
        raise CaseAuthenticationError(f"{case_id} snapshot index schema differs")
    snapshots = require_list(index.get("snapshots"), f"{case_id} snapshots")
    if require_int(index.get("snapshot_count"), f"{case_id} snapshot count") != len(snapshots):
        raise CaseAuthenticationError(f"{case_id} snapshot count differs")
    if len(snapshots) != expected_count:
        raise CaseAuthenticationError(f"{case_id} assembled snapshot count differs")
    validated: list[tuple[dict[str, object], bool]] = [
        validate_snapshot_group(case_id, value, position)
        for position, value in enumerate(snapshots)
    ]
    complete = [item for item, is_complete in validated if is_complete]
    if len(complete) != expected_complete_count:
        raise CaseAuthenticationError(f"{case_id} complete snapshot count differs")
    times = [require_finite(item.get("time"), f"{case_id} snapshot time") for item, _ in validated]
    if any(right <= left for left, right in zip(times, times[1:])):
        raise CaseAuthenticationError(f"{case_id} snapshot times are not strictly increasing")
    selected = complete[-1:] if snapshot_policy == "latest" and complete else complete
    audits = [
        inspect_snapshot_group(case_id, group, snapshots.index(group))
        for group in selected
    ]
    summary = {
        "binding": index_binding,
        "snapshot_count": len(snapshots),
        "complete_snapshot_count": len(complete),
        "incomplete_snapshot_count": len(snapshots) - len(complete),
        "selected_snapshot_count": len(selected),
        "authenticated_snapshot_count": len(audits),
        "authenticated_rank_file_count": sum(
            int(record["rank_count"]) for record in audits
        ),
        "authenticated_bytes": sum(int(record["authenticated_bytes"]) for record in audits),
        "time_first": times[0] if times else None,
        "time_last": times[-1] if times else None,
    }
    return summary, audits


def empty_native_restart_record(
    required_times: list[float], reason: str, threshold: float | None = None
) -> dict[str, object]:
    return {
        "status": "partial",
        "result": "inconclusive",
        "numerical_result": "inconclusive",
        "coverage_complete": False,
        "campaign_authority_eligible": False,
        "release_authorizing": False,
        "required_state_times": required_times,
        "audited_state_times": [],
        "missing_required_state_times": required_times,
        "restart_count": 0,
        "sampled_meshblocks": 0,
        "state_meshblocks": 0,
        "sampled_meshblock_fraction": None,
        "meshblock_coverage_complete": False,
        "normalized_ct_divb_lt": threshold,
        "maximum_normalized_ct_divb": None,
        "ct_evidence_available": False,
        "ct_claim_supported": False,
        "native_parser": (
            "cgl_lf_stage_i_scientific_acceptance.audit_restart_file plus "
            "validate_ct_state_records(require_complete=True)"
        ),
        "public_audit_ct_divb_without_inventory": {
            "can_consume_rank_file_bindings": True,
            "can_establish_complete_multirank_state_coverage": False,
            "reason": (
                "free-form audit-ct-divb has no expected rank-count authority; this "
                "adapter supplies the authenticated selected-lineage rank count directly "
                "to the same reviewed native parser and coverage validator"
            ),
        },
        "authority_limitation": DIRECT_FAST_AUTHORITY_LIMITATION,
        "restart_inventory_limitation": RESTART_INDEX_LIMITATION,
        "authority_chain": {
            "status": "not_applicable",
            "campaign_authority_eligible": False,
            "release_authorizing": False,
            "blocker": None,
        },
        "authority_blockers": [],
        "authenticated_fast_segments": [],
        "authenticated_accepted_bundles": [],
        "discovered_state_times": [],
        "discovered_complete_group_count": 0,
        "incomplete_groups": [],
        "duplicate_required_state_candidates": [],
        "selected_state_sources": [],
        "state_audits": [],
        "blockers": [reason],
        "errors": [],
    }


def missing_case_record(
    case_id: str, required_times: list[float], threshold: float
) -> dict[str, object]:
    return {
        "case_id": case_id,
        "assembled_status": "not_present",
        "audit_status": "partial_inconclusive",
        "provenance_authenticated": False,
        "snapshot_content_authenticated": False,
        "ct_result": "inconclusive",
        "ct_evidence_available": False,
        "ct_claim_supported": False,
        "reason": "case is not present in the partial assembled inventory",
        "errors": [],
        "snapshot_summary": {
            "snapshot_count": 0,
            "complete_snapshot_count": 0,
            "incomplete_snapshot_count": 0,
            "selected_snapshot_count": 0,
            "authenticated_snapshot_count": 0,
            "authenticated_rank_file_count": 0,
            "authenticated_bytes": 0,
            "time_first": None,
            "time_last": None,
        },
        "snapshot_audits": [],
        "native_restart_ct": empty_native_restart_record(
            required_times,
            "case is not present in the partial assembled inventory",
            threshold,
        ),
    }


def audit_case(
    acceptance: object,
    policy: dict[str, object],
    case_id: str,
    value: object,
    output_root: Path,
    matrix_sha256: str,
    snapshot_policy: str,
    required_times: list[float],
    threshold: float,
) -> dict[str, object]:
    case = require_dict(value, f"{case_id} assembled case")
    base: dict[str, object] = {
        "case_id": case_id,
        "case_name": case.get("case_name"),
        "assembled_status": case.get("status"),
        "ct_result": "inconclusive",
        "ct_evidence_available": False,
        "ct_claim_supported": False,
        "reason": CT_ABSENCE_REASON,
        "errors": [],
        "snapshot_audits": [],
        "native_restart_ct": empty_native_restart_record(
            required_times,
            "assembled-case provenance has not yet been authenticated",
            threshold,
        ),
    }
    try:
        case_dir = output_root / "cases" / case_id
        lineage_binding, verified = verify_case_lineage(
            case_id, case, case_dir, matrix_sha256
        )
        snapshot_record = require_dict(case.get("snapshots"), f"{case_id} snapshots record")
        snapshot_path = Path(
            require_text(snapshot_record.get("path"), f"{case_id} snapshot index path")
        )
        expected_snapshot_path = case_dir / "snapshots.json"
        summary, audits = audit_snapshot_index(
            case_id,
            snapshot_path,
            expected_snapshot_path,
            require_int(snapshot_record.get("snapshot_count"), f"{case_id} snapshot count"),
            require_int(
                snapshot_record.get("complete_snapshot_count"),
                f"{case_id} complete snapshot count",
            ),
            snapshot_policy,
        )
        native = audit_selected_lineage_restarts(
            acceptance,
            policy,
            case_id,
            case,
            matrix_sha256,
            required_times,
            threshold,
        )
        native_status = str(native["status"])
        ct_result = str(native["result"])
        native_errors = [
            str(error) for error in require_list(native.get("errors"), "native CT errors")
        ]
        reason = (
            (
                "exact authenticated selected-lineage native restart states t=9,10 are "
                "below the reviewed sampled-state CT-divB threshold"
            )
            if ct_result == "pass"
            else (
                "an authenticated selected-lineage native restart state exceeds the "
                "reviewed CT-divB threshold"
            )
            if ct_result == "fail"
            else (
                "available authenticated native restart states provide partial sampled "
                "CT-divB evidence, but exact t=9,10 coverage is incomplete"
            )
            if native.get("ct_evidence_available")
            else (
                CT_ABSENCE_REASON
                if audits else
                "no complete retained snapshot or exact required native restart state "
                "was available for authenticated inspection"
            )
        )
        base.update({
            "audit_status": (
                "authentication_failed"
                if native_status == "authentication_failed"
                else "authenticated_complete"
                if ct_result in ("pass", "fail")
                else "authenticated_inconclusive"
                if audits or native.get("ct_evidence_available")
                else "partial_inconclusive"
            ),
            "provenance_authenticated": native_status != "authentication_failed",
            "snapshot_content_authenticated": bool(audits),
            "lineage_binding": lineage_binding,
            "verified_provenance_artifacts": verified,
            "snapshot_summary": summary,
            "snapshot_audits": audits,
            "native_restart_ct": native,
            "ct_result": ct_result,
            "ct_evidence_available": bool(native.get("ct_evidence_available")),
            "ct_claim_supported": bool(native.get("ct_claim_supported")),
            "errors": native_errors,
            "reason": reason,
        })
    except (FastCtAuditError, OSError) as error:
        base.update({
            "audit_status": "authentication_failed",
            "provenance_authenticated": False,
            "snapshot_content_authenticated": False,
            "errors": [str(error)],
            "reason": "assembled-case provenance or retained snapshot authentication failed",
            "native_restart_ct": empty_native_restart_record(
                required_times,
                "assembled-case provenance or retained snapshot authentication failed",
                threshold,
            ),
            "snapshot_summary": {
                "snapshot_count": 0,
                "complete_snapshot_count": 0,
                "incomplete_snapshot_count": 0,
                "selected_snapshot_count": 0,
                "authenticated_snapshot_count": 0,
                "authenticated_rank_file_count": 0,
                "authenticated_bytes": 0,
                "time_first": None,
                "time_last": None,
            },
        })
    return base


def build_audit(
    inventory_path: Path,
    expected_inventory_sha256: str,
    requested_cases: Iterable[str],
    snapshot_policy: str,
) -> dict[str, object]:
    """Build one deterministic accelerated-lineage CT audit."""

    expected_sha = require_sha256(expected_inventory_sha256, "expected inventory SHA-256")
    acceptance, policy, required_times, threshold, native_sources = (
        load_native_ct_authority()
    )
    inventory, inventory_binding = read_bound_json(
        inventory_path, "direct-fast assembled inventory", expected_sha256=expected_sha
    )
    if inventory.get("schema_version") != 1:
        raise FastCtAuditError("direct-fast assembled inventory schema differs")
    output_root = Path(require_text(inventory.get("output"), "assembled output")).absolute()
    if inventory_path.expanduser().absolute().parent != output_root:
        raise FastCtAuditError("assembled inventory path differs from declared output root")
    matrix_binding = verify_declared_binding(inventory.get("matrix"), "assembled matrix")
    reporter_binding = verify_declared_binding(
        inventory.get("adapter"), "direct-fast reporter adapter"
    )
    cases = require_dict(inventory.get("cases"), "assembled cases")
    selected = expand_cases(requested_cases, cases)
    records: dict[str, object] = {}
    for case_id in selected:
        records[case_id] = (
            audit_case(
                acceptance,
                policy,
                case_id,
                cases[case_id],
                output_root,
                str(matrix_binding["sha256"]),
                snapshot_policy,
                required_times,
                threshold,
            )
            if case_id in cases
            else missing_case_record(case_id, required_times, threshold)
        )
    authenticated = sum(
        bool(require_dict(value, "case audit").get("provenance_authenticated"))
        for value in records.values()
    )
    failed = sum(
        require_dict(value, "case audit").get("audit_status") == "authentication_failed"
        for value in records.values()
    )
    snapshots = sum(
        int(require_dict(require_dict(value, "case audit")["snapshot_summary"], "summary")
            ["authenticated_snapshot_count"])
        for value in records.values()
    )
    evidence_cases = sum(
        bool(require_dict(value, "case audit").get("ct_evidence_available"))
        for value in records.values()
    )
    pass_cases = sum(
        require_dict(value, "case audit").get("ct_result") == "pass"
        for value in records.values()
    )
    fail_cases = sum(
        require_dict(value, "case audit").get("ct_result") == "fail"
        for value in records.values()
    )
    complete_cases = sum(
        bool(
            require_dict(
                require_dict(value, "case audit").get("native_restart_ct"),
                "native restart CT",
            ).get("coverage_complete")
        )
        for value in records.values()
    )
    overall = (
        "authentication_failed"
        if failed else "fail"
        if fail_cases else "pass"
        if records and pass_cases == len(records) else "inconclusive"
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "result": overall,
        "inventory": inventory_binding,
        "source_bindings": {
            "matrix": matrix_binding,
            "direct_fast_reporter": reporter_binding,
            "native_restart_ct": native_sources,
        },
        "selection": {
            "cases": selected,
            "snapshot_policy": snapshot_policy,
            "required_native_restart_state_times": required_times,
        },
        "claim_boundary": {
            "snapshot_ct_claim_supported": False,
            "authenticated_selected_lineage_native_restart_ct_claim_supported": True,
            "campaign_authority_eligible": False,
            "release_authorizing": False,
            "cell_centered_divergence_is_ct_evidence": False,
            "required_for_sampled_ct_divb": (
                "exact complete rank-local native restart groups at t=9 and t=10 "
                "from authenticated selected direct-fast roots or the exact selected "
                "accepted-R02 bundle and accepted segment declarations"
            ),
            "snapshot_reason": CT_ABSENCE_REASON,
            "authority_limitation": (
                "Sampled numerical CT evidence from this adapter is non-authorizing; "
                "any reviewed F118/F116 authority-chain blocker remains explicit."
            ),
        },
        "summary": {
            "requested_case_count": len(selected),
            "present_case_count": sum(case_id in cases for case_id in selected),
            "provenance_authenticated_case_count": authenticated,
            "authentication_failed_case_count": failed,
            "authenticated_snapshot_count": snapshots,
            "ct_evidence_case_count": evidence_cases,
            "ct_complete_coverage_case_count": complete_cases,
            "ct_pass_case_count": pass_cases,
            "ct_fail_case_count": fail_cases,
        },
        "cases": records,
    }
    return seal_evidence(result)


CSV_FIELDS = (
    "case_id",
    "assembled_status",
    "audit_status",
    "provenance_authenticated",
    "snapshot_content_authenticated",
    "snapshot_count",
    "complete_snapshot_count",
    "selected_snapshot_count",
    "authenticated_snapshot_count",
    "authenticated_rank_file_count",
    "time_last",
    "native_restart_ct_status",
    "native_restart_ct_coverage_complete",
    "audited_native_restart_state_times",
    "missing_required_native_restart_state_times",
    "maximum_normalized_ct_divb",
    "normalized_ct_divb_lt",
    "campaign_authority_eligible",
    "authority_chain_status",
    "authority_blockers",
    "ct_result",
    "ct_evidence_available",
    "ct_claim_supported",
    "reason",
    "errors",
)


def render_csv(audit: dict[str, object]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS, lineterminator="\n")
    writer.writeheader()
    cases = require_dict(audit.get("cases"), "audit cases")
    for case_id in sorted(cases):
        case = require_dict(cases[case_id], f"{case_id} audit")
        summary = require_dict(case.get("snapshot_summary"), f"{case_id} snapshot summary")
        native = require_dict(case.get("native_restart_ct"), f"{case_id} native restart CT")
        writer.writerow({
            "case_id": case_id,
            "assembled_status": case.get("assembled_status"),
            "audit_status": case.get("audit_status"),
            "provenance_authenticated": str(
                bool(case.get("provenance_authenticated"))
            ).lower(),
            "snapshot_content_authenticated": str(
                bool(case.get("snapshot_content_authenticated"))
            ).lower(),
            "snapshot_count": summary.get("snapshot_count"),
            "complete_snapshot_count": summary.get("complete_snapshot_count"),
            "selected_snapshot_count": summary.get("selected_snapshot_count"),
            "authenticated_snapshot_count": summary.get("authenticated_snapshot_count"),
            "authenticated_rank_file_count": summary.get("authenticated_rank_file_count"),
            "time_last": summary.get("time_last"),
            "native_restart_ct_status": native.get("status"),
            "native_restart_ct_coverage_complete": str(
                bool(native.get("coverage_complete"))
            ).lower(),
            "audited_native_restart_state_times": " ".join(
                format(float(value), ".17g")
                for value in native.get("audited_state_times", [])
            ),
            "missing_required_native_restart_state_times": " ".join(
                format(float(value), ".17g")
                for value in native.get("missing_required_state_times", [])
            ),
            "maximum_normalized_ct_divb": native.get("maximum_normalized_ct_divb"),
            "normalized_ct_divb_lt": native.get("normalized_ct_divb_lt"),
            "campaign_authority_eligible": str(
                bool(native.get("campaign_authority_eligible"))
            ).lower(),
            "authority_chain_status": require_dict(
                native.get("authority_chain"), f"{case_id} authority chain"
            ).get("status"),
            "authority_blockers": " | ".join(
                str(value) for value in native.get("authority_blockers", [])
            ),
            "ct_result": case.get("ct_result"),
            "ct_evidence_available": str(bool(case.get("ct_evidence_available"))).lower(),
            "ct_claim_supported": str(bool(case.get("ct_claim_supported"))).lower(),
            "reason": case.get("reason"),
            "errors": " | ".join(str(value) for value in case.get("errors", [])),
        })
    return stream.getvalue()


def render_markdown(audit: dict[str, object]) -> str:
    summary = require_dict(audit.get("summary"), "audit summary")
    lines = [
        "# Accelerated Stage I CT Snapshot Audit",
        "",
        f"Overall result: **{audit['result']}**.",
        "",
        "This audit authenticates assembled-case provenance and retained snapshot bytes. "
        "It does not treat cell-centered magnetic-field output as discrete constrained-"
        "transport divergence evidence. Exact complete native restart groups beneath "
        "authenticated selected direct-fast roots or declared by the exact selected "
        "accepted-R02 bundle are evaluated with the reviewed face-centered restart "
        "parser when available.",
        "",
        "## Claim Boundary",
        "",
        f"- {CT_ABSENCE_REASON}",
        f"- {DIRECT_FAST_AUTHORITY_LIMITATION}",
        f"- {ACCEPTED_R02_AUTHORITY_LIMITATION}",
        "- A CT pass is never inferred from snapshot availability, finite cell-centered "
        "magnetic values, or incomplete/non-exact-time restart groups.",
        "",
        "## Summary",
        "",
        f"- Requested cases: {summary['requested_case_count']}",
        f"- Present cases: {summary['present_case_count']}",
        f"- Provenance-authenticated cases: {summary['provenance_authenticated_case_count']}",
        f"- Authentication failures: {summary['authentication_failed_case_count']}",
        f"- Authenticated snapshots: {summary['authenticated_snapshot_count']}",
        f"- Cases with direct CT evidence: {summary['ct_evidence_case_count']}",
        f"- Cases with complete exact t=9,10 restart coverage: "
        f"{summary['ct_complete_coverage_case_count']}",
        f"- Sampled numerical CT passes: {summary['ct_pass_case_count']}",
        f"- Sampled numerical CT failures: {summary['ct_fail_case_count']}",
        "",
        "## Cases",
        "",
        (
            "| Case | Assembled | Audit | Snapshots | Native states | "
            "Max normalized CT-divB | CT result |"
        ),
        "|---|---|---|---:|---|---:|---|",
    ]
    cases = require_dict(audit.get("cases"), "audit cases")
    for case_id in sorted(cases):
        case = require_dict(cases[case_id], f"{case_id} audit")
        snapshot = require_dict(case.get("snapshot_summary"), f"{case_id} summary")
        native = require_dict(case.get("native_restart_ct"), f"{case_id} native CT")
        audited = ",".join(
            format(float(value), ".17g")
            for value in native.get("audited_state_times", [])
        )
        maximum = native.get("maximum_normalized_ct_divb")
        lines.append(
            f"| {case_id} | {case.get('assembled_status')} | {case.get('audit_status')} | "
            f"{snapshot.get('authenticated_snapshot_count')} | "
            f"{audited} | "
            f"{'' if maximum is None else format(float(maximum), '.17g')} | "
            f"{case.get('ct_result')} |"
        )
        for error in case.get("errors", []):
            lines.append(f"\n- **{case_id} authentication error:** {error}")
        for blocker in native.get("blockers", []):
            lines.append(f"\n- **{case_id} CT blocker:** {blocker}")
        for blocker in native.get("authority_blockers", []):
            lines.append(f"\n- **{case_id} authority-chain blocker:** {blocker}")
    return "\n".join(lines) + "\n"


def write_outputs(output: Path, audit: dict[str, object]) -> None:
    atomic_write(
        output / "ct_audit.json",
        json.dumps(audit, indent=2, sort_keys=True, allow_nan=False).encode("utf-8")
        + b"\n",
    )
    atomic_write(output / "ct_audit.csv", render_csv(audit).encode("utf-8"))
    atomic_write(output / "ct_audit.md", render_markdown(audit).encode("utf-8"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--expected-inventory-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--cases",
        action="append",
        default=[],
        help="case ID, comma list, or inclusive range; default is every present case",
    )
    parser.add_argument(
        "--snapshot-policy",
        choices=("all", "latest"),
        default="all",
        help="authenticate and inspect every complete snapshot or only the latest per case",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        inventory_root = args.inventory.expanduser().absolute().parent
        if safe_relative_to(args.output, inventory_root):
            raise FastCtAuditError("output beneath the assembled input root is forbidden")
        audit = build_audit(
            args.inventory,
            args.expected_inventory_sha256,
            args.cases,
            args.snapshot_policy,
        )
        write_outputs(args.output, audit)
    except (FastCtAuditError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    print(json.dumps({
        "result": audit["result"],
        "output": str(args.output.expanduser().absolute()),
        "cases": audit["selection"]["cases"],
        "authenticated_snapshots": audit["summary"]["authenticated_snapshot_count"],
    }, sort_keys=True))
    return 1 if audit["result"] == "authentication_failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
