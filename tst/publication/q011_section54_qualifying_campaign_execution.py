#!/usr/bin/env python3
"""Materialize one immutable Q-011 Section 5.4 qualifying campaign plan.

This source-local tool creates a review artifact only.  It validates a human
pressure-selection receipt, binds one frozen clean candidate, and writes
deterministic launch-prohibited handoff contracts.  It never mutates live
policy, invokes a scheduler, submits work, infers pressure, or authorizes a
launch.
"""

from __future__ import annotations

import argparse
import ctypes
import errno
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
import sys
import types
from typing import Any, Mapping, Sequence
import uuid

if __package__:
    from . import immutable_orion_tree
    from . import q011_section54_model as model
    from . import q011_section54_pressure_pilot_execution as pressure_pilot_execution
    from . import q011_section54_pressure_selection as pressure_selection
    from . import q011_section54_restart as restart
else:
    import immutable_orion_tree
    import q011_section54_model as model
    import q011_section54_pressure_pilot_execution as pressure_pilot_execution
    import q011_section54_pressure_selection as pressure_selection
    import q011_section54_restart as restart


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"
MATERIALIZER_SOURCE = Path(__file__).resolve()
AUTHORIZED_ORION_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_SOURCE_ROOT = Path("/ccs/home/dfielding/athenak-pic")
REVIEWED_ENVIRONMENT_PROFILE_SOURCE = pressure_pilot_execution.ENVIRONMENT_PROFILE_SOURCE
CONTROL_PLANE_COMMON_SOURCE = (
    REPO_ROOT / "tst/publication/frontier_control_plane/control_plane_common.py"
)
CONTROL_PLANE_OPERATOR_ATTESTATION_SOURCE = (
    REPO_ROOT / "tst/publication/frontier_control_plane/operator_attestation.py"
)
QUALIFYING_PREREGISTRATION = (
    READINESS_ROOT
    / "q011_section54_qualifying_campaign_preregistration_successor_v2_2026-06-01.json"
)
RESTART_PREREGISTRATION = (
    READINESS_ROOT / "q011_section54_restart_continuation_preregistration_2026-06-01.json"
)
PAPER_DECK = (
    REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
QUALIFYING_PREREGISTRATION_SHA256 = (
    "6fd9ebbc247b6cace69f0ff61553cf198241577b457410d57d26afcbf27cdc35"
)
RESTART_PREREGISTRATION_SHA256 = (
    "c3360694dc90d391c5ccf7a0620ae576733e87beea3fa974c69602c82dd866ab"
)
PAPER_DECK_SHA256 = "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
CAMPAIGN_ID = "Q011-SECTION54-QUALIFYING-CAMPAIGN"
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
ARTIFACT_ROLE = "q011_section54_source_local_immutable_qualifying_campaign_plan"
QUALIFICATION_EFFECT = "plan_only_no_execution_authorization_no_claim_closure"
PLAN_RECORD_TYPE = "q011_section54_qualifying_campaign_execution_plan"
POLICY_FRAGMENT_RECORD_TYPE = "q011_section54_qualifying_campaign_nonauthorizing_policy_fragment"
RECOMPUTE_RECORD_TYPE = "q011_section54_independent_raw_artifact_recompute_plan"
RESTART_CARRIER_RECORD_TYPE = "q011_section54_amr_restart_continuation_carrier"
ATTEMPT_RECORD_TYPE = "q011_section54_baseline_attempt_descriptor"
LAUNCH_CONTRACT_RECORD_TYPE = "q011_section54_launch_prohibited_handoff_contract"
PLANNER_RETENTION_ROLE = "q011_section54_deterministic_retained_attempt"
MATERIALIZATION_RECEIPT_RECORD_TYPE = (
    "q011_section54_qualifying_campaign_plan_materialization_receipt"
)
MATERIALIZATION_RECEIPT_NAME = "materialization_receipt.json"
EXPECTED_BASELINE_ATTEMPTS = 24
_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_SAFE_SEGMENT = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_RENAME_NOREPLACE = 1
_STAGING_ROOT_NAME = "publishable"
_MATERIALIZED_MEMBER_INVENTORY_ALGORITHM = (
    "sha256 of '<file_sha256>  <root-relative-path>\\n' entries ordered "
    "lexically by root-relative path"
)
_PLAN_FREEZE_RECEIPT = {
    "schema_version": 1,
    "artifact_role": ARTIFACT_ROLE,
    "qualification_effect": QUALIFICATION_EFFECT,
    "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}


CANONICAL_VARIANT_IDS = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
QUALIFYING_SEEDS = (
    23050101,
    23050102,
    23050103,
    23050104,
    23050105,
    23050106,
    23050107,
    23050108,
)
_FIXED_HELPER_SOURCES = (
    "tst/publication/q011_section54_model.py",
    "tst/publication/q011_section54_pressure_pilot_execution.py",
    "tst/publication/q011_section54_pressure_selection.py",
    "tst/publication/q011_section54_restart.py",
    "tst/publication/analyze_q011_section54_outputs.py",
    "tst/publication/analyze_q011_section54_campaign.py",
    "tst/publication/analyze_q011_section54_numerical_qualification.py",
    "tst/publication/q011_section54_particles.py",
    "tst/publication/q011_section54_spatial.py",
    "tst/publication/q011_section54_artifacts.py",
    "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
    "tst/publication/analyze_q011_section54_pressure_pilot.py",
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    "tst/publication/frontier_f1_structured_artifacts.py",
    "tst/publication/q011_section54_attempt_manifest_materializer.py",
    "tst/publication/publish_q011_section54_campaign_attempt.py",
    "tst/publication/immutable_orion_tree.py",
    "tst/publication/pvtk_particles.py",
    "tst/publication/q011_parallel_shock_storage_estimator.py",
    "tst/publication/frontier_control_plane/control_plane_common.py",
    "tst/publication/frontier_control_plane/ledger.py",
    "tst/publication/frontier_control_plane/operator_attestation.py",
    "tst/publication/q011_section54_qualifying_campaign_execution.py",
)
class CampaignPlanError(ValueError):
    """Raised when campaign-plan materialization fails closed."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise CampaignPlanError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise CampaignPlanError("campaign plan contains noncanonical JSON values") from error


def _digest_value(value: object) -> str:
    try:
        payload = json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise CampaignPlanError("campaign plan digest input is not canonical JSON") from error
    return _sha256_bytes(payload)


def _decode_json(payload: bytes, *, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise CampaignPlanError(f"{label} contains forbidden JSON constant {value}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise CampaignPlanError(f"{label} is not valid UTF-8 JSON") from error


def _canonical_existing_directory(path: Path, *, label: str) -> Path:
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute(), f"{label} must be absolute")
    try:
        resolved = lexical.resolve(strict=True)
        mode = os.lstat(lexical).st_mode
    except OSError as error:
        raise CampaignPlanError(f"{label} is unavailable: {lexical}") from error
    _require(resolved == lexical, f"{label} must not use a symlink or path alias")
    _require(stat.S_ISDIR(mode), f"{label} must be a real directory")
    return lexical


def _stable_regular_bytes(
    path: Path,
    *,
    label: str,
    require_read_only: bool = False,
    require_executable: bool = False,
) -> tuple[Path, bytes]:
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute(), f"{label} must be absolute")
    try:
        resolved = lexical.resolve(strict=True)
    except OSError as error:
        raise CampaignPlanError(f"{label} is unavailable: {lexical}") from error
    _require(resolved == lexical, f"{label} must not use a symlink or path alias")
    try:
        descriptor = os.open(lexical, _FILE_FLAGS)
    except OSError as error:
        raise CampaignPlanError(f"{label} is not an openable regular file: {lexical}") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} must be a regular file")
        _require(before.st_nlink == 1, f"{label} must not be hard linked")
        if require_read_only:
            _require(not before.st_mode & _WRITE_BITS, f"{label} must be read-only")
        if require_executable:
            _require(bool(before.st_mode & 0o111), f"{label} must be executable")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(lexical, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        _require(identity(before) == identity(after), f"{label} changed while reading")
        _require(
            (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino),
            f"{label} path changed while reading",
        )
        return lexical, bytes(payload)
    finally:
        os.close(descriptor)


def _stable_regular_bytes_below(
    path: Path,
    root: Path,
    *,
    label: str,
    require_read_only: bool = False,
    require_executable: bool = False,
) -> tuple[Path, bytes, int]:
    """Read one regular file through pinned canonical ancestor descriptors."""
    lexical_root = _canonical_existing_directory(root, label=f"{label} root")
    lexical = Path(os.path.abspath(path))
    _require(path.is_absolute() and lexical == path, f"{label} must use a canonical absolute path")
    try:
        relative = lexical.relative_to(lexical_root)
    except ValueError as error:
        raise CampaignPlanError(f"{label} escaped its authorized root") from error
    _require(
        relative.parts
        and all(part not in {"", ".", ".."} for part in relative.parts),
        f"{label} path is unsafe",
    )
    root_descriptor = os.open(lexical_root, _DIRECTORY_FLAGS)
    directory_descriptor = os.dup(root_descriptor)
    descriptor = -1
    try:
        for part in relative.parts[:-1]:
            child = os.open(part, _DIRECTORY_FLAGS, dir_fd=directory_descriptor)
            os.close(directory_descriptor)
            directory_descriptor = child
        descriptor = os.open(relative.parts[-1], _FILE_FLAGS, dir_fd=directory_descriptor)
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} must be a regular file")
        _require(before.st_nlink == 1, f"{label} must not be hard linked")
        if require_read_only:
            _require(not before.st_mode & _WRITE_BITS, f"{label} must be read-only")
        if require_executable:
            _require(
                stat.S_IMODE(before.st_mode) & 0o111 == 0o111,
                f"{label} must retain all execute bits",
            )
        payload = bytearray()
        offset = 0
        while chunk := os.pread(descriptor, 1024 * 1024, offset):
            payload.extend(chunk)
            offset += len(chunk)
        after = os.fstat(descriptor)
        current = os.stat(relative.parts[-1], dir_fd=directory_descriptor, follow_symlinks=False)
        opened_root = os.fstat(root_descriptor)
        current_root = os.stat(lexical_root, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        _require(identity(before) == identity(after), f"{label} changed while reading")
        _require(
            (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino),
            f"{label} path changed while reading",
        )
        _require(
            (opened_root.st_dev, opened_root.st_ino)
            == (current_root.st_dev, current_root.st_ino),
            f"{label} root changed while reading",
        )
        return lexical, bytes(payload), before.st_mode
    except OSError as error:
        raise CampaignPlanError(f"{label} is unavailable below its authorized root") from error
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(directory_descriptor)
        os.close(root_descriptor)


def _load_control_plane_common() -> types.ModuleType:
    """Load reviewed shared validation bytes without path-based import races."""
    _, operator_payload = _stable_regular_bytes(
        CONTROL_PLANE_OPERATOR_ATTESTATION_SOURCE,
        label="control-plane operator-attestation helper",
    )
    _, common_payload = _stable_regular_bytes(
        CONTROL_PLANE_COMMON_SOURCE,
        label="control-plane shared validator",
    )
    operator_module = types.ModuleType("operator_attestation")
    operator_module.__file__ = str(CONTROL_PLANE_OPERATOR_ATTESTATION_SOURCE)
    exec(
        compile(
            operator_payload,
            str(CONTROL_PLANE_OPERATOR_ATTESTATION_SOURCE),
            "exec",
            dont_inherit=True,
        ),
        operator_module.__dict__,
    )
    previous = sys.modules.get("operator_attestation")
    sys.modules["operator_attestation"] = operator_module
    try:
        common_module = types.ModuleType("_q011_section54_control_plane_common")
        common_module.__file__ = str(CONTROL_PLANE_COMMON_SOURCE)
        exec(
            compile(
                common_payload,
                str(CONTROL_PLANE_COMMON_SOURCE),
                "exec",
                dont_inherit=True,
            ),
            common_module.__dict__,
        )
        return common_module
    finally:
        if previous is None:
            del sys.modules["operator_attestation"]
        else:
            sys.modules["operator_attestation"] = previous


def _installed_environment_profile(
    path: Path, *, authorized_root: Path
) -> tuple[Path, bytes, str, bytes]:
    lexical = Path(os.path.abspath(path))
    _require(
        path.is_absolute() and lexical == path,
        "installed environment profile must use a canonical absolute path",
    )
    try:
        relative = lexical.relative_to(authorized_root)
    except ValueError as error:
        raise CampaignPlanError("installed environment profile escaped the authorized Orion root") from error
    _require(
        len(relative.parts) == 3
        and relative.parts[0] == "control_plane"
        and _SHA256.fullmatch(relative.parts[1]) is not None
        and relative.parts[2] == "frontier_pic_environment.sh",
        "installed environment profile must use the canonical control_plane/<version> path",
    )
    environment_path, environment_payload, _ = _stable_regular_bytes_below(
        lexical,
        authorized_root,
        label="installed environment profile",
        require_read_only=True,
    )
    _, reviewed_payload = _stable_regular_bytes(
        REVIEWED_ENVIRONMENT_PROFILE_SOURCE,
        label="reviewed pressure-pilot environment profile source",
    )
    _require(
        _sha256_bytes(environment_payload) == _sha256_bytes(reviewed_payload),
        "installed environment profile differs from reviewed pressure-pilot source bytes",
    )
    return environment_path, environment_payload, relative.parts[1], reviewed_payload


def _text(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value.strip()), f"{label} must be nonempty text")
    return value


def _git_commit(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    _require(_GIT_COMMIT.fullmatch(text) is not None, f"{label} must be 40 lowercase hex digits")
    return text


def _exact_int(value: object, *, label: str) -> int:
    _require(type(value) is int, f"{label} must be an exact integer")
    return value


def _exact_float(value: object, *, label: str) -> float:
    _require(type(value) is float and math.isfinite(value), f"{label} must be a finite JSON float")
    return value


def _relative_path(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    path = PurePosixPath(text)
    _require(
        not path.is_absolute()
        and text == path.as_posix()
        and text != "."
        and all(part not in {"", ".", ".."} for part in path.parts)
        and not any(character.isspace() for character in text),
        f"{label} must be a canonical relative path",
    )
    return text


def _uuid(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    try:
        parsed = uuid.UUID(text)
    except ValueError as error:
        raise CampaignPlanError(f"{label} must be a canonical UUID") from error
    _require(str(parsed) == text, f"{label} must be a canonical UUID")
    return text


def _object(value: object, expected: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label} must be an object")
    _require(set(value) == expected, f"{label} schema drifted")
    return value


def _list(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label} must be an array")
    return value


def _binding(path: str, payload: bytes) -> dict[str, str]:
    return {"path": path, "sha256": _sha256_bytes(payload)}


def _repo_binding(relative: str) -> dict[str, str]:
    path = REPO_ROOT / relative
    _, payload = _stable_regular_bytes(path, label=f"helper source {relative}")
    return _binding(relative, payload)


def validate_campaign_matrix(value: object) -> dict[str, object]:
    """Require the exact preregistered three-by-eight baseline matrix."""
    matrix = _object(
        value,
        {
            "physical_mode",
            "grid_variants",
            "qualifying_seeds",
            "expected_baseline_attempts",
            "paired_seed_rule",
        },
        label="campaign matrix",
    )
    variants = [
        _text(item, label=f"campaign matrix/grid_variants[{index}]")
        for index, item in enumerate(_list(matrix["grid_variants"], label="campaign matrix/grid_variants"))
    ]
    _require(len(variants) == len(set(variants)), "campaign matrix contains duplicate variants")
    seeds = [
        _exact_int(item, label=f"campaign matrix/qualifying_seeds[{index}]")
        for index, item in enumerate(_list(matrix["qualifying_seeds"], label="campaign matrix/qualifying_seeds"))
    ]
    _require(len(seeds) == len(set(seeds)), "campaign matrix contains duplicate seeds")
    expected_count = _exact_int(
        matrix["expected_baseline_attempts"],
        label="campaign matrix/expected_baseline_attempts",
    )
    _require(
        matrix["physical_mode"] == PHYSICAL_MODE
        and tuple(variants) == CANONICAL_VARIANT_IDS
        and tuple(seeds) == QUALIFYING_SEEDS
        and expected_count == EXPECTED_BASELINE_ATTEMPTS
        and expected_count == len(variants) * len(seeds),
        "campaign matrix drifted from the fixed three-by-eight preregistration",
    )
    _require(
        matrix["paired_seed_rule"]
        == "Use the same qualifying seed for coarse-uniform, AMR and fine-uniform variants.",
        "campaign matrix paired-seed rule drifted",
    )
    return {
        "physical_mode": PHYSICAL_MODE,
        "grid_variants": variants,
        "qualifying_seeds": seeds,
        "expected_baseline_attempts": expected_count,
        "paired_seed_rule": matrix["paired_seed_rule"],
    }


def _load_qualifying_preregistration(path: Path) -> tuple[bytes, dict[str, Any], dict[str, object]]:
    _, payload = _stable_regular_bytes(path, label="qualifying campaign preregistration")
    _require(
        _sha256_bytes(payload) == QUALIFYING_PREREGISTRATION_SHA256,
        "qualifying campaign preregistration SHA-256 drifted",
    )
    preregistration = _decode_json(payload, label="qualifying campaign preregistration")
    _require(type(preregistration) is dict, "qualifying campaign preregistration must be an object")
    _require(
        preregistration.get("record_type") == "q011_section54_qualifying_campaign_preregistration"
        and type(preregistration.get("schema_version")) is int
        and preregistration["schema_version"] == 1,
        "qualifying campaign preregistration identity drifted",
    )
    criteria = preregistration.get("athenak_selected_release_criteria")
    _require(type(criteria) is dict, "qualifying campaign release criteria drifted")
    matrix = validate_campaign_matrix(criteria.get("campaign_matrix"))
    execution = preregistration.get("qualifying_execution_bindings")
    _require(type(execution) is dict, "qualifying execution boundary drifted")
    _require(
        execution.get("frontier_execution_authorized_by_this_record") is False,
        "qualifying preregistration unexpectedly authorizes Frontier execution",
    )
    return payload, preregistration, matrix


def _load_restart_preregistration(path: Path) -> tuple[bytes, dict[str, object]]:
    _, payload = _stable_regular_bytes(path, label="restart-continuation preregistration")
    _require(
        _sha256_bytes(payload) == RESTART_PREREGISTRATION_SHA256,
        "restart-continuation preregistration SHA-256 drifted",
    )
    try:
        preregistration = restart.decode_preregistration(payload.decode("utf-8"))
        restart.validate_preregistration(preregistration)
    except (UnicodeDecodeError, restart.RestartPolicyError) as error:
        raise CampaignPlanError("restart-continuation preregistration drifted") from error
    return payload, preregistration


def _load_pressure_receipt(
    path: Path, *, authorized_pic_root: Path
) -> tuple[bytes, dict[str, object]]:
    _, payload = _stable_regular_bytes(
        path, label="human pressure-selection receipt", require_read_only=True
    )
    try:
        receipt = pressure_selection.validate_pressure_selection_receipt_bytes(
            payload,
            authorized_pic_root=authorized_pic_root,
        )
    except pressure_selection.PressureSelectionReceiptError as error:
        raise CampaignPlanError(f"human pressure-selection receipt is invalid: {error}") from error
    return payload, receipt


def _candidate_member_path(candidate_root: Path, value: object, suffix: str, *, label: str) -> Path:
    path = Path(_text(value, label=label))
    expected = candidate_root / suffix
    _require(path == expected, f"{label} must use the frozen clean-candidate path {expected}")
    return path


def _candidate_binding(
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    paper_deck_sha256: str,
    admission_analyzer_sha256: str,
    output_primitives_sha256: str,
) -> tuple[bytes, bytes, dict[str, object]]:
    authorized_root = _canonical_existing_directory(
        AUTHORIZED_ORION_ROOT, label="authorized Orion PIC root"
    )
    common = _load_control_plane_common()
    manifest_path, manifest_payload, _ = _stable_regular_bytes_below(
        clean_candidate_manifest,
        authorized_root,
        label="clean-candidate manifest",
        require_read_only=True,
    )
    manifest = _decode_json(manifest_payload, label="clean-candidate manifest")
    _require(type(manifest) is dict, "clean-candidate manifest must be an object")
    _require(
        type(manifest.get("schema_version")) is int and manifest["schema_version"] == 4,
        "clean-candidate manifest schema version drifted",
    )
    freeze_id = _uuid(manifest.get("freeze_id"), label="clean candidate/freeze_id")
    candidate_root = authorized_root / "clean_candidates" / freeze_id
    _require(
        manifest_path == candidate_root / "clean_candidate_manifest.json",
        "clean-candidate manifest path escaped its authorized Orion freeze root",
    )
    _canonical_existing_directory(candidate_root, label="clean-candidate freeze root")
    source = manifest.get("source")
    build = manifest.get("build")
    _require(type(source) is dict, "clean-candidate source attestation must be an object")
    _require(type(build) is dict, "clean-candidate build attestation must be an object")
    git_commit = _git_commit(source.get("git_commit"), label="clean candidate/source/git_commit")
    archive = _candidate_member_path(candidate_root, source["archive_path"], "source.tar", label="clean candidate/source/archive_path")
    commit = _candidate_member_path(candidate_root, source["commit_path"], "source.commit", label="clean candidate/source/commit_path")
    _, archive_payload, _ = _stable_regular_bytes_below(archive, candidate_root, label="clean-candidate source archive", require_read_only=True)
    _, commit_payload, _ = _stable_regular_bytes_below(commit, candidate_root, label="clean-candidate source commit", require_read_only=True)
    raw_submodules = source.get("submodules")
    _require(type(raw_submodules) is list, "clean-candidate submodules must be an array")
    submodule_archives = []
    submodule_commits = []
    for index, raw in enumerate(raw_submodules):
        _require(type(raw) is dict, f"clean-candidate submodule {index} must be an object")
        relative = _relative_path(raw.get("path"), label=f"clean candidate/source/submodules[{index}]/path")
        submodule_archive = _candidate_member_path(
            candidate_root,
            raw.get("archive_path"),
            f"submodules/{index:04d}.tar",
            label=f"clean candidate/source/submodules[{index}]/archive_path",
        )
        submodule_commit = _candidate_member_path(
            candidate_root,
            raw.get("commit_path"),
            f"submodules/{index:04d}.commit",
            label=f"clean candidate/source/submodules[{index}]/commit_path",
        )
        _, archive_bytes, _ = _stable_regular_bytes_below(
            submodule_archive,
            candidate_root,
            label=f"clean-candidate submodule archive {relative}",
            require_read_only=True,
        )
        _, commit_bytes, _ = _stable_regular_bytes_below(
            submodule_commit,
            candidate_root,
            label=f"clean-candidate submodule commit {relative}",
            require_read_only=True,
        )
        submodule_archives.append(archive_bytes)
        submodule_commits.append(commit_bytes)
    profile = _candidate_member_path(candidate_root, build["profile_path"], "build_profile.json", label="clean candidate/build/profile_path")
    profile_receipt = _candidate_member_path(candidate_root, build["profile_receipt_path"], "profile_receipt.json", label="clean candidate/build/profile_receipt_path")
    _, profile_payload, _ = _stable_regular_bytes_below(profile, candidate_root, label="clean-candidate build profile", require_read_only=True)
    _, profile_receipt_payload, _ = _stable_regular_bytes_below(profile_receipt, candidate_root, label="clean-candidate profile receipt", require_read_only=True)
    executable_path = _candidate_member_path(candidate_root, build["executable_path"], "athena", label="clean candidate/build/executable_path")
    supplied_executable = Path(os.path.abspath(executable))
    _require(executable.is_absolute() and supplied_executable == executable_path, "supplied executable differs from clean-candidate binding")
    _, executable_payload, executable_mode = _stable_regular_bytes_below(
        executable_path,
        candidate_root,
        label="clean-candidate executable",
        require_read_only=True,
        require_executable=True,
    )
    executable_sha256 = _sha256_bytes(executable_payload)
    provenance_root = candidate_root / "build_provenance"
    _canonical_existing_directory(provenance_root, label="clean-candidate build provenance root")
    build_provenance = {
        provenance_label: _stable_regular_bytes_below(
            provenance_root / filename,
            candidate_root,
            label=f"clean-candidate build provenance {provenance_label}",
            require_read_only=True,
        )[1]
        for provenance_label, filename in common.BUILD_PROVENANCE_FILENAMES.items()
    }
    environment_path, environment_payload, control_plane_version, reviewed_environment_payload = (
        _installed_environment_profile(environment_profile, authorized_root=authorized_root)
    )
    try:
        installed_inventory = common.verify_installed_control_plane(
            environment_path.parent,
            authorized_pic_root=authorized_root,
        )
        _require(
            installed_inventory["version"] == control_plane_version,
            "installed environment profile control-plane inventory drifted",
        )
        validated_submodules = common.validate_clean_candidate_bundle(
            manifest,
            source_archive=archive_payload,
            source_commit=commit_payload,
            submodule_archives=submodule_archives,
            submodule_commits=submodule_commits,
            build_profile=profile_payload,
            build_profile_receipt=profile_receipt_payload,
            build_provenance=build_provenance,
            executable_sha256=executable_sha256,
            expected_control_plane_version=control_plane_version,
            authorized_pic_root=authorized_root,
            authorized_source_root=AUTHORIZED_SOURCE_ROOT,
        )
        sealed_executable = immutable_orion_tree._sealed_memfd(
            executable_payload,
            executable_mode,
            name="q011-section54-clean-candidate-athena",
        )
        try:
            immutable_orion_tree.validate_executable_elf(
                Path(f"/proc/self/fd/{sealed_executable}"),
                executable_sha256,
                error_type=CampaignPlanError,
                label="clean-candidate executable",
            )
        finally:
            os.close(sealed_executable)
    except CampaignPlanError:
        raise
    except (OSError, ValueError) as error:
        raise CampaignPlanError(f"clean-candidate bundle failed shared validation: {error}") from error
    prepared = manifest["prepared_artifacts"]
    measured_prepared = {
        record["path"]: record["sha256"]
        for record in [*prepared["paper_decks"], *prepared["analyzers"]]
    }
    for path, digest in {
        model.BASE_DECK_PATH: paper_deck_sha256,
        "tst/publication/analyze_q011_section54_campaign.py": admission_analyzer_sha256,
        "tst/publication/analyze_q011_section54_outputs.py": output_primitives_sha256,
    }.items():
        _require(
            measured_prepared.get(path) == digest,
            f"clean-candidate prepared-artifact binding drifted: {path}",
        )
    profile_sha256 = _sha256_bytes(profile_payload)
    profile_receipt_sha256 = _sha256_bytes(profile_receipt_payload)
    return manifest_payload, environment_payload, {
        "clean_candidate_manifest": {
            "path": str(manifest_path),
            "sha256": _sha256_bytes(manifest_payload),
        },
        "freeze_id": freeze_id,
        "git_commit": git_commit,
        "git_tree": source["git_tree"],
        "source_archive_sha256": source["archive_sha256"],
        "source_commit_sha256": source["commit_sha256"],
        "source_bundle_sha256": source["source_bundle_sha256"],
        "prepared_artifact_inventory_sha256": prepared["inventory_sha256"],
        "validated_submodules": validated_submodules,
        "build_profile": {"path": str(profile), "sha256": profile_sha256},
        "build_profile_receipt": {
            "path": str(profile_receipt),
            "sha256": profile_receipt_sha256,
        },
        "build_invocations_sha256": build["build_invocations_sha256"],
        "executable": {"path": str(executable_path), "sha256": executable_sha256},
        "environment_profile": {
            "path": str(environment_path),
            "sha256": _sha256_bytes(environment_payload),
            "control_plane_version": control_plane_version,
            "reviewed_source": _binding(
                "tst/publication/frontier_control_plane/frontier_pic_environment.sh",
                reviewed_environment_payload,
            ),
        },
    }


def _attempt_id(index: int, variant: str, seed: int) -> str:
    identifier = f"baseline-{index:03d}-{variant}-seed-{seed}"
    _require(_SAFE_SEGMENT.fullmatch(identifier) is not None, "generated attempt ID is unsafe")
    return identifier


def _restart_carrier_id(seed: int) -> str:
    return f"amr-restart-continuation-seed-{seed}"


def _selected_pressure_token(value: float) -> str:
    _exact_float(value, label="selected problem/ps_p0")
    return repr(value)


def _variant_binding(variant_id: str) -> model.VariantBinding:
    try:
        binding = model.variant_binding(variant_id)
    except model.ModelContractError as error:
        raise CampaignPlanError(f"campaign variant binding drifted: {error}") from error
    _require(
        binding.deck_path == model.BASE_DECK_PATH,
        f"campaign variant deck binding drifted: {variant_id}",
    )
    return binding


def _baseline_launch_contract(
    *,
    attempt_id: str,
    variant: model.VariantBinding,
    seed: int,
    selected_ps_p0: float,
    candidate: Mapping[str, object],
    paper_deck_binding: Mapping[str, str],
    artifact_root: Path,
) -> dict[str, object]:
    seed_overrides = (
        f"particles/pic_random_seed={seed}",
        f"problem/ps_inject_seed={seed}",
        f"problem/ps_seed_noise_seed={seed}",
    )
    return {
        "record_type": LAUNCH_CONTRACT_RECORD_TYPE,
        "schema_version": 1,
        "contract_role": "source_local_review_handoff_only",
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "live_policy_mutation_authorized": False,
        "attempt_id": attempt_id,
        "variant": variant.variant,
        "qualifying_seed": seed,
        "selected_problem_ps_p0": selected_ps_p0,
        "executable": candidate["executable"],
        "environment_profile": candidate["environment_profile"],
        "paper_deck": paper_deck_binding,
        "authorized_orion_attempt_root": str(artifact_root),
        "argv": [
            "-i",
            "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
            "-d",
            str(artifact_root / "raw"),
            f"job/basename={attempt_id}",
            f"problem/ps_p0={_selected_pressure_token(selected_ps_p0)}",
            *seed_overrides,
            *variant.model_launch_overrides,
        ],
        "required_separate_boundary": (
            "review_and_promote_a_registered_frontier_submission_policy_then_use_"
            "the_installed_control_plane_wrapper"
        ),
    }


def _restart_launch_contract(
    *,
    carrier_id: str,
    source_attempt: Mapping[str, object],
    restart_preregistration: Mapping[str, object],
    candidate: Mapping[str, object],
    paper_deck_binding: Mapping[str, str],
    artifact_root: Path,
) -> dict[str, object]:
    continuation = restart_preregistration["continuation_contract"]
    checkpoint = continuation["checkpoint_time_omega0_inverse"]
    return {
        "record_type": LAUNCH_CONTRACT_RECORD_TYPE,
        "schema_version": 1,
        "contract_role": "source_local_restart_review_handoff_only",
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "live_policy_mutation_authorized": False,
        "carrier_id": carrier_id,
        "source_baseline_attempt_id": source_attempt["attempt_id"],
        "variant": source_attempt["variant"],
        "qualifying_seed": source_attempt["qualifying_seed"],
        "selected_problem_ps_p0": source_attempt["selected_problem_ps_p0"],
        "executable": candidate["executable"],
        "environment_profile": candidate["environment_profile"],
        "paper_deck": paper_deck_binding,
        "authorized_orion_attempt_root": str(artifact_root),
        "checkpoint_time_omega0_inverse": checkpoint,
        "checkpoint_input": (
            f"retain_from_{source_attempt['attempt_id']}_at_t{int(checkpoint)}_"
            "then_bind_exact_checksum_before_any_separately_authorized_continuation"
        ),
        "argv_template": [
            "-r",
            "<exact-retained-checkpoint-path-bound-by-separate-reviewed-successor>",
            "-d",
            str(artifact_root / "raw"),
        ],
        "required_separate_boundary": (
            "materialize_a_checkpoint_checksum_bound_registered_restart_successor_"
            "then_review_and_promote_policy_before_using_the_installed_control_plane_wrapper"
        ),
    }


def _attempt_descriptor(
    *,
    index: int,
    attempt_id: str,
    variant: model.VariantBinding,
    seed: int,
    selected_ps_p0: float,
    candidate: Mapping[str, object],
    artifact_root: Path,
    contract_path: str,
    contract_payload: bytes,
) -> dict[str, object]:
    return {
        "record_type": ATTEMPT_RECORD_TYPE,
        "schema_version": 1,
        "attempt_index": index,
        "attempt_id": attempt_id,
        "status": "planned_not_authorized",
        "variant": variant.variant,
        "qualifying_seed": seed,
        "physical_mode": PHYSICAL_MODE,
        "selected_problem_ps_p0": selected_ps_p0,
        "candidate_binding": {
            "git_commit": candidate["git_commit"],
            "source_bundle_sha256": candidate["source_bundle_sha256"],
            "executable_sha256": candidate["executable"]["sha256"],
            "environment_profile_sha256": candidate["environment_profile"]["sha256"],
        },
        "authorized_orion_attempt_root": str(artifact_root),
        "launch_contract": _binding(contract_path, contract_payload),
        "artifact_retention": (
            "retain_every_emitted_raw_artifact_for_every_attempt_including_failed_attempts"
        ),
    }


def _independent_recompute_plan(
    *, plan_id: str, campaign_root: Path, qualifying_preregistration_binding: Mapping[str, str]
) -> dict[str, object]:
    return {
        "record_type": RECOMPUTE_RECORD_TYPE,
        "schema_version": 1,
        "plan_id": plan_id,
        "status": "plan_frozen_independent_implementation_and_review_artifacts_open",
        "qualifying_preregistration": dict(qualifying_preregistration_binding),
        "authorized_orion_campaign_root": str(campaign_root),
        "implementation_rule": (
            "Recompute every primary metric from archived raw artifacts using a "
            "reviewer-owned script or an independently implemented analyzer. The "
            "recompute script must not import production metric-extraction functions "
            "or the local analyzer and helper sources bound by this plan."
        ),
        "raw_input_policy": (
            "consume_root_relative_sha256_inventories_and_archived_raw_bin_pvtk_rst_"
            "bytes_for_every_attempt_including_failed_attempts"
        ),
        "required_records": [
            "independent_script_path_or_archive_locator",
            "independent_script_sha256",
            "independent_environment_lock",
            "attempt_inventory_sha256_values",
            "input_artifact_checksums",
            "metric_comparison_table",
            "reviewer_identity",
            "reviewer_disposition",
        ],
        "required_metric_table_columns": [
            "attempt_id",
            "variant",
            "qualifying_seed",
            "observable",
            "production_metric",
            "independent_metric",
            "absolute_difference",
            "relative_difference",
            "declared_tolerance",
            "disposition",
        ],
        "production_helper_imports_authorized": False,
        "claim_closure_authorized": False,
        "frontier_execution_authorized": False,
    }


def _policy_fragment(
    *,
    plan_id: str,
    campaign_root: Path,
    candidate: Mapping[str, object],
    pressure_receipt_binding: Mapping[str, str],
    contract_bindings: Sequence[Mapping[str, str]],
    restart_contract_binding: Mapping[str, str],
) -> dict[str, object]:
    return {
        "record_type": POLICY_FRAGMENT_RECORD_TYPE,
        "schema_version": 1,
        "plan_id": plan_id,
        "status": "review_fragment_only_not_live_policy",
        "qualification_effect": QUALIFICATION_EFFECT,
        "authorized_orion_root": str(AUTHORIZED_ORION_ROOT),
        "authorized_orion_campaign_root": str(campaign_root),
        "selected_pressure_receipt": dict(pressure_receipt_binding),
        "candidate_binding": dict(candidate),
        "baseline_launch_contracts": [dict(binding) for binding in contract_bindings],
        "restart_continuation_launch_contract": dict(restart_contract_binding),
        "integration_policy": (
            "requires_separate_reviewed_registered_frontier_submission_policy_"
            "successor_and_installed_control_plane_promotion"
        ),
        "mutates_live_policy": False,
        "scheduler_calls_authorized": False,
        "scheduler_submission_authorized": False,
        "frontier_execution_authorized": False,
        "launch_authorized": False,
        "claim_closure_authorized": False,
    }


def _require_same_directory(path: Path, descriptor: int, *, label: str) -> None:
    expected = os.fstat(descriptor)
    try:
        actual = os.stat(path, follow_symlinks=False)
    except OSError as error:
        raise CampaignPlanError(f"{label} changed during publication") from error
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (expected.st_dev, expected.st_ino) == (actual.st_dev, actual.st_ino),
        f"{label} changed during publication",
    )


def _require_absent_at(parent_descriptor: int, name: str, *, label: str) -> None:
    try:
        os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return
    except OSError as error:
        raise CampaignPlanError(f"cannot inspect {label}") from error
    raise CampaignPlanError(f"{label} already exists")


def _require_same_directory_at(
    parent_descriptor: int, name: str, descriptor: int, *, label: str
) -> None:
    _require(
        "/" not in name,
        "descriptor-relative campaign-plan directory check received a nested path",
    )
    try:
        actual = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    except OSError as error:
        raise CampaignPlanError(f"{label} changed during publication") from error
    expected = os.fstat(descriptor)
    _require(
        stat.S_ISDIR(actual.st_mode)
        and (expected.st_dev, expected.st_ino) == (actual.st_dev, actual.st_ino),
        f"{label} changed during publication",
    )


def _rename_no_replace_at(
    source_parent_descriptor: int,
    source_name: str,
    destination_parent_descriptor: int,
    destination_name: str,
) -> None:
    _require(
        "/" not in source_name and "/" not in destination_name,
        "descriptor-relative campaign-plan rename received a nested path",
    )
    libc = ctypes.CDLL(None, use_errno=True)
    renameat2 = getattr(libc, "renameat2", None)
    error_number = errno.ENOSYS
    if renameat2 is not None:
        renameat2.argtypes = [
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        ]
        renameat2.restype = ctypes.c_int
        if renameat2(
            source_parent_descriptor,
            os.fsencode(source_name),
            destination_parent_descriptor,
            os.fsencode(destination_name),
            _RENAME_NOREPLACE,
        ) == 0:
            return
        error_number = ctypes.get_errno()
    if error_number == errno.EEXIST:
        raise CampaignPlanError(
            f"deterministic campaign-plan output root already exists: {destination_name}"
        )
    unsupported = {
        errno.EINVAL,
        errno.ENOSYS,
        getattr(errno, "ENOTSUP", errno.EINVAL),
        getattr(errno, "EOPNOTSUPP", errno.EINVAL),
    }
    if error_number not in unsupported:
        raise OSError(error_number, os.strerror(error_number), destination_name)
    raise CampaignPlanError(
        "campaign-plan publication requires atomic no-replace rename support"
    )


def _remove_anchored_tree_at(
    parent_descriptor: int, name: str, descriptor: int, *, label: str
) -> None:
    """Remove one hidden tree without reopening its ancestor path."""
    _require_same_directory_at(parent_descriptor, name, descriptor, label=label)

    def remove_members(directory_descriptor: int) -> None:
        try:
            status = os.fstat(directory_descriptor)
            os.fchmod(
                directory_descriptor,
                stat.S_IMODE(status.st_mode)
                | stat.S_IRUSR
                | stat.S_IWUSR
                | stat.S_IXUSR,
            )
            names = os.listdir(directory_descriptor)
        except OSError as error:
            raise CampaignPlanError(f"cannot prepare {label} for removal") from error
        for member_name in names:
            try:
                observed = os.stat(
                    member_name,
                    dir_fd=directory_descriptor,
                    follow_symlinks=False,
                )
                if stat.S_ISDIR(observed.st_mode):
                    child_descriptor = os.open(
                        member_name, _DIRECTORY_FLAGS, dir_fd=directory_descriptor
                    )
                    try:
                        opened = os.fstat(child_descriptor)
                        _require(
                            (observed.st_dev, observed.st_ino)
                            == (opened.st_dev, opened.st_ino),
                            f"{label} changed during descriptor-relative removal",
                        )
                        remove_members(child_descriptor)
                        current = os.stat(
                            member_name,
                            dir_fd=directory_descriptor,
                            follow_symlinks=False,
                        )
                        _require(
                            (current.st_dev, current.st_ino)
                            == (opened.st_dev, opened.st_ino),
                            f"{label} changed during descriptor-relative removal",
                        )
                    finally:
                        os.close(child_descriptor)
                    os.rmdir(member_name, dir_fd=directory_descriptor)
                else:
                    os.unlink(member_name, dir_fd=directory_descriptor)
            except CampaignPlanError:
                raise
            except OSError as error:
                raise CampaignPlanError(f"cannot remove member from {label}") from error

    remove_members(descriptor)
    _require_same_directory_at(parent_descriptor, name, descriptor, label=label)
    try:
        os.rmdir(name, dir_fd=parent_descriptor)
    except OSError as error:
        raise CampaignPlanError(f"cannot remove {label}") from error


def _rollback_published_destination(
    parent_descriptor: int, destination_name: str, descriptor: int
) -> None:
    """Withdraw only the pinned invalid plan, never a path replacement."""
    rollback_name = f".{destination_name}.rollback-{uuid.uuid4()}"

    def quarantine_moved_original() -> None:
        expected = os.fstat(descriptor)
        candidates = []
        for name in os.listdir(parent_descriptor):
            try:
                observed = os.stat(
                    name, dir_fd=parent_descriptor, follow_symlinks=False
                )
            except FileNotFoundError:
                continue
            if (
                stat.S_ISDIR(observed.st_mode)
                and (observed.st_dev, observed.st_ino)
                == (expected.st_dev, expected.st_ino)
            ):
                candidates.append(name)
        if len(candidates) != 1:
            raise CampaignPlanError(
                "cannot locate raced campaign-plan publication for quarantine"
            )
        _rename_no_replace_at(
            parent_descriptor, candidates[0], parent_descriptor, rollback_name
        )
        os.fsync(parent_descriptor)
        _require_same_directory_at(
            parent_descriptor,
            rollback_name,
            descriptor,
            label="campaign-plan rollback tree",
        )

    raced_replacement = False
    try:
        _rename_no_replace_at(
            parent_descriptor, destination_name, parent_descriptor, rollback_name
        )
    except FileNotFoundError:
        raced_replacement = True
        os.fsync(parent_descriptor)
        quarantine_moved_original()
    else:
        os.fsync(parent_descriptor)
        try:
            _require_same_directory_at(
                parent_descriptor,
                rollback_name,
                descriptor,
                label="campaign-plan rollback tree",
            )
        except CampaignPlanError:
            raced_replacement = True
            _rename_no_replace_at(
                parent_descriptor, rollback_name, parent_descriptor, destination_name
            )
            os.fsync(parent_descriptor)
            quarantine_moved_original()
    if not raced_replacement:
        _require_absent_at(
            parent_descriptor,
            destination_name,
            label="invalid published campaign-plan output root after rollback",
        )
    rollback_descriptor = os.open(
        rollback_name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor
    )
    try:
        _remove_anchored_tree_at(
            parent_descriptor,
            rollback_name,
            rollback_descriptor,
            label="campaign-plan rollback tree",
        )
    finally:
        os.close(rollback_descriptor)
    os.fsync(parent_descriptor)
    if raced_replacement:
        raise CampaignPlanError(
            "campaign-plan rollback rejected a substituted public destination"
        )


def _cleanup_private_container(
    parent_descriptor: int, name: str, descriptor: int
) -> None:
    """Best-effort removal of one empty pinned private staging container."""
    try:
        _require_same_directory_at(
            parent_descriptor,
            name,
            descriptor,
            label="campaign-plan private staging container",
        )
        os.rmdir(name, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
    except (CampaignPlanError, OSError):
        return


def _write_new_file(root: Path, root_descriptor: int, relative: str, payload: bytes) -> None:
    path = PurePosixPath(relative)
    _require(
        not path.is_absolute()
        and path.parts
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"generated output member path is unsafe: {relative!r}",
    )
    descriptor = os.dup(root_descriptor)
    try:
        for part in path.parts[:-1]:
            try:
                os.mkdir(part, mode=0o700, dir_fd=descriptor)
            except FileExistsError:
                pass
            else:
                os.fsync(descriptor)
            child = os.open(part, _DIRECTORY_FLAGS, dir_fd=descriptor)
            os.close(descriptor)
            descriptor = child
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
        try:
            output = os.open(path.parts[-1], flags, 0o600, dir_fd=descriptor)
        except OSError as error:
            raise CampaignPlanError(f"refusing to overwrite output file: {root / relative}") from error
        try:
            offset = 0
            while offset < len(payload):
                written = os.write(output, payload[offset:])
                _require(written > 0, f"short write while publishing output file: {root / relative}")
                offset += written
            os.fsync(output)
        finally:
            os.close(output)
        os.fsync(descriptor)
        _require_same_directory(root, root_descriptor, label="campaign-plan output root")
    except OSError as error:
        raise CampaignPlanError(f"cannot publish generated output member: {root / relative}") from error
    finally:
        os.close(descriptor)


def _member_inventory_payload(members: Mapping[str, bytes]) -> bytes:
    return "".join(
        f"{_sha256_bytes(payload)}  {relative}\n"
        for relative, payload in sorted(members.items())
    ).encode("utf-8")


def _expected_directories(members: Mapping[str, bytes]) -> set[str]:
    directories: set[str] = set()
    for relative in members:
        parent = PurePosixPath(relative).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _validate_staged_inventory(
    root_descriptor: int,
    expected_members: Mapping[str, bytes],
    *,
    label: str,
) -> None:
    snapshot = immutable_orion_tree._scan_anchored_tree(
        root_descriptor,
        hash_regular=True,
        error_type=CampaignPlanError,
        label=label,
    )
    measured_files = {
        entry.relative: entry.sha256
        for entry in snapshot.entries
        if entry.entry_type == "file"
    }
    expected_files = {
        relative: _sha256_bytes(payload) for relative, payload in expected_members.items()
    }
    _require(
        measured_files == expected_files,
        f"{label} file whitelist or inventory drifted",
    )
    measured_directories = {
        entry.relative
        for entry in snapshot.entries
        if entry.entry_type == "directory"
    }
    _require(
        measured_directories == _expected_directories(expected_members),
        f"{label} directory whitelist drifted",
    )


def _reserve_staging_root(
    output_parent: Path, plan_id: str
) -> tuple[Path, Path, Path, int, int, int]:
    parent = _canonical_existing_directory(output_parent, label="authorized output parent")
    name = f"q011-section54-qualifying-campaign-plan-{plan_id}"
    destination = parent / name
    private_container = parent / f".{name}.staging-{uuid.uuid4()}"
    staging = private_container / _STAGING_ROOT_NAME
    parent_descriptor = os.open(parent, _DIRECTORY_FLAGS)
    private_descriptor = -1
    root_descriptor = -1
    try:
        _require_same_directory(
            parent, parent_descriptor, label="campaign-plan output parent"
        )
        _require_absent_at(
            parent_descriptor, destination.name, label="deterministic campaign-plan output root"
        )
        os.mkdir(private_container.name, mode=0o700, dir_fd=parent_descriptor)
        os.fsync(parent_descriptor)
        private_descriptor = os.open(
            private_container.name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor
        )
        os.mkdir(staging.name, mode=0o700, dir_fd=private_descriptor)
        os.fsync(private_descriptor)
        root_descriptor = os.open(
            staging.name, _DIRECTORY_FLAGS, dir_fd=private_descriptor
        )
    except OSError as error:
        if root_descriptor >= 0:
            os.close(root_descriptor)
        if private_descriptor >= 0:
            try:
                os.rmdir(staging.name, dir_fd=private_descriptor)
            except OSError:
                pass
            os.close(private_descriptor)
        try:
            os.rmdir(private_container.name, dir_fd=parent_descriptor)
        except OSError:
            pass
        os.close(parent_descriptor)
        raise CampaignPlanError("cannot reserve hidden campaign-plan staging root") from error
    except BaseException:
        if root_descriptor >= 0:
            os.close(root_descriptor)
        if private_descriptor >= 0:
            try:
                os.rmdir(staging.name, dir_fd=private_descriptor)
            except OSError:
                pass
            os.close(private_descriptor)
        try:
            os.rmdir(private_container.name, dir_fd=parent_descriptor)
        except OSError:
            pass
        os.close(parent_descriptor)
        raise
    try:
        os.fsync(parent_descriptor)
        _require_same_directory(
            staging, root_descriptor, label="campaign-plan staging root"
        )
    except BaseException:
        try:
            _remove_anchored_tree_at(
                private_descriptor,
                staging.name,
                root_descriptor,
                label="campaign-plan staging root",
            )
        except BaseException:
            pass
        os.close(root_descriptor)
        _cleanup_private_container(
            parent_descriptor, private_container.name, private_descriptor
        )
        os.close(private_descriptor)
        os.close(parent_descriptor)
        raise
    return (
        destination,
        private_container,
        staging,
        parent_descriptor,
        private_descriptor,
        root_descriptor,
    )


def _helper_source_closure() -> list[dict[str, str]]:
    records = [_repo_binding(relative) for relative in _FIXED_HELPER_SOURCES]
    _require(
        [record["path"] for record in records] == list(_FIXED_HELPER_SOURCES),
        "helper source closure order drifted",
    )
    _require(
        len({record["path"] for record in records}) == len(records),
        "helper source closure contains duplicate paths",
    )
    return records


def _retained_binding(value: object, *, label: str) -> dict[str, str]:
    binding = _object(value, {"path", "sha256"}, label=f"{label} binding")
    relative = _relative_path(binding["path"], label=f"{label} binding/path")
    digest = _text(binding["sha256"], label=f"{label} binding/sha256")
    _require(
        _SHA256.fullmatch(digest) is not None,
        f"{label} binding/sha256 must be 64 lowercase hex digits",
    )
    return {"path": relative, "sha256": digest}


def _retained_member_payload(
    snapshot: Any, value: object, *, label: str
) -> bytes:
    binding = _retained_binding(value, label=label)
    try:
        payload = snapshot.member_path(binding["path"]).read_bytes()
    except (OSError, ValueError) as error:
        raise CampaignPlanError(
            f"{label} is unavailable in the immutable qualifying planner tree"
        ) from error
    _require(
        _sha256_bytes(payload) == binding["sha256"],
        f"{label} SHA-256 drifted",
    )
    return payload


def _retained_json_member(
    snapshot: Any, value: object, *, label: str
) -> dict[str, Any]:
    payload = _retained_member_payload(snapshot, value, label=label)
    decoded = _decode_json(payload, label=label)
    _require(type(decoded) is dict, f"{label} must be an object")
    return decoded


def _validate_retained_helper_source_closure(
    snapshot: Any, plan: Mapping[str, Any], *, plan_id: str
) -> list[dict[str, str]]:
    binding = _retained_binding(
        plan["helper_source_closure"], label="helper source closure"
    )
    _require(
        binding["path"] == "helper_source_closure.json",
        "helper source closure path drifted",
    )
    closure = _retained_json_member(
        snapshot, binding, label="helper source closure"
    )
    _object(
        closure,
        {"record_type", "schema_version", "plan_id", "sources"},
        label="helper source closure",
    )
    _require(
        closure["record_type"] == "q011_section54_helper_source_closure"
        and _exact_int(
            closure["schema_version"], label="helper source closure/schema_version"
        )
        == 1
        and closure["plan_id"] == plan_id,
        "helper source closure identity drifted",
    )
    expected_sources = _helper_source_closure()
    _require(
        _list(closure["sources"], label="helper source closure/sources")
        == expected_sources,
        "helper source closure drifted from reviewed source bytes",
    )
    return expected_sources


def materialize_planner_retention(
    *,
    planner_root: str | Path,
    planner_inventory_sha256: str,
    attempt_id: str,
    authorized_pic_root: str | Path = AUTHORIZED_ORION_ROOT,
) -> dict[str, object]:
    """Emit one source-bound, launch-prohibited pre-submit retention overlay."""
    authorized_root = _canonical_existing_directory(
        Path(authorized_pic_root), label="authorized Orion PIC root"
    )
    inventory_sha256 = _text(
        planner_inventory_sha256, label="qualifying planner inventory SHA-256"
    )
    _require(
        _SHA256.fullmatch(inventory_sha256) is not None,
        "qualifying planner inventory SHA-256 must be 64 lowercase hex digits",
    )
    selected_attempt_id = _text(attempt_id, label="baseline attempt ID")
    _require(
        _SAFE_SEGMENT.fullmatch(selected_attempt_id) is not None,
        "baseline attempt ID is unsafe",
    )
    retained_planner_root = Path(planner_root)
    _require(
        retained_planner_root.is_absolute()
        and retained_planner_root.parent == authorized_root / "plans",
        "qualifying planner root must use the dedicated retained plans namespace",
    )
    try:
        with immutable_orion_tree.staged_verified_frozen_tree(
            retained_planner_root,
            inventory_sha256,
            authorized_root=authorized_root,
            error_type=CampaignPlanError,
            label="Q-011 immutable qualifying campaign plan",
        ) as (_planner_report, snapshot):
            try:
                plan_payload = snapshot.member_path("campaign_plan.json").read_bytes()
            except (OSError, ValueError) as error:
                raise CampaignPlanError(
                    "campaign plan is unavailable in the immutable "
                    "qualifying planner tree"
                ) from error
            plan = _decode_json(plan_payload, label="campaign plan")
            plan = _object(
                plan,
                {
                    "record_type",
                    "schema_version",
                    "plan_id",
                    "artifact_role",
                    "qualification_effect",
                    "status",
                    "authorized_orion_root",
                    "authorized_orion_campaign_root",
                    "selected_pressure",
                    "candidate_binding",
                    "source_bindings",
                    "helper_source_closure",
                    "campaign_matrix",
                    "baseline_attempt_count",
                    "baseline_attempt_descriptors",
                    "restart_continuation_carrier",
                    "independent_raw_artifact_recompute_plan",
                    "nonauthorizing_policy_fragment",
                    "execution_boundary",
                    "preregistration_execution_boundary",
                },
                label="campaign plan",
            )
            plan_id = _text(plan["plan_id"], label="campaign plan/plan_id")
            _require(
                _SHA256.fullmatch(plan_id) is not None,
                "campaign plan/plan_id must be 64 lowercase hex digits",
            )
            _require(
                plan["record_type"] == PLAN_RECORD_TYPE
                and _exact_int(
                    plan["schema_version"], label="campaign plan/schema_version"
                )
                == 1
                and plan["artifact_role"] == ARTIFACT_ROLE
                and plan["qualification_effect"] == QUALIFICATION_EFFECT
                and plan["status"] == "source_local_immutable_review_plan_only"
                and plan["authorized_orion_root"] == str(authorized_root),
                "campaign plan identity or authorization boundary drifted",
            )
            campaign_root = (
                authorized_root / "campaigns" / f"q011-section54-{plan_id}"
            )
            _require(
                plan["authorized_orion_campaign_root"] == str(campaign_root),
                "campaign plan authorized Orion campaign root drifted",
            )
            execution_boundary = _object(
                plan["execution_boundary"],
                {
                    "mutates_live_policy",
                    "scheduler_calls",
                    "submits_jobs",
                    "infers_pressure_selection",
                    "launch_authorized",
                    "frontier_execution_authorized",
                    "claim_closure_authorized",
                },
                label="campaign plan/execution_boundary",
            )
            _require(
                all(value is False for value in execution_boundary.values()),
                "campaign plan must remain launch-prohibited",
            )
            matrix = validate_campaign_matrix(plan["campaign_matrix"])
            selected_pressure = _object(
                plan["selected_pressure"],
                {"selection_method", "selected_case", "receipt"},
                label="campaign plan/selected_pressure",
            )
            selected_case = _object(
                selected_pressure["selected_case"],
                {"case_id", "problem_ps_p0"},
                label="campaign plan/selected_pressure/selected_case",
            )
            selected_ps_p0 = _exact_float(
                selected_case["problem_ps_p0"],
                label="campaign plan/selected_pressure/selected_case/problem_ps_p0",
            )
            source_bindings = _object(
                plan["source_bindings"],
                {
                    "pressure_selection_receipt",
                    "clean_candidate_manifest",
                    "environment_profile",
                    "qualifying_preregistration",
                    "restart_preregistration",
                    "paper_deck",
                },
                label="campaign plan/source_bindings",
            )
            normalized_source_bindings = {
                name: _retained_binding(binding, label=f"source binding {name}")
                for name, binding in source_bindings.items()
            }
            _require(
                selected_pressure["receipt"]
                == normalized_source_bindings["pressure_selection_receipt"],
                "campaign plan selected-pressure receipt binding drifted",
            )
            helper_sources = _validate_retained_helper_source_closure(
                snapshot, plan, plan_id=plan_id
            )
            candidate = plan["candidate_binding"]
            _require(
                type(candidate) is dict,
                "campaign plan candidate binding must be an object",
            )
            expected_plan_id = _digest_value(
                {
                    "record_type": PLAN_RECORD_TYPE,
                    "schema_version": 1,
                    "pressure_selection_receipt_sha256": normalized_source_bindings[
                        "pressure_selection_receipt"
                    ]["sha256"],
                    "selected_case": selected_case,
                    "candidate_binding": candidate,
                    "source_binding_sha256": {
                        name: binding["sha256"]
                        for name, binding in normalized_source_bindings.items()
                    },
                    "helper_source_closure": helper_sources,
                    "campaign_matrix": matrix,
                    "authorized_orion_root": str(authorized_root),
                }
            )
            _require(
                plan_id == expected_plan_id,
                "campaign plan ID drifted from its reviewed source-bound basis",
            )
            descriptor_bindings = _list(
                plan["baseline_attempt_descriptors"],
                label="campaign plan/baseline_attempt_descriptors",
            )
            _require(
                _exact_int(
                    plan["baseline_attempt_count"],
                    label="campaign plan/baseline_attempt_count",
                )
                == EXPECTED_BASELINE_ATTEMPTS
                and len(descriptor_bindings) == EXPECTED_BASELINE_ATTEMPTS,
                "campaign plan baseline attempt count drifted",
            )
            expected_attempts = []
            index = 0
            for variant_id in matrix["grid_variants"]:
                variant = _variant_binding(variant_id)
                for seed in matrix["qualifying_seeds"]:
                    index += 1
                    generated_attempt_id = _attempt_id(index, variant.variant, seed)
                    expected_attempts.append(
                        (index, generated_attempt_id, variant, seed)
                    )
            normalized_descriptor_bindings = [
                _retained_binding(binding, label="baseline attempt descriptor")
                for binding in descriptor_bindings
            ]
            _require(
                [
                    binding["path"]
                    for binding in normalized_descriptor_bindings
                ]
                == [
                    f"attempts/baseline/{generated_attempt_id}.json"
                    for _, generated_attempt_id, _, _ in expected_attempts
                ],
                "campaign plan baseline descriptor path ordering drifted",
            )
            selected = [
                (binding, expected)
                for binding, expected in zip(
                    normalized_descriptor_bindings, expected_attempts
                )
                if expected[1] == selected_attempt_id
            ]
            _require(
                len(selected) == 1,
                "attempt ID does not select exactly one baseline planner descriptor",
            )
            descriptor_binding, (
                attempt_index,
                _generated_attempt_id,
                variant,
                qualifying_seed,
            ) = selected[0]
            descriptor = _retained_json_member(
                snapshot,
                descriptor_binding,
                label="baseline attempt descriptor",
            )
            attempt_root = campaign_root / "baseline" / selected_attempt_id
            contract_path = f"launch_contracts/baseline/{selected_attempt_id}.json"
            contract_payload = _retained_member_payload(
                snapshot,
                descriptor.get("launch_contract"),
                label="baseline attempt launch contract",
            )
            contract = _decode_json(
                contract_payload, label="baseline attempt launch contract"
            )
            _require(
                type(contract) is dict,
                "baseline attempt launch contract must be an object",
            )
            try:
                expected_contract = _baseline_launch_contract(
                    attempt_id=selected_attempt_id,
                    variant=variant,
                    seed=qualifying_seed,
                    selected_ps_p0=selected_ps_p0,
                    candidate=candidate,
                    paper_deck_binding=normalized_source_bindings["paper_deck"],
                    artifact_root=attempt_root,
                )
            except (KeyError, TypeError) as error:
                raise CampaignPlanError(
                    "campaign plan candidate binding drifted"
                ) from error
            immutable_orion_tree.require_exact_primitive_types(
                contract,
                expected_contract,
                error_type=CampaignPlanError,
                label="baseline attempt launch contract",
            )
            _require(
                contract == expected_contract,
                "baseline attempt launch contract drifted",
            )
            expected_descriptor = _attempt_descriptor(
                index=attempt_index,
                attempt_id=selected_attempt_id,
                variant=variant,
                seed=qualifying_seed,
                selected_ps_p0=selected_ps_p0,
                candidate=candidate,
                artifact_root=attempt_root,
                contract_path=contract_path,
                contract_payload=contract_payload,
            )
            immutable_orion_tree.require_exact_primitive_types(
                descriptor,
                expected_descriptor,
                error_type=CampaignPlanError,
                label="baseline attempt descriptor",
            )
            _require(
                descriptor == expected_descriptor,
                "baseline attempt descriptor drifted",
            )
            materialization_receipt_payload = _retained_member_payload(
                snapshot,
                {
                    "path": MATERIALIZATION_RECEIPT_NAME,
                    "sha256": _sha256_bytes(
                        snapshot.member_path(MATERIALIZATION_RECEIPT_NAME).read_bytes()
                    ),
                },
                label="campaign-plan materialization receipt",
            )
            materialization_receipt = _decode_json(
                materialization_receipt_payload,
                label="campaign-plan materialization receipt",
            )
            _require(
                type(materialization_receipt) is dict
                and materialization_receipt.get("record_type")
                == MATERIALIZATION_RECEIPT_RECORD_TYPE
                and materialization_receipt.get("schema_version") == 1
                and materialization_receipt.get("plan_id") == plan_id,
                "campaign-plan materialization receipt identity drifted",
            )
            return {
                "schema_version": 1,
                "retention_role": PLANNER_RETENTION_ROLE,
                "planner_root": str(Path(planner_root)),
                "planner_inventory_sha256": inventory_sha256,
                "planner_plan_id": plan_id,
                "planner_materialization_receipt": {
                    "path": MATERIALIZATION_RECEIPT_NAME,
                    "sha256": _sha256_bytes(materialization_receipt_payload),
                },
                "attempt_id": selected_attempt_id,
                "authorized_orion_attempt_root": str(attempt_root),
                "authorized_orion_raw_root": str(attempt_root / "raw"),
                "argv": list(contract["argv"]),
            }
    except CampaignPlanError:
        raise
    except (OSError, TypeError, ValueError) as error:
        raise CampaignPlanError(
            f"immutable qualifying planner retention materialization failed: {error}"
        ) from error


def materialize_qualifying_campaign_plan(
    *,
    output_parent: Path,
    pressure_selection_receipt: Path,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    qualifying_preregistration: Path = QUALIFYING_PREREGISTRATION,
    restart_preregistration: Path = RESTART_PREREGISTRATION,
    paper_deck: Path = PAPER_DECK,
) -> dict[str, object]:
    """Create one deterministic recursively read-only source-local plan tree."""
    pressure_payload, pressure_receipt = _load_pressure_receipt(
        pressure_selection_receipt,
        authorized_pic_root=AUTHORIZED_ORION_ROOT,
    )
    qualifying_payload, qualifying_policy, matrix = _load_qualifying_preregistration(
        qualifying_preregistration
    )
    restart_payload, restart_policy = _load_restart_preregistration(restart_preregistration)
    _, deck_payload = _stable_regular_bytes(paper_deck, label="Section 5.4 paper deck")
    _require(_sha256_bytes(deck_payload) == PAPER_DECK_SHA256, "Section 5.4 paper deck SHA-256 drifted")
    helper_sources = _helper_source_closure()
    helper_by_path = {record["path"]: record["sha256"] for record in helper_sources}
    manifest_payload, environment_payload, candidate = _candidate_binding(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
        paper_deck_sha256=PAPER_DECK_SHA256,
        admission_analyzer_sha256=helper_by_path["tst/publication/analyze_q011_section54_campaign.py"],
        output_primitives_sha256=helper_by_path["tst/publication/analyze_q011_section54_outputs.py"],
    )
    selected_case = pressure_receipt["selected_case"]
    selected_ps_p0 = _exact_float(
        selected_case["problem_ps_p0"], label="selected pressure receipt/problem_ps_p0"
    )
    source_bindings = {
        "pressure_selection_receipt": _binding(
            "bindings/human_pressure_selection_receipt.json", pressure_payload
        ),
        "clean_candidate_manifest": _binding(
            "bindings/clean_candidate_manifest.json", manifest_payload
        ),
        "environment_profile": _binding(
            "bindings/environment_profile.sh", environment_payload
        ),
        "qualifying_preregistration": _binding(
            "bindings/q011_section54_qualifying_campaign_preregistration.json",
            qualifying_payload,
        ),
        "restart_preregistration": _binding(
            "bindings/q011_section54_restart_continuation_preregistration.json",
            restart_payload,
        ),
        "paper_deck": _binding(
            "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
            deck_payload,
        ),
    }
    plan_basis = {
        "record_type": PLAN_RECORD_TYPE,
        "schema_version": 1,
        "pressure_selection_receipt_sha256": source_bindings["pressure_selection_receipt"]["sha256"],
        "selected_case": selected_case,
        "candidate_binding": candidate,
        "source_binding_sha256": {
            name: binding["sha256"] for name, binding in source_bindings.items()
        },
        "helper_source_closure": helper_sources,
        "campaign_matrix": matrix,
        "authorized_orion_root": str(AUTHORIZED_ORION_ROOT),
    }
    plan_id = _digest_value(plan_basis)
    campaign_root = AUTHORIZED_ORION_ROOT / "campaigns" / f"q011-section54-{plan_id}"
    (
        destination,
        private_container,
        root,
        parent_descriptor,
        private_descriptor,
        root_descriptor,
    ) = _reserve_staging_root(output_parent, plan_id)
    expected_members: dict[str, bytes] = {}
    renamed = False

    def write(relative: str, payload: bytes) -> None:
        _require(
            relative not in expected_members,
            f"generated output member path is duplicated: {relative}",
        )
        _write_new_file(root, root_descriptor, relative, payload)
        expected_members[relative] = payload

    try:
        write(source_bindings["pressure_selection_receipt"]["path"], pressure_payload)
        write(source_bindings["clean_candidate_manifest"]["path"], manifest_payload)
        write(source_bindings["environment_profile"]["path"], environment_payload)
        write(source_bindings["qualifying_preregistration"]["path"], qualifying_payload)
        write(source_bindings["restart_preregistration"]["path"], restart_payload)
        write(source_bindings["paper_deck"]["path"], deck_payload)

        contract_bindings = []
        descriptor_bindings = []
        descriptors = []
        index = 0
        for variant_id in matrix["grid_variants"]:
            variant = _variant_binding(variant_id)
            for seed in matrix["qualifying_seeds"]:
                index += 1
                attempt_id = _attempt_id(index, variant.variant, seed)
                artifact_root = campaign_root / "baseline" / attempt_id
                contract_path = f"launch_contracts/baseline/{attempt_id}.json"
                contract = _baseline_launch_contract(
                    attempt_id=attempt_id,
                    variant=variant,
                    seed=seed,
                    selected_ps_p0=selected_ps_p0,
                    candidate=candidate,
                    paper_deck_binding=source_bindings["paper_deck"],
                    artifact_root=artifact_root,
                )
                contract_payload = _json_bytes(contract)
                write(contract_path, contract_payload)
                contract_binding = _binding(contract_path, contract_payload)
                contract_bindings.append(contract_binding)
                descriptor = _attempt_descriptor(
                    index=index,
                    attempt_id=attempt_id,
                    variant=variant,
                    seed=seed,
                    selected_ps_p0=selected_ps_p0,
                    candidate=candidate,
                    artifact_root=artifact_root,
                    contract_path=contract_path,
                    contract_payload=contract_payload,
                )
                descriptor_path = f"attempts/baseline/{attempt_id}.json"
                descriptor_payload = _json_bytes(descriptor)
                write(descriptor_path, descriptor_payload)
                descriptor_bindings.append(_binding(descriptor_path, descriptor_payload))
                descriptors.append(descriptor)
        _require(index == EXPECTED_BASELINE_ATTEMPTS, "materialized baseline attempt count drifted")

        source_attempt = next(
            descriptor
            for descriptor in descriptors
            if descriptor["variant"] == "three_level_amr_root_dx12_finest_dx3"
            and descriptor["qualifying_seed"] == matrix["qualifying_seeds"][0]
        )
        carrier_id = _restart_carrier_id(source_attempt["qualifying_seed"])
        restart_artifact_root = campaign_root / "restart_continuation" / carrier_id
        restart_contract_path = f"launch_contracts/restart_continuation/{carrier_id}.json"
        restart_contract = _restart_launch_contract(
            carrier_id=carrier_id,
            source_attempt=source_attempt,
            restart_preregistration=restart_policy,
            candidate=candidate,
            paper_deck_binding=source_bindings["paper_deck"],
            artifact_root=restart_artifact_root,
        )
        restart_contract_payload = _json_bytes(restart_contract)
        write(restart_contract_path, restart_contract_payload)
        restart_contract_binding = _binding(restart_contract_path, restart_contract_payload)
        restart_carrier = {
            "record_type": RESTART_CARRIER_RECORD_TYPE,
            "schema_version": 1,
            "carrier_id": carrier_id,
            "status": "planned_not_authorized",
            "source_baseline_attempt_id": source_attempt["attempt_id"],
            "variant": source_attempt["variant"],
            "qualifying_seed": source_attempt["qualifying_seed"],
            "selected_problem_ps_p0": selected_ps_p0,
            "authorized_orion_attempt_root": str(restart_artifact_root),
            "restart_preregistration": source_bindings["restart_preregistration"],
            "checkpoint_time_omega0_inverse": restart_policy["continuation_contract"][
                "checkpoint_time_omega0_inverse"
            ],
            "retained_output_schedule_after_checkpoint_omega0_inverse": restart_policy[
                "continuation_contract"
            ]["retained_output_schedule_after_checkpoint_omega0_inverse"],
            "comparison_tolerances_max_absolute_difference": restart_policy[
                "continuation_contract"
            ]["comparison_tolerances_max_absolute_difference"],
            "launch_contract": restart_contract_binding,
        }
        restart_carrier_payload = _json_bytes(restart_carrier)
        restart_carrier_path = "restart_continuation/amr_restart_continuation_carrier.json"
        write(restart_carrier_path, restart_carrier_payload)
        restart_carrier_binding = _binding(restart_carrier_path, restart_carrier_payload)

        helper_closure = {
            "record_type": "q011_section54_helper_source_closure",
            "schema_version": 1,
            "plan_id": plan_id,
            "sources": helper_sources,
        }
        helper_closure_payload = _json_bytes(helper_closure)
        helper_closure_path = "helper_source_closure.json"
        write(helper_closure_path, helper_closure_payload)
        helper_closure_binding = _binding(helper_closure_path, helper_closure_payload)

        recompute = _independent_recompute_plan(
            plan_id=plan_id,
            campaign_root=campaign_root,
            qualifying_preregistration_binding=source_bindings["qualifying_preregistration"],
        )
        recompute_payload = _json_bytes(recompute)
        recompute_path = "independent_raw_artifact_recompute_plan.json"
        write(recompute_path, recompute_payload)
        recompute_binding = _binding(recompute_path, recompute_payload)

        fragment = _policy_fragment(
            plan_id=plan_id,
            campaign_root=campaign_root,
            candidate=candidate,
            pressure_receipt_binding=source_bindings["pressure_selection_receipt"],
            contract_bindings=contract_bindings,
            restart_contract_binding=restart_contract_binding,
        )
        fragment_payload = _json_bytes(fragment)
        fragment_path = "nonauthorizing_policy_fragment.json"
        write(fragment_path, fragment_payload)
        fragment_binding = _binding(fragment_path, fragment_payload)

        campaign_plan = {
            "record_type": PLAN_RECORD_TYPE,
            "schema_version": 1,
            "plan_id": plan_id,
            "artifact_role": ARTIFACT_ROLE,
            "qualification_effect": QUALIFICATION_EFFECT,
            "status": "source_local_immutable_review_plan_only",
            "authorized_orion_root": str(AUTHORIZED_ORION_ROOT),
            "authorized_orion_campaign_root": str(campaign_root),
            "selected_pressure": {
                "selection_method": pressure_receipt["selection_method"],
                "selected_case": selected_case,
                "receipt": source_bindings["pressure_selection_receipt"],
            },
            "candidate_binding": candidate,
            "source_bindings": source_bindings,
            "helper_source_closure": helper_closure_binding,
            "campaign_matrix": matrix,
            "baseline_attempt_count": len(descriptor_bindings),
            "baseline_attempt_descriptors": descriptor_bindings,
            "restart_continuation_carrier": restart_carrier_binding,
            "independent_raw_artifact_recompute_plan": recompute_binding,
            "nonauthorizing_policy_fragment": fragment_binding,
            "execution_boundary": {
                "mutates_live_policy": False,
                "scheduler_calls": False,
                "submits_jobs": False,
                "infers_pressure_selection": False,
                "launch_authorized": False,
                "frontier_execution_authorized": False,
                "claim_closure_authorized": False,
            },
            "preregistration_execution_boundary": qualifying_policy[
                "qualifying_execution_bindings"
            ],
        }
        campaign_plan_payload = _json_bytes(campaign_plan)
        write("campaign_plan.json", campaign_plan_payload)
        materialized_inventory_payload = _member_inventory_payload(expected_members)
        materialization_receipt = {
            "record_type": MATERIALIZATION_RECEIPT_RECORD_TYPE,
            "schema_version": 1,
            "plan_id": plan_id,
            "campaign_plan": _binding("campaign_plan.json", campaign_plan_payload),
            "helper_source_closure": helper_closure_binding,
            "tree_inventory": {
                "algorithm": _MATERIALIZED_MEMBER_INVENTORY_ALGORITHM,
                "scope": (
                    "all materialized campaign-plan members before this receipt "
                    "and recursive-freeze metadata"
                ),
                "excludes": [
                    MATERIALIZATION_RECEIPT_NAME,
                    immutable_orion_tree.FREEZE_RECEIPT_NAME,
                    immutable_orion_tree.INVENTORY_NAME,
                ],
                "sha256": _sha256_bytes(materialized_inventory_payload),
                "inventoried_file_count": len(expected_members),
            },
        }
        materialization_receipt_payload = _json_bytes(materialization_receipt)
        write(MATERIALIZATION_RECEIPT_NAME, materialization_receipt_payload)
        os.fsync(root_descriptor)
        _require_same_directory(root, root_descriptor, label="campaign-plan staging root")
        _validate_staged_inventory(
            root_descriptor,
            expected_members,
            label="campaign-plan materialized staging tree",
        )
        freeze = immutable_orion_tree.freeze_tree_anchored(
            root,
            root_descriptor,
            _PLAN_FREEZE_RECEIPT,
            authorized_root=_canonical_existing_directory(
                output_parent, label="authorized output parent"
            ),
            error_type=CampaignPlanError,
            label="Q-011 Section 5.4 immutable campaign plan",
        )
        frozen_inventory_payload = _member_inventory_payload(
            {
                **expected_members,
                immutable_orion_tree.FREEZE_RECEIPT_NAME: _json_bytes(
                    _PLAN_FREEZE_RECEIPT
                ),
            }
        )
        _require(
            freeze["inventory_sha256"] == _sha256_bytes(frozen_inventory_payload),
            "campaign-plan frozen inventory drifted from expected materialized members",
        )
        frozen_members = {
            **expected_members,
            immutable_orion_tree.FREEZE_RECEIPT_NAME: _json_bytes(
                _PLAN_FREEZE_RECEIPT
            ),
            immutable_orion_tree.INVENTORY_NAME: frozen_inventory_payload,
        }
        _require_same_directory(root, root_descriptor, label="campaign-plan staging root")
        _validate_staged_inventory(
            root_descriptor,
            frozen_members,
            label="campaign-plan frozen staging tree",
        )
        immutable_orion_tree._verify_frozen_tree_anchored(
            root,
            root_descriptor,
            freeze["inventory_sha256"],
            authorized_root=_canonical_existing_directory(
                output_parent, label="authorized output parent"
            ),
            error_type=CampaignPlanError,
            label="Q-011 Section 5.4 immutable campaign plan",
        )
        _require_same_directory(
            output_parent,
            parent_descriptor,
            label="campaign-plan output parent",
        )
        _require_same_directory_at(
            private_descriptor,
            root.name,
            root_descriptor,
            label="campaign-plan staging root",
        )
        _require_absent_at(
            parent_descriptor,
            destination.name,
            label="deterministic campaign-plan output root",
        )
        # Orion rejects cross-parent rename of a read-only directory. Descendants
        # remain frozen; make the root owner-write-only and non-traversable for
        # the rename, then restore its exact frozen mode through the pinned fd.
        root_mode = stat.S_IMODE(os.fstat(root_descriptor).st_mode)
        os.fchmod(root_descriptor, stat.S_IWUSR)
        os.fsync(root_descriptor)
        _rename_no_replace_at(
            private_descriptor, root.name, parent_descriptor, destination.name
        )
        renamed = True
        os.fchmod(root_descriptor, root_mode)
        os.fsync(root_descriptor)
        os.fsync(private_descriptor)
        os.fsync(parent_descriptor)
        _require_same_directory_at(
            parent_descriptor,
            destination.name,
            root_descriptor,
            label="published campaign-plan output root",
        )
        verified = immutable_orion_tree._verify_frozen_tree_anchored(
            destination,
            root_descriptor,
            freeze["inventory_sha256"],
            authorized_root=_canonical_existing_directory(
                output_parent, label="authorized output parent"
            ),
            error_type=CampaignPlanError,
            label="Q-011 Section 5.4 immutable campaign plan",
        )
        _validate_staged_inventory(
            root_descriptor,
            frozen_members,
            label="published campaign-plan output tree",
        )
        os.close(root_descriptor)
        root_descriptor = -1
        _cleanup_private_container(
            parent_descriptor, private_container.name, private_descriptor
        )
        os.close(private_descriptor)
        private_descriptor = -1
        os.close(parent_descriptor)
        parent_descriptor = -1
        return {
            "plan_root": str(destination),
            "plan_id": plan_id,
            "campaign_plan_sha256": _sha256_bytes(campaign_plan_payload),
            "materialization_receipt": _binding(
                MATERIALIZATION_RECEIPT_NAME, materialization_receipt_payload
            ),
            "materialized_member_inventory_sha256": _sha256_bytes(
                materialized_inventory_payload
            ),
            "inventory_sha256": freeze["inventory_sha256"],
            "inventoried_file_count": freeze["inventoried_file_count"],
            "baseline_attempt_count": len(descriptor_bindings),
            "restart_continuation_carrier_count": 1,
            "recursively_read_only": verified["recursively_read_only"],
        }
    except BaseException:
        rollback_error: BaseException | None = None
        try:
            if root_descriptor >= 0 and renamed:
                try:
                    _rollback_published_destination(
                        parent_descriptor, destination.name, root_descriptor
                    )
                except BaseException as error:
                    rollback_error = error
            elif root_descriptor >= 0:
                try:
                    _remove_anchored_tree_at(
                        private_descriptor,
                        root.name,
                        root_descriptor,
                        label="campaign-plan staging root",
                    )
                except BaseException:
                    pass
        finally:
            if root_descriptor >= 0:
                os.close(root_descriptor)
            if private_descriptor >= 0:
                _cleanup_private_container(
                    parent_descriptor, private_container.name, private_descriptor
                )
                os.close(private_descriptor)
            if parent_descriptor >= 0:
                os.close(parent_descriptor)
        if rollback_error is not None:
            raise CampaignPlanError(
                "cannot remove invalid published campaign-plan output root"
            ) from rollback_error
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-parent", type=Path)
    parser.add_argument("--pressure-selection-receipt", type=Path)
    parser.add_argument("--clean-candidate-manifest", type=Path)
    parser.add_argument("--executable", type=Path)
    parser.add_argument("--environment-profile", type=Path)
    parser.add_argument("--materialize-planner-retention", action="store_true")
    parser.add_argument("--planner-root", type=Path)
    parser.add_argument("--planner-inventory-sha256")
    parser.add_argument("--attempt-id")
    parser.add_argument("--authorized-pic-root", type=Path)
    arguments = parser.parse_args()
    if arguments.materialize_planner_retention:
        _require(
            arguments.planner_root is not None
            and arguments.planner_inventory_sha256 is not None
            and arguments.attempt_id is not None
            and arguments.authorized_pic_root is not None,
            "planner-retention validation mode requires its complete immutable binding",
        )
        _require(
            arguments.output_parent is None
            and arguments.pressure_selection_receipt is None
            and arguments.clean_candidate_manifest is None
            and arguments.executable is None
            and arguments.environment_profile is None,
            "planner-retention validation mode cannot create a campaign plan",
        )
        result = materialize_planner_retention(
            planner_root=arguments.planner_root,
            planner_inventory_sha256=arguments.planner_inventory_sha256,
            attempt_id=arguments.attempt_id,
            authorized_pic_root=arguments.authorized_pic_root,
        )
    else:
        _require(
            arguments.output_parent is not None
            and arguments.pressure_selection_receipt is not None
            and arguments.clean_candidate_manifest is not None
            and arguments.executable is not None
            and arguments.environment_profile is not None,
            "campaign-plan materialization requires its complete source binding",
        )
        _require(
            arguments.planner_root is None
            and arguments.planner_inventory_sha256 is None
            and arguments.attempt_id is None
            and arguments.authorized_pic_root is None,
            "campaign-plan materialization cannot consume a planner-retention overlay",
        )
        result = materialize_qualifying_campaign_plan(
            output_parent=arguments.output_parent,
            pressure_selection_receipt=arguments.pressure_selection_receipt,
            clean_candidate_manifest=arguments.clean_candidate_manifest,
            executable=arguments.executable,
            environment_profile=arguments.environment_profile,
        )
    print(_json_bytes(result).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
