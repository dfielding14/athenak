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
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
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
EXPECTED_BASELINE_ATTEMPTS = 24
_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_SAFE_SEGMENT = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_FILE_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH


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
    "tst/publication/publish_q011_section54_campaign_attempt.py",
    "tst/publication/immutable_orion_tree.py",
    "tst/publication/pvtk_particles.py",
    "tst/publication/q011_parallel_shock_storage_estimator.py",
    "tst/publication/frontier_control_plane/control_plane_common.py",
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


def _cleanup_created_tree(root: Path) -> None:
    if not os.path.lexists(root):
        return
    for directory, names, filenames in os.walk(root, topdown=False, followlinks=False):
        base = Path(directory)
        for name in filenames:
            path = base / name
            if not path.is_symlink():
                os.chmod(path, 0o600, follow_symlinks=False)
        for name in names:
            path = base / name
            if not path.is_symlink():
                os.chmod(path, 0o700, follow_symlinks=False)
        if not base.is_symlink():
            os.chmod(base, 0o700, follow_symlinks=False)
    shutil.rmtree(root)


def _reserve_output_root(output_parent: Path, plan_id: str) -> tuple[Path, int, int]:
    parent = _canonical_existing_directory(output_parent, label="authorized output parent")
    name = f"q011-section54-qualifying-campaign-plan-{plan_id}"
    root = parent / name
    parent_descriptor = os.open(parent, _DIRECTORY_FLAGS)
    created = False
    try:
        os.mkdir(name, mode=0o700, dir_fd=parent_descriptor)
        created = True
        root_descriptor = os.open(name, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
    except OSError as error:
        if created:
            try:
                os.rmdir(name, dir_fd=parent_descriptor)
            except OSError:
                pass
        os.close(parent_descriptor)
        if not created and os.path.lexists(root):
            raise CampaignPlanError("deterministic campaign-plan output root already exists") from error
        raise CampaignPlanError("cannot reserve deterministic campaign-plan output root") from error
    try:
        os.fsync(parent_descriptor)
        _require_same_directory(root, root_descriptor, label="campaign-plan output root")
    except BaseException:
        os.close(root_descriptor)
        os.close(parent_descriptor)
        _cleanup_created_tree(root)
        raise
    return root, parent_descriptor, root_descriptor


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
    root, parent_descriptor, root_descriptor = _reserve_output_root(output_parent, plan_id)
    try:
        _write_new_file(root, root_descriptor, source_bindings["pressure_selection_receipt"]["path"], pressure_payload)
        _write_new_file(root, root_descriptor, source_bindings["clean_candidate_manifest"]["path"], manifest_payload)
        _write_new_file(root, root_descriptor, source_bindings["environment_profile"]["path"], environment_payload)
        _write_new_file(root, root_descriptor, source_bindings["qualifying_preregistration"]["path"], qualifying_payload)
        _write_new_file(root, root_descriptor, source_bindings["restart_preregistration"]["path"], restart_payload)
        _write_new_file(root, root_descriptor, source_bindings["paper_deck"]["path"], deck_payload)

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
                _write_new_file(root, root_descriptor, contract_path, contract_payload)
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
                _write_new_file(root, root_descriptor, descriptor_path, descriptor_payload)
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
        _write_new_file(root, root_descriptor, restart_contract_path, restart_contract_payload)
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
        _write_new_file(root, root_descriptor, restart_carrier_path, restart_carrier_payload)
        restart_carrier_binding = _binding(restart_carrier_path, restart_carrier_payload)

        helper_closure = {
            "record_type": "q011_section54_helper_source_closure",
            "schema_version": 1,
            "plan_id": plan_id,
            "sources": helper_sources,
        }
        helper_closure_payload = _json_bytes(helper_closure)
        helper_closure_path = "helper_source_closure.json"
        _write_new_file(root, root_descriptor, helper_closure_path, helper_closure_payload)
        helper_closure_binding = _binding(helper_closure_path, helper_closure_payload)

        recompute = _independent_recompute_plan(
            plan_id=plan_id,
            campaign_root=campaign_root,
            qualifying_preregistration_binding=source_bindings["qualifying_preregistration"],
        )
        recompute_payload = _json_bytes(recompute)
        recompute_path = "independent_raw_artifact_recompute_plan.json"
        _write_new_file(root, root_descriptor, recompute_path, recompute_payload)
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
        _write_new_file(root, root_descriptor, fragment_path, fragment_payload)
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
        _write_new_file(root, root_descriptor, "campaign_plan.json", campaign_plan_payload)
        os.fsync(root_descriptor)
        _require_same_directory(root, root_descriptor, label="campaign-plan output root")
        os.close(root_descriptor)
        root_descriptor = -1
        os.close(parent_descriptor)
        parent_descriptor = -1
        freeze = immutable_orion_tree.freeze_tree(
            root,
            {
                "schema_version": 1,
                "artifact_role": ARTIFACT_ROLE,
                "qualification_effect": QUALIFICATION_EFFECT,
                "inventory_excludes": immutable_orion_tree.INVENTORY_NAME,
                "freeze_policy": "remove all owner, group and other write bits recursively",
            },
            authorized_root=_canonical_existing_directory(
                output_parent, label="authorized output parent"
            ),
            error_type=CampaignPlanError,
            label="Q-011 Section 5.4 immutable campaign plan",
        )
        return {
            "plan_root": str(root),
            "plan_id": plan_id,
            "campaign_plan_sha256": _sha256_bytes(campaign_plan_payload),
            "inventory_sha256": freeze["inventory_sha256"],
            "inventoried_file_count": freeze["inventoried_file_count"],
            "baseline_attempt_count": len(descriptor_bindings),
            "restart_continuation_carrier_count": 1,
            "recursively_read_only": freeze["recursively_read_only"],
        }
    except BaseException:
        if root_descriptor >= 0:
            os.close(root_descriptor)
        if parent_descriptor >= 0:
            os.close(parent_descriptor)
        _cleanup_created_tree(root)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-parent", required=True, type=Path)
    parser.add_argument("--pressure-selection-receipt", required=True, type=Path)
    parser.add_argument("--clean-candidate-manifest", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--environment-profile", required=True, type=Path)
    arguments = parser.parse_args()
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
