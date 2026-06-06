#!/usr/bin/env python3
"""Non-authorizing Q-011 Section 5.4 production-science admission successor.

This module verifies immutable execution evidence and raw artifact bytes before
building attempt admissions and a phase-aware diagnostic work graph.  It does
not mutate evidence, launch work, promote policy, authorize qualification, or
close a scientific claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
from pathlib import PurePosixPath
import re
import stat
import struct
from typing import Any, Mapping, Sequence

if __package__:
    from . import q011_section54_model as model
else:
    import q011_section54_model as model


SCHEMA_VERSION = 2
SUCCESSOR_ID = (
    "q011_section54_production_science_admission_orchestration_successor_v1"
)
ATTEMPT_RECORD_TYPE = (
    "q011_section54_production_science_attempt_admission_successor_v1"
)
CAMPAIGN_RECORD_TYPE = (
    "q011_section54_production_science_campaign_orchestration_successor_v1"
)
CAMPAIGN_ID = "Q011-SECTION54-PRODUCTION-SCIENCE-SUCCESSOR-V1"
CLAIM_ID = "CLAIM-PAPER-SHOCK-001"
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"

GRID_VARIANTS = (
    "coarse_uniform_dx12",
    "three_level_amr_root_dx12_finest_dx3",
    "fine_uniform_dx3",
)
AMR_VARIANT = "three_level_amr_root_dx12_finest_dx3"
FINE_VARIANT = "fine_uniform_dx3"
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
EXPECTED_MATRIX_CELLS = tuple(
    (variant, seed) for variant in GRID_VARIANTS for seed in QUALIFYING_SEEDS
)
CORE_QUALIFYING_SEEDS = QUALIFYING_SEEDS[:3]
CAMPAIGN_PHASES = {
    f"phase_{count}_paired_triads": QUALIFYING_SEEDS[:count]
    for count in range(3, len(QUALIFYING_SEEDS) + 1)
}
REQUIRED_NOMINAL_SLOTS = tuple(float(value) for value in range(0, 1300, 100))
SCIENCE_SLOTS = (500.0, 1200.0)

SNAPSHOT_PRODUCT_KINDS = (
    "mhd_w_bcc",
    "prtcl_rho",
    "prtcl_jx",
    "prtcl_jy",
    "prtcl_jz",
    "mhd_j2",
    "prtcl_all",
)
RUN_PRODUCT_KINDS = ("hst", "stdout")
RESTART_PRODUCT_KIND = "rst"
ALL_RAW_PRODUCT_KINDS = frozenset(
    (*SNAPSHOT_PRODUCT_KINDS, *RUN_PRODUCT_KINDS, RESTART_PRODUCT_KIND)
)

PREDECESSOR_PREREGISTRATION = (
    "tst/publication/readiness/"
    "q011_section54_qualifying_campaign_preregistration_successor_v3_2026-06-06.json"
)
SUCCESSOR_DECK = (
    "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)
SUCCESSOR_DECK_LAUNCH_PATH = (
    "bindings/pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)
THIS_REDUCER = (
    "tst/publication/"
    "q011_section54_production_science_admission_orchestration_successor_v1.py"
)
REQUIRED_REDUCER_PATHS = frozenset(
    {
        THIS_REDUCER,
        "tst/publication/q011_section54_production_science_successor_v1.py",
        "tst/publication/q011_section54_particles.py",
        "tst/publication/analyze_q011_section54_outputs.py",
        "tst/publication/q011_section54_model.py",
    }
)
REQUIRED_SOURCE_PATHS = frozenset(
    {
        "src/outputs/derived_variables.cpp",
        "src/pgen/tests/pic_parallel_shock.cpp",
        "tst/publication/frontier_control_plane/frontier_pic_environment.sh",
    }
)
SEED_OVERRIDE_NAMES = (
    "particles/pic_random_seed",
    "problem/ps_inject_seed",
    "problem/ps_seed_noise_seed",
)
RESOURCE_TELEMETRY_FIELDS = (
    "scheduler_accounting",
    "runtime_telemetry",
    "active_cell_history",
    "meshblock_history",
    "particle_count_and_updates",
    "memory_high_water_mark",
    "output_timers",
    "artifact_byte_counts",
)

_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_ID = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_KIND_PATH_PATTERNS = {
    "mhd_w_bcc": re.compile(r"(?:[^/]+/)*[^/]+\.mhd_w_bcc\.[^/]+\.bin"),
    "prtcl_rho": re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_rho\.[^/]+\.bin"),
    "prtcl_jx": re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_jx\.[^/]+\.bin"),
    "prtcl_jy": re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_jy\.[^/]+\.bin"),
    "prtcl_jz": re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_jz\.[^/]+\.bin"),
    "mhd_j2": re.compile(r"(?:[^/]+/)*[^/]+\.(?:mhd_j2|j2)\.[^/]+\.bin"),
    "prtcl_all": re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_all\.[^/]+\.part\.vtk"),
    "hst": re.compile(r"(?:[^/]+/)*[^/]+\.hst"),
    "stdout": re.compile(r"stdout\.txt"),
    "rst": re.compile(
        r"rst/(?:[^/]+/)*[^/]+\.rst(?:\.complete|\.manifest(?:\.complete)?)?"
    ),
}
_RESTART_CYCLE = re.compile(rb"(?:^|\n)cycle=([0-9]+)(?:\n|$)")
_RESTART_MARKER = re.compile(
    rb"ATHENAK_RESTART_COMPLETE_V1\n"
    rb"size=(0|[1-9][0-9]*)\n"
    rb"fnv1a64=([0-9a-f]{16})\n"
)
_TERMINAL_STATE = re.compile(
    r"^time=([^\s]+) cycle=(0|[1-9][0-9]*)$", re.MULTILINE
)
_TERMINAL_LIMIT = re.compile(r"^tlim=([^\s]+) nlim=([^\s]+)$", re.MULTILINE)


class ProductionScienceAdmissionError(ValueError):
    """Reject an incomplete, mixed-lineage, or authority-bearing record."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ProductionScienceAdmissionError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _array(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected array")
    return value


def _text(value: object, *, label: str) -> str:
    _require(
        type(value) is str and bool(value) and not any(char.isspace() for char in value),
        f"{label}: expected nonempty whitespace-free text",
    )
    return value


def _relative_path(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    path = PurePosixPath(text)
    _require(
        not path.is_absolute()
        and path.as_posix() == text
        and text != "."
        and all(part not in ("", ".", "..") for part in path.parts),
        f"{label}: unsafe root-relative path",
    )
    return text


def _absolute_path(value: object, *, label: str) -> str:
    text = _text(value, label=label)
    path = PurePosixPath(text)
    _require(
        path.is_absolute()
        and path.as_posix() == text
        and all(part not in ("", ".", "..") for part in path.parts),
        f"{label}: expected canonical absolute path",
    )
    return text


def _sha256(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _SHA256.fullmatch(value) is not None,
        f"{label}: malformed SHA-256",
    )
    return value


def _git_commit(value: object, *, label: str) -> str:
    _require(
        type(value) is str and _GIT_COMMIT.fullmatch(value) is not None,
        f"{label}: malformed Git commit",
    )
    return value


def _finite_float(value: object, *, label: str, minimum: float = 0.0) -> float:
    _require(
        type(value) is float and math.isfinite(value) and value >= minimum,
        f"{label}: expected finite float >= {minimum}",
    )
    return value


def _nonnegative_integer(value: object, *, label: str) -> int:
    _require(
        type(value) is int and value >= 0,
        f"{label}: expected nonnegative integer",
    )
    return value


def _positive_integer(value: object, *, label: str) -> int:
    _require(type(value) is int and value > 0, f"{label}: expected positive integer")
    return value


def _canonical_json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(
                value,
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise ProductionScienceAdmissionError(
            "record is not canonical finite JSON"
        ) from error


def canonical_sha256(value: object) -> str:
    """Return the deterministic SHA-256 of one finite JSON record."""
    return hashlib.sha256(_canonical_json_bytes(value)).hexdigest()


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


def _binding(
    value: object, *, label: str, absolute_path: bool = False
) -> dict[str, str]:
    record = _object(value, {"path", "sha256"}, label=label)
    path = (
        _absolute_path(record["path"], label=f"{label}/path")
        if absolute_path
        else _relative_path(record["path"], label=f"{label}/path")
    )
    return {"path": path, "sha256": _sha256(record["sha256"], label=f"{label}/sha256")}


def _read_verified_regular_file(
    binding: Mapping[str, str], *, label: str, require_executable: bool = False
) -> bytes:
    path = Path(binding["path"])
    descriptor = None
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
        before = os.fstat(descriptor)
        chunks = []
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
    except OSError as error:
        raise ProductionScienceAdmissionError(f"{label}: bound file is unavailable") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
    _require(
        stat.S_ISREG(before.st_mode)
        and not path.is_symlink()
        and (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns),
        f"{label}: bound file is not one stable regular file",
    )
    _require(
        (current.st_dev, current.st_ino, current.st_size, current.st_mtime_ns)
        == (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns),
        f"{label}: bound file pathname changed during verification",
    )
    _require(
        not require_executable or bool(before.st_mode & 0o111),
        f"{label}: bound executable lacks an execute bit",
    )
    _require(
        hashlib.sha256(payload).hexdigest() == binding["sha256"],
        f"{label}: bound file SHA-256 drifted",
    )
    return payload


def _json_payload(payload: bytes, *, label: str) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
        value = json.loads(
            text,
            object_pairs_hook=lambda pairs: _pairs_to_unique_object(pairs, label=label),
            parse_constant=lambda value: (_ for _ in ()).throw(
                ProductionScienceAdmissionError(
                    f"{label}: non-finite JSON constant {value}"
                )
            ),
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ProductionScienceAdmissionError(f"{label}: invalid UTF-8 JSON") from error
    _require(type(value) is dict, f"{label}: expected JSON object")
    return value


def _pairs_to_unique_object(
    pairs: Sequence[tuple[str, object]], *, label: str
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        _require(key not in result, f"{label}: duplicate JSON key {key!r}")
        result[key] = value
    return result


def _source_member_path(source_root: Path, relative: str, *, label: str) -> Path:
    _require(source_root.is_absolute(), f"{label}: source root must be absolute")
    try:
        resolved_root = source_root.resolve(strict=True)
        _require(resolved_root == source_root, f"{label}: root uses a symlink alias")
        path = (resolved_root / relative).resolve(strict=True)
        path.relative_to(resolved_root)
    except (OSError, ValueError) as error:
        raise ProductionScienceAdmissionError(
            f"{label}: source member escaped or is unavailable"
        ) from error
    _require(path == resolved_root / relative, f"{label}: source member uses a symlink alias")
    return path


def _verified_source_binding(
    value: object, *, source_root: Path, label: str
) -> tuple[dict[str, str], bytes]:
    binding = _binding(value, label=label)
    absolute = _source_member_path(source_root, binding["path"], label=label)
    absolute_binding = {"path": str(absolute), "sha256": binding["sha256"]}
    return binding, _read_verified_regular_file(absolute_binding, label=label)


def _closure(
    value: object,
    *,
    source_root: Path,
    label: str,
    required_paths: frozenset[str] = frozenset(),
) -> dict[str, object]:
    record = _object(value, {"members", "sha256"}, label=label)
    members = [
        _binding(item, label=f"{label}/members[{index}]")
        for index, item in enumerate(_array(record["members"], label=f"{label}/members"))
    ]
    _require(bool(members), f"{label}: closure must not be empty")
    paths = [item["path"] for item in members]
    _require(paths == sorted(paths), f"{label}: members must be ordered by path")
    _require(len(paths) == len(set(paths)), f"{label}: duplicate member path")
    missing = sorted(required_paths - set(paths))
    _require(not missing, f"{label}: required closure paths missing: {missing}")
    for member in members:
        _verified_source_binding(member, source_root=source_root, label=f"{label}/{member['path']}")
    digest = _sha256(record["sha256"], label=f"{label}/sha256")
    _require(
        digest == canonical_sha256(members),
        f"{label}: closure SHA-256 does not bind its members",
    )
    return {"members": members, "sha256": digest}


def _attempt_identity(value: object) -> dict[str, object]:
    record = _object(
        value,
        {"campaign_id", "attempt_id", "variant", "seed", "physical_mode"},
        label="attempt_identity",
    )
    _require(record["campaign_id"] == CAMPAIGN_ID, "attempt campaign ID drifted")
    attempt_id = _text(record["attempt_id"], label="attempt_identity/attempt_id")
    _require(
        _ATTEMPT_ID.fullmatch(attempt_id) is not None,
        "attempt_identity/attempt_id: noncanonical attempt ID",
    )
    variant = _text(record["variant"], label="attempt_identity/variant")
    _require(variant in GRID_VARIANTS, "attempt_identity/variant drifted")
    seed = _positive_integer(record["seed"], label="attempt_identity/seed")
    _require(seed in QUALIFYING_SEEDS, "attempt_identity/seed is not preregistered")
    _require(
        record["physical_mode"] == PHYSICAL_MODE,
        "attempt_identity/physical_mode drifted",
    )
    return {
        "campaign_id": CAMPAIGN_ID,
        "attempt_id": attempt_id,
        "variant": variant,
        "seed": seed,
        "physical_mode": PHYSICAL_MODE,
    }


def _candidate_binding(value: object, *, source_root: Path) -> dict[str, object]:
    record = _object(
        value,
        {
            "git_commit",
            "source_bundle_sha256",
            "clean_candidate_manifest",
            "executable",
            "environment_profile",
            "deck",
            "source_closure",
            "reducer_closure",
        },
        label="candidate_binding",
    )
    deck, deck_payload = _verified_source_binding(
        record["deck"], source_root=source_root, label="candidate_binding/deck"
    )
    _require(deck["path"] == SUCCESSOR_DECK, "candidate_binding/deck path drifted")
    git_commit = _git_commit(record["git_commit"], label="candidate_binding/git_commit")
    source_bundle_sha256 = _sha256(
        record["source_bundle_sha256"],
        label="candidate_binding/source_bundle_sha256",
    )
    clean_manifest = _binding(
        record["clean_candidate_manifest"],
        label="candidate_binding/clean_candidate_manifest",
        absolute_path=True,
    )
    executable = _binding(
        record["executable"], label="candidate_binding/executable", absolute_path=True
    )
    environment_record = _object(
        record["environment_profile"],
        {"path", "sha256", "control_plane_version", "reviewed_source"},
        label="candidate_binding/environment_profile",
    )
    reviewed_environment, _ = _verified_source_binding(
        environment_record["reviewed_source"],
        source_root=source_root,
        label="candidate_binding/environment_profile/reviewed_source",
    )
    environment = {
        **_binding(
            {"path": environment_record["path"], "sha256": environment_record["sha256"]},
            label="candidate_binding/environment_profile",
            absolute_path=True,
        ),
        "control_plane_version": _sha256(
            environment_record["control_plane_version"],
            label="candidate_binding/environment_profile/control_plane_version",
        ),
        "reviewed_source": reviewed_environment,
    }
    manifest_value = _json_payload(
        _read_verified_regular_file(clean_manifest, label="clean-candidate manifest"),
        label="clean-candidate manifest",
    )
    _require(
        manifest_value.get("schema_version") == 4
        and type(manifest_value.get("source")) is dict
        and type(manifest_value.get("build")) is dict,
        "clean-candidate manifest identity drifted",
    )
    _require(
        manifest_value["source"].get("git_commit") == git_commit
        and manifest_value["source"].get("source_bundle_sha256") == source_bundle_sha256
        and manifest_value["build"].get("executable_path") == executable["path"],
        "clean-candidate manifest candidate cross-link drifted",
    )
    _read_verified_regular_file(
        executable, label="candidate executable", require_executable=True
    )
    _read_verified_regular_file(environment, label="candidate environment profile")
    _require(
        environment["sha256"] == reviewed_environment["sha256"],
        "candidate environment profile differs from reviewed source bytes",
    )
    try:
        deck_contract = model.parse_deck_contract(deck_payload.decode("utf-8"))
        _require(
            deck_contract.variant == AMR_VARIANT,
            "candidate_binding/deck: successor base deck is not the AMR production layout",
        )
    except (UnicodeDecodeError, model.ModelContractError) as error:
        raise ProductionScienceAdmissionError(
            f"candidate_binding/deck: frozen Section 5.4 deck is invalid: {error}"
        ) from error
    return {
        "git_commit": git_commit,
        "source_bundle_sha256": source_bundle_sha256,
        "clean_candidate_manifest": clean_manifest,
        "executable": executable,
        "environment_profile": environment,
        "deck": deck,
        "source_closure": _closure(
            record["source_closure"],
            source_root=source_root,
            label="candidate_binding/source_closure",
            required_paths=REQUIRED_SOURCE_PATHS,
        ),
        "reducer_closure": _closure(
            record["reducer_closure"],
            source_root=source_root,
            label="candidate_binding/reducer_closure",
            required_paths=REQUIRED_REDUCER_PATHS,
        ),
    }


def _seed_overrides(seed: int) -> dict[str, int]:
    return {name: seed for name in SEED_OVERRIDE_NAMES}


def _expected_argv(
    identity: Mapping[str, object], raw_output_root: str
) -> list[str]:
    try:
        variant = model.variant_binding(identity["variant"])
    except model.ModelContractError as error:
        raise ProductionScienceAdmissionError(f"variant binding drifted: {error}") from error
    return [
        "-i",
        SUCCESSOR_DECK_LAUNCH_PATH,
        "-d",
        raw_output_root,
        f"job/basename={identity['attempt_id']}",
        "problem/ps_p0=1.0",
        *(f"{name}={identity['seed']}" for name in SEED_OVERRIDE_NAMES),
        *variant.model_launch_overrides,
    ]


def _document_binding(
    value: object, *, label: str
) -> tuple[dict[str, str], dict[str, Any]]:
    binding = _binding(value, label=label, absolute_path=True)
    return binding, _json_payload(_read_verified_regular_file(binding, label=label), label=label)


def _execution_binding(
    value: object,
    *,
    identity: Mapping[str, object],
    candidate: Mapping[str, object],
) -> dict[str, object]:
    record = _object(
        value,
        {
            "raw_output_root",
            "artifact_dir",
            "attempt_contract",
            "registered_execution_receipt",
            "selected_pressure_receipt",
        },
        label="execution_binding",
    )
    raw_root = _absolute_path(
        record["raw_output_root"], label="execution_binding/raw_output_root"
    )
    artifact_dir = _absolute_path(
        record["artifact_dir"], label="execution_binding/artifact_dir"
    )
    _require(
        raw_root != artifact_dir,
        "execution_binding: raw root and artifact directory coincide",
    )
    pressure_binding, pressure = _document_binding(
        record["selected_pressure_receipt"], label="selected pressure receipt"
    )
    _require(
        pressure.get("schema_version") == 3
        and pressure.get("record_type") == "q011_section54_pressure_selection_receipt"
        and pressure.get("selection_method") == "human_review_only"
        and type(pressure.get("selected_case")) is dict
        and pressure["selected_case"].get("case_id") == "ps_p0_1p00"
        and pressure["selected_case"].get("problem_ps_p0") == 1.0,
        "selected pressure receipt does not bind ps_p0_1p00 / problem_ps_p0=1.0",
    )
    contract_binding, contract = _document_binding(
        record["attempt_contract"], label="attempt contract"
    )
    receipt_binding, receipt = _document_binding(
        record["registered_execution_receipt"], label="registered execution receipt"
    )
    expected_argv = _expected_argv(identity, raw_root)
    expected_variant_overrides = list(
        model.variant_binding(identity["variant"]).model_launch_overrides
    )
    _require(
        contract.get("attempt_id") == identity["attempt_id"]
        and contract.get("variant") == identity["variant"]
        and contract.get("qualifying_seed") == identity["seed"]
        and contract.get("selected_problem_ps_p0") == 1.0
        and contract.get("launch_authorized") is False
        and contract.get("scheduler_submission_authorized") is False
        and contract.get("argv") == expected_argv
        and contract.get("executable") == candidate["executable"]
        and contract.get("environment_profile") == candidate["environment_profile"]
        and contract.get("authorized_orion_attempt_root") == str(Path(raw_root).parent)
        and Path(raw_root).name == "raw"
        and type(contract.get("paper_deck")) is dict
        and contract["paper_deck"].get("sha256") == candidate["deck"]["sha256"],
        "attempt contract candidate, pressure, variant, seed, or argv cross-link drifted",
    )
    _require(
        receipt.get("record_type") == "q011_section54_reconciled_registered_execution_receipt"
        and receipt.get("schema_version") == 1
        and receipt.get("receipt_role") == "immutable_reconciled_registered_execution"
        and receipt.get("registration_scope") == "registered_science"
        and receipt.get("reconciled") is True
        and type(receipt.get("reservation_id")) is str
        and _UUID.fullmatch(receipt["reservation_id"]) is not None
        and type(receipt.get("submission_id")) is str
        and _UUID.fullmatch(receipt["submission_id"]) is not None
        and type(receipt.get("reconciliation_event_sha256")) is str
        and _SHA256.fullmatch(receipt["reconciliation_event_sha256"]) is not None
        and receipt.get("attempt_id") == identity["attempt_id"]
        and receipt.get("source_commit") == candidate["git_commit"]
        and receipt.get("executable_sha256") == candidate["executable"]["sha256"]
        and receipt.get("deck_sha256") == candidate["deck"]["sha256"]
        and receipt.get("environment_sha256") == candidate["environment_profile"]["sha256"]
        and receipt.get("control_plane_version")
        == candidate["environment_profile"]["control_plane_version"]
        and receipt.get("argv") == expected_argv
        and receipt.get("slurm_terminal_state") == "COMPLETED"
        and type(receipt.get("slurm_job_id")) is str
        and re.fullmatch(r"[1-9][0-9]*", receipt["slurm_job_id"]) is not None
        and receipt.get("raw_output_root") == raw_root
        and receipt.get("artifact_dir") == artifact_dir
        and type(receipt.get("pre_submit_manifest_sha256")) is str
        and _SHA256.fullmatch(receipt["pre_submit_manifest_sha256"]) is not None,
        "registered execution receipt immutable cross-link drifted",
    )
    retention = receipt.get("planner_retention")
    _require(
        type(retention) is dict
        and retention.get("attempt_id") == identity["attempt_id"]
        and retention.get("authorized_orion_raw_root") == raw_root
        and retention.get("argv") == expected_argv,
        "registered execution receipt planner-retention cross-link drifted",
    )
    return {
        "raw_output_root": raw_root,
        "artifact_dir": artifact_dir,
        "attempt_contract": contract_binding,
        "registered_execution_receipt": receipt_binding,
        "selected_pressure_receipt": pressure_binding,
        "model_launch_overrides": expected_variant_overrides,
        "seed_overrides": _seed_overrides(identity["seed"]),
        "launch_argv": expected_argv,
        "registered_execution_identity": {
            "reservation_id": receipt["reservation_id"],
            "submission_id": receipt["submission_id"],
            "reconciliation_event_sha256": receipt["reconciliation_event_sha256"],
            "slurm_job_id": receipt["slurm_job_id"],
            "pre_submit_manifest_sha256": receipt["pre_submit_manifest_sha256"],
        },
    }


def _raw_artifact(
    value: object,
    *,
    identity: Mapping[str, object],
    raw_output_root: Path,
    index: int,
) -> tuple[dict[str, object], bytes]:
    label = f"raw_artifacts[{index}]"
    record = _object(
        value,
        {
            "path",
            "kind",
            "sha256",
            "byte_count",
            "attempt_id",
            "variant",
            "seed",
            "nominal_slot_time",
            "cycle",
            "observed_committed_time",
        },
        label=label,
    )
    path = _relative_path(record["path"], label=f"{label}/path")
    kind = _text(record["kind"], label=f"{label}/kind")
    _require(kind in ALL_RAW_PRODUCT_KINDS, f"{label}: unknown raw artifact kind")
    _require(
        _KIND_PATH_PATTERNS[kind].fullmatch(path) is not None,
        f"{label}: path does not match raw artifact kind {kind}",
    )
    _require(
        record["attempt_id"] == identity["attempt_id"]
        and record["variant"] == identity["variant"]
        and type(record["seed"]) is int
        and record["seed"] == identity["seed"],
        f"{label}: mixed attempt lineage",
    )
    artifact = {
        "path": path,
        "kind": kind,
        "sha256": _sha256(record["sha256"], label=f"{label}/sha256"),
        "byte_count": _positive_integer(
            record["byte_count"], label=f"{label}/byte_count"
        ),
        "attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "seed": identity["seed"],
        "nominal_slot_time": record["nominal_slot_time"],
        "cycle": record["cycle"],
        "observed_committed_time": record["observed_committed_time"],
    }
    if kind in RUN_PRODUCT_KINDS:
        _require(
            artifact["nominal_slot_time"] is None
            and artifact["cycle"] is None
            and artifact["observed_committed_time"] is None,
            f"{label}: run-level artifact must not claim snapshot lineage",
        )
    else:
        artifact["nominal_slot_time"] = _finite_float(
            artifact["nominal_slot_time"], label=f"{label}/nominal_slot_time"
        )
        artifact["cycle"] = _nonnegative_integer(
            artifact["cycle"], label=f"{label}/cycle"
        )
        artifact["observed_committed_time"] = _finite_float(
            artifact["observed_committed_time"],
            label=f"{label}/observed_committed_time",
        )
        _require(
            artifact["nominal_slot_time"] in REQUIRED_NOMINAL_SLOTS,
            f"{label}: unregistered nominal slot",
        )
    absolute = _source_member_path(raw_output_root, path, label=label)
    payload = _read_verified_regular_file(
        {"path": str(absolute), "sha256": artifact["sha256"]}, label=label
    )
    _require(len(payload) == artifact["byte_count"], f"{label}: byte count drifted")
    return artifact, payload


def _artifact_binding(artifact: Mapping[str, object]) -> dict[str, object]:
    return {
        "path": artifact["path"],
        "sha256": artifact["sha256"],
        "byte_count": artifact["byte_count"],
    }


def _fnv1a64(payload: bytes) -> int:
    value = 0xCBF29CE484222325
    for byte in payload:
        value ^= byte
        value = (value * 0x100000001B3) & 0xFFFFFFFFFFFFFFFF
    return value


def _validate_restart_marker(marker: bytes, payload: bytes, *, label: str) -> None:
    match = _RESTART_MARKER.fullmatch(marker)
    _require(match is not None, f"{label}: malformed restart completion marker")
    _require(
        int(match.group(1)) == len(payload)
        and int(match.group(2), 16) == _fnv1a64(payload),
        f"{label}: restart completion marker does not bind payload bytes",
    )


def _restart_cycle(payload: bytes, *, label: str) -> int:
    header, marker, body = payload.partition(b"<par_end>\n")
    _require(
        marker == b"<par_end>\n" and bool(header) and bool(body),
        f"{label}: invalid restart payload",
    )
    match = _RESTART_CYCLE.search(header)
    _require(match is not None, f"{label}: restart header lacks cycle")
    return int(match.group(1))


def _validate_restart_publication(
    artifacts: Sequence[Mapping[str, object]],
    payload_bytes: Mapping[str, bytes],
    *,
    slot: float,
    cycle: int,
) -> None:
    paths = {str(item["path"]) for item in artifacts}
    payload_paths = {path for path in paths if path.endswith(".rst")}
    payload_completions = {
        path.removesuffix(".complete")
        for path in paths
        if path.endswith(".rst.complete")
    }
    manifests = {path for path in paths if path.endswith(".rst.manifest")}
    manifest_completions = {
        path.removesuffix(".complete")
        for path in paths
        if path.endswith(".rst.manifest.complete")
    }
    _require(
        bool(payload_paths),
        f"nominal slot {slot}: expected at least one restart payload",
    )
    _require(
        payload_completions == payload_paths,
        f"nominal slot {slot}: restart payload completion inventory is incomplete",
    )
    _require(
        len(manifests) == 1 and manifest_completions == manifests,
        f"nominal slot {slot}: restart manifest publication is incomplete",
    )
    _require(
        paths
        == payload_paths | {f"{path}.complete" for path in payload_paths}
        | manifests
        | {f"{path}.complete" for path in manifests},
        f"nominal slot {slot}: unsupported restart publication member",
    )
    for path in sorted(payload_bytes):
        if path.endswith(".rst"):
            _require(
                _restart_cycle(payload_bytes[path], label=path) == cycle,
                f"nominal slot {slot}: restart header cycle drifted",
            )
            _validate_restart_marker(
                payload_bytes[f"{path}.complete"],
                payload_bytes[path],
                label=f"{path}.complete",
            )
    manifest_path = next(iter(manifests))
    manifest = _json_payload(payload_bytes[manifest_path], label=manifest_path)
    _require(
        set(manifest) == {"schema", "members"}
        and manifest["schema"] == "ATHENAK_RESTART_MANIFEST_V1"
        and type(manifest["members"]) is list
        and bool(manifest["members"]),
        f"nominal slot {slot}: invalid restart manifest",
    )
    parsed_members = {}
    for index, raw_member in enumerate(manifest["members"]):
        label = f"{manifest_path}/members[{index}]"
        member = _object(raw_member, {"path", "size", "fnv1a64"}, label=label)
        path = _relative_path(member["path"], label=f"{label}/path")
        _require(path not in parsed_members, f"{label}: duplicate restart member")
        size = _nonnegative_integer(member["size"], label=f"{label}/size")
        digest = _text(member["fnv1a64"], label=f"{label}/fnv1a64")
        _require(
            re.fullmatch(r"[0-9a-f]{16}", digest) is not None,
            f"{label}: invalid FNV-1a digest",
        )
        parsed_members[path] = (size, digest)
    _require(
        set(parsed_members) == payload_paths,
        f"nominal slot {slot}: restart manifest does not bind every payload",
    )
    for path, (size, digest) in parsed_members.items():
        _require(
            size == len(payload_bytes[path])
            and int(digest, 16) == _fnv1a64(payload_bytes[path]),
            f"nominal slot {slot}: restart manifest member digest drifted",
        )
    _validate_restart_marker(
        payload_bytes[f"{manifest_path}.complete"],
        payload_bytes[manifest_path],
        label=f"{manifest_path}.complete",
    )


def _float32(value: float, *, label: str) -> float:
    try:
        projected = struct.unpack("=f", struct.pack("=f", value))[0]
    except (OverflowError, struct.error) as error:
        raise ProductionScienceAdmissionError(
            f"{label}: cannot be represented as binary32"
        ) from error
    _require(math.isfinite(projected), f"{label}: non-finite binary32 projection")
    return projected


def _validate_observed_time_for_slot(slot: float, observed: float) -> None:
    index = REQUIRED_NOMINAL_SLOTS.index(slot)
    if index == 0:
        valid = observed == 0.0
    elif index == len(REQUIRED_NOMINAL_SLOTS) - 1:
        valid = observed == 1200.0
    else:
        next_slot = REQUIRED_NOMINAL_SLOTS[index + 1]
        valid = (
            _float32(observed, label=f"slot {slot} observed time")
            >= _float32(slot, label=f"slot {slot}")
            and _float32(observed, label=f"slot {slot} observed time")
            < _float32(next_slot, label=f"slot {next_slot}")
        )
    _require(valid, f"nominal slot {slot}: frozen float32 slot assignment failed")
    if slot == 500.0:
        _require(
            observed <= 500.1,
            "t=500 nominal slot exceeds the frozen maximum observed-time lateness",
        )


def _validate_stdout_terminal(payload: bytes, *, terminal_cycle: int) -> dict[str, object]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ProductionScienceAdmissionError("stdout: invalid UTF-8") from error
    termination = [line for line in text.splitlines() if line.startswith("Terminating on ")]
    states = list(_TERMINAL_STATE.finditer(text))
    limits = list(_TERMINAL_LIMIT.finditer(text))
    _require(
        termination == ["Terminating on time limit"]
        and len(states) == 1
        and len(limits) == 1,
        "stdout: terminal time-limit evidence is missing or ambiguous",
    )
    try:
        time = float(states[0].group(1))
        cycle = int(states[0].group(2))
        tlim = float(limits[0].group(1))
    except ValueError as error:
        raise ProductionScienceAdmissionError("stdout: terminal evidence is non-numeric") from error
    _require(
        math.isfinite(time)
        and math.isfinite(tlim)
        and time == 1200.0
        and tlim == 1200.0
        and cycle == terminal_cycle,
        "stdout: terminal time, tlim, or cycle differs from terminal snapshot",
    )
    return {
        "termination_reason": termination[0],
        "time_omega0_inverse": time,
        "cycle": cycle,
        "tlim_omega0_inverse": tlim,
    }


def _snapshot_bindings(
    artifacts: Sequence[Mapping[str, object]], payloads: Mapping[str, bytes]
) -> list[dict[str, object]]:
    snapshots = []
    previous_cycle: int | None = None
    previous_time: float | None = None
    for slot in REQUIRED_NOMINAL_SLOTS:
        selected = [item for item in artifacts if item["nominal_slot_time"] == slot]
        by_kind = {
            kind: [item for item in selected if item["kind"] == kind]
            for kind in ALL_RAW_PRODUCT_KINDS
        }
        for kind in SNAPSHOT_PRODUCT_KINDS:
            _require(
                len(by_kind[kind]) == 1,
                f"nominal slot {slot}: expected exactly one {kind} artifact",
            )
        _require(
            len(by_kind[RESTART_PRODUCT_KIND]) >= 1,
            f"nominal slot {slot}: expected at least one rst artifact",
        )
        _require(
            not by_kind["hst"],
            f"nominal slot {slot}: hst must remain run-level",
        )
        _require(
            not by_kind["stdout"],
            f"nominal slot {slot}: stdout must remain run-level",
        )
        cycles = {item["cycle"] for item in selected}
        times = {item["observed_committed_time"] for item in selected}
        _require(
            len(cycles) == 1 and len(times) == 1,
            f"nominal slot {slot}: raw artifacts do not share one exact cycle/time",
        )
        cycle = next(iter(cycles))
        observed_time = next(iter(times))
        _require(type(cycle) is int, f"nominal slot {slot}: invalid cycle")
        _require(
            type(observed_time) is float,
            f"nominal slot {slot}: invalid observed committed time",
        )
        _validate_restart_publication(
            by_kind[RESTART_PRODUCT_KIND],
            {
                str(item["path"]): payloads[str(item["path"])]
                for item in by_kind[RESTART_PRODUCT_KIND]
            },
            slot=slot,
            cycle=cycle,
        )
        _validate_observed_time_for_slot(slot, observed_time)
        if previous_cycle is not None and previous_time is not None:
            _require(
                cycle > previous_cycle and observed_time > previous_time,
                f"nominal slot {slot}: cycle/time sequence did not increase",
            )
        previous_cycle = cycle
        previous_time = observed_time
        snapshots.append(
            {
                "nominal_slot_time": slot,
                "cycle": cycle,
                "observed_committed_time": observed_time,
                "artifact_bindings": {
                    **{
                        kind: _artifact_binding(by_kind[kind][0])
                        for kind in SNAPSHOT_PRODUCT_KINDS
                    },
                    "rst": [
                        _artifact_binding(item)
                        for item in sorted(
                            by_kind[RESTART_PRODUCT_KIND],
                            key=lambda item: str(item["path"]),
                        )
                    ],
                },
            }
        )
    return snapshots


def analysis_policy() -> dict[str, object]:
    """Return the fixed successor policy without granting policy authority."""
    return {
        "shock_front_preregistration_supersession": {
            "predecessor_path": PREDECESSOR_PREREGISTRATION,
            "superseded_erroneous_detector": "unique_strongest_positive_density_gradient",
            "successor_detector": "unique_strongest_negative_density_gradient",
            "source_field": "dens_from_mhd_w_bcc",
            "profile": "y_area_weighted_mean",
            "search_window_relative_to_ideal_surface_c_over_omega_pi": [
                -1200.0,
                1200.0,
            ],
            "maximum_absolute_offset_from_ideal_surface_c_over_omega_pi": 600.0,
            "orientation_rationale": (
                "The reflecting wall and high-density downstream are at lower x1; "
                "density therefore decreases across the outward-moving shock."
            ),
            "promotion_status": "defined_not_registered_not_authoritative",
        },
        "magnetic_amplification": {
            "future_qualifying_metric_after_separate_policy_promotion": {
                "reference_surface": "ideal_injection_surface",
                "window_relative_to_ideal_surface_c_over_omega_pi": [120.0, 1200.0],
                "snapshot_nominal_slot_omega0_inverse": 500.0,
                "source": "mhd_w_bcc_derived_abs_b",
                "observable": "area_weighted_mean_abs_b_over_b0",
                "reference_b0": 1.0,
            },
            "supplemental_nonqualifying_metric": {
                "reference_surface": "detected_negative_gradient_front",
                "window_relative_to_detected_front_c_over_omega_pi": [120.0, 1200.0],
                "snapshot_nominal_slot_omega0_inverse": 500.0,
                "source": "mhd_w_bcc_derived_abs_b",
                "observables": [
                    "area_weighted_mean_abs_b_over_b0",
                    "area_weighted_rms_abs_b_over_b0",
                    "local_max_abs_b_over_b0",
                    "area_weighted_q50_q90_q99_abs_b_over_b0",
                    "area_fraction_abs_b_over_b0_ge_2_3_4",
                ],
                "qualification_use": "forbidden",
            },
        },
        "downstream_spectra": {
            "nominal_slots_omega0_inverse": [500.0, 1200.0],
            "raw_source": "prtcl_all",
            "downstream_classifier": "x1_less_than_ideal_injection_surface",
            "particle_filter": "cr_source==1_and_birth_time>=45",
            "reducer": "tst/publication/q011_section54_particles.py:reduce_particle_snapshot",
            "late_tail_slope_slot_omega0_inverse": 1200.0,
            "qualification_status": "path_defined_not_authorized",
        },
        "current_profiles": {
            "raw_sources": [
                "mhd_w_bcc",
                "prtcl_rho",
                "prtcl_jx",
                "prtcl_jy",
                "prtcl_jz",
                "mhd_j2",
            ],
            "deposited_representation": "J_CR_over_c",
            "required_frames": [
                "lab_frame_deposited_J_CR_over_c",
                (
                    "gas_frame_deposited_J_CR_over_c_equals_lab_frame_deposited_"
                    "J_CR_over_c_minus_prtcl_rho_times_gas_velocity"
                ),
            ],
            "mhd_j2_label": "mhd_current_magnitude_squared_not_particle_current",
            "qualification_status": "path_defined_not_authorized",
        },
        "acceleration_histories": {
            "raw_source": "prtcl_all_all_nominal_slots",
            "primary_tail_observables": ["q990_chi", "q999_chi"],
            "supplemental_observable": "maximum_chi",
            "minimum_fit_slots": 4,
            "float32_projection_uncertainty_required": True,
            "qualification_status": "path_defined_not_authorized",
        },
        "energy_partitions": {
            "raw_sources": ["hst", "mhd_w_bcc", "prtcl_all", "rst"],
            "instantaneous_partition_status": "diagnostic_not_conservation_closure",
            "qualification_status": "path_defined_not_authorized",
        },
        "float32_uncertainty": {
            "raw_source": "prtcl_all",
            "required_observables": [
                "particle_chi_interval",
                "q990_chi_interval",
                "q999_chi_interval",
                "late_tail_slope_interval",
            ],
            "qualification_status": "path_defined_not_authorized",
        },
        "morphology": {
            "nominal_slot_omega0_inverse": 500.0,
            "raw_sources": ["mhd_w_bcc", "prtcl_all"],
            "required_products": [
                "energy_weighted_cr_spatial_distribution",
                "refinement_overlay",
                "publication_bound_particle_current_b_field_profiles",
            ],
            "qualification_status": "hook_defined_not_authorized",
        },
        "conserved_energy_closure": {
            "raw_sources": ["hst", "rst", "prtcl_all"],
            "required_ledger_terms": [
                "initial_mhd_plus_cr_state",
                "boundary_transport",
                "injected_cr_ledger",
                "removed_cr_ledger",
                "gas_subtraction_transaction",
                "current_mhd_plus_cr_state",
            ],
            "required_normalizations": [
                "initial_total_energy",
                "cumulative_absolute_energy_exchange",
                "current_total_energy",
            ],
            "qualification_status": "hook_defined_not_authorized",
        },
        "paired_grid_amr_residuals": {
            "pairing_rule": "same_seed_three_level_amr_with_fine_uniform",
            "common_grid_rule": (
                "conservatively_restrict_spatial_fields_to_dx12_then_y_area_average"
            ),
            "spectrum_rule": (
                "compare_same_fixed_chi_bins_after_identical_ideal_surface_downstream_filter"
            ),
            "required_nominal_slots_omega0_inverse": [500.0, 1200.0],
            "required_observables": [
                "shock_front_position_at_t500",
                "ideal_surface_relative_upstream_magnetic_amplification_at_t500",
                "rho_y_average_at_t500",
                "bmag_y_average_at_t500",
                "normalized_downstream_chi_f_chi_at_t500",
                "normalized_downstream_chi_f_chi_at_t1200",
            ],
            "qualification_status": "path_defined_not_authorized",
        },
    }


def _authorization_boundary() -> dict[str, object]:
    return {
        "launch_authorized": False,
        "scheduler_submission_authorized": False,
        "policy_mutation_authorized": False,
        "qualifying_output_inspection_authorized": False,
        "qualification_authorized": False,
        "claim_closure_authorized": False,
        "required_next_boundary": (
            "separate_reviewed_registration_and_policy_promotion_of_this_exact_successor"
        ),
    }


def build_attempt_admission(
    *,
    attempt_identity: Mapping[str, object],
    candidate_binding: Mapping[str, object],
    execution_binding: Mapping[str, object],
    raw_artifacts: Sequence[Mapping[str, object]],
    source_root: str | Path,
) -> dict[str, object]:
    """Verify and build one non-authorizing expanded-raw-set attempt admission."""
    identity = _attempt_identity(attempt_identity)
    source = Path(source_root)
    candidate = _candidate_binding(candidate_binding, source_root=source)
    execution = _execution_binding(
        execution_binding, identity=identity, candidate=candidate
    )
    raw_root = Path(execution["raw_output_root"])
    _require(raw_root.is_absolute(), "raw output root must be absolute")
    _require(type(raw_artifacts) is list, "raw_artifacts: expected array")
    verified = [
        _raw_artifact(
            item, identity=identity, raw_output_root=raw_root, index=index
        )
        for index, item in enumerate(raw_artifacts)
    ]
    artifacts = [item[0] for item in verified]
    payloads = {str(item[0]["path"]): item[1] for item in verified}
    artifacts.sort(key=lambda item: str(item["path"]))
    paths = [str(item["path"]) for item in artifacts]
    _require(len(paths) == len(set(paths)), "raw_artifacts: duplicate path")
    hst = [item for item in artifacts if item["kind"] == "hst"]
    stdout = [item for item in artifacts if item["kind"] == "stdout"]
    _require(len(hst) == 1, "raw_artifacts: expected exactly one hst artifact")
    _require(len(stdout) == 1, "raw_artifacts: expected exactly one stdout artifact")
    snapshots = _snapshot_bindings(artifacts, payloads)
    terminal_completion = _validate_stdout_terminal(
        payloads[str(stdout[0]["path"])], terminal_cycle=snapshots[-1]["cycle"]
    )
    observed_kinds = {str(item["kind"]) for item in artifacts}
    _require(
        observed_kinds == ALL_RAW_PRODUCT_KINDS,
        "raw_artifacts: expanded raw product set is incomplete",
    )
    policy = analysis_policy()
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": ATTEMPT_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": CAMPAIGN_ID,
        "claim_id": CLAIM_ID,
        "status": "admitted_to_non_authorizing_production_science_orchestration_only",
        "qualification_effect": (
            "raw_lineage_and_orchestration_admission_only_no_launch_no_qualification"
        ),
        "attempt_identity": identity,
        "candidate_binding": candidate,
        "candidate_binding_sha256": canonical_sha256(candidate),
        "execution_binding": execution,
        "execution_binding_sha256": canonical_sha256(execution),
        "raw_inventory_sha256": canonical_sha256(artifacts),
        "raw_artifacts": artifacts,
        "run_bindings": {
            "hst": _artifact_binding(hst[0]),
            "stdout": _artifact_binding(stdout[0]),
            "terminal_completion": terminal_completion,
        },
        "snapshot_bindings": snapshots,
        "analysis_policy": policy,
        "analysis_policy_sha256": canonical_sha256(policy),
        "authorization": _authorization_boundary(),
    }


def validate_attempt_admission(
    value: object, *, source_root: str | Path
) -> dict[str, object]:
    """Rebuild and exactly validate one attempt admission record."""
    record = _object(
        value,
        {
            "schema_version",
            "record_type",
            "successor_id",
            "campaign_id",
            "claim_id",
            "status",
            "qualification_effect",
            "attempt_identity",
            "candidate_binding",
            "candidate_binding_sha256",
            "execution_binding",
            "execution_binding_sha256",
            "raw_inventory_sha256",
            "raw_artifacts",
            "run_bindings",
            "snapshot_bindings",
            "analysis_policy",
            "analysis_policy_sha256",
            "authorization",
        },
        label="attempt_admission",
    )
    rebuilt = build_attempt_admission(
        attempt_identity=record["attempt_identity"],
        candidate_binding=record["candidate_binding"],
        execution_binding={
            key: record["execution_binding"][key]
            for key in (
                "raw_output_root",
                "artifact_dir",
                "attempt_contract",
                "registered_execution_receipt",
                "selected_pressure_receipt",
            )
        },
        raw_artifacts=record["raw_artifacts"],
        source_root=source_root,
    )
    _require(
        _strict_equal(record, rebuilt),
        "attempt_admission: derived fields or authorization boundary drifted",
    )
    return rebuilt


def _snapshot_for_slot(
    admission: Mapping[str, object], slot: float
) -> Mapping[str, object]:
    matches = [
        item
        for item in admission["snapshot_bindings"]
        if item["nominal_slot_time"] == slot
    ]
    _require(len(matches) == 1, f"attempt admission lacks exact slot {slot}")
    return matches[0]


def _orchestration_binding(admission: Mapping[str, object]) -> dict[str, object]:
    identity = admission["attempt_identity"]
    return {
        "attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "seed": identity["seed"],
        "attempt_admission_sha256": canonical_sha256(admission),
        "raw_inventory_sha256": admission["raw_inventory_sha256"],
        "candidate_binding_sha256": admission["candidate_binding_sha256"],
        "execution_binding_sha256": admission["execution_binding_sha256"],
        "analysis_policy_sha256": admission["analysis_policy_sha256"],
    }


def _raw_binding_for_slot(
    admission: Mapping[str, object], slot: float, kind: str
) -> dict[str, object]:
    snapshot = _snapshot_for_slot(admission, slot)
    return {
        "cycle": snapshot["cycle"],
        "observed_committed_time": snapshot["observed_committed_time"],
        "artifact": snapshot["artifact_bindings"][kind],
    }


def _resource_phase_binding(phase_id: str, value: object) -> dict[str, object]:
    _require(phase_id in CAMPAIGN_PHASES, "campaign phase is not preregistered")
    selected_seeds = CAMPAIGN_PHASES[phase_id]
    if selected_seeds == CORE_QUALIFYING_SEEDS:
        _require(value is None, "core phase must not claim an expansion decision")
        return {
            "phase_id": phase_id,
            "selected_seeds": list(selected_seeds),
            "decision": "preregistered_core_three_paired_triads",
            "expansion_evidence": None,
        }
    record = _object(
        value,
        {
            "record_type",
            "selected_phase_id",
            "previous_phase_id",
            "decision_basis_fields",
            "qualifying_output_inspected",
            "resource_model_sha256",
        },
        label="resource_expansion_evidence",
    )
    count = len(selected_seeds)
    previous_phase = f"phase_{count - 1}_paired_triads"
    _require(
        record["record_type"] == "q011_resource_telemetry_expansion_decision_v1"
        and record["selected_phase_id"] == phase_id
        and record["previous_phase_id"] == previous_phase
        and record["qualifying_output_inspected"] is False
        and record["decision_basis_fields"] == list(RESOURCE_TELEMETRY_FIELDS),
        "resource expansion is not the exact telemetry-only preregistered phase transition",
    )
    evidence = {
        **record,
        "resource_model_sha256": _sha256(
            record["resource_model_sha256"],
            label="resource_expansion_evidence/resource_model_sha256",
        ),
    }
    return {
        "phase_id": phase_id,
        "selected_seeds": list(selected_seeds),
        "decision": "cumulative_next_paired_triad_from_resource_telemetry_only",
        "expansion_evidence": evidence,
    }


def _attempt_job_identity(admission: Mapping[str, object]) -> dict[str, object]:
    identity = admission["attempt_identity"]
    return {
        "attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "seed": identity["seed"],
        "attempt_admission_sha256": canonical_sha256(admission),
    }


def _all_slot_bindings(
    admission: Mapping[str, object], kinds: Sequence[str]
) -> list[dict[str, object]]:
    return [
        {
            "nominal_slot_time": slot,
            "bindings": {
                kind: _raw_binding_for_slot(admission, slot, kind) for kind in kinds
            },
        }
        for slot in REQUIRED_NOMINAL_SLOTS
    ]


def _complete_work_graph(
    ordered: Sequence[Mapping[str, object]],
    by_cell: Mapping[tuple[str, int], Mapping[str, object]],
    selected_seeds: Sequence[int],
) -> dict[str, object]:
    graph: dict[str, object] = {
        "shock_front_jobs": [],
        "magnetic_amplification_jobs": [],
        "current_profile_jobs": [],
        "downstream_spectrum_jobs": [],
        "acceleration_history_jobs": [],
        "energy_partition_jobs": [],
        "float32_uncertainty_jobs": [],
        "morphology_jobs": [],
        "conserved_energy_closure_jobs": [],
        "paired_grid_amr_residual_jobs": [],
    }
    for admission in ordered:
        job_identity = _attempt_job_identity(admission)
        graph["shock_front_jobs"].append(
            {
                **job_identity,
                "nominal_slot_time": 500.0,
                "mhd_w_bcc_binding": _raw_binding_for_slot(admission, 500.0, "mhd_w_bcc"),
                "detector": "unique_strongest_negative_density_gradient",
            }
        )
        graph["magnetic_amplification_jobs"].append(
            {
                **job_identity,
                "nominal_slot_time": 500.0,
                "mhd_w_bcc_binding": _raw_binding_for_slot(admission, 500.0, "mhd_w_bcc"),
                "qualifying_reference_surface": "ideal_injection_surface",
                "detected_front_relative_role": "supplemental_nonqualifying",
            }
        )
        for slot in SCIENCE_SLOTS:
            graph["current_profile_jobs"].append(
                {
                    **job_identity,
                    "nominal_slot_time": slot,
                    "raw_bindings": {
                        kind: _raw_binding_for_slot(admission, slot, kind)
                        for kind in (
                            "mhd_w_bcc",
                            "prtcl_rho",
                            "prtcl_jx",
                            "prtcl_jy",
                            "prtcl_jz",
                            "mhd_j2",
                        )
                    },
                    "deposited_current_representation": "J_CR_over_c",
                    "frames": [
                        "lab_frame_deposited_J_CR_over_c",
                        "gas_frame_deposited_J_CR_over_c",
                    ],
                }
            )
            graph["downstream_spectrum_jobs"].append(
                {
                    **job_identity,
                    "nominal_slot_time": slot,
                    "prtcl_all_binding": _raw_binding_for_slot(admission, slot, "prtcl_all"),
                    "downstream_classifier": "ideal_injection_surface",
                    "evaluate_late_slope": slot == 1200.0,
                }
            )
            graph["float32_uncertainty_jobs"].append(
                {
                    **job_identity,
                    "nominal_slot_time": slot,
                    "prtcl_all_binding": _raw_binding_for_slot(admission, slot, "prtcl_all"),
                    "required_intervals": [
                        "particle_chi",
                        "q990_chi",
                        "q999_chi",
                        "late_tail_slope_if_t1200",
                    ],
                }
            )
        graph["acceleration_history_jobs"].append(
            {
                **job_identity,
                "slot_bindings": _all_slot_bindings(admission, ("prtcl_all",)),
                "primary_observables": ["q990_chi", "q999_chi"],
                "supplemental_observable": "maximum_chi",
                "minimum_fit_slots": 4,
            }
        )
        graph["energy_partition_jobs"].append(
            {
                **job_identity,
                "hst_binding": admission["run_bindings"]["hst"],
                "slot_bindings": _all_slot_bindings(
                    admission, ("mhd_w_bcc", "prtcl_all")
                ),
                "role": "instantaneous_partition_diagnostic_not_conservation_closure",
            }
        )
        graph["morphology_jobs"].append(
            {
                **job_identity,
                "nominal_slot_time": 500.0,
                "raw_bindings": {
                    kind: _raw_binding_for_slot(admission, 500.0, kind)
                    for kind in (
                        "mhd_w_bcc",
                        "prtcl_all",
                        "prtcl_rho",
                        "prtcl_jx",
                        "prtcl_jy",
                        "prtcl_jz",
                    )
                },
                "required_products": [
                    "energy_weighted_cr_spatial_distribution",
                    "refinement_overlay",
                    "publication_bound_profiles",
                ],
            }
        )
        graph["conserved_energy_closure_jobs"].append(
            {
                **job_identity,
                "hst_binding": admission["run_bindings"]["hst"],
                "stdout_binding": admission["run_bindings"]["stdout"],
                "slot_bindings": _all_slot_bindings(admission, ("prtcl_all", "rst")),
                "required_terms": analysis_policy()["conserved_energy_closure"][
                    "required_ledger_terms"
                ],
                "status": "hook_requires_exact_boundary_transport_instrumentation",
            }
        )
    for seed in selected_seeds:
        amr = by_cell[(AMR_VARIANT, seed)]
        fine = by_cell[(FINE_VARIANT, seed)]
        graph["paired_grid_amr_residual_jobs"].append(
            {
                "seed": seed,
                "amr_attempt_admission_sha256": canonical_sha256(amr),
                "fine_attempt_admission_sha256": canonical_sha256(fine),
                "slot_bindings": [
                    {
                        "nominal_slot_time": slot,
                        "amr": {
                            kind: _raw_binding_for_slot(amr, slot, kind)
                            for kind in ("mhd_w_bcc", "prtcl_all")
                        },
                        "fine": {
                            kind: _raw_binding_for_slot(fine, slot, kind)
                            for kind in ("mhd_w_bcc", "prtcl_all")
                        },
                    }
                    for slot in SCIENCE_SLOTS
                ],
                "common_grid_rule": (
                    "conservatively_restrict_spatial_fields_to_dx12_then_y_area_average"
                ),
            }
        )
    return graph


def build_campaign_orchestration(
    attempt_admissions: Sequence[Mapping[str, object]],
    *,
    source_root: str | Path,
    phase_id: str = "phase_3_paired_triads",
    resource_expansion_evidence: object = None,
) -> dict[str, object]:
    """Build one cumulative, telemetry-only-phase, non-authorizing work graph."""
    _require(type(attempt_admissions) is list, "attempt_admissions: expected array")
    phase = _resource_phase_binding(phase_id, resource_expansion_evidence)
    selected_seeds = tuple(phase["selected_seeds"])
    expected_cells = tuple(
        (variant, seed) for variant in GRID_VARIANTS for seed in selected_seeds
    )
    _require(
        len(attempt_admissions) == len(expected_cells),
        f"attempt_admissions: expected complete {len(expected_cells)}-attempt phase matrix",
    )
    validated = [
        validate_attempt_admission(item, source_root=source_root)
        for item in attempt_admissions
    ]
    by_cell: dict[tuple[str, int], dict[str, object]] = {}
    attempt_ids = set()
    admission_sha256_values = set()
    candidate_sha256_values = set()
    policy_sha256_values = set()
    raw_roots = set()
    artifact_dirs = set()
    receipt_sha256_values = set()
    contract_sha256_values = set()
    reservation_ids = set()
    submission_ids = set()
    reconciliation_event_sha256_values = set()
    slurm_job_ids = set()
    pre_submit_manifest_sha256_values = set()
    for admission in validated:
        identity = admission["attempt_identity"]
        cell = (identity["variant"], identity["seed"])
        _require(cell not in by_cell, f"attempt_admissions: duplicate matrix cell {cell}")
        digest = canonical_sha256(admission)
        execution = admission["execution_binding"]
        registered_identity = execution["registered_execution_identity"]
        _require(
            identity["attempt_id"] not in attempt_ids
            and digest not in admission_sha256_values
            and execution["raw_output_root"] not in raw_roots
            and execution["artifact_dir"] not in artifact_dirs
            and execution["registered_execution_receipt"]["sha256"]
            not in receipt_sha256_values
            and execution["attempt_contract"]["sha256"] not in contract_sha256_values
            and registered_identity["reservation_id"] not in reservation_ids
            and registered_identity["submission_id"] not in submission_ids
            and registered_identity["reconciliation_event_sha256"]
            not in reconciliation_event_sha256_values
            and registered_identity["slurm_job_id"] not in slurm_job_ids
            and registered_identity["pre_submit_manifest_sha256"]
            not in pre_submit_manifest_sha256_values,
            "attempt_admissions: relabeled or reused execution evidence",
        )
        by_cell[cell] = admission
        attempt_ids.add(identity["attempt_id"])
        admission_sha256_values.add(digest)
        raw_roots.add(execution["raw_output_root"])
        artifact_dirs.add(execution["artifact_dir"])
        receipt_sha256_values.add(execution["registered_execution_receipt"]["sha256"])
        contract_sha256_values.add(execution["attempt_contract"]["sha256"])
        reservation_ids.add(registered_identity["reservation_id"])
        submission_ids.add(registered_identity["submission_id"])
        reconciliation_event_sha256_values.add(
            registered_identity["reconciliation_event_sha256"]
        )
        slurm_job_ids.add(registered_identity["slurm_job_id"])
        pre_submit_manifest_sha256_values.add(
            registered_identity["pre_submit_manifest_sha256"]
        )
        candidate_sha256_values.add(admission["candidate_binding_sha256"])
        policy_sha256_values.add(admission["analysis_policy_sha256"])
    _require(
        set(by_cell) == set(expected_cells),
        "attempt_admissions: phase matrix is incomplete or drifted",
    )
    _require(
        len(candidate_sha256_values) == 1,
        "attempt_admissions: mixed candidate/deck/source/reducer lineage",
    )
    _require(
        len(policy_sha256_values) == 1,
        "attempt_admissions: mixed analysis policy lineage",
    )
    ordered = [by_cell[cell] for cell in expected_cells]
    policy = analysis_policy()
    graph = _complete_work_graph(ordered, by_cell, selected_seeds)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": CAMPAIGN_RECORD_TYPE,
        "successor_id": SUCCESSOR_ID,
        "campaign_id": CAMPAIGN_ID,
        "claim_id": CLAIM_ID,
        "status": "complete_non_authorizing_production_science_work_graph",
        "qualification_effect": "orchestration_only_no_launch_no_policy_promotion_no_qualification",
        "campaign_phase": phase,
        "campaign_matrix": {
            "physical_mode": PHYSICAL_MODE,
            "grid_variants": list(GRID_VARIANTS),
            "qualifying_seeds": list(selected_seeds),
            "expected_attempt_count": len(expected_cells),
            "paired_seed_rule": "same_seed_three_level_amr_with_fine_uniform",
        },
        "attempt_admission_bindings": [_orchestration_binding(admission) for admission in ordered],
        "common_candidate_binding_sha256": next(iter(candidate_sha256_values)),
        "analysis_policy": policy,
        "analysis_policy_sha256": canonical_sha256(policy),
        "work_graph": graph,
        "work_graph_sha256": canonical_sha256(graph),
        "authorization": _authorization_boundary(),
    }


def validate_campaign_orchestration(
    value: object,
    attempt_admissions: Sequence[Mapping[str, object]],
    *,
    source_root: str | Path,
) -> dict[str, object]:
    """Rebuild and exactly validate one persisted campaign work graph."""
    _require(type(value) is dict, "campaign orchestration: expected object")
    phase = value.get("campaign_phase")
    _require(type(phase) is dict, "campaign orchestration: missing campaign phase")
    rebuilt = build_campaign_orchestration(
        attempt_admissions,
        source_root=source_root,
        phase_id=phase.get("phase_id"),
        resource_expansion_evidence=phase.get("expansion_evidence"),
    )
    _require(
        _strict_equal(value, rebuilt),
        "campaign orchestration: persisted work graph or authorization boundary drifted",
    )
    return rebuilt


__all__ = [
    "ALL_RAW_PRODUCT_KINDS",
    "AMR_VARIANT",
    "ATTEMPT_RECORD_TYPE",
    "CAMPAIGN_ID",
    "CAMPAIGN_PHASES",
    "CAMPAIGN_RECORD_TYPE",
    "CORE_QUALIFYING_SEEDS",
    "EXPECTED_MATRIX_CELLS",
    "FINE_VARIANT",
    "GRID_VARIANTS",
    "ProductionScienceAdmissionError",
    "QUALIFYING_SEEDS",
    "REQUIRED_NOMINAL_SLOTS",
    "REQUIRED_REDUCER_PATHS",
    "REQUIRED_SOURCE_PATHS",
    "RESOURCE_TELEMETRY_FIELDS",
    "SCIENCE_SLOTS",
    "SNAPSHOT_PRODUCT_KINDS",
    "SUCCESSOR_DECK",
    "SUCCESSOR_DECK_LAUNCH_PATH",
    "SUCCESSOR_ID",
    "analysis_policy",
    "build_attempt_admission",
    "build_campaign_orchestration",
    "canonical_sha256",
    "validate_attempt_admission",
    "validate_campaign_orchestration",
]
