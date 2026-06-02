#!/usr/bin/env python3
"""Fail-closed admission gate for one Q-011 Section 5.4 campaign run.

This campaign-bound admission slice validates immutable provenance and retained
raw-artifact structure.  It intentionally does not claim that the later
numerical or qualitative Section 5.4 gates have passed.
"""

from __future__ import annotations

import argparse
from datetime import datetime
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import re
from typing import Any, Mapping, Sequence
import uuid

import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from .frontier_control_plane.control_plane_common import (
        BUILD_PROVENANCE_FILENAMES,
    )
    from .frontier_control_plane.control_plane_common import (
        read_json_bytes as read_control_plane_json_bytes,
    )
    from .frontier_control_plane.control_plane_common import (
        read_stable_regular_file_below,
    )
    from .frontier_control_plane.control_plane_common import (
        validate_clean_candidate_bundle,
    )
    from .immutable_orion_tree import authorized_tree_root
    from .immutable_orion_tree import loads_json_reject_duplicate_keys
    from .immutable_orion_tree import require_exact_primitive_types
    from .immutable_orion_tree import staged_verified_frozen_tree
    from .immutable_orion_tree import validate_executable_elf
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    import analyze_q011_section54_outputs as output_primitives
    from frontier_control_plane.control_plane_common import BUILD_PROVENANCE_FILENAMES
    from frontier_control_plane.control_plane_common import (
        read_json_bytes as read_control_plane_json_bytes,
    )
    from frontier_control_plane.control_plane_common import read_stable_regular_file_below
    from frontier_control_plane.control_plane_common import validate_clean_candidate_bundle
    from immutable_orion_tree import authorized_tree_root
    from immutable_orion_tree import loads_json_reject_duplicate_keys
    from immutable_orion_tree import require_exact_primitive_types
    from immutable_orion_tree import staged_verified_frozen_tree
    from immutable_orion_tree import validate_executable_elf
    from pvtk_particles import ParticleVTKData, read_particle_vtk


REPO_ROOT = Path(__file__).resolve().parents[2]
ANALYZER_PATH = Path(__file__).resolve()
PREREGISTRATION_PATH = (
    REPO_ROOT
    / "tst/publication/readiness/"
    "q011_section54_qualifying_campaign_preregistration_successor_v2_2026-06-01.json"
)
EXPECTED_PREREGISTRATION_SHA256 = (
    "6fd9ebbc247b6cace69f0ff61553cf198241577b457410d57d26afcbf27cdc35"
)
ORION_BULK_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
ACTIVE_DECK_SOURCE_PATH = (
    "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
CAMPAIGN_ANALYZER_SOURCE_PATH = (
    "tst/publication/analyze_q011_section54_campaign.py"
)
PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH = (
    "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json"
)
CAMPAIGN_ID = "Q011-SECTION54-QUALIFYING-CAMPAIGN"
MANIFEST_NAME = "campaign_manifest.json"
INVENTORY_NAME = "artifact_inventory.sha256"
FREEZE_RECEIPT_NAME = "freeze_receipt.json"
ARTIFACT_ROLE = "q011_section54_qualifying_campaign_attempt"
QUALIFICATION_SCOPE = "artifact_admission_only_no_final_claim_closure"
RESULT_RECORD_TYPE = "q011_section54_campaign_admission_result"
_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
_COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_PATTERN = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_RFC3339_UTC_PATTERN = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"(?:\.[0-9]{1,6})?Z"
)
_ENDPOINT_PRODUCT_KINDS = ("rho", "bmag", "prtcl_jx", "j2", "prtcl_all")
_RESTART_PRODUCT_KINDS = (
    "restart",
    "restart_manifest",
    "restart_complete",
    "restart_manifest_complete",
)
_RUN_PRODUCT_KINDS = ("stdout",)
_SNAPSHOT_PRODUCT_KINDS = frozenset((*_ENDPOINT_PRODUCT_KINDS, *_RESTART_PRODUCT_KINDS))
_ALL_PRODUCT_KINDS = frozenset((*_ENDPOINT_PRODUCT_KINDS, *_RUN_PRODUCT_KINDS))
_ALL_PRODUCT_KINDS = frozenset((*_ALL_PRODUCT_KINDS, *_RESTART_PRODUCT_KINDS))
_BINDING_NAMES = (
    "clean_candidate_manifest",
    "executable",
    "deck",
    "analyzer",
    "preregistration",
)
_PVTK_SCALARS = frozenset(
    {
        "gid",
        "ptag",
        "species",
        "cr_source",
        "macro_weight",
        "birth_time",
        "deltaf_f0",
        "deltaf_weight",
    }
)
_PVTK_INTEGER_SCALARS = frozenset({"gid", "ptag", "species", "cr_source"})
_PVTK_REAL_SCALARS = _PVTK_SCALARS - _PVTK_INTEGER_SCALARS
_EXPECTED_BIN_FIELDS = {
    "rho": ("dens",),
    "bmag": ("bmag",),
    "prtcl_jx": ("prtcl_jx",),
    "j2": ("j2",),
}
_PVTK_EXECUTION_PATTERN = re.compile(
    rb"^# vtk DataFile Version 2\.0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  "
    rb"nranks= (0|[1-9][0-9]*)  cycle=(0|[1-9][0-9]*)  variables=([^\n]+)\n"
)
_RESTART_MARKER_PATTERN = re.compile(
    rb"ATHENAK_RESTART_COMPLETE_V1\n"
    rb"size=(0|[1-9][0-9]*)\n"
    rb"fnv1a64=([0-9a-f]{16})\n"
)
_Q017_TELEMETRY_PATTERN = re.compile(
    r"^q017\.telemetry\.([A-Za-z0-9_.]+)=([^\s]+)$", re.MULTILINE
)
_Q017_REQUIRED_NAMES = frozenset(
    {
        "schema_version",
        "mpi.ranks",
        "timer.driver.seconds_rank_max",
        "timer.task_lists.seconds_rank_max",
        "timer.task_lists.seconds_rank_mean",
        "timer.task_lists.calls_rank_max",
        "timer.task_list.before_timeintegrator.seconds_rank_max",
        "timer.task_list.before_stagen.seconds_rank_max",
        "timer.task_list.stagen.seconds_rank_max",
        "timer.task_list.after_stagen.seconds_rank_max",
        "timer.task_list.after_timeintegrator.seconds_rank_max",
        "timer.output_publication.seconds_rank_max",
        "timer.output_publication.calls_rank_max",
        "timer.amr_load_balance.seconds_rank_max",
        "timer.amr_load_balance.calls_rank_max",
        "timer.particle.adaptive_deltaf.seconds_rank_max",
        "timer.particle.adaptive_deltaf.calls_rank_max",
        "timer.particle.push.seconds_rank_max",
        "timer.particle.push.calls_rank_max",
        "timer.particle.deposition.seconds_rank_max",
        "timer.particle.deposition.calls_rank_max",
        "timer.particle.migration.seconds_rank_max",
        "timer.particle.migration.calls_rank_max",
        "cycles",
        "meshblocks.total",
        "meshblocks.rank_min",
        "meshblocks.rank_max",
        "mesh.cells_per_meshblock",
        "mesh.active_cells",
        "particles.total",
        "particles.rank_min",
        "particles.rank_max",
        "updates.meshblock_cycles",
        "updates.particle_updates",
        "throughput.zone_cycles_per_second",
        "throughput.particle_updates_per_second",
        "load.meshblock_efficiency",
        "load.particle_efficiency",
        "load.cost.total",
        "load.cost.rank_min",
        "load.cost.rank_max",
        "load.cost.efficiency",
        "load.cost.invalid_meshblocks",
        "amr.enabled",
        "amr.meshblocks_created",
        "amr.meshblocks_deleted",
        "amr.meshblocks_communicated",
        "particle_memory.sync_kernel_timers",
        "particle_memory.record_bytes",
        "particle_memory.root_level",
        "particle_memory.max_level",
        "particle_memory.resident_records.bytes_total",
        "particle_memory.direct_views.allocated_snapshot_bytes_total",
        "particle_memory.direct_views.allocated_snapshot_bytes_rank_max",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_snapshot_bytes_total",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_snapshot_bytes_rank_max",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_high_water_bytes_rank_sum",
        "particle_memory.athenak_owned_tracked_kokkos_views.allocated_high_water_bytes_rank_max",
        "particle_memory.paper_smooth_host_transport.allocated_snapshot_bytes_total",
        "particle_memory.paper_smooth_host_transport.allocated_snapshot_bytes_rank_max",
        "particle_memory.paper_smooth_host_transport.allocated_high_water_bytes_rank_sum",
        "particle_memory.paper_smooth_host_transport.allocated_high_water_bytes_rank_max",
        "particle_memory.invalid_records",
        "particle_memory.species.0.count",
        "particle_memory.species.0.resident_bytes",
        "particle_memory.species.1.count",
        "particle_memory.species.1.resident_bytes",
        "particle_memory.level.1.count",
        "particle_memory.level.1.resident_bytes",
        "particle_memory.level.2.count",
        "particle_memory.level.2.resident_bytes",
        "particle_memory.species.0.level.1.count",
        "particle_memory.species.0.level.2.count",
        "particle_memory.species.1.level.1.count",
        "particle_memory.species.1.level.2.count",
    }
)
_EXPECTED_FREEZE_RECEIPT = {
    "schema_version": 1,
    "artifact_role": ARTIFACT_ROLE,
    "qualification_effect": "retained_qualifying_campaign_attempt",
    "inventory_excludes": INVENTORY_NAME,
    "freeze_policy": "remove all owner, group and other write bits recursively",
}
_EXPECTED_POLICY_PROJECTION = {
    "campaign_matrix": {
        "physical_mode": "paper_mhd_pic_vl2_tsc",
        "grid_variants": [
            "coarse_uniform_dx12",
            "three_level_amr_root_dx12_finest_dx3",
            "fine_uniform_dx3",
        ],
        "qualifying_seeds": [
            23050101,
            23050102,
            23050103,
            23050104,
            23050105,
            23050106,
            23050107,
            23050108,
        ],
        "expected_baseline_attempts": 24,
        "paired_seed_rule": (
            "Use the same qualifying seed for coarse-uniform, AMR and "
            "fine-uniform variants."
        ),
    },
    "snapshot_selection": {
        "time_unit": "omega0_inverse",
        "required_times": [500.0, 1200.0],
        "absolute_match_tolerance": 1.0e-06,
        "missing_or_ambiguous_snapshot": "fail_endpoint",
    },
    "raw_bin_ids": ["rho", "bmag", "prtcl_jx", "j2"],
    "raw_pvtk_ids": ["prtcl_all"],
    "required_output_times_omega0_inverse": [
        0.0,
        100.0,
        200.0,
        300.0,
        400.0,
        500.0,
        600.0,
        700.0,
        800.0,
        900.0,
        1000.0,
        1100.0,
        1200.0,
    ],
    "raw_rst_policy": "retain_each_emitted_restart_checkpoint_and_completion_metadata",
}


class QualificationError(ValueError):
    """Raised when one campaign attempt cannot enter numerical qualification."""

    def __init__(self, message: str, *, code: str = "invalid_campaign_tree"):
        super().__init__(message)
        self.code = code


def _fail(code: str, message: str) -> None:
    raise QualificationError(message, code=code)


def _require(condition: bool, code: str, message: str) -> None:
    if not condition:
        _fail(code, message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _sha256_path(path: Path) -> str:
    try:
        return _sha256_bytes(path.read_bytes())
    except OSError as error:
        _fail("unreadable_bound_file", f"unable to read bound file {path}: {error}")


def _object(value: object, expected: set[str], label: str) -> dict[str, Any]:
    _require(type(value) is dict, "schema_type_error", f"{label}: expected object")
    mapping = value
    _require(set(mapping) == expected, "schema_key_error", f"{label}: keys drifted")
    return mapping


def _list(value: object, label: str) -> list[Any]:
    _require(type(value) is list, "schema_type_error", f"{label}: expected list")
    return value


def _text(value: object, label: str) -> str:
    _require(
        type(value) is str and bool(value),
        "schema_type_error",
        f"{label}: expected text",
    )
    return value


def _exact_int(value: object, label: str) -> int:
    _require(type(value) is int, "schema_type_error", f"{label}: expected integer")
    return value


def _finite_float(value: object, label: str) -> float:
    _require(type(value) is float, "schema_type_error", f"{label}: expected float")
    _require(math.isfinite(value), "schema_type_error", f"{label}: expected finite float")
    return value


def _sha256(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _SHA256_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected 64 lowercase hexadecimal digits",
    )
    return text


def _relative_path(value: object, label: str) -> str:
    text = _text(value, label)
    candidate = PurePosixPath(text)
    _require(
        not candidate.is_absolute()
        and text == candidate.as_posix()
        and text != "."
        and not any(character.isspace() for character in text)
        and all(part not in ("", ".", "..") for part in candidate.parts),
        "unsafe_relative_path",
        f"{label}: unsafe root-relative path {text!r}",
    )
    return text


def _git_sha1(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _COMMIT_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected 40 lowercase hexadecimal digits",
    )
    return text


def _uuid(value: object, label: str) -> str:
    text = _text(value, label)
    try:
        parsed = uuid.UUID(text)
    except ValueError:
        _fail("schema_type_error", f"{label}: expected canonical UUID")
    _require(
        str(parsed) == text, "schema_type_error", f"{label}: expected canonical UUID"
    )
    return text


def _utc_timestamp(value: object, label: str) -> str:
    text = _text(value, label)
    _require(
        _RFC3339_UTC_PATTERN.fullmatch(text) is not None,
        "schema_type_error",
        f"{label}: expected canonical RFC-3339 UTC timestamp",
    )
    try:
        datetime.fromisoformat(text.removesuffix("Z") + "+00:00")
    except ValueError:
        _fail("schema_type_error", f"{label}: expected valid UTC timestamp")
    return text


def _frozen_candidate_path(value: object, freeze_id: str, suffix: str, label: str) -> str:
    text = _text(value, label)
    expected = (ORION_BULK_ROOT / "clean_candidates" / freeze_id / suffix).as_posix()
    _require(
        text == expected,
        "clean_candidate_schema_error",
        f"{label}: expected frozen clean-candidate path {expected!r}",
    )
    return text


def _load_json_bytes(payload: bytes, label: str) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_json", f"{label}: JSON is not UTF-8: {error}")
    decoded = loads_json_reject_duplicate_keys(
        text,
        error_type=QualificationError,
        label=label,
    )
    _require(type(decoded) is dict, "schema_type_error", f"{label}: expected object")
    return decoded


def _exact_match(actual: object, expected: object, label: str) -> None:
    require_exact_primitive_types(
        actual,
        expected,
        error_type=QualificationError,
        label=label,
    )
    _require(actual == expected, "policy_drift", f"{label}: value drifted")


def _policy_projection(policy: Mapping[str, Any]) -> dict[str, Any]:
    criteria = _object(
        policy.get("athenak_selected_release_criteria"),
        {
            "criteria_provenance",
            "campaign_matrix",
            "snapshot_selection",
            "ideal_injection_surface_classifier",
            "particle_filter",
            "spectrum",
            "shock_front",
            "upstream_magnetic_amplification",
            "amr_vs_fine_uniform_residuals",
            "qualitative_figure_requirements",
        },
        "preregistration/athenak_selected_release_criteria",
    )
    retention = _object(
        policy.get("artifact_retention_policy"),
        {
            "retention_scope",
            "required_output_times_omega0_inverse",
            "raw_bin_ids",
            "raw_pvtk_ids",
            "raw_rst_policy",
            "inventory_policy",
            "derived_artifact_policy",
        },
        "preregistration/artifact_retention_policy",
    )
    return {
        "campaign_matrix": criteria["campaign_matrix"],
        "snapshot_selection": criteria["snapshot_selection"],
        "raw_bin_ids": retention["raw_bin_ids"],
        "raw_pvtk_ids": retention["raw_pvtk_ids"],
        "required_output_times_omega0_inverse": retention[
            "required_output_times_omega0_inverse"
        ],
        "raw_rst_policy": retention["raw_rst_policy"],
    }


def _load_bound_policy(snapshot: Any, binding: Mapping[str, str]) -> dict[str, Any]:
    retained_path = snapshot.member_path(binding["path"])
    retained_payload = retained_path.read_bytes()
    _require(
        _sha256_bytes(retained_payload) == binding["sha256"],
        "hash_drift",
        "preregistration SHA-256 drifted",
    )
    _require(
        binding["sha256"] == EXPECTED_PREREGISTRATION_SHA256,
        "hash_drift",
        "retained preregistration does not match the frozen complete-policy SHA-256",
    )
    _require(
        _sha256_path(PREREGISTRATION_PATH) == binding["sha256"],
        "hash_drift",
        "invoked preregistration source SHA-256 drifted",
    )
    _require(
        retained_payload == PREREGISTRATION_PATH.read_bytes(),
        "hash_drift",
        "retained preregistration bytes differ from invoked frozen policy",
    )
    policy = _load_json_bytes(retained_payload, "preregistration")
    _require(
        policy.get("record_type") == "q011_section54_qualifying_campaign_preregistration"
        and type(policy.get("schema_version")) is int
        and policy["schema_version"] == 1,
        "policy_drift",
        "preregistration identity drifted",
    )
    projection = _policy_projection(policy)
    _exact_match(projection, _EXPECTED_POLICY_PROJECTION, "preregistration projection")
    return policy


def _validate_binding(value: object, label: str) -> dict[str, str]:
    binding = _object(value, {"path", "sha256"}, label)
    return {
        "path": _relative_path(binding["path"], f"{label}/path"),
        "sha256": _sha256(binding["sha256"], f"{label}/sha256"),
    }


def _validate_manifest_schema(manifest: object, root: Path) -> dict[str, Any]:
    item = _object(
        manifest,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "qualification_scope",
            "authorized_orion_campaign_root",
            "run_identity",
            "candidate_binding",
            "artifact_bindings",
            "products",
        },
        "campaign manifest",
    )
    _require(
        _exact_int(item["schema_version"], "campaign manifest/schema_version") == 1,
        "schema_value_error",
        "campaign manifest/schema_version: expected 1",
    )
    _require(
        _text(item["record_type"], "campaign manifest/record_type")
        == "q011_section54_campaign_run_manifest",
        "schema_value_error",
        "campaign manifest/record_type drifted",
    )
    _require(
        _text(item["campaign_id"], "campaign manifest/campaign_id") == CAMPAIGN_ID,
        "schema_value_error",
        "campaign manifest/campaign_id drifted",
    )
    _require(
        _text(item["qualification_scope"], "campaign manifest/qualification_scope")
        == QUALIFICATION_SCOPE,
        "schema_value_error",
        "campaign manifest/qualification_scope drifted",
    )
    _require(
        _text(
            item["authorized_orion_campaign_root"],
            "campaign manifest/authorized_orion_campaign_root",
        )
        == str(root),
        "orion_root_binding_drift",
        "campaign manifest authorized Orion root binding drifted",
    )

    identity = _object(
        item["run_identity"],
        {"variant", "seed", "physical_mode", "attempt_id"},
        "campaign manifest/run_identity",
    )
    attempt_id = _text(
        identity["attempt_id"], "campaign manifest/run_identity/attempt_id"
    )
    _require(
        _ATTEMPT_PATTERN.fullmatch(attempt_id) is not None,
        "schema_value_error",
        "campaign manifest/run_identity/attempt_id is noncanonical",
    )
    parsed_identity = {
        "variant": _text(identity["variant"], "campaign manifest/run_identity/variant"),
        "seed": _exact_int(identity["seed"], "campaign manifest/run_identity/seed"),
        "physical_mode": _text(
            identity["physical_mode"], "campaign manifest/run_identity/physical_mode"
        ),
        "attempt_id": attempt_id,
    }

    candidate = _object(
        item["candidate_binding"],
        {"git_commit", "source_bundle_sha256", "clean_candidate_manifest"},
        "campaign manifest/candidate_binding",
    )
    commit = _text(
        candidate["git_commit"], "campaign manifest/candidate_binding/git_commit"
    )
    _require(
        _COMMIT_PATTERN.fullmatch(commit) is not None,
        "schema_type_error",
        "campaign manifest/candidate_binding/git_commit must be 40 lowercase hex digits",
    )
    parsed_candidate = {
        "git_commit": commit,
        "source_bundle_sha256": _sha256(
            candidate["source_bundle_sha256"],
            "campaign manifest/candidate_binding/source_bundle_sha256",
        ),
        "clean_candidate_manifest": _validate_binding(
            candidate["clean_candidate_manifest"],
            "campaign manifest/candidate_binding/clean_candidate_manifest",
        ),
    }

    raw_bindings = _object(
        item["artifact_bindings"],
        set(_BINDING_NAMES[1:]),
        "campaign manifest/artifact_bindings",
    )
    bindings = {
        name: _validate_binding(
            raw_bindings[name], f"campaign manifest/artifact_bindings/{name}"
        )
        for name in _BINDING_NAMES[1:]
    }
    products = [
        _validate_product(product, index)
        for index, product in enumerate(
            _list(item["products"], "campaign manifest/products")
        )
    ]
    _require(bool(products), "missing_product", "campaign manifest/products is empty")
    return {
        "run_identity": parsed_identity,
        "candidate_binding": parsed_candidate,
        "artifact_bindings": bindings,
        "products": products,
    }


def _validate_product(value: object, index: int) -> dict[str, Any]:
    label = f"campaign manifest/products[{index}]"
    product = _object(value, {"kind", "path", "sha256", "snapshot_time"}, label)
    kind = _text(product["kind"], f"{label}/kind")
    _require(
        kind in _ALL_PRODUCT_KINDS,
        "unknown_product",
        f"{label}: unknown kind {kind!r}",
    )
    snapshot_time = product["snapshot_time"]
    if kind in _SNAPSHOT_PRODUCT_KINDS:
        snapshot_time = _finite_float(snapshot_time, f"{label}/snapshot_time")
    else:
        _require(
            snapshot_time is None,
            "schema_type_error",
            f"{label}/snapshot_time: run-level products require null",
        )
    return {
        "kind": kind,
        "path": _relative_path(product["path"], f"{label}/path"),
        "sha256": _sha256(product["sha256"], f"{label}/sha256"),
        "snapshot_time": snapshot_time,
    }


def _validate_prepared_records(value: object, label: str) -> list[dict[str, str]]:
    raw_records = _list(value, label)
    _require(
        1 <= len(raw_records) <= 4096,
        "clean_candidate_schema_error",
        f"{label}: expected between 1 and 4096 records",
    )
    records = [
        _validate_binding(record, f"{label}[{index}]")
        for index, record in enumerate(raw_records)
    ]
    paths = [record["path"] for record in records]
    _require(
        paths == sorted(set(paths)),
        "clean_candidate_schema_error",
        f"{label}: paths must be unique and canonically ordered",
    )
    return records


def _validate_clean_candidate_submodules(
    value: object, freeze_id: str
) -> list[dict[str, str]]:
    raw_records = _list(value, "clean candidate/source/submodules")
    records: list[dict[str, str]] = []
    for index, raw_record in enumerate(raw_records):
        label = f"clean candidate/source/submodules[{index}]"
        record = _object(
            raw_record,
            {
                "path",
                "archive_path",
                "archive_sha256",
                "commit_path",
                "commit_sha256",
                "git_commit",
                "git_tree",
                "worktree_status",
            },
            label,
        )
        _require(
            _text(record["worktree_status"], f"{label}/worktree_status") == "clean",
            "clean_candidate_schema_error",
            f"{label}/worktree_status: expected 'clean'",
        )
        records.append(
            {
                "path": _relative_path(record["path"], f"{label}/path"),
                "archive_path": _frozen_candidate_path(
                    record["archive_path"],
                    freeze_id,
                    f"submodules/{index:04d}.tar",
                    f"{label}/archive_path",
                ),
                "archive_sha256": _sha256(
                    record["archive_sha256"], f"{label}/archive_sha256"
                ),
                "commit_path": _frozen_candidate_path(
                    record["commit_path"],
                    freeze_id,
                    f"submodules/{index:04d}.commit",
                    f"{label}/commit_path",
                ),
                "commit_sha256": _sha256(
                    record["commit_sha256"], f"{label}/commit_sha256"
                ),
                "git_commit": _git_sha1(record["git_commit"], f"{label}/git_commit"),
                "git_tree": _git_sha1(record["git_tree"], f"{label}/git_tree"),
            }
        )
    paths = [record["path"] for record in records]
    _require(
        paths == sorted(set(paths)),
        "clean_candidate_schema_error",
        "clean candidate/source/submodules: paths must be unique and canonically ordered",
    )
    return records


def _validate_clean_candidate_manifest(payload: bytes) -> dict[str, Any]:
    label = "clean candidate"
    item = _object(
        _load_json_bytes(payload, label),
        {
            "schema_version",
            "freeze_id",
            "created_utc",
            "prepared_artifacts",
            "source",
            "build",
        },
        label,
    )
    _require(
        _exact_int(item["schema_version"], f"{label}/schema_version") == 4,
        "clean_candidate_schema_error",
        f"{label}/schema_version: expected 4",
    )
    freeze_id = _uuid(item["freeze_id"], f"{label}/freeze_id")
    created_utc = _utc_timestamp(item["created_utc"], f"{label}/created_utc")

    prepared = _object(
        item["prepared_artifacts"],
        {"inventory_path", "inventory_sha256", "paper_decks", "analyzers"},
        f"{label}/prepared_artifacts",
    )
    inventory_path = _relative_path(
        prepared["inventory_path"], f"{label}/prepared_artifacts/inventory_path"
    )
    _require(
        inventory_path == PREPARED_ARTIFACT_INVENTORY_SOURCE_PATH,
        "clean_candidate_schema_error",
        "clean candidate prepared-artifact inventory path drifted",
    )
    paper_decks = _validate_prepared_records(
        prepared["paper_decks"], f"{label}/prepared_artifacts/paper_decks"
    )
    analyzers = _validate_prepared_records(
        prepared["analyzers"], f"{label}/prepared_artifacts/analyzers"
    )
    prepared_paths = [record["path"] for record in (*paper_decks, *analyzers)]
    _require(
        inventory_path not in prepared_paths
        and len(prepared_paths) == len(set(prepared_paths)),
        "clean_candidate_schema_error",
        "clean candidate prepared-artifact paths and inventory must be distinct",
    )

    source = _object(
        item["source"],
        {
            "archive_path",
            "archive_sha256",
            "commit_path",
            "commit_sha256",
            "source_bundle_sha256",
            "git_commit",
            "git_tree",
            "worktree_status",
            "submodule_status",
            "submodules",
        },
        f"{label}/source",
    )
    _require(
        _text(source["worktree_status"], f"{label}/source/worktree_status") == "clean",
        "clean_candidate_schema_error",
        f"{label}/source/worktree_status: expected 'clean'",
    )
    submodule_status = _text(
        source["submodule_status"], f"{label}/source/submodule_status"
    )
    submodules = _validate_clean_candidate_submodules(source["submodules"], freeze_id)
    expected_submodule_status = "clean_pinned_archived" if submodules else "absent"
    _require(
        submodule_status == expected_submodule_status,
        "clean_candidate_schema_error",
        f"{label}/source/submodule_status: does not match retained submodules",
    )
    parsed_source = {
        "archive_path": _frozen_candidate_path(
            source["archive_path"],
            freeze_id,
            "source.tar",
            f"{label}/source/archive_path",
        ),
        "archive_sha256": _sha256(
            source["archive_sha256"], f"{label}/source/archive_sha256"
        ),
        "commit_path": _frozen_candidate_path(
            source["commit_path"],
            freeze_id,
            "source.commit",
            f"{label}/source/commit_path",
        ),
        "commit_sha256": _sha256(
            source["commit_sha256"], f"{label}/source/commit_sha256"
        ),
        "source_bundle_sha256": _sha256(
            source["source_bundle_sha256"], f"{label}/source/source_bundle_sha256"
        ),
        "git_commit": _git_sha1(source["git_commit"], f"{label}/source/git_commit"),
        "git_tree": _git_sha1(source["git_tree"], f"{label}/source/git_tree"),
        "submodules": submodules,
    }

    build = _object(
        item["build"],
        {
            "profile_id",
            "profile_path",
            "profile_sha256",
            "profile_receipt_path",
            "profile_receipt_sha256",
            "source_archive_sha256",
            "source_commit_sha256",
            "source_bundle_sha256",
            "toolchain",
            "build_invocations_sha256",
            "executable_path",
            "executable_sha256",
        },
        f"{label}/build",
    )
    profile_id = _text(build["profile_id"], f"{label}/build/profile_id")
    _require(
        profile_id == "hip-mpi-release-paper-pic",
        "clean_candidate_schema_error",
        f"{label}/build/profile_id: expected production PIC build profile",
    )
    toolchain = _text(build["toolchain"], f"{label}/build/toolchain")
    _require(
        bool(toolchain.strip()),
        "clean_candidate_schema_error",
        f"{label}/build/toolchain: expected nonblank text",
    )
    parsed_build = {
        "profile_id": profile_id,
        "profile_path": _frozen_candidate_path(
            build["profile_path"],
            freeze_id,
            "build_profile.json",
            f"{label}/build/profile_path",
        ),
        "profile_sha256": _sha256(
            build["profile_sha256"], f"{label}/build/profile_sha256"
        ),
        "profile_receipt_path": _frozen_candidate_path(
            build["profile_receipt_path"],
            freeze_id,
            "profile_receipt.json",
            f"{label}/build/profile_receipt_path",
        ),
        "profile_receipt_sha256": _sha256(
            build["profile_receipt_sha256"], f"{label}/build/profile_receipt_sha256"
        ),
        "source_archive_sha256": _sha256(
            build["source_archive_sha256"], f"{label}/build/source_archive_sha256"
        ),
        "source_commit_sha256": _sha256(
            build["source_commit_sha256"], f"{label}/build/source_commit_sha256"
        ),
        "source_bundle_sha256": _sha256(
            build["source_bundle_sha256"], f"{label}/build/source_bundle_sha256"
        ),
        "toolchain": toolchain,
        "build_invocations_sha256": _sha256(
            build["build_invocations_sha256"],
            f"{label}/build/build_invocations_sha256",
        ),
        "executable_path": _frozen_candidate_path(
            build["executable_path"],
            freeze_id,
            "athena",
            f"{label}/build/executable_path",
        ),
        "executable_sha256": _sha256(
            build["executable_sha256"], f"{label}/build/executable_sha256"
        ),
    }
    _require(
        parsed_build["source_archive_sha256"] == parsed_source["archive_sha256"]
        and parsed_build["source_commit_sha256"] == parsed_source["commit_sha256"]
        and parsed_build["source_bundle_sha256"]
        == parsed_source["source_bundle_sha256"],
        "clean_candidate_binding_drift",
        "clean candidate build is not bound to its frozen source",
    )
    return {
        "schema_version": 4,
        "freeze_id": freeze_id,
        "created_utc": created_utc,
        "prepared_artifacts": {
            "inventory_path": inventory_path,
            "inventory_sha256": _sha256(
                prepared["inventory_sha256"],
                f"{label}/prepared_artifacts/inventory_sha256",
            ),
            "paper_decks": paper_decks,
            "analyzers": analyzers,
        },
        "source": parsed_source,
        "build": parsed_build,
    }


def _validate_identity(identity: Mapping[str, Any], policy: Mapping[str, Any]) -> None:
    matrix = policy["athenak_selected_release_criteria"]["campaign_matrix"]
    _require(
        identity["variant"] in matrix["grid_variants"],
        "invalid_run_identity",
        "run identity variant is not preregistered",
    )
    _require(
        identity["seed"] in matrix["qualifying_seeds"],
        "invalid_run_identity",
        "run identity seed is not preregistered",
    )
    _require(
        identity["physical_mode"] == matrix["physical_mode"],
        "invalid_run_identity",
        "run identity physical mode is not preregistered",
    )


def _member_payload(snapshot: Any, binding: Mapping[str, str], label: str) -> bytes:
    try:
        path = snapshot.member_path(binding["path"])
        payload = path.read_bytes()
    except (OSError, ValueError) as error:
        _fail("missing_bound_file", f"{label}: retained member is unavailable: {error}")
    measured = _sha256_bytes(payload)
    _require(measured == binding["sha256"], "hash_drift", f"{label}: SHA-256 drifted")
    _require(bool(payload), "empty_bound_file", f"{label}: retained member is empty")
    return payload


def _member_sha256(snapshot: Any, binding: Mapping[str, str], label: str) -> str:
    return _sha256_bytes(_member_payload(snapshot, binding, label))


def _validate_prepared_binding(
    records: Sequence[Mapping[str, str]],
    required_path: str,
    binding: Mapping[str, str],
    label: str,
) -> None:
    matches = [record for record in records if record["path"] == required_path]
    _require(
        len(matches) == 1,
        "clean_candidate_binding_drift",
        f"{label}: frozen prepared-artifact path is missing or ambiguous",
    )
    _require(
        matches[0]["sha256"] == binding["sha256"],
        "clean_candidate_binding_drift",
        f"{label}: SHA-256 differs from frozen prepared artifact",
    )


def _cross_bind_clean_candidate(
    parsed: Mapping[str, Any], frozen: Mapping[str, Any]
) -> None:
    candidate = parsed["candidate_binding"]
    bindings = parsed["artifact_bindings"]
    source = frozen["source"]
    build = frozen["build"]
    prepared = frozen["prepared_artifacts"]
    _require(
        candidate["git_commit"] == source["git_commit"],
        "clean_candidate_binding_drift",
        "campaign Git commit differs from frozen clean candidate",
    )
    _require(
        candidate["source_bundle_sha256"] == source["source_bundle_sha256"]
        and candidate["source_bundle_sha256"] == build["source_bundle_sha256"],
        "clean_candidate_binding_drift",
        "campaign source bundle differs from frozen clean candidate",
    )
    _require(
        bindings["executable"]["sha256"] == build["executable_sha256"],
        "clean_candidate_binding_drift",
        "campaign executable differs from frozen clean candidate",
    )
    _validate_prepared_binding(
        prepared["paper_decks"],
        ACTIVE_DECK_SOURCE_PATH,
        bindings["deck"],
        "campaign active VL2 deck",
    )
    _validate_prepared_binding(
        prepared["analyzers"],
        CAMPAIGN_ANALYZER_SOURCE_PATH,
        bindings["analyzer"],
        "campaign admission analyzer",
    )


def _validate_elf_identity(
    path: Path, payload: bytes, expected_sha256: str, label: str
) -> dict[str, str | bool]:
    try:
        report = validate_executable_elf(
            path,
            expected_sha256,
            error_type=QualificationError,
            label=label,
        )
    except ValueError as error:
        _fail("invalid_executable_elf", f"{label}: {error}")
    _require(
        len(payload) >= 64
        and payload[:4] == b"\x7fELF"
        and payload[4] in (1, 2)
        and payload[5] in (1, 2)
        and payload[6] == 1,
        "invalid_executable_elf",
        f"{label}: expected retained ELF executable identity",
    )
    return report


def _read_external_candidate_member(path: str, label: str) -> bytes:
    try:
        return read_stable_regular_file_below(
            Path(path),
            ORION_BULK_ROOT,
            require_read_only_mode=True,
        )
    except (OSError, ValueError) as error:
        _fail(
            "external_clean_candidate_closure_error",
            f"{label}: unable to read frozen clean-candidate member: {error}",
        )


def _validate_external_clean_candidate_closure(
    frozen: Mapping[str, Any],
    retained_manifest_payload: bytes,
    executable_sha256: str,
) -> dict[str, Any]:
    """Revalidate the referenced external source/build closure from frozen bytes."""
    candidate_root = ORION_BULK_ROOT / "clean_candidates" / frozen["freeze_id"]
    external_manifest_payload = _read_external_candidate_member(
        str(candidate_root / "clean_candidate_manifest.json"),
        "external clean candidate manifest",
    )
    _require(
        external_manifest_payload == retained_manifest_payload,
        "external_clean_candidate_closure_error",
        "retained clean-candidate manifest differs from referenced external freeze",
    )
    try:
        candidate = read_control_plane_json_bytes(
            external_manifest_payload,
            label="external clean candidate manifest",
        )
        source = frozen["source"]
        build = frozen["build"]
        build_profile = _read_external_candidate_member(
            build["profile_path"], "external clean-candidate build profile"
        )
        profile = read_control_plane_json_bytes(
            build_profile,
            label="external clean-candidate build profile",
        )
        provenance_inputs = profile.get("provenance_inputs")
        _require(
            type(provenance_inputs) is dict
            and set(provenance_inputs) == set(BUILD_PROVENANCE_FILENAMES),
            "external_clean_candidate_closure_error",
            "external clean-candidate build provenance input set is incomplete",
        )
        build_provenance: dict[str, bytes] = {}
        for name in BUILD_PROVENANCE_FILENAMES:
            record = provenance_inputs[name]
            _require(
                type(record) is dict
                and set(record) == {"path", "sha256"}
                and type(record["path"]) is str,
                "external_clean_candidate_closure_error",
                f"external clean-candidate build provenance record is malformed: {name}",
            )
            build_provenance[name] = _read_external_candidate_member(
                record["path"], f"external clean-candidate build provenance {name}"
            )
        external_executable = _read_external_candidate_member(
            build["executable_path"], "external clean-candidate executable"
        )
        executable_report = _validate_elf_identity(
            Path(build["executable_path"]),
            external_executable,
            executable_sha256,
            "external clean-candidate executable",
        )
        _require(
            _sha256_bytes(external_executable) == executable_sha256,
            "external_clean_candidate_closure_error",
            "external clean-candidate executable differs from retained executable",
        )
        submodules = source["submodules"]
        validated_submodules = validate_clean_candidate_bundle(
            candidate,
            source_archive=_read_external_candidate_member(
                source["archive_path"], "external clean-candidate source archive"
            ),
            source_commit=_read_external_candidate_member(
                source["commit_path"], "external clean-candidate source commit"
            ),
            submodule_archives=[
                _read_external_candidate_member(
                    record["archive_path"],
                    f"external clean-candidate submodule archive {index}",
                )
                for index, record in enumerate(submodules)
            ],
            submodule_commits=[
                _read_external_candidate_member(
                    record["commit_path"],
                    f"external clean-candidate submodule commit {index}",
                )
                for index, record in enumerate(submodules)
            ],
            build_profile=build_profile,
            build_profile_receipt=_read_external_candidate_member(
                build["profile_receipt_path"],
                "external clean-candidate build-profile receipt",
            ),
            build_provenance=build_provenance,
            executable_sha256=executable_sha256,
            authorized_pic_root=ORION_BULK_ROOT,
        )
    except QualificationError:
        raise
    except (OSError, ValueError) as error:
        _fail(
            "external_clean_candidate_closure_error",
            f"external clean-candidate closure validation failed: {error}",
        )
    return {
        "freeze_id": frozen["freeze_id"],
        "referenced_manifest_sha256": _sha256_bytes(external_manifest_payload),
        "retained_executable_sha256": executable_sha256,
        "executable": executable_report,
        "validated_submodule_count": len(validated_submodules),
        "validation": "control_plane_validate_clean_candidate_bundle",
    }


def _validate_bound_files(
    parsed: Mapping[str, Any], snapshot: Any
) -> tuple[dict[str, Any], dict[str, Any]]:
    candidate = parsed["candidate_binding"]["clean_candidate_manifest"]
    retained_candidate_payload = _member_payload(
        snapshot, candidate, "clean candidate manifest"
    )
    frozen = _validate_clean_candidate_manifest(retained_candidate_payload)
    retained_artifacts = {}
    for name, binding in parsed["artifact_bindings"].items():
        retained_artifacts[name] = _member_payload(snapshot, binding, name)
    _cross_bind_clean_candidate(parsed, frozen)
    _validate_elf_identity(
        snapshot.member_path(parsed["artifact_bindings"]["executable"]["path"]),
        retained_artifacts["executable"],
        parsed["artifact_bindings"]["executable"]["sha256"],
        "retained campaign executable",
    )
    external_closure = _validate_external_clean_candidate_closure(
        frozen,
        retained_candidate_payload,
        parsed["artifact_bindings"]["executable"]["sha256"],
    )
    _require(
        parsed["artifact_bindings"]["analyzer"]["sha256"] == _sha256_path(ANALYZER_PATH),
        "hash_drift",
        "invoked analyzer source SHA-256 drifted",
    )
    return frozen, external_closure


def _required_retained_times(policy: Mapping[str, Any]) -> list[float]:
    return policy["artifact_retention_policy"]["required_output_times_omega0_inverse"]


def _snapshot_tolerance(policy: Mapping[str, Any]) -> float:
    snapshot_policy = policy["athenak_selected_release_criteria"]["snapshot_selection"]
    return snapshot_policy["absolute_match_tolerance"]


def _required_time_key(snapshot_time: float, policy: Mapping[str, Any]) -> str:
    matches = [
        required_time
        for required_time in _required_retained_times(policy)
        if abs(snapshot_time - required_time) <= _snapshot_tolerance(policy)
    ]
    _require(
        len(matches) == 1,
        "snapshot_cadence_drift",
        f"snapshot time {snapshot_time!r} does not identify exactly one retained cadence time",
    )
    return f"{matches[0]:.1f}"


def _select_retained_snapshot_products(
    products: Sequence[Mapping[str, Any]], policy: Mapping[str, Any]
) -> dict[str, dict[str, dict[str, Any]]]:
    retained: dict[str, dict[str, dict[str, Any]]] = {}
    raw_products = [product for product in products if product["kind"] in _ENDPOINT_PRODUCT_KINDS]
    for product in raw_products:
        _required_time_key(product["snapshot_time"], policy)
    for required_time in _required_retained_times(policy):
        snapshot: dict[str, dict[str, Any]] = {}
        for kind in _ENDPOINT_PRODUCT_KINDS:
            matches = [
                product
                for product in raw_products
                if product["kind"] == kind
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            _require(
                len(matches) == 1,
                "ambiguous_snapshot_cadence" if matches else "missing_snapshot_cadence",
                f"{kind} snapshot at t={required_time:.1f} has {len(matches)} matches",
            )
            snapshot[kind] = dict(matches[0])
        retained[f"{required_time:.1f}"] = snapshot
    return retained


def _select_products(
    retained: Mapping[str, Mapping[str, Mapping[str, Any]]],
    policy: Mapping[str, Any],
) -> dict[str, dict[str, dict[str, Any]]]:
    snapshot_policy = policy["athenak_selected_release_criteria"]["snapshot_selection"]
    selected: dict[str, dict[str, dict[str, Any]]] = {}
    for required_time in snapshot_policy["required_times"]:
        key = f"{required_time:.1f}"
        _require(key in retained, "missing_endpoint", f"missing endpoint snapshot t={key}")
        selected[key] = {kind: dict(product) for kind, product in retained[key].items()}
    return selected


def _select_stdout_product(products: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    matches = [product for product in products if product["kind"] == "stdout"]
    _require(
        len(matches) == 1,
        "duplicate_product" if matches else "missing_product",
        f"stdout product has {len(matches)} matches",
    )
    return dict(matches[0])


def _fnv1a64(payload: bytes) -> int:
    digest = 14695981039346656037
    for value in payload:
        digest ^= value
        digest = (digest * 1099511628211) & ((1 << 64) - 1)
    return digest


def _validate_restart_marker(marker_payload: bytes, artifact_payload: bytes, label: str) -> None:
    match = _RESTART_MARKER_PATTERN.fullmatch(marker_payload)
    _require(
        match is not None,
        "restart_marker_corruption",
        f"{label}: malformed restart completion marker",
    )
    size = int(match.group(1))
    digest = int(match.group(2), 16)
    _require(
        size == len(artifact_payload) and digest == _fnv1a64(artifact_payload),
        "restart_marker_corruption",
        f"{label}: completion marker does not bind retained artifact bytes",
    )


def _validate_restart_payload_identity(payload: bytes, label: str) -> None:
    marker = b"<par_end>\n"
    location = payload.find(marker)
    _require(
        0 < location <= 160 * 1024 and location + len(marker) < len(payload),
        "invalid_restart_payload",
        f"{label}: retained restart lacks serialized parameter header and binary body",
    )
    try:
        header = payload[:location].decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_restart_payload", f"{label}: parameter header is not UTF-8: {error}")
    _require(
        header.startswith("<") and "=" in header,
        "invalid_restart_payload",
        f"{label}: serialized parameter header is malformed",
    )


def _parse_restart_manifest(payload: bytes, manifest_path: str) -> list[dict[str, Any]]:
    item = _object(
        _load_json_bytes(payload, f"restart manifest {manifest_path}"),
        {"schema", "members"},
        f"restart manifest {manifest_path}",
    )
    _require(
        item["schema"] == "ATHENAK_RESTART_MANIFEST_V1",
        "invalid_restart_manifest",
        f"{manifest_path}: unsupported restart manifest schema",
    )
    members = []
    for index, raw_member in enumerate(_list(item["members"], f"{manifest_path}/members")):
        label = f"{manifest_path}/members[{index}]"
        member = _object(raw_member, {"path", "size", "fnv1a64"}, label)
        path = _relative_path(member["path"], f"{label}/path")
        size = _exact_int(member["size"], f"{label}/size")
        digest = _text(member["fnv1a64"], f"{label}/fnv1a64")
        _require(
            size >= 0 and re.fullmatch(r"[0-9a-f]{16}", digest) is not None,
            "invalid_restart_manifest",
            f"{label}: invalid restart member digest",
        )
        members.append({"path": path, "size": size, "fnv1a64": digest})
    _require(
        bool(members) and len({member["path"] for member in members}) == len(members),
        "invalid_restart_manifest",
        f"{manifest_path}: restart members are empty or duplicated",
    )
    return members


def _validate_restart_manifest_layout(
    manifest_path: str, members: Sequence[Mapping[str, Any]]
) -> str:
    manifest = PurePosixPath(manifest_path)
    _require(
        manifest.parent == PurePosixPath("rst") and manifest.name.endswith(".rst.manifest"),
        "restart_publication_layout_drift",
        f"{manifest_path}: restart manifest must be directly below rst/",
    )
    payload_name = manifest.name.removesuffix(".manifest")
    ranked = []
    for member in members:
        path = PurePosixPath(member["path"])
        _require(
            path.name == payload_name and path.parts[0] == "rst",
            "restart_publication_layout_drift",
            f"{manifest_path}: restart member filename or rst/ root drifted",
        )
        parent = path.parent.name
        match = re.fullmatch(r"rank_([0-9]{8})", parent)
        _require(
            path.parent == PurePosixPath("rst")
            or (match is not None and path.parent.parent == PurePosixPath("rst")),
            "restart_publication_layout_drift",
            f"{manifest_path}: restart member directory layout drifted",
        )
        ranked.append(None if match is None else int(match.group(1)))
    _require(
        all(rank is None for rank in ranked) or all(rank is not None for rank in ranked),
        "restart_publication_layout_drift",
        f"{manifest_path}: restart manifest mixes shared and per-rank members",
    )
    if ranked[0] is None:
        _require(
            len(ranked) == 1,
            "restart_publication_layout_drift",
            f"{manifest_path}: shared restart manifest must contain one member",
        )
        return "shared_mpi_io"
    _require(
        sorted(ranked) == list(range(len(ranked))),
        "restart_publication_layout_drift",
        f"{manifest_path}: per-rank restart members are not contiguous",
    )
    return "per_rank_shards"


def _validate_restart_publications(
    products: Sequence[Mapping[str, Any]], snapshot: Any, policy: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    restart_products = [
        product for product in products if product["kind"] in _RESTART_PRODUCT_KINDS
    ]
    for product in restart_products:
        _required_time_key(product["snapshot_time"], policy)
    actual = {(product["kind"], product["path"]) for product in restart_products}
    expected: set[tuple[str, str]] = set()
    report: dict[str, dict[str, Any]] = {}
    for required_time in _required_retained_times(policy):
        key = f"{required_time:.1f}"
        manifest_matches = [
            product
            for product in restart_products
            if product["kind"] == "restart_manifest"
            and abs(product["snapshot_time"] - required_time) <= _snapshot_tolerance(policy)
        ]
        _require(
            len(manifest_matches) == 1,
            "ambiguous_restart_publication" if manifest_matches else "missing_restart_publication",
            f"restart manifest at t={key} has {len(manifest_matches)} matches",
        )
        manifest_product = manifest_matches[0]
        manifest_path = manifest_product["path"]
        manifest_payload = _member_payload(
            snapshot, manifest_product, f"restart manifest at t={key}"
        )
        members = _parse_restart_manifest(manifest_payload, manifest_path)
        layout = _validate_restart_manifest_layout(manifest_path, members)
        manifest_marker_path = manifest_path + ".complete"
        expected.add(("restart_manifest", manifest_path))
        expected.add(("restart_manifest_complete", manifest_marker_path))
        marker_matches = [
            product
            for product in restart_products
            if product["kind"] == "restart_manifest_complete"
            and product["path"] == manifest_marker_path
            and abs(product["snapshot_time"] - required_time) <= _snapshot_tolerance(policy)
        ]
        _require(
            len(marker_matches) == 1,
            "restart_publication_layout_drift",
            f"{manifest_path}: manifest completion marker is missing or ambiguous",
        )
        _validate_restart_marker(
            _member_payload(snapshot, marker_matches[0], f"{manifest_path} marker"),
            manifest_payload,
            f"{manifest_path} marker",
        )
        for member in members:
            member_path = member["path"]
            marker_path = member_path + ".complete"
            expected.add(("restart", member_path))
            expected.add(("restart_complete", marker_path))
            payload_matches = [
                product
                for product in restart_products
                if product["kind"] == "restart"
                and product["path"] == member_path
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            marker_matches = [
                product
                for product in restart_products
                if product["kind"] == "restart_complete"
                and product["path"] == marker_path
                and abs(product["snapshot_time"] - required_time)
                <= _snapshot_tolerance(policy)
            ]
            _require(
                len(payload_matches) == 1 and len(marker_matches) == 1,
                "restart_publication_layout_drift",
                f"{member_path}: payload or completion marker is missing or ambiguous",
            )
            member_payload = _member_payload(
                snapshot, payload_matches[0], f"restart payload {member_path}"
            )
            _validate_restart_payload_identity(member_payload, member_path)
            _validate_restart_marker(
                _member_payload(snapshot, marker_matches[0], f"{member_path} marker"),
                member_payload,
                f"{member_path} marker",
            )
            _require(
                member["size"] == len(member_payload)
                and int(member["fnv1a64"], 16) == _fnv1a64(member_payload),
                "invalid_restart_manifest",
                f"{manifest_path}: member digest drifted for {member_path}",
            )
        report[key] = {
            "manifest_path": manifest_path,
            "layout": layout,
            "member_count": len(members),
        }
    _require(
        actual == expected,
        "restart_publication_layout_drift",
        "restart product inventory contains missing or unexpected publication members",
    )
    return report


def _validate_declared_tree_members(parsed: Mapping[str, Any], snapshot: Any) -> None:
    declared = {MANIFEST_NAME, INVENTORY_NAME, FREEZE_RECEIPT_NAME}
    declared.add(parsed["candidate_binding"]["clean_candidate_manifest"]["path"])
    declared.update(binding["path"] for binding in parsed["artifact_bindings"].values())
    declared.update(product["path"] for product in parsed["products"])
    _require(
        len(declared)
        == 3 + 1 + len(parsed["artifact_bindings"]) + len(parsed["products"]),
        "duplicate_declared_path",
        "campaign manifest contains duplicate retained paths",
    )
    measured = snapshot.relative_files()
    _require(
        measured == declared,
        "artifact_inventory_drift",
        "campaign tree membership does not match declared retained files",
    )


def _validate_product_hashes(
    products: Sequence[Mapping[str, Any]], snapshot: Any
) -> None:
    for product in products:
        _member_sha256(snapshot, product, f"{product['kind']} product {product['path']}")


def _parse_pvtk_execution_header(payload: bytes, label: str) -> dict[str, Any]:
    match = _PVTK_EXECUTION_PATTERN.match(payload[:4096])
    _require(match is not None, "invalid_prtcl_all", f"{label}: malformed execution header")
    try:
        time = float(match.group(1))
        nranks = int(match.group(2))
        cycle = int(match.group(3))
        variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        _fail("invalid_prtcl_all", f"{label}: invalid execution header: {error}")
    _require(
        math.isfinite(time) and nranks > 0 and cycle >= 0 and variables == "prtcl_all",
        "invalid_prtcl_all",
        f"{label}: execution header binding drifted",
    )
    return {"time": time, "nranks": nranks, "cycle": cycle, "variables": variables}


def _validate_particle_snapshot(path: Path, payload: bytes, label: str) -> dict[str, Any]:
    execution = _parse_pvtk_execution_header(payload, label)
    try:
        data: ParticleVTKData = read_particle_vtk(path)
    except (OSError, ValueError) as error:
        _fail("invalid_prtcl_all", f"{label}: particle VTK decode failed: {error}")
    _require(
        set(data.scalars) == _PVTK_SCALARS,
        "invalid_prtcl_all",
        f"{label}: particle scalar inventory drifted",
    )
    _require(
        set(data.vectors) == {"vel"},
        "invalid_prtcl_all",
        f"{label}: particle vector inventory drifted",
    )
    _require(
        data.points.shape == data.vectors["vel"].shape
        and data.points.ndim == 2
        and data.points.shape[1] == 3,
        "invalid_prtcl_all",
        f"{label}: particle point/vector shape drifted",
    )
    _require(
        data.points.shape[0] > 0,
        "invalid_prtcl_all",
        f"{label}: particle payload is empty",
    )
    _require(
        all(
            np.issubdtype(data.scalars[name].dtype, np.integer)
            for name in _PVTK_INTEGER_SCALARS
        )
        and all(
            np.issubdtype(data.scalars[name].dtype, np.floating)
            for name in _PVTK_REAL_SCALARS
        ),
        "invalid_prtcl_all",
        f"{label}: particle scalar integer/real typing drifted",
    )
    _require(
        np.all(np.isfinite(data.points))
        and np.all(np.isfinite(data.vectors["vel"]))
        and all(np.all(np.isfinite(data.scalars[name])) for name in _PVTK_REAL_SCALARS),
        "invalid_prtcl_all",
        f"{label}: particle payload contains non-finite values",
    )
    _require(
        np.all(data.scalars["macro_weight"] >= 0.0)
        and np.all(data.scalars["gid"] >= 0)
        and np.all(data.scalars["ptag"] >= 0)
        and np.all(data.scalars["species"] >= 0)
        and np.all(np.isin(data.scalars["cr_source"], (0, 1))),
        "invalid_prtcl_all",
        f"{label}: particle provenance or macro weight is invalid",
    )
    _require(
        np.unique(data.scalars["ptag"]).size == data.points.shape[0],
        "invalid_prtcl_all",
        f"{label}: particle ptag values are not unique",
    )
    return {
        "particle_count": int(data.points.shape[0]),
        "execution_header": execution,
        "storage_layout": "single_mpi_io_pvtk_file_with_nranks_header",
    }


def _validate_snapshot_payloads(
    retained: Mapping[str, Mapping[str, Mapping[str, Any]]],
    snapshot: Any,
    policy: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    report = {}
    tolerance = _snapshot_tolerance(policy)
    for time, products in retained.items():
        expected_time = float(time)
        cycles: dict[str, int] = {}
        mesh_report = {}
        for kind in ("rho", "bmag", "prtcl_jx", "j2"):
            product = products[kind]
            payload = _member_payload(snapshot, product, f"{kind} snapshot at t={time}")
            try:
                dataset = output_primitives.parse_athenak_binary_bytes(
                    payload, source=product["path"]
                )
            except ValueError as error:
                _fail("invalid_mesh_bin", f"{kind} snapshot at t={time}: {error}")
            _require(
                abs(dataset.time - expected_time) <= tolerance
                and abs(dataset.time - product["snapshot_time"]) <= tolerance,
                "snapshot_metadata_drift",
                f"{kind} snapshot at t={time}: embedded time is {dataset.time!r}",
            )
            _require(
                dataset.variable_names == _EXPECTED_BIN_FIELDS[kind],
                "invalid_mesh_bin",
                f"{kind} snapshot at t={time}: binary field inventory drifted",
            )
            cycles[kind] = dataset.cycle
            mesh_report[kind] = {
                "cycle": dataset.cycle,
                "time": dataset.time,
                "meshblock_count": len(dataset.blocks),
            }
        particle = products["prtcl_all"]
        particle_payload = _member_payload(
            snapshot, particle, f"prtcl_all snapshot at t={time}"
        )
        particle_report = _validate_particle_snapshot(
            snapshot.member_path(particle["path"]),
            particle_payload,
            f"prtcl_all snapshot at t={time}",
        )
        header = particle_report["execution_header"]
        _require(
            abs(header["time"] - expected_time) <= tolerance
            and abs(header["time"] - particle["snapshot_time"]) <= tolerance,
            "snapshot_metadata_drift",
            f"prtcl_all snapshot at t={time}: embedded time is {header['time']!r}",
        )
        cycles["prtcl_all"] = header["cycle"]
        _require(
            len(set(cycles.values())) == 1,
            "snapshot_metadata_drift",
            f"snapshot at t={time}: embedded cycles differ across retained raw products",
        )
        report[time] = {
            "cycle": next(iter(cycles.values())),
            "mesh_bins": mesh_report,
            "particles": particle_report,
        }
    return report


def _validate_stdout_telemetry(payload: bytes) -> dict[str, Any]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        _fail("invalid_q017_telemetry", f"stdout telemetry is not UTF-8: {error}")
    telemetry = {}
    for match in _Q017_TELEMETRY_PATTERN.finditer(text):
        name = match.group(1)
        _require(
            name not in telemetry,
            "invalid_q017_telemetry",
            f"duplicate Q-017 stdout telemetry name: {name}",
        )
        try:
            telemetry[name] = float(match.group(2))
        except ValueError:
            _fail("invalid_q017_telemetry", f"non-numeric Q-017 telemetry value: {name}")
    missing = sorted(_Q017_REQUIRED_NAMES - set(telemetry))
    _require(
        not missing,
        "invalid_q017_telemetry",
        f"stdout is missing required Q-017 telemetry: {', '.join(missing)}",
    )
    _require(
        all(math.isfinite(value) and value >= 0.0 for value in telemetry.values())
        and telemetry["schema_version"] == 2.0,
        "invalid_q017_telemetry",
        "stdout Q-017 telemetry contains invalid values or schema version",
    )
    return {
        "schema_version": telemetry["schema_version"],
        "mpi_ranks": telemetry["mpi.ranks"],
        "retained_name_count": len(telemetry),
    }


def _pairing_schema(
    identity: Mapping[str, Any], policy: Mapping[str, Any]
) -> dict[str, Any]:
    matrix = policy["athenak_selected_release_criteria"]["campaign_matrix"]
    variant = identity["variant"]
    seed = identity["seed"]
    if variant == "three_level_amr_root_dx12_finest_dx3":
        counterpart = {"variant": "fine_uniform_dx3", "seed": seed}
    elif variant == "fine_uniform_dx3":
        counterpart = {"variant": "three_level_amr_root_dx12_finest_dx3", "seed": seed}
    else:
        counterpart = None
    return {
        "paired_seed_rule": matrix["paired_seed_rule"],
        "pair_key": {"seed": seed},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": "schema_wired_not_evaluated_by_artifact_admission_slice",
    }


def spectrum_wiring(policy: Mapping[str, Any]) -> dict[str, Any]:
    """Expose the frozen weighted-spectrum arithmetic contract without inventing data."""
    spectrum = policy["athenak_selected_release_criteria"]["spectrum"]
    return {
        "histogram_primitive": "analyze_q011_section54_outputs.fixed_histogram",
        "slope_fit_primitive": (
            "analyze_q011_section54_outputs.fit_fixed_bin_loglog_slope"
        ),
        "bin_edges": spectrum["bin_edges"],
        "histogram_weight": spectrum["histogram_weight"],
        "slope_fit_window_chi": spectrum["slope_fit_window_chi"],
        "minimum_positive_fit_bins": spectrum["minimum_positive_fit_bins"],
        "runtime_projection_status": (
            "wired_requires_reviewed_chi_values_not_inferred_from_pvtk_velocity"
        ),
    }


def weighted_spectrum_from_chi(
    chi: object, macro_weight: object, policy: Mapping[str, Any]
) -> dict[str, Any]:
    """Apply the preregistered fixed weighted histogram to caller-reviewed chi values."""
    spectrum = policy["athenak_selected_release_criteria"]["spectrum"]
    histogram = output_primitives.fixed_histogram(
        chi,
        spectrum["bin_edges"],
        weights=macro_weight,
    )
    widths = np.diff(histogram.bin_edges)
    f_chi = histogram.weighted_counts / widths
    macro_weight_in_bins = float(np.sum(histogram.weighted_counts))
    admitted_weight = float(
        macro_weight_in_bins + histogram.underflow_weight + histogram.overflow_weight
    )
    plotted = np.zeros_like(f_chi)
    if admitted_weight > 0.0:
        plotted = histogram.geometric_bin_centers * f_chi / admitted_weight
    overflow_fraction = (
        histogram.overflow_weight / admitted_weight if admitted_weight > 0.0 else 0.0
    )
    return {
        "bin_edges": histogram.bin_edges.tolist(),
        "counts": histogram.counts.tolist(),
        "weighted_counts": histogram.weighted_counts.tolist(),
        "f_chi": f_chi.tolist(),
        "normalized_chi_f_chi": plotted.tolist(),
        "underflow_count": histogram.underflow_count,
        "overflow_count": histogram.overflow_count,
        "underflow_weight": histogram.underflow_weight,
        "overflow_weight": histogram.overflow_weight,
        "macro_weight_in_bins": macro_weight_in_bins,
        "admitted_weight": admitted_weight,
        "overflow_macro_weight_fraction": overflow_fraction,
        "overflow_gate": {
            "maximum_macro_weight_fraction": spectrum[
                "max_overflow_macro_weight_fraction"
            ],
            "passed": bool(
                overflow_fraction <= spectrum["max_overflow_macro_weight_fraction"]
            ),
        },
    }


def _validate_freeze_receipt_semantics(receipt: object) -> None:
    require_exact_primitive_types(
        receipt,
        _EXPECTED_FREEZE_RECEIPT,
        error_type=QualificationError,
        label="campaign freeze receipt",
    )
    _require(
        receipt == _EXPECTED_FREEZE_RECEIPT,
        "freeze_receipt_semantics_drift",
        "campaign freeze receipt role or effect semantics drifted",
    )


def _admit_campaign(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_orion_root: Path,
) -> dict[str, Any]:
    root = authorized_tree_root(
        campaign_root,
        authorized_root=authorized_orion_root,
        error_type=QualificationError,
        label="Q-011 Section 5.4 campaign root",
    )
    with staged_verified_frozen_tree(
        root,
        expected_inventory_sha256,
        authorized_root=authorized_orion_root,
        error_type=QualificationError,
        label="Q-011 Section 5.4 campaign tree",
    ) as (tree_report, snapshot):
        _validate_freeze_receipt_semantics(tree_report["freeze_receipt"])
        manifest = _load_json_bytes(
            snapshot.member_path(MANIFEST_NAME).read_bytes(),
            "campaign manifest",
        )
        parsed = _validate_manifest_schema(manifest, root)
        _validate_declared_tree_members(parsed, snapshot)
        frozen_candidate, external_candidate_closure = _validate_bound_files(parsed, snapshot)
        policy = _load_bound_policy(
            snapshot, parsed["artifact_bindings"]["preregistration"]
        )
        _validate_identity(parsed["run_identity"], policy)
        retained_snapshot_products = _select_retained_snapshot_products(
            parsed["products"], policy
        )
        endpoint_products = _select_products(retained_snapshot_products, policy)
        stdout_product = _select_stdout_product(parsed["products"])
        _validate_product_hashes(parsed["products"], snapshot)
        snapshot_payloads = _validate_snapshot_payloads(
            retained_snapshot_products, snapshot, policy
        )
        restart_publications = _validate_restart_publications(
            parsed["products"], snapshot, policy
        )
        stdout_telemetry = _validate_stdout_telemetry(
            _member_payload(snapshot, stdout_product, "stdout product")
        )
        particle_endpoints = {
            time: snapshot_payloads[time]["particles"] for time in endpoint_products
        }
        return {
            "run_identity": parsed["run_identity"],
            "candidate_binding": parsed["candidate_binding"],
            "frozen_clean_candidate": frozen_candidate,
            "external_clean_candidate_closure": external_candidate_closure,
            "artifact_bindings": parsed["artifact_bindings"],
            "preregistration_binding": {
                "sha256": parsed["artifact_bindings"]["preregistration"]["sha256"],
                "expected_sha256": EXPECTED_PREREGISTRATION_SHA256,
                "binding_scope": "complete_retained_bytes_equal_invoked_frozen_policy",
            },
            "retained_snapshot_products": retained_snapshot_products,
            "endpoint_products": endpoint_products,
            "run_products": {"stdout": stdout_product},
            "snapshot_payloads": snapshot_payloads,
            "restart_publications": restart_publications,
            "stdout_telemetry": stdout_telemetry,
            "particle_endpoints": particle_endpoints,
            "immutable_tree": {
                "inventory_sha256": tree_report["inventory_sha256"],
                "inventoried_file_count": tree_report["inventoried_file_count"],
                "recursively_read_only": tree_report["recursively_read_only"],
            },
            "weighted_spectrum": spectrum_wiring(policy),
            "amr_pairing": _pairing_schema(parsed["run_identity"], policy),
            "numerical_qualification_status": (
                "not_evaluated_by_artifact_admission_slice"
            ),
        }


def qualify_campaign(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_orion_root: Path = ORION_BULK_ROOT,
) -> dict[str, Any]:
    """Return one deterministic admission result without partial success."""
    result = {
        "schema_version": 1,
        "record_type": RESULT_RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "qualification_scope": QUALIFICATION_SCOPE,
        "admitted_for_follow_on_numerical_qualification": False,
        "final_claim_closure": False,
        "status": "rejected",
        "failure_reasons": [],
        "admission": None,
    }
    try:
        admission = _admit_campaign(
            campaign_root,
            expected_inventory_sha256,
            authorized_orion_root=authorized_orion_root,
        )
    except (QualificationError, OSError, ValueError) as error:
        code = getattr(error, "code", "invalid_campaign")
        result["failure_reasons"] = [{"code": code, "message": str(error)}]
        return result
    result["admitted_for_follow_on_numerical_qualification"] = True
    result["status"] = "admitted_for_follow_on_numerical_qualification"
    result["admission"] = admission
    return result


def _write_json(value: Mapping[str, Any]) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("campaign_root")
    parser.add_argument("expected_inventory_sha256")
    args = parser.parse_args(argv)
    result = qualify_campaign(args.campaign_root, args.expected_inventory_sha256)
    _write_json(result)
    return 0 if result["admitted_for_follow_on_numerical_qualification"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
