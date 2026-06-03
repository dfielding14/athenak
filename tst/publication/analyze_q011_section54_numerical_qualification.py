#!/usr/bin/env python3
"""Fail-closed source-local Q-011 Section 5.4 numerical aggregation.

This module revalidates retained attempt trees, runs raw reducers over sealed
snapshot members, and consumes independently retained pair/restart recompute
receipts.  It does not discover artifacts, execute AthenaK, publish evidence,
automate human review, or close a manuscript claim.
"""

from __future__ import annotations

import copy
from contextlib import contextmanager
import ctypes
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Iterator, Mapping, Sequence

import numpy as np

if __package__:
    from . import analyze_q011_section54_campaign as admission
    from . import frontier_f1_structured_artifacts as structured_artifacts
    from . import immutable_orion_tree
    from . import q011_section54_artifacts as artifacts
    from . import q011_section54_particles as particles
    from . import q011_section54_restart as restart
    from . import q011_section54_spatial as spatial
else:
    import analyze_q011_section54_campaign as admission
    import frontier_f1_structured_artifacts as structured_artifacts
    import immutable_orion_tree
    import q011_section54_artifacts as artifacts
    import q011_section54_particles as particles
    import q011_section54_restart as restart
    import q011_section54_spatial as spatial


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_source_local_numerical_qualification"
ATTEMPT_RECORD_TYPE = "q011_section54_source_local_numerical_baseline_attempt"
PAIR_RECOMPUTE_RECORD_TYPE = "q011_section54_retained_paired_amr_fine_recompute"
RESTART_RECOMPUTE_RECORD_TYPE = "q011_section54_retained_restart_parity_recompute"
PAIR_RECOMPUTE_MEMBER = "paired_amr_fine_result.json"
RESTART_RECOMPUTE_MEMBER = "restart_parity_binding.json"
PAIR_RECOMPUTE_ANALYZER = (
    "tst/publication/analyze_q011_section54_numerical_qualification.py"
)
ARTIFACTS_RECOMPUTE_ANALYZER = "tst/publication/q011_section54_artifacts.py"
CAMPAIGN_RECOMPUTE_ANALYZER = "tst/publication/analyze_q011_section54_campaign.py"
IMMUTABLE_TREE_RECOMPUTE_ANALYZER = "tst/publication/immutable_orion_tree.py"
PARTICLES_RECOMPUTE_ANALYZER = "tst/publication/q011_section54_particles.py"
PLANNER_RECOMPUTE_ANALYZER = (
    "tst/publication/q011_section54_qualifying_campaign_execution.py"
)
MODEL_RECOMPUTE_ANALYZER = "tst/publication/q011_section54_model.py"
PRESSURE_SELECTION_RECOMPUTE_ANALYZER = (
    "tst/publication/q011_section54_pressure_selection.py"
)
PRESSURE_PILOT_EXECUTION_RECOMPUTE_ANALYZER = (
    "tst/publication/q011_section54_pressure_pilot_execution.py"
)
PRESSURE_PILOT_PUBLISHER_RECOMPUTE_ANALYZER = (
    "tst/publication/publish_q011_section54_pressure_pilot_bundle.py"
)
PRESSURE_PILOT_RECOMPUTE_ANALYZER = (
    "tst/publication/analyze_q011_section54_pressure_pilot.py"
)
PRESSURE_PILOT_CASE_RECOMPUTE_ANALYZER = (
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py"
)
STRUCTURED_ARTIFACTS_RECOMPUTE_ANALYZER = (
    "tst/publication/frontier_f1_structured_artifacts.py"
)
RESTART_RECOMPUTE_ANALYZER = "tst/publication/q011_section54_restart.py"
SPATIAL_RECOMPUTE_ANALYZER = "tst/publication/q011_section54_spatial.py"
OUTPUT_RECOMPUTE_ANALYZER = "tst/publication/analyze_q011_section54_outputs.py"
PVTK_RECOMPUTE_ANALYZER = "tst/publication/pvtk_particles.py"
CONTROL_PLANE_COMMON_RECOMPUTE_ANALYZER = (
    "tst/publication/frontier_control_plane/control_plane_common.py"
)
CONTROL_PLANE_LEDGER_RECOMPUTE_ANALYZER = (
    "tst/publication/frontier_control_plane/ledger.py"
)
CONTROL_PLANE_ATTESTATION_RECOMPUTE_ANALYZER = (
    "tst/publication/frontier_control_plane/operator_attestation.py"
)
PAIR_RECOMPUTE_ANALYZERS = (
    PAIR_RECOMPUTE_ANALYZER,
    ARTIFACTS_RECOMPUTE_ANALYZER,
    IMMUTABLE_TREE_RECOMPUTE_ANALYZER,
    PARTICLES_RECOMPUTE_ANALYZER,
    SPATIAL_RECOMPUTE_ANALYZER,
    OUTPUT_RECOMPUTE_ANALYZER,
    MODEL_RECOMPUTE_ANALYZER,
)
RESTART_RECOMPUTE_ANALYZERS = (
    *PAIR_RECOMPUTE_ANALYZERS,
    CAMPAIGN_RECOMPUTE_ANALYZER,
    PLANNER_RECOMPUTE_ANALYZER,
    MODEL_RECOMPUTE_ANALYZER,
    PRESSURE_SELECTION_RECOMPUTE_ANALYZER,
    PRESSURE_PILOT_EXECUTION_RECOMPUTE_ANALYZER,
    PRESSURE_PILOT_PUBLISHER_RECOMPUTE_ANALYZER,
    PRESSURE_PILOT_RECOMPUTE_ANALYZER,
    PRESSURE_PILOT_CASE_RECOMPUTE_ANALYZER,
    STRUCTURED_ARTIFACTS_RECOMPUTE_ANALYZER,
    RESTART_RECOMPUTE_ANALYZER,
    PVTK_RECOMPUTE_ANALYZER,
    CONTROL_PLANE_COMMON_RECOMPUTE_ANALYZER,
    CONTROL_PLANE_LEDGER_RECOMPUTE_ANALYZER,
    CONTROL_PLANE_ATTESTATION_RECOMPUTE_ANALYZER,
)
RESTART_SCREEN_SCOPE = (
    "bounded_selected_channel_restart_continuation_screen_not_full_state_equivalence"
)
PLANNED_RESTART_CARRIER_PATH = (
    "restart_continuation/amr_restart_continuation_carrier.json"
)
QUALIFICATION_SCOPE = (
    "bounded_source_local_numerical_aggregation_only_external_review_required_"
    "no_claim_closure"
)
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
EXPECTED_BASELINE_ATTEMPTS = 24
EXPECTED_MATRIX_CELLS = tuple(
    (variant, seed) for variant in GRID_VARIANTS for seed in QUALIFYING_SEEDS
)
PAIRED_SEED_RULE = (
    "Use the same qualifying seed for coarse-uniform, AMR and fine-uniform variants."
)
PAIRED_RESIDUAL_THRESHOLDS = (
    ("shock_front_position_at_t500", 240.0, None, None),
    ("upstream_magnetic_amplification_at_t500", 0.35, None, None),
    ("rho_y_average_at_t500", None, 0.2, 0.3),
    ("bmag_y_average_at_t500", None, 0.25, 0.35),
    ("normalized_downstream_chi_f_chi_at_t500", None, 0.2, 0.3),
    ("normalized_downstream_chi_f_chi_at_t1200", None, 0.2, 0.3),
)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_ATTEMPT_ID = re.compile(r"[a-z0-9][a-z0-9._-]{0,127}")
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_REGULAR_FLAGS = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
_WRITE_BITS = stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH
_INOTIFY_MUTATION_MASK = (
    0x00000004  # IN_ATTRIB
    | 0x00000100  # IN_CREATE
    | 0x00000200  # IN_DELETE
    | 0x00000400  # IN_DELETE_SELF
    | 0x00000040  # IN_MOVED_FROM
    | 0x00000080  # IN_MOVED_TO
    | 0x00000800  # IN_MOVE_SELF
)
_PRODUCTION_BINDING_KEY = object()
_REPO_ROOT = Path(__file__).resolve().parents[2]
_PAIR_PROFILE_DX = 12.0
_RESTART_MESH_MEMBERS = {
    "rho_bin": "dens",
    "bmag_bin": "bmag",
    "prtcl_jx_bin": "prtcl_jx",
    "j2_bin": "j2",
}
_RESTART_PARTICLE_MEMBER = "prtcl_all_pvtk"
_RESTART_PARTICLE_INTEGER_SCALARS = ("gid", "ptag", "species", "cr_source")
_RESTART_PARTICLE_REAL_SCALARS = (
    "macro_weight",
    "birth_time",
    "deltaf_f0",
    "deltaf_weight",
)


class NumericalQualificationError(ValueError):
    """Raised when source-local numerical aggregation cannot proceed."""


class _RegisteredRunStructuredArtifactTree(
    structured_artifacts.StructuredArtifactTree
):
    """Permit the reconciler receipt while retaining exact raw-tree closure."""

    def require_analysis_identity(self) -> None:
        expected = os.fstat(self.analysis_fd)
        actual = os.stat("analysis", dir_fd=self.root_fd, follow_symlinks=False)
        if (
            not stat.S_ISDIR(expected.st_mode)
            or not stat.S_ISDIR(actual.st_mode)
            or stat.S_IMODE(expected.st_mode) != 0o700
            or stat.S_IMODE(actual.st_mode) != 0o700
            or (expected.st_dev, expected.st_ino) != (actual.st_dev, actual.st_ino)
        ):
            raise ValueError("Structured artifact analysis directory changed during analysis")
        allowed = {
            "analysis.json",
            "offline_analysis_receipt.json",
            admission.REGISTERED_EXECUTION_RECEIPT_NAME,
        }
        for name in os.listdir(self.analysis_fd):
            if name not in allowed:
                raise ValueError("Structured artifact analysis directory has an unexpected entry")
            metadata = os.stat(name, dir_fd=self.analysis_fd, follow_symlinks=False)
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_mode & _WRITE_BITS:
                raise ValueError("Structured artifact analysis result is not read-only")


@dataclass(frozen=True)
class _RetainedAttemptBinding:
    wrapper_payload: bytes
    _key: object

    def __post_init__(self) -> None:
        _require(self._key is _PRODUCTION_BINDING_KEY, "invalid retained attempt binding")
        _require(type(self.wrapper_payload) is bytes, "invalid retained attempt payload")


@dataclass(frozen=True)
class _RetainedPairResultBinding:
    record_payload: bytes
    _key: object

    def __post_init__(self) -> None:
        _require(self._key is _PRODUCTION_BINDING_KEY, "invalid retained pair binding")
        _require(type(self.record_payload) is bytes, "invalid retained pair payload")


@dataclass(frozen=True)
class _RetainedRestartParityBinding:
    record_payload: bytes
    _key: object

    def __post_init__(self) -> None:
        _require(self._key is _PRODUCTION_BINDING_KEY, "invalid retained restart binding")
        _require(type(self.record_payload) is bytes, "invalid retained restart payload")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise NumericalQualificationError(message)


def _object(value: object, keys: set[str], *, label: str) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(set(value) == keys, f"{label}: keys drifted")
    return value


def _required_object(
    value: object, keys: set[str], *, label: str
) -> dict[str, Any]:
    _require(type(value) is dict, f"{label}: expected object")
    _require(keys <= set(value), f"{label}: required keys drifted")
    return value


def _list(value: object, *, label: str) -> list[Any]:
    _require(type(value) is list, f"{label}: expected list")
    return value


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


def _canonical_sha256(value: object) -> str:
    try:
        return artifacts.sha256_bytes(artifacts.canonical_json_bytes(value))
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("record is not canonical finite JSON") from error


def _bind_unit_only_canonical_attempt(
    attempt: Mapping[str, object],
) -> dict[str, object]:
    """Bind one synthetic attempt for unit tests only; never production evidence."""
    _require(isinstance(attempt, Mapping), "attempt binding requires a mapping")
    retained = copy.deepcopy(dict(attempt))
    return {
        "attempt_sha256": _canonical_sha256(retained),
        "attempt": retained,
    }


def _json_bytes(value: object) -> bytes:
    try:
        return artifacts.canonical_json_bytes(value)
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("record is not canonical finite JSON") from error


def _load_json_bytes(payload: bytes, *, label: str) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result = {}
        for key, value in pairs:
            _require(key not in result, f"{label}: duplicate JSON key {key!r}")
            result[key] = value
        return result

    def reject_constant(value: str) -> None:
        raise NumericalQualificationError(f"{label}: forbidden JSON constant {value}")

    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise NumericalQualificationError(f"{label}: JSON is not UTF-8") from error
    try:
        return json.loads(
            text,
            object_pairs_hook=reject_duplicates,
            parse_constant=reject_constant,
        )
    except json.JSONDecodeError as error:
        raise NumericalQualificationError(f"{label}: malformed JSON") from error


def _relative_path(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected path text")
    path = PurePosixPath(value)
    _require(
        not path.is_absolute()
        and path.as_posix() == value
        and value != "."
        and all(part not in {"", ".", ".."} for part in path.parts),
        f"{label}: unsafe relative path",
    )
    return value


def _regular_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _read_anchored_regular(
    root_descriptor: int,
    relative: str,
    *,
    label: str,
    expected_sha256: str | None = None,
) -> bytes:
    path = PurePosixPath(_relative_path(relative, label=label))
    parent_descriptor = os.dup(root_descriptor)
    descriptor: int | None = None
    try:
        for part in path.parts[:-1]:
            child_descriptor = os.open(part, _DIRECTORY_FLAGS, dir_fd=parent_descriptor)
            os.close(parent_descriptor)
            parent_descriptor = child_descriptor
        descriptor = os.open(path.parts[-1], _REGULAR_FLAGS, dir_fd=parent_descriptor)
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & _WRITE_BITS,
            f"{label}: expected one read-only regular file",
        )
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        _require(
            _regular_identity(before) == _regular_identity(after),
            f"{label}: retained file changed while reading",
        )
        measured = hashlib.sha256(payload).hexdigest()
        if expected_sha256 is not None:
            _require(measured == expected_sha256, f"{label}: SHA-256 drifted")
        return bytes(payload)
    except OSError as error:
        raise NumericalQualificationError(f"{label}: cannot read retained file: {error}") from error
    finally:
        if descriptor is not None:
            os.close(descriptor)
        os.close(parent_descriptor)


def _scan_anchored_tree(root_descriptor: int) -> tuple[object, set[str], set[str]]:
    snapshot = immutable_orion_tree._scan_anchored_tree(
        root_descriptor,
        hash_regular=False,
        error_type=NumericalQualificationError,
        label="retained numerical recompute bundle",
    )
    _require(
        not snapshot.root.mode & _WRITE_BITS,
        "retained recompute bundle root is writable",
    )
    _require(
        all(not entry.mode & _WRITE_BITS for entry in snapshot.entries),
        "retained recompute bundle contains a writable node",
    )
    files = {
        entry.relative for entry in snapshot.entries if entry.entry_type == "file"
    }
    directories = {
        entry.relative
        for entry in snapshot.entries
        if entry.entry_type == "directory"
    }
    return snapshot, files, directories


def _required_parent_directories(paths: set[str]) -> set[str]:
    directories = set()
    for relative in paths:
        parent = PurePosixPath(relative).parent
        while parent.as_posix() != ".":
            directories.add(parent.as_posix())
            parent = parent.parent
    return directories


def _load_retained_recompute_record(
    bundle_root: str | Path,
    expected_inventory_sha256: str,
    *,
    member: str,
    authorized_orion_root: Path,
) -> tuple[dict[str, object], dict[str, object], dict[str, bytes]]:
    root = immutable_orion_tree.authorized_tree_root(
        bundle_root,
        authorized_root=authorized_orion_root,
        error_type=NumericalQualificationError,
        label="retained numerical recompute bundle",
    )
    inventory_sha256 = _sha256(
        expected_inventory_sha256, label="retained recompute inventory SHA-256"
    )
    descriptor = os.open(root, _DIRECTORY_FLAGS)
    try:
        immutable_orion_tree._require_root_binding(
            root,
            descriptor,
            authorized_root=authorized_orion_root,
            error_type=NumericalQualificationError,
            label="retained numerical recompute bundle",
        )
        preflight, _, _ = _scan_anchored_tree(descriptor)
        inventory_payload = _read_anchored_regular(
            descriptor,
            artifacts.INVENTORY_NAME,
            label=artifacts.INVENTORY_NAME,
            expected_sha256=inventory_sha256,
        )
        inventory = _object(
            _load_json_bytes(inventory_payload, label=artifacts.INVENTORY_NAME),
            {"record_type", "schema_version", "members"},
            label=artifacts.INVENTORY_NAME,
        )
        _require(
            inventory["record_type"] == "q011_section54_derived_artifact_inventory"
            and type(inventory["schema_version"]) is int
            and inventory["schema_version"] == 1,
            "retained recompute inventory identity drifted",
        )
        members = _list(inventory["members"], label=f"{artifacts.INVENTORY_NAME}/members")
        declared = {artifacts.INVENTORY_NAME}
        payloads = {}
        for index, value in enumerate(members):
            label = f"{artifacts.INVENTORY_NAME}/members[{index}]"
            record = _object(value, {"path", "size", "sha256"}, label=label)
            relative = _relative_path(record["path"], label=f"{label}/path")
            _require(relative not in declared, f"{label}: duplicate retained member")
            size = record["size"]
            _require(type(size) is int and size >= 0, f"{label}: size drifted")
            digest = _sha256(record["sha256"], label=f"{label}/sha256")
            payload = _read_anchored_regular(
                descriptor, relative, label=relative, expected_sha256=digest
            )
            _require(len(payload) == size, f"{label}: retained member size drifted")
            declared.add(relative)
            payloads[relative] = payload
        boundary, files, directories = _scan_anchored_tree(descriptor)
        immutable_orion_tree._require_same_snapshot(
            preflight,
            boundary,
            error_type=NumericalQualificationError,
            label="retained numerical recompute bundle",
            phase="recompute bundle closure census",
        )
        _require(
            files == declared
            and directories == _required_parent_directories(declared),
            "retained recompute bundle tree closure drifted",
        )
        immutable_orion_tree._require_root_binding(
            root,
            descriptor,
            authorized_root=authorized_orion_root,
            error_type=NumericalQualificationError,
            label="retained numerical recompute bundle",
        )
    finally:
        os.close(descriptor)
    _require(artifacts.MANIFEST_NAME in payloads, "retained recompute manifest is missing")
    manifest_value = _load_json_bytes(
        payloads[artifacts.MANIFEST_NAME], label=artifacts.MANIFEST_NAME
    )
    try:
        manifest = artifacts.validate_derived_manifest(manifest_value)
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("retained recompute manifest is invalid") from error
    _require(member in payloads, f"retained recompute member is missing: {member}")
    record_payload = payloads[member]
    record = _load_json_bytes(record_payload, label=member)
    _require(_json_bytes(record) == record_payload, f"{member}: JSON is not canonical")
    return manifest, record, payloads


def _require_recompute_manifest(
    manifest: Mapping[str, object],
    *,
    raw_inventory_sha256: str,
    analyzer_paths: Sequence[str],
) -> None:
    _require(manifest["campaign_id"] == admission.CAMPAIGN_ID, "recompute campaign ID drifted")
    _require(
        manifest["raw_artifact_inventory_sha256"] == raw_inventory_sha256,
        "recompute raw inventory binding drifted",
    )
    analyzer_bindings = manifest["analyzer_bindings"]
    _require(type(analyzer_bindings) is dict, "recompute analyzer bindings drifted")
    _require(
        set(analyzer_bindings) == set(analyzer_paths),
        "recompute analyzer binding set drifted",
    )
    for analyzer_path in analyzer_paths:
        source_path = _REPO_ROOT / analyzer_path
        _require(source_path.is_file(), f"required recompute analyzer is unavailable: {analyzer_path}")
        expected = hashlib.sha256(source_path.read_bytes()).hexdigest()
        _require(
            analyzer_bindings.get(analyzer_path) == expected,
            f"recompute analyzer binding drifted: {analyzer_path}",
        )


def _retained_source_checkpoint_lineage(
    admitted: Mapping[str, object],
    snapshot: immutable_orion_tree.VerifiedFrozenTree,
) -> dict[str, object]:
    """Project the admitted shared-MPI restart publication at t=500."""
    publications = admitted["restart_publications"]
    _require(type(publications) is dict, "admitted restart publications are unavailable")
    publication = _object(
        publications["500.0"],
        {"manifest_path", "layout", "member_count"},
        label="admitted restart publication t=500",
    )
    _require(
        publication["layout"] == "shared_mpi_io"
        and publication["member_count"] == 1,
        "admitted t=500 restart publication is not one shared-MPI checkpoint",
    )
    manifest_path = _relative_path(
        publication["manifest_path"],
        label="admitted restart publication t=500/manifest_path",
    )
    try:
        members = admission._parse_restart_manifest(
            snapshot.member_path(manifest_path).read_bytes(),
            manifest_path,
        )
        _require(
            len(members) == 1,
            "admitted t=500 restart manifest does not bind exactly one checkpoint",
        )
        member_path = _relative_path(
            members[0]["path"],
            label="admitted restart publication t=500/member_path",
        )
        checkpoint_payload = snapshot.member_path(member_path).read_bytes()
        semantics = admitted["retained_attempt_semantics"]
        _require(
            type(semantics) is dict,
            "admitted retained attempt semantics are unavailable",
        )
        contract = semantics["attempt_contract"]
        _require(
            type(contract) is dict,
            "admitted baseline attempt contract is unavailable",
        )
        attempt_root = Path(
            _text(
                contract["authorized_orion_attempt_root"],
                label="admitted baseline attempt contract/root",
            )
        )
    except (KeyError, OSError, ValueError) as error:
        if isinstance(error, NumericalQualificationError):
            raise
        raise NumericalQualificationError(
            f"cannot bind admitted t=500 restart publication: {error}"
        ) from error
    _require(
        attempt_root.is_absolute(),
        "admitted baseline attempt root is not absolute",
    )
    return {
        "snapshot_time_omega0_inverse": 500.0,
        "retained_attempt_id": admitted["run_identity"]["attempt_id"],
        "restart_manifest_path": manifest_path,
        "restart_member_path": member_path,
        "retained_restart_member_absolute_path": str(attempt_root / member_path),
        "restart_member_sha256": hashlib.sha256(checkpoint_payload).hexdigest(),
    }


def _reduce_retained_attempt(
    result: Mapping[str, object],
    snapshot: immutable_orion_tree.VerifiedFrozenTree,
) -> dict[str, object]:
    admitted = result["admission"]
    _require(type(admitted) is dict, "admitted attempt payload is unavailable")
    retained = admitted["retained_snapshot_products"]
    _require(type(retained) is dict, "admitted retained snapshot products are unavailable")
    particle_reductions = {}
    try:
        for time, evaluate_late_slope in (("500.0", False), ("1200.0", True)):
            product = retained[time]["prtcl_all"]
            decoded = admission.read_particle_vtk(snapshot.member_path(product["path"]))
            particle_reductions[f"t{time.removesuffix('.0')}"] = (
                particles.reduce_particle_snapshot(
                    snapshot_time=float(time),
                    points=decoded.points,
                    cr_source=decoded.scalars["cr_source"],
                    birth_time=decoded.scalars["birth_time"],
                    velocity=decoded.vectors["vel"],
                    macro_weight=decoded.scalars["macro_weight"],
                    evaluate_late_slope=evaluate_late_slope,
                )
            )
        datasets = {}
        for quantity in spatial.REQUIRED_MESH_QUANTITIES:
            product = retained["500.0"][quantity]
            datasets[quantity] = admission.output_primitives.parse_athenak_binary_bytes(
                snapshot.member_path(product["path"]).read_bytes(),
                source=product["path"],
            )
        spatial_reduction = spatial.reduce_t500_spatial_snapshot(
            datasets,
            x_ideal_c_over_omega_pi=particles.ideal_surface_x1(500.0),
        )
    except (KeyError, OSError, ValueError) as error:
        raise NumericalQualificationError(
            f"retained raw attempt recompute failed: {error}"
        ) from error
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": ATTEMPT_RECORD_TYPE,
        "run_identity": copy.deepcopy(admitted["run_identity"]),
        "raw_inventory_sha256": admitted["immutable_tree"]["inventory_sha256"],
        "source_checkpoint_lineage": _retained_source_checkpoint_lineage(
            admitted, snapshot
        ),
        "admission_result": copy.deepcopy(dict(result)),
        "particle_reductions": particle_reductions,
        "spatial_reduction": spatial_reduction,
    }


def bind_retained_attempt_tree(
    campaign_root: str | Path,
    expected_inventory_sha256: str,
    *,
    authorized_orion_root: Path = admission.ORION_BULK_ROOT,
) -> _RetainedAttemptBinding:
    """Bind one admitted immutable raw tree and recompute its raw numerical reducers."""
    result = admission.qualify_campaign(
        campaign_root,
        expected_inventory_sha256,
        authorized_orion_root=authorized_orion_root,
    )
    _require(
        result["admitted_for_follow_on_numerical_qualification"] is True,
        "retained attempt was not admitted for numerical qualification",
    )
    try:
        with immutable_orion_tree.staged_verified_frozen_tree(
            campaign_root,
            expected_inventory_sha256,
            authorized_root=authorized_orion_root,
            error_type=NumericalQualificationError,
            label="Q-011 retained numerical attempt tree",
        ) as (_, snapshot):
            attempt = _reduce_retained_attempt(result, snapshot)
    except (OSError, ValueError) as error:
        if isinstance(error, NumericalQualificationError):
            raise
        raise NumericalQualificationError(f"retained attempt tree binding failed: {error}") from error
    retained = copy.deepcopy(attempt)
    wrapper = {
        "attempt_sha256": _canonical_sha256(retained),
        "attempt": retained,
    }
    return _RetainedAttemptBinding(
        wrapper_payload=_json_bytes(wrapper),
        _key=_PRODUCTION_BINDING_KEY,
    )


def _retained_attempt_payload(
    value: _RetainedAttemptBinding, *, index: int
) -> tuple[dict[str, object], dict[str, object]]:
    parsed = _validate_attempt_binding(value, index=index, allow_unit_only=False)
    wrapper = _object(
        _load_json_bytes(value.wrapper_payload, label=f"attempts[{index}]/retained_binding"),
        {"attempt_sha256", "attempt"},
        label=f"attempts[{index}]",
    )
    _require(
        wrapper["attempt_sha256"] == parsed["attempt_sha256"],
        f"attempts[{index}]: retained attempt payload drifted",
    )
    _require(type(wrapper["attempt"]) is dict, f"attempts[{index}]: retained attempt is unavailable")
    return parsed, wrapper["attempt"]


def _finite_array(value: object, *, label: str, ndim: int = 1) -> np.ndarray:
    try:
        result = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise NumericalQualificationError(f"{label}: expected numeric array") from error
    _require(result.ndim == ndim and result.size > 0, f"{label}: array shape drifted")
    _require(np.all(np.isfinite(result)), f"{label}: array contains non-finite values")
    return result


def _dx12_profile(
    value: object, *, quantity: str, label: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    record = _required_object(
        value,
        {
            "quantity",
            "source_field",
            "time_omega0_inverse",
            "x1_centers_c_over_omega_pi",
            "values_x",
            "column_areas",
        },
        label=label,
    )
    _require(
        record["quantity"] == quantity
        and record["source_field"] == spatial.MESH_QUANTITY_FIELDS[quantity]
        and record["time_omega0_inverse"] == spatial.T500_OMEGA0_INVERSE,
        f"{label}: profile identity drifted",
    )
    centers = _finite_array(record["x1_centers_c_over_omega_pi"], label=f"{label}/x1")
    profile = _finite_array(record["values_x"], label=f"{label}/values")
    areas = _finite_array(record["column_areas"], label=f"{label}/areas")
    _require(
        centers.size == profile.size == areas.size and centers.size >= 2,
        f"{label}: profile lengths drifted",
    )
    _require(np.all(areas > 0.0), f"{label}: profile areas must be positive")
    spacing = np.diff(centers)
    _require(np.all(spacing > 0.0), f"{label}: profile centers must increase")
    source_dx = float(spacing[0])
    _require(
        np.allclose(spacing, source_dx, rtol=0.0, atol=1.0e-10),
        f"{label}: profile spacing is not uniform",
    )
    ratio = _PAIR_PROFILE_DX / source_dx
    group = round(ratio)
    _require(
        group >= 1 and math.isclose(ratio, group, rel_tol=0.0, abs_tol=1.0e-10),
        f"{label}: profile spacing cannot be conservatively restricted to dx=12",
    )
    _require(profile.size % group == 0, f"{label}: dx=12 profile grouping drifted")
    grouped_areas = np.sum(areas.reshape(-1, group), axis=1)
    grouped_profile = np.sum((profile * areas).reshape(-1, group), axis=1) / grouped_areas
    grouped_centers = np.mean(centers.reshape(-1, group), axis=1)
    return grouped_centers, grouped_profile, grouped_areas


def _matched_residual(
    reference: object,
    candidate: object,
    *,
    label: str,
    cell_areas: object | None = None,
) -> dict[str, object]:
    try:
        residual = admission.output_primitives.matched_grid_residual_metrics(
            reference, candidate, cell_areas=cell_areas
        )
    except ValueError as error:
        raise NumericalQualificationError(f"{label}: residual recompute failed: {error}") from error
    return {
        "observable": label,
        "maximum_absolute_difference": residual.maximum_absolute,
        "relative_mean_absolute": residual.relative_mean_absolute,
        "relative_root_mean_square": residual.relative_root_mean_square,
    }


def _spectrum_window(
    value: object, *, upper: float, label: str
) -> np.ndarray:
    record = _required_object(
        value, {"bin_edges", "normalized_chi_f_chi"}, label=label
    )
    edges = _finite_array(record["bin_edges"], label=f"{label}/bin_edges")
    spectrum = _finite_array(
        record["normalized_chi_f_chi"], label=f"{label}/normalized_chi_f_chi"
    )
    _require(
        edges.size == spectrum.size + 1 and np.all(edges > 0.0) and np.all(np.diff(edges) > 0.0),
        f"{label}: spectrum bin layout drifted",
    )
    _require(
        np.array_equal(edges, np.asarray(particles.CHI_BIN_EDGES, dtype=np.float64)),
        f"{label}: spectrum bins differ from the frozen chi grid",
    )
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= 10.0) & (centers <= upper)
    _require(np.count_nonzero(selected) > 0, f"{label}: comparison window is empty")
    return spectrum[selected]


def _recompute_pair_residuals(
    amr_attempt: Mapping[str, object],
    fine_attempt: Mapping[str, object],
) -> list[dict[str, object]]:
    try:
        amr_spatial = amr_attempt["spatial_reduction"]
        fine_spatial = fine_attempt["spatial_reduction"]
        amr_front = amr_spatial["detected_front"]["x_front_c_over_omega_pi"]
        fine_front = fine_spatial["detected_front"]["x_front_c_over_omega_pi"]
        amr_amplification = amr_spatial["upstream_b_amplification"]["amplification_over_b0"]
        fine_amplification = fine_spatial["upstream_b_amplification"]["amplification_over_b0"]
        residuals = [
            {
                "observable": "shock_front_position_at_t500",
                "maximum_absolute_difference": float(abs(amr_front - fine_front)),
                "relative_mean_absolute": None,
                "relative_root_mean_square": None,
            },
            {
                "observable": "upstream_magnetic_amplification_at_t500",
                "maximum_absolute_difference": float(abs(amr_amplification - fine_amplification)),
                "relative_mean_absolute": None,
                "relative_root_mean_square": None,
            },
        ]
        for quantity in ("rho", "bmag"):
            observable = f"{quantity}_y_average_at_t500"
            amr_centers, amr_profile, amr_areas = _dx12_profile(
                amr_spatial["y_area_weighted_profiles"][quantity],
                quantity=quantity,
                label=f"AMR/{observable}",
            )
            fine_centers, fine_profile, fine_areas = _dx12_profile(
                fine_spatial["y_area_weighted_profiles"][quantity],
                quantity=quantity,
                label=f"fine_uniform/{observable}",
            )
            _require(
                np.allclose(amr_centers, fine_centers, rtol=0.0, atol=1.0e-10)
                and np.allclose(amr_areas, fine_areas, rtol=0.0, atol=1.0e-10),
                f"{observable}: restricted profile areas drifted",
            )
            residuals.append(
                _matched_residual(
                    amr_profile,
                    fine_profile,
                    label=observable,
                    cell_areas=amr_areas,
                )
            )
        for time, upper in (("t500", 256.0), ("t1200", 512.0)):
            observable = f"normalized_downstream_chi_f_chi_at_{time}"
            amr_spectrum = _spectrum_window(
                amr_attempt["particle_reductions"][time]["weighted_spectrum"],
                upper=upper,
                label=f"AMR/{observable}",
            )
            fine_spectrum = _spectrum_window(
                fine_attempt["particle_reductions"][time]["weighted_spectrum"],
                upper=upper,
                label=f"fine_uniform/{observable}",
            )
            residuals.append(_matched_residual(amr_spectrum, fine_spectrum, label=observable))
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, NumericalQualificationError):
            raise
        raise NumericalQualificationError(f"retained pair residual recompute failed: {error}") from error
    return residuals


def bind_retained_pair_result(
    bundle_root: str | Path,
    expected_inventory_sha256: str,
    *,
    amr_attempt: _RetainedAttemptBinding,
    fine_uniform_attempt: _RetainedAttemptBinding,
    authorized_orion_root: Path,
) -> _RetainedPairResultBinding:
    """Recompute and bind one immutable retained AMR/fine residual receipt."""
    amr, amr_payload = _retained_attempt_payload(amr_attempt, index=0)
    fine, fine_payload = _retained_attempt_payload(fine_uniform_attempt, index=1)
    _require(
        amr["identity"]["variant"] == AMR_VARIANT
        and fine["identity"]["variant"] == FINE_VARIANT
        and amr["identity"]["seed"] == fine["identity"]["seed"],
        "retained pair attempts are not one matched AMR/fine seed",
    )
    manifest, value, _ = _load_retained_recompute_record(
        bundle_root,
        expected_inventory_sha256,
        member=PAIR_RECOMPUTE_MEMBER,
        authorized_orion_root=authorized_orion_root,
    )
    _require_recompute_manifest(
        manifest,
        raw_inventory_sha256=amr["raw_inventory_sha256"],
        analyzer_paths=PAIR_RECOMPUTE_ANALYZERS,
    )
    record = _object(
        value,
        {
            "schema_version",
            "record_type",
            "seed",
            "amr_attempt_sha256",
            "amr_raw_inventory_sha256",
            "fine_uniform_attempt_sha256",
            "fine_uniform_raw_inventory_sha256",
            "residuals",
        },
        label=PAIR_RECOMPUTE_MEMBER,
    )
    _require(
        type(record["schema_version"]) is int
        and record["schema_version"] == SCHEMA_VERSION
        and record["record_type"] == PAIR_RECOMPUTE_RECORD_TYPE,
        "retained pair recompute identity drifted",
    )
    _require(
        record["seed"] == amr["identity"]["seed"]
        and record["amr_attempt_sha256"] == amr["attempt_sha256"]
        and record["amr_raw_inventory_sha256"] == amr["raw_inventory_sha256"]
        and record["fine_uniform_attempt_sha256"] == fine["attempt_sha256"]
        and record["fine_uniform_raw_inventory_sha256"] == fine["raw_inventory_sha256"],
        "retained pair recompute raw attempt binding drifted",
    )
    recomputed = _recompute_pair_residuals(amr_payload, fine_payload)
    _require(
        record["residuals"] == recomputed,
        "retained pair residual values differ from trusted raw-attempt recompute",
    )
    record["residuals"] = recomputed
    return _RetainedPairResultBinding(
        record_payload=_json_bytes(copy.deepcopy(record)),
        _key=_PRODUCTION_BINDING_KEY,
    )


def _bundle_member(payloads: Mapping[str, bytes], value: object, *, label: str) -> bytes:
    relative = _relative_path(value, label=label)
    _require(relative in payloads, f"{label}: retained bundle member is missing")
    return payloads[relative]


def _text(value: object, *, label: str) -> str:
    _require(type(value) is str and bool(value), f"{label}: expected non-empty text")
    return value


def _digest_path_binding(value: object, *, label: str) -> dict[str, str]:
    item = _object(value, {"path", "sha256"}, label=label)
    return {
        "path": _text(item["path"], label=f"{label}/path"),
        "sha256": _sha256(item["sha256"], label=f"{label}/sha256"),
    }


def _relative_digest_binding(value: object, *, label: str) -> dict[str, str]:
    binding = _digest_path_binding(value, label=label)
    binding["path"] = _relative_path(binding["path"], label=f"{label}/path")
    return binding


def _tree_json_member(
    payloads: Mapping[str, bytes],
    value: object,
    *,
    label: str,
) -> dict[str, object]:
    binding = _relative_digest_binding(value, label=label)
    _require(binding["path"] in payloads, f"{label}: planner tree member is missing")
    payload = payloads[binding["path"]]
    _require(
        hashlib.sha256(payload).hexdigest() == binding["sha256"],
        f"{label}: planner tree member SHA-256 drifted",
    )
    item = _load_json_bytes(payload, label=label)
    _require(type(item) is dict, f"{label}: expected object")
    return item


def _validate_planned_restart_carrier(
    record: Mapping[str, object],
    payloads: Mapping[str, bytes],
    *,
    source: Mapping[str, object],
    source_attempt_payload: Mapping[str, object],
    authorized_orion_root: Path,
) -> dict[str, str]:
    """Derive the restart carrier from the frozen source campaign plan."""
    receipt_member = _relative_path(
        record["planner_materialization_receipt_member"],
        label=f"{RESTART_RECOMPUTE_MEMBER}/planner_materialization_receipt_member",
    )
    campaign_plan_member = _relative_path(
        record["retained_campaign_plan_member"],
        label=f"{RESTART_RECOMPUTE_MEMBER}/retained_campaign_plan_member",
    )
    receipt_payload = _bundle_member(
        payloads,
        receipt_member,
        label=f"{RESTART_RECOMPUTE_MEMBER}/planner_materialization_receipt_member",
    )
    campaign_plan_payload = _bundle_member(
        payloads,
        campaign_plan_member,
        label=f"{RESTART_RECOMPUTE_MEMBER}/retained_campaign_plan_member",
    )
    campaign_plan_sha256 = hashlib.sha256(campaign_plan_payload).hexdigest()
    try:
        admitted_bindings = source_attempt_payload["admission_result"]["admission"][
            "artifact_bindings"
        ]
        _require(
            campaign_plan_sha256 == admitted_bindings["campaign_plan"]["sha256"]
            and hashlib.sha256(receipt_payload).hexdigest()
            == admitted_bindings["planner_materialization_receipt"]["sha256"],
            "retained restart planner bytes differ from admitted source attempt",
        )
    except (KeyError, TypeError) as error:
        raise NumericalQualificationError(
            "admitted source-attempt planner bindings are unavailable"
        ) from error
    try:
        planner_receipt = admission._validate_planner_materialization_receipt(
            receipt_payload,
            retained_campaign_plan_payload=campaign_plan_payload,
            retained_campaign_plan_sha256=campaign_plan_sha256,
            authorized_pic_root=authorized_orion_root,
        )
    except ValueError as error:
        raise NumericalQualificationError(
            "retained restart planner materialization receipt failed validation"
        ) from error
    tree_payloads = planner_receipt.get("_tree_member_payloads")
    _require(
        type(tree_payloads) is dict,
        "retained restart planner tree payloads are unavailable",
    )
    campaign_plan = _required_object(
        _load_json_bytes(campaign_plan_payload, label=campaign_plan_member),
        {
            "authorized_orion_campaign_root",
            "baseline_attempt_descriptors",
            "campaign_matrix",
            "candidate_binding",
            "selected_pressure",
            "source_bindings",
            "restart_continuation_carrier",
        },
        label=campaign_plan_member,
    )
    identity = source["identity"]
    _require(type(identity) is dict, "restart source attempt identity is unavailable")
    matrix = _required_object(
        campaign_plan["campaign_matrix"],
        {"qualifying_seeds"},
        label=f"{campaign_plan_member}/campaign_matrix",
    )
    seeds = _list(
        matrix["qualifying_seeds"],
        label=f"{campaign_plan_member}/campaign_matrix/qualifying_seeds",
    )
    _require(
        seeds
        and identity["variant"] == AMR_VARIANT
        and identity["seed"] == seeds[0],
        "restart source attempt is not the canonical source-plan AMR baseline",
    )
    selected_pressure = _required_object(
        campaign_plan["selected_pressure"],
        {"selected_case"},
        label=f"{campaign_plan_member}/selected_pressure",
    )
    selected_case = _required_object(
        selected_pressure["selected_case"],
        {"problem_ps_p0"},
        label=f"{campaign_plan_member}/selected_pressure/selected_case",
    )
    descriptor_matches = []
    for index, raw_binding in enumerate(
        _list(
            campaign_plan["baseline_attempt_descriptors"],
            label=f"{campaign_plan_member}/baseline_attempt_descriptors",
        )
    ):
        descriptor = _tree_json_member(
            tree_payloads,
            raw_binding,
            label=f"source plan baseline descriptor[{index}]",
        )
        if descriptor.get("attempt_id") == identity["attempt_id"]:
            descriptor_matches.append(descriptor)
    _require(
        len(descriptor_matches) == 1,
        "canonical source-plan baseline descriptor is missing or ambiguous",
    )
    source_descriptor = _required_object(
        descriptor_matches[0],
        {
            "attempt_id",
            "variant",
            "qualifying_seed",
            "selected_problem_ps_p0",
        },
        label="canonical source-plan baseline descriptor",
    )
    _require(
        source_descriptor["attempt_id"] == identity["attempt_id"]
        and source_descriptor["variant"] == identity["variant"] == AMR_VARIANT
        and source_descriptor["qualifying_seed"] == identity["seed"] == seeds[0]
        and source_descriptor["selected_problem_ps_p0"]
        == selected_case["problem_ps_p0"],
        "canonical source-plan baseline descriptor identity drifted",
    )
    source_bindings = _required_object(
        campaign_plan["source_bindings"],
        {"restart_preregistration", "paper_deck"},
        label=f"{campaign_plan_member}/source_bindings",
    )
    restart_preregistration = _relative_digest_binding(
        source_bindings["restart_preregistration"],
        label=f"{campaign_plan_member}/source_bindings/restart_preregistration",
    )
    _require(
        restart_preregistration["sha256"]
        == admission.EXPECTED_RESTART_PREREGISTRATION_SHA256,
        "retained restart preregistration binding drifted",
    )
    paper_deck = _relative_digest_binding(
        source_bindings["paper_deck"],
        label=f"{campaign_plan_member}/source_bindings/paper_deck",
    )
    candidate = _required_object(
        campaign_plan["candidate_binding"],
        {"clean_candidate_manifest", "executable", "environment_profile", "git_commit"},
        label=f"{campaign_plan_member}/candidate_binding",
    )
    clean_candidate = _digest_path_binding(
        candidate["clean_candidate_manifest"],
        label=f"{campaign_plan_member}/candidate_binding/clean_candidate_manifest",
    )
    executable = _digest_path_binding(
        candidate["executable"],
        label=f"{campaign_plan_member}/candidate_binding/executable",
    )
    environment = _required_object(
        candidate["environment_profile"],
        {"path", "sha256", "control_plane_version", "reviewed_source"},
        label=f"{campaign_plan_member}/candidate_binding/environment_profile",
    )
    environment_sha256 = _sha256(
        environment["sha256"],
        label=f"{campaign_plan_member}/candidate_binding/environment_profile/sha256",
    )
    control_plane_version = _sha256(
        environment["control_plane_version"],
        label=(
            f"{campaign_plan_member}/candidate_binding/"
            "environment_profile/control_plane_version"
        ),
    )
    git_commit = _git_commit(
        candidate["git_commit"],
        label=f"{campaign_plan_member}/candidate_binding/git_commit",
    )
    campaign_root = Path(
        _text(
            campaign_plan["authorized_orion_campaign_root"],
            label=f"{campaign_plan_member}/authorized_orion_campaign_root",
        )
    )
    _require(campaign_root.is_absolute(), "source plan campaign root is not absolute")
    carrier_id = admission.campaign_planner._restart_carrier_id(identity["seed"])
    restart_artifact_root = campaign_root / "restart_continuation" / carrier_id
    try:
        policy = restart.load_preregistration()
        expected_contract = admission.campaign_planner._restart_launch_contract(
            carrier_id=carrier_id,
            source_attempt=source_descriptor,
            restart_preregistration=policy,
            candidate=candidate,
            paper_deck_binding=paper_deck,
            artifact_root=restart_artifact_root,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise NumericalQualificationError(
            "retained restart planned carrier reconstruction failed"
        ) from error
    expected_contract_payload = admission.campaign_planner._json_bytes(expected_contract)
    launch_binding = {
        "path": f"launch_contracts/restart_continuation/{carrier_id}.json",
        "sha256": hashlib.sha256(expected_contract_payload).hexdigest(),
    }
    launch_contract = _tree_json_member(
        tree_payloads, launch_binding, label="planned restart launch contract"
    )
    _require(
        launch_contract == expected_contract,
        "planned restart launch contract drifted",
    )
    continuation = policy["continuation_contract"]
    expected_carrier = {
        "record_type": admission.campaign_planner.RESTART_CARRIER_RECORD_TYPE,
        "schema_version": 1,
        "carrier_id": carrier_id,
        "status": "planned_not_authorized",
        "source_baseline_attempt_id": identity["attempt_id"],
        "variant": identity["variant"],
        "qualifying_seed": identity["seed"],
        "selected_problem_ps_p0": selected_case["problem_ps_p0"],
        "authorized_orion_attempt_root": str(restart_artifact_root),
        "restart_preregistration": restart_preregistration,
        "checkpoint_time_omega0_inverse": continuation[
            "checkpoint_time_omega0_inverse"
        ],
        "retained_output_schedule_after_checkpoint_omega0_inverse": continuation[
            "retained_output_schedule_after_checkpoint_omega0_inverse"
        ],
        "comparison_tolerances_max_absolute_difference": continuation[
            "comparison_tolerances_max_absolute_difference"
        ],
        "launch_contract": launch_binding,
    }
    expected_carrier_payload = admission.campaign_planner._json_bytes(
        expected_carrier
    )
    expected_planned_binding = {
        "path": PLANNED_RESTART_CARRIER_PATH,
        "sha256": hashlib.sha256(expected_carrier_payload).hexdigest(),
    }
    planned_binding = _relative_digest_binding(
        record["planned_restart_carrier"],
        label=f"{RESTART_RECOMPUTE_MEMBER}/planned_restart_carrier",
    )
    _require(
        planned_binding == expected_planned_binding
        and campaign_plan["restart_continuation_carrier"] == expected_planned_binding,
        "retained restart planned carrier binding drifted",
    )
    carrier = _tree_json_member(
        tree_payloads,
        expected_planned_binding,
        label="planned restart carrier",
    )
    _require(carrier == expected_carrier, "planned restart carrier content drifted")
    return {
        "planned_restart_carrier_sha256": planned_binding["sha256"],
        "clean_candidate_manifest_sha256": clean_candidate["sha256"],
        "executable_sha256": executable["sha256"],
        "paper_deck_sha256": paper_deck["sha256"],
        "environment_profile_sha256": environment_sha256,
        "control_plane_version": control_plane_version,
        "git_commit": git_commit,
        "restart_carrier_id": carrier_id,
        "restart_authorized_orion_attempt_root": str(restart_artifact_root),
    }


def _directory_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _require_same_open_directory(
    path: Path,
    descriptor: int,
    expected: tuple[int, ...],
    *,
    label: str,
) -> None:
    lexical = os.open(path, _DIRECTORY_FLAGS)
    try:
        retained = os.fstat(descriptor)
        actual = os.fstat(lexical)
        _require(
            _directory_identity(retained)
            == _directory_identity(actual)
            == expected,
            f"{label}: registered execution receipt ancestry changed",
        )
    finally:
        os.close(lexical)


def _watch_open_directories(descriptors: Sequence[int], *, label: str) -> int:
    libc = ctypes.CDLL(None, use_errno=True)
    watch_descriptor = libc.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
    if watch_descriptor < 0:
        error = ctypes.get_errno()
        raise NumericalQualificationError(
            f"{label}: cannot watch registered execution receipt ancestry: "
            f"{os.strerror(error)}"
        )
    try:
        for descriptor in descriptors:
            result = libc.inotify_add_watch(
                watch_descriptor,
                f"/proc/self/fd/{descriptor}".encode("ascii"),
                _INOTIFY_MUTATION_MASK,
            )
            if result < 0:
                error = ctypes.get_errno()
                raise NumericalQualificationError(
                    f"{label}: cannot watch registered execution receipt ancestry: "
                    f"{os.strerror(error)}"
                )
        return watch_descriptor
    except BaseException:
        os.close(watch_descriptor)
        raise


def _require_no_ancestry_events(descriptor: int, *, label: str) -> None:
    try:
        payload = os.read(descriptor, 1024 * 1024)
    except BlockingIOError:
        return
    _require(
        not payload,
        f"{label}: registered execution receipt ancestry changed while pinned",
    )


@contextmanager
def _pinned_authoritative_execution_receipt(
    value: object,
    *,
    authorized_orion_root: Path,
    label: str,
) -> Iterator[tuple[Path, bytes, int]]:
    supplied = Path(_text(value, label=label))
    root = Path(os.path.abspath(authorized_orion_root))
    _require(
        supplied.is_absolute()
        and Path(os.path.abspath(supplied)) == supplied
        and supplied.name == admission.REGISTERED_EXECUTION_RECEIPT_NAME
        and supplied.parent.name == "analysis",
        f"{label}: registered execution receipt fixed path drifted",
    )
    try:
        relative = supplied.relative_to(root)
    except ValueError as error:
        raise NumericalQualificationError(
            f"{label}: registered execution receipt fixed path is unavailable"
        ) from error
    _require(
        relative.parts[:1] == ("runs",) and len(relative.parts) >= 5,
        f"{label}: registered execution receipt path is outside the run namespace",
    )
    ancestry = [os.open(root, _DIRECTORY_FLAGS)]
    descriptor: int | None = None
    try:
        for part in relative.parts[:-1]:
            ancestry.append(os.open(part, _DIRECTORY_FLAGS, dir_fd=ancestry[-1]))
        descriptor = os.open(relative.parts[-1], _REGULAR_FLAGS, dir_fd=ancestry[-1])
    except OSError as error:
        for ancestor in reversed(ancestry):
            os.close(ancestor)
        raise NumericalQualificationError(
            f"{label}: cannot open registered execution receipt authority"
        ) from error
    ancestry_identities = [_directory_identity(os.fstat(ancestor)) for ancestor in ancestry]
    watch_descriptor = _watch_open_directories(ancestry, label=label)

    def require_same() -> None:
        _require_no_ancestry_events(watch_descriptor, label=label)
        current = os.stat(relative.parts[-1], dir_fd=ancestry[-1], follow_symlinks=False)
        opened = os.fstat(descriptor)
        _require(
            (opened.st_dev, opened.st_ino) == (current.st_dev, current.st_ino),
            f"{label}: registered execution receipt authority changed while pinned",
        )
        _require_same_open_directory(
            root, ancestry[0], ancestry_identities[0], label=label
        )
        current_path = root
        for part, ancestor, identity in zip(
            relative.parts[:-1], ancestry[1:], ancestry_identities[1:]
        ):
            current_path /= part
            _require_same_open_directory(
                current_path, ancestor, identity, label=label
            )

    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & _WRITE_BITS,
            f"{label}: registered execution receipt authority is not immutable",
        )
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        _require(
            _regular_identity(before) == _regular_identity(after),
            f"{label}: registered execution receipt authority changed while reading",
        )
        require_same()
        artifact_root_descriptor = ancestry[-2]
        yield supplied, bytes(payload), artifact_root_descriptor
        require_same()
    finally:
        os.close(watch_descriptor)
        os.close(descriptor)
        for ancestor in reversed(ancestry):
            os.close(ancestor)


def _stable_authoritative_execution_receipt(
    value: object,
    *,
    authorized_orion_root: Path,
    label: str,
) -> tuple[Path, bytes]:
    with _pinned_authoritative_execution_receipt(
        value,
        authorized_orion_root=authorized_orion_root,
        label=label,
    ) as (path, payload, _):
        return path, payload


@contextmanager
def _validated_authoritative_restart_execution_receipt(
    branch: Mapping[str, object],
    payloads: Mapping[str, bytes],
    *,
    branch_role: str,
    source: Mapping[str, object],
    source_attempt_payload: Mapping[str, object],
    source_checkpoint_lineage: Mapping[str, object],
    planned: Mapping[str, str],
    authorized_orion_root: Path,
    label: str,
) -> Iterator[tuple[dict[str, object], int]]:
    member = _relative_path(
        branch["execution_receipt_member"], label=f"{label}/execution_receipt_member"
    )
    retained_payload = _bundle_member(
        payloads, member, label=f"{label}/execution_receipt_member"
    )
    admitted = source_attempt_payload["admission_result"]["admission"]
    semantics = admitted["retained_attempt_semantics"]
    baseline_receipt = semantics["registered_execution_receipt"]
    baseline_contract = semantics["attempt_contract"]
    artifact_bindings = admitted["artifact_bindings"]
    if branch_role == "uninterrupted_baseline":
        expected_attempt_id = source["identity"]["attempt_id"]
        expected_argv = baseline_contract["argv"]
        expected_raw_root = baseline_receipt["raw_output_root"]
        expected_artifact_dir = baseline_receipt["artifact_dir"]
        expected_authority = (
            Path(expected_artifact_dir)
            / "analysis"
            / admission.REGISTERED_EXECUTION_RECEIPT_NAME
        )
    else:
        expected_attempt_id = planned["restart_carrier_id"]
        expected_raw_root = f"{planned['restart_authorized_orion_attempt_root']}/raw"
        expected_argv = [
            "-r",
            source_checkpoint_lineage["retained_restart_member_absolute_path"],
            "-d",
            expected_raw_root,
        ]
        expected_artifact_dir = None
        expected_authority = None
    authority = _pinned_authoritative_execution_receipt(
        branch["execution_receipt_authority_path"],
        authorized_orion_root=authorized_orion_root,
        label=f"{label}/execution_receipt_authority_path",
    )
    with authority as (authority_path, authority_payload, artifact_root_descriptor):
        _require(
            retained_payload == authority_payload,
            f"{label}: retained receipt differs from immutable fixed-path authority",
        )
        if expected_authority is not None:
            _require(
                authority_path == expected_authority,
                f"{label}: baseline receipt authority differs from admitted source attempt",
            )
        try:
            receipt = admission._validate_registered_execution_receipt(
                authority_payload,
                expected_attempt_id=expected_attempt_id,
                expected_source_commit=planned["git_commit"],
                expected_executable_sha256=planned["executable_sha256"],
                expected_deck_sha256=planned["paper_deck_sha256"],
                expected_environment_sha256=planned["environment_profile_sha256"],
                expected_control_plane_version=planned["control_plane_version"],
                expected_argv=expected_argv,
                expected_raw_output_root=expected_raw_root,
                expected_artifact_dir=expected_artifact_dir,
            )
        except ValueError as error:
            raise NumericalQualificationError(
                f"{label}: authoritative registered execution receipt failed validation"
            ) from error
        _require(
            authority_path
            == Path(receipt["artifact_dir"])
            / "analysis"
            / admission.REGISTERED_EXECUTION_RECEIPT_NAME,
            f"{label}: authoritative receipt is not at its immutable fixed path",
        )
        try:
            with admission._validated_registered_execution_receipt_ledger_snapshot(
                receipt,
                authorized_pic_root=authorized_orion_root,
            ) as ledger_binding:
                if branch_role == "uninterrupted_baseline":
                    _require(
                        receipt == baseline_receipt
                        and ledger_binding
                        == semantics["registered_execution_ledger_binding"]
                        and artifact_bindings["executable"]["sha256"]
                        == planned["executable_sha256"]
                        and artifact_bindings["deck"]["sha256"]
                        == planned["paper_deck_sha256"],
                        f"{label}: baseline receipt differs from admitted source authority",
                    )
                yield receipt, artifact_root_descriptor
        except ValueError as error:
            if isinstance(error, NumericalQualificationError):
                raise
            raise NumericalQualificationError(
                f"{label}: authoritative registered execution receipt failed validation"
            ) from error


def _validate_authoritative_restart_execution_receipt(
    branch: Mapping[str, object],
    payloads: Mapping[str, bytes],
    **kwargs: object,
) -> dict[str, object]:
    with _validated_authoritative_restart_execution_receipt(
        branch, payloads, **kwargs
    ) as (receipt, _):
        return receipt


def _restart_branch_member_paths(value: object, *, label: str) -> set[str]:
    branch = _object(
        value,
        {
            "branch_role",
            "execution_receipt_authority_path",
            "execution_receipt_member",
            "raw_output_root",
            "structured_artifact_inventory_sha256",
            "outputs_after_checkpoint",
        },
        label=label,
    )
    paths = {
        _relative_path(
            branch["execution_receipt_member"],
            label=f"{label}/execution_receipt_member",
        )
    }
    outputs = _list(
        branch["outputs_after_checkpoint"],
        label=f"{label}/outputs_after_checkpoint",
    )
    _text(branch["raw_output_root"], label=f"{label}/raw_output_root")
    _sha256(
        branch["structured_artifact_inventory_sha256"],
        label=f"{label}/structured_artifact_inventory_sha256",
    )
    authority_paths = set()
    for index, output in enumerate(outputs):
        output_label = f"{label}/outputs_after_checkpoint[{index}]"
        item = _object(output, {"time_omega0_inverse", "members"}, label=output_label)
        members = _object(
            item["members"],
            {*_RESTART_MESH_MEMBERS, _RESTART_PARTICLE_MEMBER},
            label=f"{output_label}/members",
        )
        for name, member in members.items():
            member_label = f"{output_label}/members/{name}"
            binding = _object(
                member,
                {"bundle_member", "structured_artifact_member", "sha256"},
                label=member_label,
            )
            relative = _relative_path(
                binding["bundle_member"], label=f"{member_label}/bundle_member"
            )
            authority = _relative_path(
                binding["structured_artifact_member"],
                label=f"{member_label}/structured_artifact_member",
            )
            _sha256(binding["sha256"], label=f"{member_label}/sha256")
            _require(
                relative not in paths,
                f"{label}: retained member aliases another branch member",
            )
            _require(
                authority not in authority_paths,
                f"{label}: structured raw authority member aliases another output",
            )
            paths.add(relative)
            authority_paths.add(authority)
    return paths


def _structured_restart_member(
    payloads: Mapping[str, bytes],
    tree: structured_artifacts.StructuredArtifactTree,
    inventory: dict[str, dict[str, object]],
    value: object,
    *,
    label: str,
) -> bytes:
    binding = _object(
        value,
        {"bundle_member", "structured_artifact_member", "sha256"},
        label=label,
    )
    bundle_member = _relative_path(
        binding["bundle_member"], label=f"{label}/bundle_member"
    )
    authority_member = _relative_path(
        binding["structured_artifact_member"],
        label=f"{label}/structured_artifact_member",
    )
    digest = _sha256(binding["sha256"], label=f"{label}/sha256")
    retained_payload = _bundle_member(payloads, bundle_member, label=label)
    try:
        authoritative_payload = structured_artifacts.read_inventory_bytes(
            tree, inventory, authority_member
        )
    except ValueError as error:
        raise NumericalQualificationError(
            f"{label}: structured raw authority validation failed"
        ) from error
    _require(
        retained_payload == authoritative_payload
        and hashlib.sha256(authoritative_payload).hexdigest() == digest,
        f"{label}: retained copy differs from structured raw authority",
    )
    return authoritative_payload


def _restart_mesh_values(payload: bytes, *, field: str, time: float, label: str) -> list[float]:
    try:
        dataset = admission.output_primitives.parse_athenak_binary_bytes(payload, source=label)
    except ValueError as error:
        raise NumericalQualificationError(f"{label}: retained mesh decode failed: {error}") from error
    _require(dataset.variable_names == (field,), f"{label}: retained mesh field inventory drifted")
    _require(dataset.time == time, f"{label}: retained mesh time drifted")
    ordered = sorted(
        dataset.blocks,
        key=lambda block: (block.level, block.logical_location, block.index_bounds),
    )
    metadata = [
        float(dataset.cycle),
        *(float(value) for value in dataset.root_grid_shape),
        *(float(value) for value in dataset.meshblock_shape),
        *(float(value) for value in dataset.domain_bounds),
    ]
    return [
        *metadata,
        *(
            float(value)
            for block in ordered
            for value in (
                block.level,
                *block.logical_location,
                *block.index_bounds,
                *block.geometry,
                *block.fields[field].reshape(-1),
            )
        ),
    ]


def _restart_particle_values(
    payload: bytes, *, time: float, label: str
) -> tuple[list[int], list[float]]:
    descriptor = os.memfd_create("q011-section54-retained-pvtk", flags=os.MFD_CLOEXEC)
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            _require(written > 0, f"{label}: short write while sealing particle payload")
            view = view[written:]
        os.lseek(descriptor, 0, os.SEEK_SET)
        decoded = admission.read_particle_vtk(Path(f"/proc/self/fd/{descriptor}"))
    except (OSError, ValueError) as error:
        raise NumericalQualificationError(f"{label}: retained particle decode failed: {error}") from error
    finally:
        os.close(descriptor)
    _require(
        set(decoded.scalars)
        == set(_RESTART_PARTICLE_INTEGER_SCALARS + _RESTART_PARTICLE_REAL_SCALARS)
        and set(decoded.vectors) == {"vel"},
        f"{label}: retained particle field inventory drifted",
    )
    try:
        header = admission._parse_pvtk_execution_header(payload, label)
    except ValueError as error:
        raise NumericalQualificationError(f"{label}: retained particle header decode failed") from error
    _require(header["time"] == time, f"{label}: retained particle time drifted")
    integer_values = [
        header["nranks"],
        header["cycle"],
        *(
            int(value)
            for name in _RESTART_PARTICLE_INTEGER_SCALARS
            for value in decoded.scalars[name].reshape(-1)
        ),
    ]
    float_arrays = [
        decoded.points.reshape(-1),
        decoded.vectors["vel"].reshape(-1),
        *(
            decoded.scalars[name].reshape(-1)
            for name in _RESTART_PARTICLE_REAL_SCALARS
        ),
    ]
    _require(
        all(np.all(np.isfinite(values)) for values in float_arrays),
        f"{label}: retained particle payload contains non-finite values",
    )
    return integer_values, [float(value) for values in float_arrays for value in values]


def _extract_retained_restart_observation(
    value: object,
    payloads: Mapping[str, bytes],
    *,
    checkpoint_payload: bytes,
    execution_receipt: Mapping[str, object],
    artifact_root_descriptor: int,
    label: str,
) -> dict[str, object]:
    source = _object(
        value,
        {
            "branch_role",
            "execution_receipt_authority_path",
            "execution_receipt_member",
            "raw_output_root",
            "structured_artifact_inventory_sha256",
            "outputs_after_checkpoint",
        },
        label=label,
    )
    try:
        policy = restart.load_preregistration()
        contract = policy["continuation_contract"]
        binding = restart.bind_checkpoint_for_continuation(
            checkpoint_payload,
            checkpoint_time_omega0_inverse=contract["checkpoint_time_omega0_inverse"],
            retained_output_schedule_after_checkpoint_omega0_inverse=contract[
                "retained_output_schedule_after_checkpoint_omega0_inverse"
            ],
            comparison_tolerances_max_absolute_difference=contract[
                "comparison_tolerances_max_absolute_difference"
            ],
            preregistration=policy,
            source="source_checkpoint_member",
        )
    except restart.RestartPolicyError as error:
        raise NumericalQualificationError(f"{label}: retained checkpoint binding failed") from error
    descriptors = _list(source["outputs_after_checkpoint"], label=f"{label}/outputs_after_checkpoint")
    schedule = binding["retained_output_schedule_after_checkpoint_omega0_inverse"]
    _require(len(descriptors) == len(schedule), f"{label}: retained output schedule length drifted")
    raw_output_root = _text(source["raw_output_root"], label=f"{label}/raw_output_root")
    _require(
        raw_output_root == execution_receipt["raw_output_root"],
        f"{label}: raw-output authority root differs from registered execution receipt",
    )
    inventory_sha256 = _sha256(
        source["structured_artifact_inventory_sha256"],
        label=f"{label}/structured_artifact_inventory_sha256",
    )
    outputs = []
    artifact_dir = Path(_text(execution_receipt["artifact_dir"], label=f"{label}/artifact_dir"))
    try:
        with _RegisteredRunStructuredArtifactTree(
            artifact_dir, inherited_root_fd=artifact_root_descriptor
        ) as tree:
            structured_artifacts.require_inventory_sha256(tree, inventory_sha256)
            inventory = structured_artifacts.load_inventory(tree)
            for index, (descriptor, expected_time) in enumerate(zip(descriptors, schedule)):
                output_label = f"{label}/outputs_after_checkpoint[{index}]"
                item = _object(
                    descriptor, {"time_omega0_inverse", "members"}, label=output_label
                )
                _require(
                    item["time_omega0_inverse"] == expected_time,
                    f"{output_label}: time drifted",
                )
                members = _object(
                    item["members"],
                    {*_RESTART_MESH_MEMBERS, _RESTART_PARTICLE_MEMBER},
                    label=f"{output_label}/members",
                )
                fields = {
                    name: _restart_mesh_values(
                        _structured_restart_member(
                            payloads,
                            tree,
                            inventory,
                            members[name],
                            label=f"{output_label}/members/{name}",
                        ),
                        field=field,
                        time=expected_time,
                        label=f"{output_label}/{name}",
                    )
                    for name, field in _RESTART_MESH_MEMBERS.items()
                }
                integers, floats = _restart_particle_values(
                    _structured_restart_member(
                        payloads,
                        tree,
                        inventory,
                        members[_RESTART_PARTICLE_MEMBER],
                        label=f"{output_label}/members/{_RESTART_PARTICLE_MEMBER}",
                    ),
                    time=expected_time,
                    label=f"{output_label}/{_RESTART_PARTICLE_MEMBER}",
                )
                fields["prtcl_all_pvtk_integer_payload"] = integers
                fields["prtcl_all_pvtk_float_payload"] = floats
                outputs.append({"time_omega0_inverse": expected_time, "fields": fields})
            tree.require_tree_closure()
    except ValueError as error:
        if isinstance(error, NumericalQualificationError):
            raise
        raise NumericalQualificationError(
            f"{label}: structured raw-output authority failed validation"
        ) from error
    return {"binding": binding, "outputs_after_checkpoint": outputs}


def bind_retained_restart_parity(
    bundle_root: str | Path,
    expected_inventory_sha256: str,
    *,
    source_attempt: _RetainedAttemptBinding,
    authorized_orion_root: Path,
) -> _RetainedRestartParityBinding:
    """Extract and bind immutable restart observations for parity comparison."""
    source, source_attempt_payload = _retained_attempt_payload(source_attempt, index=0)
    _require(
        source["identity"]["variant"] == AMR_VARIANT,
        "restart parity source attempt must be an AMR baseline",
    )
    manifest, value, payloads = _load_retained_recompute_record(
        bundle_root,
        expected_inventory_sha256,
        member=RESTART_RECOMPUTE_MEMBER,
        authorized_orion_root=authorized_orion_root,
    )
    _require_recompute_manifest(
        manifest,
        raw_inventory_sha256=source["raw_inventory_sha256"],
        analyzer_paths=RESTART_RECOMPUTE_ANALYZERS,
    )
    record = _object(
        value,
        {
            "schema_version",
            "record_type",
            "source_attempt_sha256",
            "raw_inventory_sha256",
            "planner_materialization_receipt_member",
            "retained_campaign_plan_member",
            "planned_restart_carrier",
            "source_checkpoint_member",
            "source_checkpoint_sha256",
            "source_checkpoint_lineage",
            "screen_scope",
            "uninterrupted",
            "continued",
        },
        label=RESTART_RECOMPUTE_MEMBER,
    )
    _require(
        type(record["schema_version"]) is int
        and record["schema_version"] == SCHEMA_VERSION
        and record["record_type"] == RESTART_RECOMPUTE_RECORD_TYPE,
        "retained restart recompute identity drifted",
    )
    _require(
        record["source_attempt_sha256"] == source["attempt_sha256"]
        and record["raw_inventory_sha256"] == source["raw_inventory_sha256"],
        "retained restart recompute raw attempt binding drifted",
    )
    _require(
        record["screen_scope"] == RESTART_SCREEN_SCOPE,
        "retained restart recompute screen scope drifted",
    )
    planned = _validate_planned_restart_carrier(
        record,
        payloads,
        source=source,
        source_attempt_payload=source_attempt_payload,
        authorized_orion_root=authorized_orion_root,
    )
    checkpoint_member = _relative_path(
        record["source_checkpoint_member"],
        label=f"{RESTART_RECOMPUTE_MEMBER}/source_checkpoint_member",
    )
    checkpoint_payload = _bundle_member(
        payloads,
        checkpoint_member,
        label=f"{RESTART_RECOMPUTE_MEMBER}/source_checkpoint_member",
    )
    checkpoint_sha256 = _sha256(
        record["source_checkpoint_sha256"],
        label=f"{RESTART_RECOMPUTE_MEMBER}/source_checkpoint_sha256",
    )
    _require(
        record["source_checkpoint_lineage"] == source["source_checkpoint_lineage"]
        and checkpoint_sha256
        == source["source_checkpoint_lineage"]["restart_member_sha256"]
        == hashlib.sha256(checkpoint_payload).hexdigest(),
        "retained restart source checkpoint lineage drifted",
    )
    branch_paths = {}
    branches = {}
    receipts = {}
    observations = {}
    for role, key in (
        ("uninterrupted_baseline", "uninterrupted"),
        ("checkpoint_restart_continuation", "continued"),
    ):
        branch = _object(
            record[key],
            {
                "branch_role",
                "execution_receipt_authority_path",
                "execution_receipt_member",
                "raw_output_root",
                "structured_artifact_inventory_sha256",
                "outputs_after_checkpoint",
            },
            label=key,
        )
        _require(branch["branch_role"] == role, f"{key}: branch role drifted")
        branches[key] = branch
        branch_paths[key] = _restart_branch_member_paths(branch, label=key)
    metadata_members = {
        checkpoint_member,
        _relative_path(
            record["planner_materialization_receipt_member"],
            label=f"{RESTART_RECOMPUTE_MEMBER}/planner_materialization_receipt_member",
        ),
        _relative_path(
            record["retained_campaign_plan_member"],
            label=f"{RESTART_RECOMPUTE_MEMBER}/retained_campaign_plan_member",
        ),
    }
    _require(
        branch_paths["uninterrupted"].isdisjoint(branch_paths["continued"])
        and metadata_members.isdisjoint(branch_paths["uninterrupted"])
        and metadata_members.isdisjoint(branch_paths["continued"]),
        "retained restart branches or metadata members alias",
    )
    for role, key in (
        ("uninterrupted_baseline", "uninterrupted"),
        ("checkpoint_restart_continuation", "continued"),
    ):
        with _validated_authoritative_restart_execution_receipt(
            branches[key],
            payloads,
            branch_role=role,
            source=source,
            source_attempt_payload=source_attempt_payload,
            source_checkpoint_lineage=source["source_checkpoint_lineage"],
            planned=planned,
            authorized_orion_root=authorized_orion_root,
            label=f"{key}/execution_receipt_member",
        ) as (receipt, artifact_root_descriptor):
            receipts[key] = receipt
            observations[key] = _extract_retained_restart_observation(
                branches[key],
                payloads,
                checkpoint_payload=checkpoint_payload,
                execution_receipt=receipt,
                artifact_root_descriptor=artifact_root_descriptor,
                label=key,
            )
    _require(
        receipts["uninterrupted"]["submission_id"]
        != receipts["continued"]["submission_id"]
        and receipts["uninterrupted"]["reservation_id"]
        != receipts["continued"]["reservation_id"]
        and receipts["uninterrupted"]["slurm_job_id"]
        != receipts["continued"]["slurm_job_id"]
        and receipts["uninterrupted"]["reconciliation_event_sha256"]
        != receipts["continued"]["reconciliation_event_sha256"],
        "retained restart authoritative branch execution identities alias",
    )
    retained = copy.deepcopy(record)
    retained["uninterrupted"] = observations["uninterrupted"]
    retained["continued"] = observations["continued"]
    return _RetainedRestartParityBinding(
        record_payload=_json_bytes(retained),
        _key=_PRODUCTION_BINDING_KEY,
    )


def _expected_amr_pairing(identity: Mapping[str, object]) -> dict[str, object]:
    variant = identity["variant"]
    seed = identity["seed"]
    counterpart = None
    if variant == AMR_VARIANT:
        counterpart = {"variant": FINE_VARIANT, "seed": seed}
    elif variant == FINE_VARIANT:
        counterpart = {"variant": AMR_VARIANT, "seed": seed}
    return {
        "paired_seed_rule": PAIRED_SEED_RULE,
        "pair_key": {"seed": seed},
        "amr_fine_uniform_counterpart": counterpart,
        "comparison_status": (
            "schema_wired_not_evaluated_by_artifact_admission_slice"
        ),
    }


def _validate_identity(value: object, *, label: str) -> dict[str, object]:
    identity = _object(
        value,
        {"variant", "seed", "physical_mode", "attempt_id"},
        label=label,
    )
    _require(identity["variant"] in GRID_VARIANTS, f"{label}: variant is not canonical")
    _require(
        type(identity["seed"]) is int and identity["seed"] in QUALIFYING_SEEDS,
        f"{label}: seed is not canonical",
    )
    _require(identity["physical_mode"] == PHYSICAL_MODE, f"{label}: physical mode drifted")
    _require(
        type(identity["attempt_id"]) is str
        and _ATTEMPT_ID.fullmatch(identity["attempt_id"]) is not None,
        f"{label}: attempt ID is not canonical",
    )
    return identity


def _validate_admission(
    value: object,
    *,
    identity: Mapping[str, object],
    raw_inventory_sha256: str,
    label: str,
) -> None:
    result = _object(
        value,
        {
            "schema_version",
            "record_type",
            "campaign_id",
            "qualification_scope",
            "admitted_for_follow_on_numerical_qualification",
            "final_claim_closure",
            "status",
            "failure_reasons",
            "admission",
        },
        label=label,
    )
    _require(
        type(result["schema_version"]) is int and result["schema_version"] == 1,
        f"{label}: schema version drifted",
    )
    _require(result["record_type"] == admission.RESULT_RECORD_TYPE, f"{label}: record type drifted")
    _require(result["campaign_id"] == admission.CAMPAIGN_ID, f"{label}: campaign ID drifted")
    _require(
        result["qualification_scope"] == admission.QUALIFICATION_SCOPE,
        f"{label}: qualification scope drifted",
    )
    _require(
        result["admitted_for_follow_on_numerical_qualification"] is True,
        f"{label}: attempt was not admitted",
    )
    _require(result["final_claim_closure"] is False, f"{label}: claim closure must remain false")
    _require(
        result["status"] == "admitted_for_follow_on_numerical_qualification",
        f"{label}: admission status drifted",
    )
    _require(result["failure_reasons"] == [], f"{label}: admitted result has failure reasons")
    admitted = _required_object(
        result["admission"],
        {
            "run_identity",
            "preregistration_binding",
            "immutable_tree",
            "amr_pairing",
            "numerical_qualification_status",
        },
        label=f"{label}/admission",
    )
    _require(admitted["run_identity"] == identity, f"{label}: run identity binding drifted")
    immutable_tree = _required_object(
        admitted["immutable_tree"], {"inventory_sha256"}, label=f"{label}/immutable_tree"
    )
    _require(
        immutable_tree["inventory_sha256"] == raw_inventory_sha256,
        f"{label}: raw inventory binding drifted",
    )
    preregistration = _required_object(
        admitted["preregistration_binding"],
        {"sha256", "expected_sha256"},
        label=f"{label}/preregistration_binding",
    )
    _require(
        preregistration["sha256"] == admission.EXPECTED_PREREGISTRATION_SHA256
        and preregistration["expected_sha256"] == admission.EXPECTED_PREREGISTRATION_SHA256,
        f"{label}: preregistration binding drifted",
    )
    _require(
        admitted["amr_pairing"] == _expected_amr_pairing(identity),
        f"{label}: admission AMR pairing drifted",
    )
    _require(
        admitted["numerical_qualification_status"]
        == "not_evaluated_by_artifact_admission_slice",
        f"{label}: admission slice overclaimed numerical qualification",
    )


def _validate_overflow_gate(spectrum: object, *, label: str) -> bool:
    record = _required_object(
        spectrum,
        {
            "f_chi",
            "overflow_macro_weight",
            "total_post_filter_macro_weight",
            "overflow_macro_weight_fraction",
            "overflow_gate",
        },
        label=label,
    )
    overflow = _finite_float(record["overflow_macro_weight"], label=f"{label}/overflow")
    total = _finite_float(
        record["total_post_filter_macro_weight"], label=f"{label}/total"
    )
    _require(total > 0.0, f"{label}: total admitted macro weight must be positive")
    fraction = _finite_float(
        record["overflow_macro_weight_fraction"], label=f"{label}/overflow fraction"
    )
    _require(fraction == overflow / total, f"{label}: overflow fraction is inconsistent")
    gate = _object(
        record["overflow_gate"],
        {"maximum_macro_weight_fraction", "passed"},
        label=f"{label}/overflow_gate",
    )
    _require(
        gate["maximum_macro_weight_fraction"]
        == particles.MAX_OVERFLOW_MACRO_WEIGHT_FRACTION,
        f"{label}: overflow threshold drifted",
    )
    _require(type(gate["passed"]) is bool, f"{label}: overflow gate result must be boolean")
    _require(
        gate["passed"] == (fraction <= particles.MAX_OVERFLOW_MACRO_WEIGHT_FRACTION),
        f"{label}: overflow gate result is inconsistent",
    )
    return gate["passed"]


def _validate_particle_reduction(
    value: object, *, expected_time: float, require_late_slope: bool, label: str
) -> list[bool]:
    record = _required_object(
        value,
        {"schema_version", "record_type", "snapshot_time_omega0_inverse", "weighted_spectrum"},
        label=label,
    )
    _require(
        type(record["schema_version"]) is int and record["schema_version"] == particles.SCHEMA_VERSION,
        f"{label}: schema version drifted",
    )
    _require(
        record["record_type"] == "q011_section54_particle_snapshot_reduction",
        f"{label}: record type drifted",
    )
    _require(record["snapshot_time_omega0_inverse"] == expected_time, f"{label}: time drifted")
    try:
        particles.canonical_record_bytes(record)
    except particles.ParticleReducerError as error:
        raise NumericalQualificationError(f"{label}: particle record is invalid") from error
    gates = [_validate_overflow_gate(record["weighted_spectrum"], label=f"{label}/weighted_spectrum")]
    if require_late_slope:
        _require("late_slope" in record, f"{label}: late slope record is missing")
        try:
            recomputed = particles.late_slope_record(record["weighted_spectrum"]["f_chi"])
        except particles.ParticleReducerError as error:
            raise NumericalQualificationError(f"{label}: late slope record is invalid") from error
        _require(record["late_slope"] == recomputed, f"{label}: late slope record drifted")
        gates.append(recomputed["slope_gate_passed"])
    else:
        _require("late_slope" not in record, f"{label}: unexpected early late-slope record")
    return gates


def _validate_spatial_reduction(value: object, *, label: str) -> bool:
    record = _required_object(
        value,
        {"schema_version", "record_type", "time_omega0_inverse", "upstream_b_amplification"},
        label=label,
    )
    _require(
        type(record["schema_version"]) is int and record["schema_version"] == 1,
        f"{label}: schema version drifted",
    )
    _require(
        record["record_type"] == "q011_section54_t500_spatial_reduction",
        f"{label}: record type drifted",
    )
    _require(record["time_omega0_inverse"] == spatial.T500_OMEGA0_INVERSE, f"{label}: time drifted")
    amplification = _object(
        record["upstream_b_amplification"],
        {
            "time_omega0_inverse",
            "x_ideal_c_over_omega_pi",
            "upstream_window_c_over_omega_pi",
            "selected_cell_count",
            "selected_area",
            "mean_magnetic_magnitude",
            "reference_b0",
            "amplification_over_b0",
            "acceptance_range",
            "passes_gate",
        },
        label=f"{label}/upstream_b_amplification",
    )
    try:
        parsed = spatial.UpstreamBAmplificationRecord(
            time_omega0_inverse=amplification["time_omega0_inverse"],
            x_ideal_c_over_omega_pi=amplification["x_ideal_c_over_omega_pi"],
            upstream_window_c_over_omega_pi=tuple(amplification["upstream_window_c_over_omega_pi"]),
            selected_cell_count=amplification["selected_cell_count"],
            selected_area=amplification["selected_area"],
            mean_magnetic_magnitude=amplification["mean_magnetic_magnitude"],
            reference_b0=amplification["reference_b0"],
            amplification_over_b0=amplification["amplification_over_b0"],
            acceptance_range=tuple(amplification["acceptance_range"]),
            passes_gate=amplification["passes_gate"],
        )
        normalized = spatial.upstream_b_amplification_record(parsed)
    except (TypeError, spatial.AnalysisError) as error:
        raise NumericalQualificationError(f"{label}: spatial amplification is invalid") from error
    _require(amplification == normalized, f"{label}: spatial amplification record drifted")
    return parsed.passes_gate


def _validate_source_checkpoint_lineage(
    value: object,
    *,
    identity: Mapping[str, object],
    label: str,
) -> dict[str, object]:
    lineage = _object(
        value,
        {
            "snapshot_time_omega0_inverse",
            "retained_attempt_id",
            "restart_manifest_path",
            "restart_member_path",
            "retained_restart_member_absolute_path",
            "restart_member_sha256",
        },
        label=label,
    )
    manifest_path = _relative_path(
        lineage["restart_manifest_path"], label=f"{label}/restart_manifest_path"
    )
    member_path = _relative_path(
        lineage["restart_member_path"], label=f"{label}/restart_member_path"
    )
    absolute = Path(
        _text(
            lineage["retained_restart_member_absolute_path"],
            label=f"{label}/retained_restart_member_absolute_path",
        )
    )
    _require(
        lineage["snapshot_time_omega0_inverse"] == 500.0
        and lineage["retained_attempt_id"] == identity["attempt_id"]
        and manifest_path.endswith(".rst.manifest")
        and member_path.endswith(".rst")
        and absolute.is_absolute()
        and absolute.as_posix().endswith(f"/{member_path}"),
        f"{label}: admitted t=500 restart lineage drifted",
    )
    _sha256(lineage["restart_member_sha256"], label=f"{label}/restart_member_sha256")
    return copy.deepcopy(lineage)


def _validate_attempt_wrapper(value: object, *, index: int) -> dict[str, object]:
    label = f"attempts[{index}]"
    wrapper = _object(value, {"attempt_sha256", "attempt"}, label=label)
    expected_sha256 = _sha256(wrapper["attempt_sha256"], label=f"{label}/attempt_sha256")
    _require(
        _canonical_sha256(wrapper["attempt"]) == expected_sha256,
        f"{label}: canonical attempt SHA-256 drifted",
    )
    attempt = _object(
        wrapper["attempt"],
        {
            "schema_version",
            "record_type",
            "run_identity",
            "raw_inventory_sha256",
            "source_checkpoint_lineage",
            "admission_result",
            "particle_reductions",
            "spatial_reduction",
        },
        label=f"{label}/attempt",
    )
    _require(
        type(attempt["schema_version"]) is int and attempt["schema_version"] == SCHEMA_VERSION,
        f"{label}: schema version drifted",
    )
    _require(attempt["record_type"] == ATTEMPT_RECORD_TYPE, f"{label}: record type drifted")
    identity = _validate_identity(attempt["run_identity"], label=f"{label}/run_identity")
    inventory = _sha256(attempt["raw_inventory_sha256"], label=f"{label}/raw_inventory_sha256")
    checkpoint_lineage = _validate_source_checkpoint_lineage(
        attempt["source_checkpoint_lineage"],
        identity=identity,
        label=f"{label}/source_checkpoint_lineage",
    )
    _validate_admission(
        attempt["admission_result"],
        identity=identity,
        raw_inventory_sha256=inventory,
        label=f"{label}/admission_result",
    )
    reductions = _object(
        attempt["particle_reductions"], {"t500", "t1200"}, label=f"{label}/particle_reductions"
    )
    gates = _validate_particle_reduction(
        reductions["t500"], expected_time=500.0, require_late_slope=False, label=f"{label}/t500"
    )
    gates.extend(
        _validate_particle_reduction(
            reductions["t1200"],
            expected_time=particles.LATE_SLOPE_SNAPSHOT_TIME,
            require_late_slope=True,
            label=f"{label}/t1200",
        )
    )
    gates.append(_validate_spatial_reduction(attempt["spatial_reduction"], label=f"{label}/spatial"))
    return {
        "attempt_sha256": expected_sha256,
        "identity": identity,
        "raw_inventory_sha256": inventory,
        "source_checkpoint_lineage": checkpoint_lineage,
        "gates_passed": all(gates),
    }


def _validate_attempt_binding(
    value: object, *, index: int, allow_unit_only: bool
) -> dict[str, object]:
    label = f"attempts[{index}]"
    if type(value) is _RetainedAttemptBinding:
        _require(value._key is _PRODUCTION_BINDING_KEY, f"{label}: invalid retained binding")
        wrapper = _load_json_bytes(value.wrapper_payload, label=f"{label}/retained_binding")
    else:
        _require(
            allow_unit_only and type(value) is dict,
            f"{label}: expected a retained attempt tree binding",
        )
        wrapper = value
    return _validate_attempt_wrapper(wrapper, index=index)


def _validate_residual_metric(
    value: object,
    *,
    expected: tuple[str, float | None, float | None, float | None],
    label: str,
) -> tuple[dict[str, object], bool]:
    metric = _object(
        value,
        {
            "observable",
            "maximum_absolute_difference",
            "relative_mean_absolute",
            "relative_root_mean_square",
        },
        label=label,
    )
    observable, maximum_bound, mean_bound, rms_bound = expected
    _require(metric["observable"] == observable, f"{label}: observable order drifted")
    maximum = _finite_float(metric["maximum_absolute_difference"], label=f"{label}/maximum")
    passed = maximum_bound is None or maximum <= maximum_bound
    for field, bound in (
        ("relative_mean_absolute", mean_bound),
        ("relative_root_mean_square", rms_bound),
    ):
        if bound is None:
            _require(metric[field] is None, f"{label}/{field}: expected null")
        else:
            parsed = _finite_float(metric[field], label=f"{label}/{field}")
            passed = passed and parsed <= bound
    return metric, passed


def _validate_pairs(
    value: object,
    *,
    by_sha256: Mapping[str, Mapping[str, object]],
    allow_unit_only: bool,
) -> list[dict[str, object]]:
    pairs = _list(value, label="paired_amr_fine_results")
    _require(len(pairs) == len(QUALIFYING_SEEDS), "paired AMR/fine result count drifted")
    results = []
    for index, (pair_value, seed) in enumerate(zip(pairs, QUALIFYING_SEEDS)):
        label = f"paired_amr_fine_results[{index}]"
        retained = type(pair_value) is _RetainedPairResultBinding
        if retained:
            _require(
                pair_value._key is _PRODUCTION_BINDING_KEY,
                f"{label}: invalid retained pair binding",
            )
            pair = _object(
                _load_json_bytes(pair_value.record_payload, label=f"{label}/retained_binding"),
                {
                    "schema_version",
                    "record_type",
                    "seed",
                    "amr_attempt_sha256",
                    "amr_raw_inventory_sha256",
                    "fine_uniform_attempt_sha256",
                    "fine_uniform_raw_inventory_sha256",
                    "residuals",
                },
                label=label,
            )
            _require(
                type(pair["schema_version"]) is int
                and pair["schema_version"] == SCHEMA_VERSION
                and pair["record_type"] == PAIR_RECOMPUTE_RECORD_TYPE,
                f"{label}: retained pair recompute identity drifted",
            )
        else:
            _require(
                allow_unit_only and type(pair_value) is dict,
                f"{label}: expected a retained pair recompute binding",
            )
            pair = _object(
                pair_value,
                {"seed", "amr_attempt_sha256", "fine_uniform_attempt_sha256", "residuals"},
                label=label,
            )
        _require(type(pair["seed"]) is int and pair["seed"] == seed, f"{label}: seed order drifted")
        amr_sha256 = _sha256(pair["amr_attempt_sha256"], label=f"{label}/amr_attempt_sha256")
        fine_sha256 = _sha256(
            pair["fine_uniform_attempt_sha256"], label=f"{label}/fine_uniform_attempt_sha256"
        )
        _require(amr_sha256 in by_sha256 and fine_sha256 in by_sha256, f"{label}: attempt binding is unknown")
        _require(
            (by_sha256[amr_sha256]["identity"]["variant"], by_sha256[amr_sha256]["identity"]["seed"])
            == (AMR_VARIANT, seed),
            f"{label}: AMR attempt is not the matched seed",
        )
        _require(
            (by_sha256[fine_sha256]["identity"]["variant"], by_sha256[fine_sha256]["identity"]["seed"])
            == (FINE_VARIANT, seed),
            f"{label}: fine-uniform attempt is not the matched seed",
        )
        if retained:
            _require(
                pair["amr_raw_inventory_sha256"]
                == by_sha256[amr_sha256]["raw_inventory_sha256"]
                and pair["fine_uniform_raw_inventory_sha256"]
                == by_sha256[fine_sha256]["raw_inventory_sha256"],
                f"{label}: retained pair raw inventory binding drifted",
            )
        metrics = _list(pair["residuals"], label=f"{label}/residuals")
        _require(
            len(metrics) == len(PAIRED_RESIDUAL_THRESHOLDS),
            f"{label}: residual metric count drifted",
        )
        parsed_metrics = []
        gates = []
        for metric_index, (metric, expected) in enumerate(
            zip(metrics, PAIRED_RESIDUAL_THRESHOLDS)
        ):
            parsed, passed = _validate_residual_metric(
                metric, expected=expected, label=f"{label}/residuals[{metric_index}]"
            )
            parsed_metrics.append(parsed)
            gates.append(passed)
        results.append(
            {
                "seed": seed,
                "amr_attempt_sha256": amr_sha256,
                "fine_uniform_attempt_sha256": fine_sha256,
                "residuals": parsed_metrics,
                "gates_passed": all(gates),
            }
        )
    return results


def _validate_restart_binding(
    value: object,
    *,
    by_sha256: Mapping[str, Mapping[str, object]],
    allow_unit_only: bool,
) -> dict[str, object]:
    retained = type(value) is _RetainedRestartParityBinding
    if retained:
        _require(value._key is _PRODUCTION_BINDING_KEY, "invalid retained restart binding")
        binding = _object(
            _load_json_bytes(value.record_payload, label="restart_parity_binding/retained_binding"),
            {
                "schema_version",
                "record_type",
                "source_attempt_sha256",
                "raw_inventory_sha256",
                "planner_materialization_receipt_member",
                "retained_campaign_plan_member",
                "planned_restart_carrier",
                "source_checkpoint_member",
                "source_checkpoint_sha256",
                "source_checkpoint_lineage",
                "screen_scope",
                "uninterrupted",
                "continued",
            },
            label="restart_parity_binding",
        )
        _require(
            type(binding["schema_version"]) is int
            and binding["schema_version"] == SCHEMA_VERSION
            and binding["record_type"] == RESTART_RECOMPUTE_RECORD_TYPE,
            "retained restart recompute identity drifted",
        )
    else:
        _require(
            allow_unit_only and type(value) is dict,
            "restart_parity_binding: expected a retained restart recompute binding",
        )
        binding = _object(
            value,
            {"source_attempt_sha256", "uninterrupted", "continued"},
            label="restart_parity_binding",
        )
    source_sha256 = _sha256(
        binding["source_attempt_sha256"], label="restart_parity_binding/source_attempt_sha256"
    )
    _require(source_sha256 in by_sha256, "restart parity source attempt is unknown")
    _require(
        by_sha256[source_sha256]["identity"]["variant"] == AMR_VARIANT,
        "restart parity source attempt must be an AMR baseline",
    )
    if retained:
        _require(
            binding["raw_inventory_sha256"]
            == by_sha256[source_sha256]["raw_inventory_sha256"],
            "retained restart raw inventory binding drifted",
        )
        _require(
            binding["screen_scope"] == RESTART_SCREEN_SCOPE,
            "retained restart screen scope drifted",
        )
        _require(
            binding["source_checkpoint_lineage"]
            == by_sha256[source_sha256]["source_checkpoint_lineage"],
            "retained restart checkpoint lineage binding drifted",
        )
    binding_sha256 = _canonical_sha256(binding)
    try:
        result = restart.compare_deterministic_continuation_parity(
            binding["uninterrupted"], binding["continued"]
        )
    except restart.RestartPolicyError as error:
        raise NumericalQualificationError("restart continuation parity failed") from error
    _require(
        result.get("result") == "pass_deterministic_continuation_parity",
        "restart continuation comparator did not return a passing result",
    )
    return {
        "source_attempt_sha256": source_sha256,
        "restart_parity_binding_sha256": binding_sha256,
        "screen_scope": (
            RESTART_SCREEN_SCOPE if retained else "unit_only_synthetic_restart_comparator"
        ),
        "full_state_equivalence_claimed": False,
        "result": result,
    }


def _qualify_bound_numerical_aggregate(
    *,
    attempts: Sequence[object],
    ordered_raw_inventory_sha256_values: Sequence[object],
    paired_amr_fine_results: Sequence[object],
    restart_parity_binding: object,
    reviewer_disposition: object,
    allow_unit_only: bool,
) -> dict[str, object]:
    """Aggregate one exact 3x8 bound matrix without closing external review."""
    _require(type(attempts) is list, "attempts: expected list")
    _require(
        len(attempts) == EXPECTED_BASELINE_ATTEMPTS,
        "attempts: expected exactly 24 canonical baseline attempts",
    )
    parsed_attempts = [
        _validate_attempt_binding(
            value, index=index, allow_unit_only=allow_unit_only
        )
        for index, value in enumerate(attempts)
    ]
    attempt_sha256_values = [item["attempt_sha256"] for item in parsed_attempts]
    _require(
        len(set(attempt_sha256_values)) == EXPECTED_BASELINE_ATTEMPTS,
        "attempts: duplicate canonical baseline attempt",
    )
    attempt_ids = [item["identity"]["attempt_id"] for item in parsed_attempts]
    _require(len(set(attempt_ids)) == EXPECTED_BASELINE_ATTEMPTS, "attempts: duplicate attempt ID")
    cells = [
        (item["identity"]["variant"], item["identity"]["seed"])
        for item in parsed_attempts
    ]
    _require(cells == list(EXPECTED_MATRIX_CELLS), "attempts: canonical 3x8 matrix order drifted")

    inventories = _list(
        ordered_raw_inventory_sha256_values,
        label="ordered_raw_inventory_sha256_values",
    )
    _require(
        len(inventories) == EXPECTED_BASELINE_ATTEMPTS,
        "ordered raw inventory digest count drifted",
    )
    parsed_inventories = [
        _sha256(value, label=f"ordered_raw_inventory_sha256_values[{index}]")
        for index, value in enumerate(inventories)
    ]
    _require(
        parsed_inventories
        == [item["raw_inventory_sha256"] for item in parsed_attempts],
        "ordered raw inventory digest binding drifted",
    )

    by_sha256 = {item["attempt_sha256"]: item for item in parsed_attempts}
    pair_results = _validate_pairs(
        paired_amr_fine_results,
        by_sha256=by_sha256,
        allow_unit_only=allow_unit_only,
    )
    restart_result = _validate_restart_binding(
        restart_parity_binding,
        by_sha256=by_sha256,
        allow_unit_only=allow_unit_only,
    )
    try:
        disposition = artifacts.validate_reviewer_disposition(reviewer_disposition)
    except artifacts.DerivedArtifactError as error:
        raise NumericalQualificationError("reviewer disposition is invalid") from error
    _require(
        disposition["status"] == "pending_external_review",
        "numerical aggregation must remain pending external review",
    )
    numerical_gates_passed = (
        all(item["gates_passed"] for item in parsed_attempts)
        and all(item["gates_passed"] for item in pair_results)
    )
    result = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "campaign_id": admission.CAMPAIGN_ID,
        "qualification_scope": QUALIFICATION_SCOPE,
        "final_claim_closure": False,
        "status": "pending_external_review",
        "numerical_gate_status": "passed" if numerical_gates_passed else "failed",
        "numerical_gates_passed": numerical_gates_passed,
        "campaign_matrix": {
            "physical_mode": PHYSICAL_MODE,
            "grid_variants": list(GRID_VARIANTS),
            "qualifying_seeds": list(QUALIFYING_SEEDS),
            "expected_baseline_attempts": EXPECTED_BASELINE_ATTEMPTS,
            "paired_seed_rule": PAIRED_SEED_RULE,
        },
        "attempt_count": len(parsed_attempts),
        "attempt_sha256_values": attempt_sha256_values,
        "ordered_raw_inventory_sha256_values": parsed_inventories,
        "attempt_gate_results": parsed_attempts,
        "paired_amr_fine_results": pair_results,
        "restart_parity": restart_result,
        "reviewer_disposition": disposition,
    }
    _canonical_sha256(result)
    return result


def _retained_tree_descriptor(value: object, *, label: str) -> tuple[str, str]:
    descriptor = _object(
        value,
        {"bundle_root", "expected_inventory_sha256"},
        label=label,
    )
    _require(
        type(descriptor["bundle_root"]) is str,
        f"{label}/bundle_root: expected path text",
    )
    inventory = _sha256(
        descriptor["expected_inventory_sha256"],
        label=f"{label}/expected_inventory_sha256",
    )
    return descriptor["bundle_root"], inventory


def _attempt_tree_descriptor(value: object, *, label: str) -> tuple[str, str]:
    descriptor = _object(
        value,
        {"campaign_root", "expected_inventory_sha256"},
        label=label,
    )
    _require(
        type(descriptor["campaign_root"]) is str,
        f"{label}/campaign_root: expected path text",
    )
    inventory = _sha256(
        descriptor["expected_inventory_sha256"],
        label=f"{label}/expected_inventory_sha256",
    )
    return descriptor["campaign_root"], inventory


def qualify_numerical_aggregate(
    *,
    attempts: Sequence[object],
    ordered_raw_inventory_sha256_values: Sequence[object],
    paired_amr_fine_results: Sequence[object],
    restart_parity_binding: object,
    reviewer_disposition: object,
    authorized_orion_root: Path = admission.ORION_BULK_ROOT,
) -> dict[str, object]:
    """Load retained evidence descriptors and aggregate one exact admitted matrix."""
    attempt_descriptors = _list(attempts, label="attempts")
    _require(
        len(attempt_descriptors) == EXPECTED_BASELINE_ATTEMPTS,
        "attempts: expected exactly 24 immutable retained attempt-tree descriptors",
    )
    retained_attempts = []
    for index, value in enumerate(attempt_descriptors):
        root, inventory = _attempt_tree_descriptor(
            value, label=f"attempts[{index}]"
        )
        retained_attempts.append(
            bind_retained_attempt_tree(
                root,
                inventory,
                authorized_orion_root=authorized_orion_root,
            )
        )
    parsed_attempts = [
        _validate_attempt_binding(value, index=index, allow_unit_only=False)
        for index, value in enumerate(retained_attempts)
    ]
    by_cell = {
        (item["identity"]["variant"], item["identity"]["seed"]): value
        for item, value in zip(parsed_attempts, retained_attempts)
    }
    pair_descriptors = _list(
        paired_amr_fine_results, label="paired_amr_fine_results"
    )
    _require(
        len(pair_descriptors) == len(QUALIFYING_SEEDS),
        "paired AMR/fine retained recompute descriptor count drifted",
    )
    retained_pairs = []
    for index, (value, seed) in enumerate(zip(pair_descriptors, QUALIFYING_SEEDS)):
        root, inventory = _retained_tree_descriptor(
            value, label=f"paired_amr_fine_results[{index}]"
        )
        try:
            amr = by_cell[(AMR_VARIANT, seed)]
            fine = by_cell[(FINE_VARIANT, seed)]
        except KeyError as error:
            raise NumericalQualificationError(
                f"paired_amr_fine_results[{index}]: matched retained attempts are unavailable"
            ) from error
        retained_pairs.append(
            bind_retained_pair_result(
                root,
                inventory,
                amr_attempt=amr,
                fine_uniform_attempt=fine,
                authorized_orion_root=authorized_orion_root,
            )
        )
    restart_descriptor = _object(
        restart_parity_binding,
        {"bundle_root", "expected_inventory_sha256", "source_attempt_sha256"},
        label="restart_parity_binding",
    )
    _require(
        type(restart_descriptor["bundle_root"]) is str,
        "restart_parity_binding/bundle_root: expected path text",
    )
    restart_inventory = _sha256(
        restart_descriptor["expected_inventory_sha256"],
        label="restart_parity_binding/expected_inventory_sha256",
    )
    restart_source = _sha256(
        restart_descriptor["source_attempt_sha256"],
        label="restart_parity_binding/source_attempt_sha256",
    )
    by_sha256 = {
        item["attempt_sha256"]: value
        for item, value in zip(parsed_attempts, retained_attempts)
    }
    _require(
        restart_source in by_sha256,
        "restart parity source attempt descriptor is unknown",
    )
    retained_restart = bind_retained_restart_parity(
        restart_descriptor["bundle_root"],
        restart_inventory,
        source_attempt=by_sha256[restart_source],
        authorized_orion_root=authorized_orion_root,
    )
    return _qualify_bound_numerical_aggregate(
        attempts=retained_attempts,
        ordered_raw_inventory_sha256_values=ordered_raw_inventory_sha256_values,
        paired_amr_fine_results=retained_pairs,
        restart_parity_binding=retained_restart,
        reviewer_disposition=reviewer_disposition,
        allow_unit_only=False,
    )


def _qualify_unit_only_numerical_aggregate(
    *,
    attempts: Sequence[object],
    ordered_raw_inventory_sha256_values: Sequence[object],
    paired_amr_fine_results: Sequence[object],
    restart_parity_binding: object,
    reviewer_disposition: object,
) -> dict[str, object]:
    """Exercise aggregation mechanics with synthetic mappings in unit tests only."""
    return _qualify_bound_numerical_aggregate(
        attempts=attempts,
        ordered_raw_inventory_sha256_values=ordered_raw_inventory_sha256_values,
        paired_amr_fine_results=paired_amr_fine_results,
        restart_parity_binding=restart_parity_binding,
        reviewer_disposition=reviewer_disposition,
        allow_unit_only=True,
    )


__all__ = [
    "AMR_VARIANT",
    "ATTEMPT_RECORD_TYPE",
    "EXPECTED_BASELINE_ATTEMPTS",
    "EXPECTED_MATRIX_CELLS",
    "FINE_VARIANT",
    "GRID_VARIANTS",
    "NumericalQualificationError",
    "PAIRED_RESIDUAL_THRESHOLDS",
    "PHYSICAL_MODE",
    "QUALIFYING_SEEDS",
    "QUALIFICATION_SCOPE",
    "RECORD_TYPE",
    "SCHEMA_VERSION",
    "bind_retained_attempt_tree",
    "bind_retained_pair_result",
    "bind_retained_restart_parity",
    "qualify_numerical_aggregate",
]
