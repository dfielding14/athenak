#!/usr/bin/env python3
"""Bind one completed Q-011 raw attempt to one immutable planner tree.

This source-local bridge does not launch work, call a scheduler, mutate live
policy, or claim qualification.  It adds the provenance bindings and manifest
required by the existing raw-attempt freezer and campaign admission gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Mapping, Sequence
import uuid

if __package__:
    from . import analyze_q011_section54_campaign as campaign
    from . import analyze_q011_section54_outputs as output_primitives
    from . import immutable_orion_tree
    from . import publish_q011_section54_campaign_attempt as attempt_publisher
    from . import q011_section54_model as model
    from . import q011_section54_qualifying_campaign_execution as planner
else:
    import analyze_q011_section54_campaign as campaign
    import analyze_q011_section54_outputs as output_primitives
    import immutable_orion_tree
    import publish_q011_section54_campaign_attempt as attempt_publisher
    import q011_section54_model as model
    import q011_section54_qualifying_campaign_execution as planner


MANIFEST_RECORD_TYPE = "q011_section54_campaign_run_manifest"
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_MESH_NAME = re.compile(r"(?:[^/]+/)*[^/]+\.(rho|bmag|prtcl_jx|j2)\.[^/]+\.bin")
_PARTICLE_NAME = re.compile(r"(?:[^/]+/)*[^/]+\.prtcl_all\.[^/]+\.part\.vtk")
_RESTART_NAME = re.compile(r"rst/(?:rank_[0-9]{8}/)?[^/]+\.rst")
_RESTART_MANIFEST_NAME = re.compile(r"rst/[^/]+\.rst\.manifest")
_RESTART_CYCLE = re.compile(rb"(?:^|\n)cycle=([0-9]+)(?:\n|$)")
_RESTART_TIME = re.compile(
    rb"(?:^|\n)time=([-+]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)"
    rb"(?:[eE][-+]?[0-9]+)?)(?:\n|$)"
)
_BINDING_PAYLOADS = {
    "clean_candidate_manifest": "bindings/clean_candidate_manifest.json",
    "deck": "bindings/pic_parallel_shock_section54_paper_vl2_tsc.athinput",
    "analyzer": "bindings/analyze_q011_section54_campaign.py",
    "preregistration": (
        "bindings/q011_section54_qualifying_campaign_preregistration.json"
    ),
    "campaign_plan": "bindings/q011_section54_campaign_plan.json",
    "planner_materialization_receipt": (
        "bindings/q011_section54_planner_materialization_receipt.json"
    ),
    "attempt_contract": "bindings/q011_section54_attempt_contract.json",
    "selected_pressure_receipt": (
        "bindings/q011_section54_selected_pressure_receipt.json"
    ),
    "analyzer_helper_source_closure_manifest": (
        "bindings/q011_section54_analyzer_helper_source_closure_manifest.json"
    ),
    "registered_execution_receipt": (
        "bindings/q011_section54_registered_execution_receipt.json"
    ),
}


class AttemptManifestMaterializationError(ValueError):
    """Reject an unsafe, incomplete, or cross-linked raw attempt."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AttemptManifestMaterializationError(message)


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    try:
        return (
            json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode("utf-8")
    except (TypeError, ValueError) as error:
        raise AttemptManifestMaterializationError(
            "attempt manifest contains noncanonical JSON values"
        ) from error


def _json_payload(payload: bytes, *, label: str) -> dict[str, Any]:
    try:
        value = immutable_orion_tree.loads_json_reject_duplicate_keys(
            payload.decode("utf-8"),
            error_type=AttemptManifestMaterializationError,
            label=label,
        )
    except UnicodeDecodeError as error:
        raise AttemptManifestMaterializationError(f"{label} is not UTF-8") from error
    _require(type(value) is dict, f"{label} must be an object")
    return value


def _canonical_directory(value: str | Path, *, label: str) -> Path:
    path = Path(value)
    _require(path.is_absolute(), f"{label} must be absolute")
    try:
        resolved = path.resolve(strict=True)
        mode = os.lstat(path).st_mode
    except OSError as error:
        raise AttemptManifestMaterializationError(f"{label} is unavailable") from error
    _require(resolved == path, f"{label} must be canonical without symlink aliases")
    _require(
        os.path.isdir(path) and not os.path.islink(path),
        f"{label} must be a real directory",
    )
    _require(mode != 0, f"{label} mode is unavailable")
    return path


def _require_same_directory(path: Path, descriptor: int, *, label: str) -> None:
    """Require one open source directory to retain its public pathname."""
    try:
        current = os.stat(path, follow_symlinks=False)
        opened = os.fstat(descriptor)
    except OSError as error:
        raise AttemptManifestMaterializationError(f"{label} is unavailable") from error
    _require(
        stat.S_ISDIR(current.st_mode)
        and (current.st_dev, current.st_ino) == (opened.st_dev, opened.st_ino),
        f"{label} public pathname changed",
    )


def _absent_retained_attempt_root(
    value: str | Path, *, authorized_root: Path
) -> Path:
    """Validate the planner destination lexically without requiring it to exist."""
    path = Path(value)
    campaigns_root = authorized_root / "campaigns"
    _require(
        path.is_absolute() and Path(os.path.abspath(path)) == path,
        "planner-authorized retained attempt root must be canonical",
    )
    try:
        relative = path.relative_to(campaigns_root)
    except ValueError as error:
        raise AttemptManifestMaterializationError(
            "planner-authorized retained attempt root escaped campaigns namespace"
        ) from error
    current = campaigns_root
    for part in relative.parts:
        current /= part
        _require(
            not current.is_symlink(),
            "planner-authorized retained attempt root traverses a symlink",
        )
    _require(
        not path.exists(),
        "planner-authorized retained attempt root must be absent before publication",
    )
    return path


def _binding(relative: str, payload: bytes) -> dict[str, str]:
    return {"path": relative, "sha256": _sha256_bytes(payload)}


def _member_payload(
    snapshot: Any, binding: Mapping[str, Any], *, label: str
) -> bytes:
    _require(
        type(binding) is dict and set(binding) == {"path", "sha256"},
        f"{label} binding drifted",
    )
    relative = binding["path"]
    digest = binding["sha256"]
    _require(type(relative) is str and type(digest) is str, f"{label} binding drifted")
    try:
        payload = snapshot.member_path(relative).read_bytes()
    except (OSError, ValueError) as error:
        raise AttemptManifestMaterializationError(
            f"{label} is unavailable in immutable planner tree"
        ) from error
    _require(_sha256_bytes(payload) == digest, f"{label} SHA-256 drifted")
    return payload


def _stable_bound_source_payload(
    path: Path, *, expected_sha256: str, label: str
) -> bytes:
    """Read one live helper source only while it matches its retained binding."""
    try:
        descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise AttemptManifestMaterializationError(f"{label} is unavailable") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} is not a regular file")
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        try:
            current = os.stat(path, follow_symlinks=False)
        except OSError as error:
            raise AttemptManifestMaterializationError(
                f"{label} path changed while reading"
            ) from error
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
        result = bytes(payload)
        _require(_sha256_bytes(result) == expected_sha256, f"{label} SHA-256 drifted")
        return result
    finally:
        os.close(descriptor)


def _stable_registered_execution_receipt(
    path: str | Path,
    *,
    expected_path: Path,
    authorized_root: Path,
) -> bytes:
    """Read the reconciled execution receipt only from its immutable fixed path."""
    supplied = Path(path)
    _require(supplied == expected_path, "registered execution receipt path drifted")
    try:
        supplied.relative_to(authorized_root)
        resolved = supplied.resolve(strict=True)
    except (OSError, ValueError) as error:
        raise AttemptManifestMaterializationError(
            "registered execution receipt is unavailable below authorized PIC root"
        ) from error
    _require(
        resolved == supplied,
        "registered execution receipt must be canonical without symlink aliases",
    )
    try:
        descriptor = os.open(supplied, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except OSError as error:
        raise AttemptManifestMaterializationError(
            "registered execution receipt is unavailable"
        ) from error
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and stat.S_IMODE(before.st_mode) & 0o222 == 0,
            "registered execution receipt must be one read-only regular file",
        )
        payload = bytearray()
        while chunk := os.read(descriptor, 1024 * 1024):
            payload.extend(chunk)
        after = os.fstat(descriptor)
        current = os.stat(supplied, follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        _require(
            identity(before) == identity(after)
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino),
            "registered execution receipt changed while reading",
        )
        return bytes(payload)
    finally:
        os.close(descriptor)


def _write_member(root_fd: int, relative: str, payload: bytes, mode: int = 0o644) -> None:
    try:
        attempt_publisher._write_anchored_member(root_fd, relative, payload, mode)
    except attempt_publisher.PublicationError as error:
        raise AttemptManifestMaterializationError(str(error)) from error


def _raw_snapshot(root_fd: int) -> Any:
    return immutable_orion_tree._scan_anchored_tree(
        root_fd,
        hash_regular=True,
        error_type=AttemptManifestMaterializationError,
        label="Q-011 completed raw output tree",
    )


def _read_raw_member(root_fd: int, entry: Any) -> bytes:
    try:
        payload, _ = attempt_publisher._read_anchored_source_member(
            root_fd, entry.relative, expected_sha256=entry.sha256
        )
    except attempt_publisher.PublicationError as error:
        raise AttemptManifestMaterializationError(str(error)) from error
    return payload


def _product(
    kind: str, relative: str, payload: bytes, time: float | None
) -> dict[str, Any]:
    return {
        "kind": kind,
        "path": relative,
        "sha256": _sha256_bytes(payload),
        "snapshot_time": time,
    }


def _restart_time(
    payload: bytes, cycle_times: Mapping[int, float], *, label: str
) -> float:
    header = payload.partition(b"<par_end>\n")[0]
    time_match = _RESTART_TIME.search(header)
    if time_match is not None:
        try:
            value = float(time_match.group(1))
        except ValueError as error:
            raise AttemptManifestMaterializationError(
                f"{label} embeds an invalid restart time"
            ) from error
        _require(math.isfinite(value), f"{label} embeds a non-finite restart time")
        return value
    cycle_match = _RESTART_CYCLE.search(header)
    _require(cycle_match is not None, f"{label} lacks a correlatable restart cycle")
    cycle = int(cycle_match.group(1))
    _require(cycle in cycle_times, f"{label} restart cycle has no matching mesh snapshot")
    return cycle_times[cycle]


def _classify_products(root_fd: int, snapshot: Any) -> list[dict[str, Any]]:
    snapshot_paths = {entry.relative for entry in snapshot.entries}
    _require(
        campaign.MANIFEST_NAME not in snapshot_paths
        and not any(
            path == "bindings"
            or path.startswith("bindings/")
            or path.startswith(".q011-attempt-manifest-staging-")
            for path in snapshot_paths
        ),
        "completed raw output tree already contains manifest bindings",
    )
    entries = {
        entry.relative: entry for entry in snapshot.entries if entry.entry_type == "file"
    }
    _require(entries, "completed raw output tree is empty")
    payloads = {
        relative: _read_raw_member(root_fd, entry)
        for relative, entry in entries.items()
    }
    products: dict[str, dict[str, Any]] = {}
    cycle_times: dict[int, float] = {}
    restart_times: dict[str, float] = {}

    for relative, payload in sorted(payloads.items()):
        match = _MESH_NAME.fullmatch(relative)
        if match is None:
            continue
        kind = match.group(1)
        try:
            dataset = output_primitives.parse_athenak_binary_bytes(
                payload, source=relative
            )
        except ValueError as error:
            raise AttemptManifestMaterializationError(
                f"{relative} is not a valid Athena binary output: {error}"
            ) from error
        existing = cycle_times.setdefault(dataset.cycle, dataset.time)
        _require(
            existing == dataset.time,
            f"mesh outputs disagree on cycle {dataset.cycle}",
        )
        products[relative] = _product(kind, relative, payload, dataset.time)

    for relative, payload in sorted(payloads.items()):
        if _PARTICLE_NAME.fullmatch(relative) is None:
            continue
        try:
            header = campaign._parse_pvtk_execution_header(payload, relative)
        except campaign.QualificationError as error:
            raise AttemptManifestMaterializationError(str(error)) from error
        existing = cycle_times.setdefault(header["cycle"], header["time"])
        _require(
            existing == header["time"],
            f"particle output disagrees on cycle {header['cycle']}",
        )
        products[relative] = _product("prtcl_all", relative, payload, header["time"])

    for relative, payload in sorted(payloads.items()):
        if _RESTART_NAME.fullmatch(relative) is None:
            continue
        restart_times[relative] = _restart_time(payload, cycle_times, label=relative)
        products[relative] = _product(
            "restart", relative, payload, restart_times[relative]
        )

    for relative, payload in sorted(payloads.items()):
        if not relative.endswith(".rst.complete"):
            continue
        restart_path = relative.removesuffix(".complete")
        _require(restart_path in restart_times, f"{relative} lacks its restart payload")
        products[relative] = _product(
            "restart_complete", relative, payload, restart_times[restart_path]
        )

    manifest_times: dict[str, float] = {}
    for relative, payload in sorted(payloads.items()):
        if _RESTART_MANIFEST_NAME.fullmatch(relative) is None:
            continue
        try:
            members = campaign._parse_restart_manifest(payload, relative)
        except campaign.QualificationError as error:
            raise AttemptManifestMaterializationError(str(error)) from error
        times = {restart_times.get(member["path"]) for member in members}
        _require(
            None not in times and len(times) == 1,
            f"{relative} restart members do not identify one snapshot",
        )
        time = next(iter(times))
        _require(type(time) is float, f"{relative} restart time is invalid")
        manifest_times[relative] = time
        products[relative] = _product("restart_manifest", relative, payload, time)

    for relative, payload in sorted(payloads.items()):
        if not relative.endswith(".rst.manifest.complete"):
            continue
        manifest_path = relative.removesuffix(".complete")
        _require(
            manifest_path in manifest_times,
            f"{relative} lacks its restart manifest",
        )
        products[relative] = _product(
            "restart_manifest_complete", relative, payload, manifest_times[manifest_path]
        )

    _require("stdout.txt" in payloads, "completed raw output tree lacks stdout.txt")
    products["stdout.txt"] = _product(
        "stdout", "stdout.txt", payloads["stdout.txt"], None
    )
    unknown = sorted(set(payloads) - set(products))
    _require(
        not unknown,
        f"completed raw output tree contains unsupported files: {unknown}",
    )
    return [products[path] for path in sorted(products)]


def _unlink_same_file_at(parent_fd: int, name: str, descriptor: int) -> None:
    try:
        current = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        expected = os.fstat(descriptor)
    except OSError as error:
        raise AttemptManifestMaterializationError(
            f"cannot inspect transactional member {name}: {error}"
        ) from error
    _require(
        stat.S_ISREG(current.st_mode)
        and (current.st_dev, current.st_ino) == (expected.st_dev, expected.st_ino),
        f"transactional member {name} changed before rollback",
    )
    try:
        os.unlink(name, dir_fd=parent_fd)
    except OSError as error:
        raise AttemptManifestMaterializationError(
            f"cannot rollback transactional member {name}: {error}"
        ) from error


def _remove_tree_at(parent_fd: int, name: str, descriptor: int, *, label: str) -> None:
    try:
        attempt_publisher._remove_anchored_tree_at(
            parent_fd, name, descriptor, label=label
        )
    except attempt_publisher.PublicationError as error:
        raise AttemptManifestMaterializationError(str(error)) from error


def _transactionally_write_manifest(
    raw_fd: int,
    *,
    raw_root: Path,
    bindings: Mapping[str, Mapping[str, str]],
    payloads: Mapping[str, bytes],
    manifest_payload: bytes,
) -> None:
    """Stage all bridge-owned members and withdraw the full tranche on failure."""
    staging_name = f".q011-attempt-manifest-staging-{uuid.uuid4()}"
    staging_fd: int | None = None
    bindings_fd: int | None = None
    manifest_fd: int | None = None
    published_bindings = False
    published_manifest = False
    try:
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        os.mkdir(staging_name, mode=0o700, dir_fd=raw_fd)
        staging_fd = os.open(staging_name, _DIRECTORY_FLAGS, dir_fd=raw_fd)
        for name in sorted(payloads):
            mode = 0o755 if name == "executable" else 0o644
            _write_member(staging_fd, bindings[name]["path"], payloads[name], mode)
        _write_member(staging_fd, campaign.MANIFEST_NAME, manifest_payload)
        expected_files = {
            **{
                binding["path"]: binding["sha256"]
                for binding in bindings.values()
            },
            campaign.MANIFEST_NAME: _sha256_bytes(manifest_payload),
        }
        staged = _raw_snapshot(staging_fd)
        staged_files = {
            entry.relative: entry.sha256
            for entry in staged.entries
            if entry.entry_type == "file"
        }
        staged_directories = {
            entry.relative
            for entry in staged.entries
            if entry.entry_type == "directory"
        }
        _require(
            staged_files == expected_files and staged_directories == {"bindings"},
            "transactional attempt-manifest staging tree drifted",
        )
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        bindings_fd = os.open("bindings", _DIRECTORY_FLAGS, dir_fd=staging_fd)
        manifest_fd = os.open(
            campaign.MANIFEST_NAME,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=staging_fd,
        )
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        attempt_publisher._rename_no_replace_at(
            staging_fd, "bindings", raw_fd, "bindings"
        )
        published_bindings = True
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        attempt_publisher._rename_no_replace_at(
            staging_fd, campaign.MANIFEST_NAME, raw_fd, campaign.MANIFEST_NAME
        )
        published_manifest = True
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        os.rmdir(staging_name, dir_fd=raw_fd)
        os.close(staging_fd)
        staging_fd = None
        os.fsync(raw_fd)
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
    except (
        AttemptManifestMaterializationError,
        OSError,
        attempt_publisher.PublicationError,
    ) as error:
        rollback_errors: list[str] = []
        try:
            if published_manifest and manifest_fd is not None:
                _unlink_same_file_at(raw_fd, campaign.MANIFEST_NAME, manifest_fd)
            if published_bindings and bindings_fd is not None:
                _remove_tree_at(
                    raw_fd,
                    "bindings",
                    bindings_fd,
                    label="transactional attempt-manifest bindings",
                )
            if staging_fd is not None:
                _remove_tree_at(
                    raw_fd,
                    staging_name,
                    staging_fd,
                    label="transactional attempt-manifest staging tree",
                )
                os.close(staging_fd)
                staging_fd = None
            os.fsync(raw_fd)
        except (AttemptManifestMaterializationError, OSError) as rollback_error:
            rollback_errors.append(str(rollback_error))
        message = str(error)
        if rollback_errors:
            message += f"; rollback failed: {'; '.join(rollback_errors)}"
        raise AttemptManifestMaterializationError(message) from error
    finally:
        if manifest_fd is not None:
            os.close(manifest_fd)
        if bindings_fd is not None:
            os.close(bindings_fd)
        if staging_fd is not None:
            os.close(staging_fd)


def _planner_result_receipt(
    planner_root: Path,
    planner_inventory_sha256: str,
    planner_report: Mapping[str, Any],
    plan_payload: bytes,
    internal_receipt_payload: bytes,
) -> bytes:
    internal = _json_payload(
        internal_receipt_payload, label="planner materialization receipt"
    )
    plan = _json_payload(plan_payload, label="campaign plan")
    value = {
        "plan_root": str(planner_root),
        "plan_id": plan["plan_id"],
        "campaign_plan_sha256": _sha256_bytes(plan_payload),
        "materialization_receipt": _binding(
            planner.MATERIALIZATION_RECEIPT_NAME, internal_receipt_payload
        ),
        "materialized_member_inventory_sha256": internal["tree_inventory"]["sha256"],
        "inventory_sha256": planner_inventory_sha256,
        "inventoried_file_count": planner_report["inventoried_file_count"],
        "baseline_attempt_count": plan["baseline_attempt_count"],
        "restart_continuation_carrier_count": 1,
        "recursively_read_only": planner_report["recursively_read_only"],
    }
    return _json_bytes(value)


def materialize_completed_attempt_manifest(
    *,
    planner_root: str | Path,
    planner_inventory_sha256: str,
    raw_output_root: str | Path,
    registered_execution_receipt: str | Path,
    attempt_id: str,
    authorized_pic_root: str | Path = campaign.ORION_BULK_ROOT,
) -> dict[str, Any]:
    """Add a fail-closed completed-attempt manifest without authorizing launch."""
    authorized_root = _canonical_directory(
        authorized_pic_root, label="authorized PIC root"
    )
    immutable_plan = _canonical_directory(planner_root, label="immutable planner root")
    raw_root = _canonical_directory(raw_output_root, label="completed raw output root")
    operational_artifact_dir = _canonical_directory(
        raw_root.parent, label="operational run artifact directory"
    )
    _require(
        raw_root.name == "raw"
        and operational_artifact_dir.parent.parent == authorized_root / "runs",
        "completed raw output root must be runs/<campaign>/<submission-id>/raw",
    )
    raw_fd = os.open(raw_root, _DIRECTORY_FLAGS)
    try:
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        raw_before = _raw_snapshot(raw_fd)
        products = _classify_products(raw_fd, raw_before)
        raw_after = _raw_snapshot(raw_fd)
        immutable_orion_tree._require_same_snapshot(
            raw_before,
            raw_after,
            error_type=AttemptManifestMaterializationError,
            label="Q-011 completed raw output tree",
            phase="attempt-manifest product census",
        )
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        with immutable_orion_tree.staged_verified_frozen_tree(
            immutable_plan,
            planner_inventory_sha256,
            authorized_root=authorized_root,
            error_type=AttemptManifestMaterializationError,
            label="Q-011 immutable qualifying campaign plan",
        ) as (planner_report, snapshot):
            plan_payload = snapshot.member_path("campaign_plan.json").read_bytes()
            plan = _json_payload(plan_payload, label="campaign plan")
            internal_receipt_payload = snapshot.member_path(
                planner.MATERIALIZATION_RECEIPT_NAME
            ).read_bytes()
            planner_receipt_payload = _planner_result_receipt(
                immutable_plan,
                planner_inventory_sha256,
                planner_report,
                plan_payload,
                internal_receipt_payload,
            )
            try:
                planner_receipt = campaign._validate_planner_materialization_receipt(
                    planner_receipt_payload,
                    retained_campaign_plan_payload=plan_payload,
                    retained_campaign_plan_sha256=_sha256_bytes(plan_payload),
                    authorized_pic_root=authorized_root,
                )
            except campaign.QualificationError as error:
                raise AttemptManifestMaterializationError(str(error)) from error
            source_bindings = plan["source_bindings"]
            candidate_manifest_payload = _member_payload(
                snapshot,
                source_bindings["clean_candidate_manifest"],
                label="clean candidate manifest",
            )
            deck_payload = _member_payload(
                snapshot, source_bindings["paper_deck"], label="paper deck"
            )
            preregistration_payload = _member_payload(
                snapshot,
                source_bindings["qualifying_preregistration"],
                label="qualifying preregistration",
            )
            pressure_payload = _member_payload(
                snapshot,
                source_bindings["pressure_selection_receipt"],
                label="selected pressure receipt",
            )
            helper_payload = _member_payload(
                snapshot, plan["helper_source_closure"], label="helper source closure"
            )
            try:
                helper_closure = campaign._validate_helper_source_closure(helper_payload)
            except campaign.QualificationError as error:
                raise AttemptManifestMaterializationError(str(error)) from error
            planner_projection = {
                "candidate_binding": {
                    "clean_candidate_manifest": source_bindings[
                        "clean_candidate_manifest"
                    ],
                },
                "artifact_bindings": {
                    "deck": source_bindings["paper_deck"],
                    "preregistration": source_bindings["qualifying_preregistration"],
                    "selected_pressure_receipt": source_bindings[
                        "pressure_selection_receipt"
                    ],
                },
            }
            try:
                campaign._validate_materialized_planner_graph(
                    plan,
                    planner_projection,
                    planner_receipt,
                    helper_closure,
                )
            except (
                campaign.QualificationError,
                KeyError,
                TypeError,
                ValueError,
            ) as error:
                raise AttemptManifestMaterializationError(
                    f"production planner graph validation failed: {error}"
                ) from error
            selected_descriptors = []
            for binding in plan["baseline_attempt_descriptors"]:
                descriptor_payload = _member_payload(
                    snapshot, binding, label="baseline attempt descriptor"
                )
                descriptor = _json_payload(
                    descriptor_payload, label="baseline attempt descriptor"
                )
                if descriptor.get("attempt_id") == attempt_id:
                    selected_descriptors.append(descriptor)
            _require(
                len(selected_descriptors) == 1,
                "attempt ID does not select exactly one planner descriptor",
            )
            descriptor = selected_descriptors[0]
            _require(
                descriptor.get("status") == "planned_not_authorized",
                "baseline attempt descriptor authorization boundary drifted",
            )
            contract_payload = _member_payload(
                snapshot, descriptor["launch_contract"], label="attempt contract"
            )
            contract = _json_payload(contract_payload, label="attempt contract")
            _require(
                contract.get("launch_authorized") is False
                and contract.get("scheduler_submission_authorized") is False
                and contract.get("live_policy_mutation_authorized") is False,
                "attempt contract must remain launch-prohibited",
            )
            attempt_root = Path(descriptor["authorized_orion_attempt_root"])
            _require(
                attempt_root.is_absolute(),
                "planner-authorized attempt root must be absolute",
            )
            _absent_retained_attempt_root(
                attempt_root, authorized_root=authorized_root
            )
            analyzer_records = [
                record
                for record in helper_closure["sources"]
                if record["path"] == "tst/publication/analyze_q011_section54_campaign.py"
            ]
            _require(
                len(analyzer_records) == 1,
                "helper source closure does not bind exactly one campaign analyzer",
            )
            executable_path = Path(plan["candidate_binding"]["executable"]["path"])
            try:
                executable_payload = executable_path.read_bytes()
            except OSError as error:
                raise AttemptManifestMaterializationError(
                    "frozen executable is unavailable"
                ) from error
            _require(
                _sha256_bytes(executable_payload)
                == plan["candidate_binding"]["executable"]["sha256"],
                "frozen executable SHA-256 drifted",
            )
            analyzer_payload = _stable_bound_source_payload(
                campaign.ANALYZER_PATH,
                expected_sha256=analyzer_records[0]["sha256"],
                label="campaign admission analyzer",
            )
            registered_execution_receipt_payload = (
                _stable_registered_execution_receipt(
                    registered_execution_receipt,
                    expected_path=(
                        operational_artifact_dir
                        / "analysis"
                        / campaign.REGISTERED_EXECUTION_RECEIPT_NAME
                    ),
                    authorized_root=authorized_root,
                )
            )
            try:
                validated_registered_execution_receipt = (
                    campaign._validate_registered_execution_receipt(
                        registered_execution_receipt_payload,
                        expected_attempt_id=attempt_id,
                        expected_source_commit=plan["candidate_binding"][
                            "git_commit"
                        ],
                        expected_executable_sha256=plan["candidate_binding"][
                            "executable"
                        ]["sha256"],
                        expected_deck_sha256=source_bindings["paper_deck"]["sha256"],
                        expected_environment_sha256=plan["candidate_binding"][
                            "environment_profile"
                        ]["sha256"],
                        expected_control_plane_version=plan["candidate_binding"][
                            "environment_profile"
                        ]["control_plane_version"],
                        expected_argv=contract["argv"],
                        expected_raw_output_root=f"{attempt_root}/raw",
                        expected_artifact_dir=str(operational_artifact_dir),
                    )
                )
                campaign._validate_registered_execution_receipt_ledger_binding(
                    validated_registered_execution_receipt,
                    authorized_pic_root=authorized_root,
                )
            except campaign.QualificationError as error:
                raise AttemptManifestMaterializationError(str(error)) from error

        payloads = {
            "clean_candidate_manifest": candidate_manifest_payload,
            "executable": executable_payload,
            "deck": deck_payload,
            "analyzer": analyzer_payload,
            "preregistration": preregistration_payload,
            "campaign_plan": plan_payload,
            "planner_materialization_receipt": planner_receipt_payload,
            "attempt_contract": contract_payload,
            "selected_pressure_receipt": pressure_payload,
            "analyzer_helper_source_closure_manifest": helper_payload,
            "registered_execution_receipt": registered_execution_receipt_payload,
        }
        binding_paths = dict(_BINDING_PAYLOADS)
        binding_paths["executable"] = "bindings/athena"
        bindings = {
            name: _binding(binding_paths[name], payload)
            for name, payload in payloads.items()
        }
        variant = descriptor["variant"]
        seed = descriptor["qualifying_seed"]
        try:
            model_binding = model.variant_binding(variant)
        except model.ModelContractError as error:
            raise AttemptManifestMaterializationError(str(error)) from error
        manifest = {
            "schema_version": 1,
            "record_type": MANIFEST_RECORD_TYPE,
            "campaign_id": campaign.CAMPAIGN_ID,
            "qualification_scope": campaign.QUALIFICATION_SCOPE,
            "authorized_orion_campaign_root": str(raw_root),
            "run_identity": {
                "variant": variant,
                "seed": seed,
                "physical_mode": descriptor["physical_mode"],
                "attempt_id": attempt_id,
            },
            "candidate_binding": {
                "git_commit": plan["candidate_binding"]["git_commit"],
                "source_bundle_sha256": plan["candidate_binding"]["source_bundle_sha256"],
                "clean_candidate_manifest": bindings["clean_candidate_manifest"],
            },
            "artifact_bindings": {
                name: bindings[name]
                for name in campaign._BINDING_NAMES[1:]
            },
            "attempt_identity": {
                "campaign_plan_sha256": bindings["campaign_plan"]["sha256"],
                "planner_materialization_receipt_sha256": bindings[
                    "planner_materialization_receipt"
                ]["sha256"],
                "attempt_contract_sha256": bindings["attempt_contract"]["sha256"],
                "selected_pressure_receipt_sha256": bindings[
                    "selected_pressure_receipt"
                ]["sha256"],
                "model_launch_overrides": list(model_binding.model_launch_overrides),
                "seed_overrides": {
                    name: seed for name in campaign._SEED_OVERRIDE_NAMES
                },
                "analyzer_helper_source_closure_manifest_sha256": bindings[
                    "analyzer_helper_source_closure_manifest"
                ]["sha256"],
                "registered_execution_receipt_sha256": bindings[
                    "registered_execution_receipt"
                ]["sha256"],
            },
            "products": products,
        }
        campaign._validate_manifest_schema(manifest, raw_root)
        manifest_payload = _json_bytes(manifest)
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        _transactionally_write_manifest(
            raw_fd,
            raw_root=raw_root,
            bindings=bindings,
            payloads=payloads,
            manifest_payload=manifest_payload,
        )
        _require_same_directory(raw_root, raw_fd, label="completed raw output root")
        return {
            "schema_version": 1,
            "record_type": "q011_section54_completed_attempt_manifest_materialization",
            "qualification_effect": "manifest_only_no_launch_no_claim_closure",
            "raw_output_root": str(raw_root),
            "attempt_id": attempt_id,
            "planner_root": str(immutable_plan),
            "planner_inventory_sha256": planner_inventory_sha256,
            "campaign_manifest": _binding(campaign.MANIFEST_NAME, manifest_payload),
            "product_count": len(products),
            "launch_authorized": False,
            "scheduler_submission_authorized": False,
            "claim_closure_authorized": False,
        }
    finally:
        os.close(raw_fd)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--planner-root", required=True)
    parser.add_argument("--planner-inventory-sha256", required=True)
    parser.add_argument("--raw-output-root", required=True)
    parser.add_argument("--registered-execution-receipt", required=True)
    parser.add_argument("--attempt-id", required=True)
    parser.add_argument("--authorized-pic-root", default=str(campaign.ORION_BULK_ROOT))
    arguments = parser.parse_args(argv)
    result = materialize_completed_attempt_manifest(
        planner_root=arguments.planner_root,
        planner_inventory_sha256=arguments.planner_inventory_sha256,
        raw_output_root=arguments.raw_output_root,
        registered_execution_receipt=arguments.registered_execution_receipt,
        attempt_id=arguments.attempt_id,
        authorized_pic_root=arguments.authorized_pic_root,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
