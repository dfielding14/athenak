#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Reconcile Q043 and publish deterministic registered-execution evidence."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import base64
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import PinnedDirectoryAncestry, atomic_write_bytes_at
from control_plane_common import durable_mkdir_parents, fsync_directory
from control_plane_common import launch_contract_sha256, project_home_ledger_root
from control_plane_common import read_json_bytes, read_stable_regular_file_below
from control_plane_common import record_for_role, require_canonical_path_below
from control_plane_common import require_same_directory, validate_launch_contract
from control_plane_common import verify_installed_control_plane
from ledger import require_explicit_genesis, validate_mirrored_state, validate_receipts
from reconcile_frontier_job import reconcile
from validate_and_reserve_frontier_job import _require_run_artifact_dir


SCRIPT_DIR = Path(__file__).absolute().parent
ENTRYPOINT_NAME = "reconcile_q043_registered_execution.py"
TRAMPOLINE_ENTRYPOINT = "launch_trampoline.py"
REGISTERED_CAMPAIGN = "q043_registered_execution_raw_oracle_successor_v1"
CAMPAIGN_ID = "Q043-BELL-DEPOSITED-J-OVER-C-VOLUME-AWARE"
AUTHORIZATION_PREFIX = "q043-"
EXECUTION_RECEIPT_NAME = "q043_registered_execution_receipt.json"
EXECUTION_RECEIPT_RECORD_TYPE = "q043_reconciled_registered_execution_receipt"
TERMINAL_RECEIPT_NAME = "terminal_receipt.json"
TERMINAL_RECEIPT_RECORD_TYPE = "q043_registered_execution_terminal_receipt"
PROJECT_HOME_MIRROR_NAMESPACE = Path("ledger/q043_registered_execution_receipts")
ARTIFACT_INVENTORY_NAME = "artifact_inventory.json"
TRAMPOLINE_COMPLETION_NAMESPACE = Path("ledger/trampoline_completion_receipts")
TRAMPOLINE_COMPLETION_NAME = "trampoline_completion_receipt.json"
TRAMPOLINE_COMPLETION_RECORD_TYPE = "trusted_trampoline_completion_receipt"
TRAMPOLINE_COMPLETION_LEDGER_RECORD_TYPE = (
    "trusted_trampoline_completion_ledger_binding"
)
REQUIRED_FIELDS = ("prtcl_rho", "prtcl_jx", "prtcl_jy", "prtcl_jz")
REQUIRED_CYCLES = (0, 1)
_SHA256 = re.compile(r"[0-9a-f]{64}")
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}"
)
_JOB_ID = re.compile(r"[1-9][0-9]*")
_DIRECTORY_OPEN_FLAGS = (
    os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
)


def _json_bytes(value: dict[str, object]) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _sha256(value: object, *, label: str) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise ValueError(f"{label} is not one lowercase SHA-256")
    return value


def _safe_relative(value: object, *, label: str) -> str:
    text = str(value)
    relative = PurePosixPath(text)
    if (
        not text
        or relative.is_absolute()
        or relative.as_posix() != text
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise ValueError(f"{label} is not one safe relative path")
    return text


def _read_json(
    path: Path, *, root: Path, label: str, read_only: bool = True
) -> tuple[dict[str, object], bytes]:
    payload = read_stable_regular_file_below(
        path, root, require_read_only_mode=read_only
    )
    return read_json_bytes(payload, label=label), payload


def _inventory_file_digest(inventory: dict[str, object], name: str) -> str:
    records = inventory.get("files")
    if not isinstance(records, list):
        raise ValueError("Installed control-plane inventory files are malformed")
    matches = [
        record
        for record in records
        if isinstance(record, dict) and record.get("path") == name
    ]
    if len(matches) != 1:
        raise ValueError(f"Installed control plane lacks required entrypoint: {name}")
    return _sha256(matches[0].get("sha256"), label=f"{name} digest")


def _producer_binding(inventory: dict[str, object]) -> dict[str, object]:
    return {
        "entrypoint": ENTRYPOINT_NAME,
        "entrypoint_sha256": _inventory_file_digest(inventory, ENTRYPOINT_NAME),
        "launch_trampoline_sha256": _inventory_file_digest(
            inventory, TRAMPOLINE_ENTRYPOINT
        ),
        "control_plane_version": _sha256(
            inventory.get("version"), label="Q043 producer control-plane version"
        ),
    }


def _require_bootstrapped_installed_producer(
    control_plane_dir: Path,
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path,
) -> None:
    if not getattr(_sys, "_pic_control_plane_bootstrapped", False):
        raise ValueError("Q043 reconciliation requires bootstrapped installed execution")
    captured = globals().get("_PIC_CAPTURED_CONTROL_PLANE_BINDING")
    if captured is not None:
        producer = _producer_binding(inventory)
        if (
            not isinstance(captured, dict)
            or set(captured) != {"directory", "version", "filename", "sha256"}
            or captured.get("directory") != str(Path(os.path.abspath(control_plane_dir)))
            or captured.get("directory") != str(Path(os.path.abspath(SCRIPT_DIR)))
            or captured.get("version") != inventory.get("version")
            or captured.get("filename") != ENTRYPOINT_NAME
            or captured.get("sha256") != producer["entrypoint_sha256"]
        ):
            raise ValueError(
                "Q043 captured executing entrypoint differs from verified generation"
            )
        return
    root = Path(os.path.abspath(authorized_pic_root))
    installed = require_canonical_path_below(control_plane_dir, root)
    script_dir = Path(os.path.abspath(SCRIPT_DIR))
    entrypoint = Path(os.path.abspath(__file__))
    if (
        installed != script_dir
        or installed.resolve(strict=True) != script_dir
        or entrypoint != installed / ENTRYPOINT_NAME
        or entrypoint.resolve(strict=True) != entrypoint
    ):
        raise ValueError(
            "Q043 executing entrypoint is not the exact installed control-plane generation"
        )
    payload = read_stable_regular_file_below(
        entrypoint, root, require_read_only_mode=True
    )
    if hashlib.sha256(payload).hexdigest() != _producer_binding(inventory)[
        "entrypoint_sha256"
    ]:
        raise ValueError("Q043 executing entrypoint digest differs from installed inventory")


def _read_exact_identity(
    path: Path, payload: bytes, *, authorized_root: Path
) -> tuple[int, int]:
    require_canonical_path_below(path, authorized_root)
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_mode & 0o222
        ):
            raise ValueError("Q043 reconciliation evidence is not one read-only file")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            observed = stream.read()
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_nlink,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        if (
            observed != payload
            or identity(before) != identity(after)
            or (after.st_dev, after.st_ino) != (current.st_dev, current.st_ino)
        ):
            raise ValueError("Q043 reconciliation evidence retry identity differs")
        return before.st_dev, before.st_ino
    finally:
        os.close(descriptor)


def _publish_exact(
    path: Path, payload: bytes, *, authorized_root: Path
) -> tuple[int, int]:
    durable_mkdir_parents(path.parent, mode=0o700, root=authorized_root)
    with PinnedDirectoryAncestry(
        path.parent, root=Path(os.path.abspath(authorized_root))
    ) as ancestry:
        metadata = os.fstat(ancestry.descriptor)
        if (
            not stat.S_ISDIR(metadata.st_mode)
            or stat.S_IMODE(metadata.st_mode) not in {0o500, 0o700}
        ):
            raise ValueError("Q043 evidence directory is not private")
        try:
            os.stat(path.name, dir_fd=ancestry.descriptor, follow_symlinks=False)
        except FileNotFoundError:
            atomic_write_bytes_at(
                ancestry.descriptor,
                path.name,
                payload,
                mode=0o444,
                replace=False,
                post_publish_check=ancestry.require_same,
            )
        ancestry.require_same()
        identity = _read_exact_identity(
            path, payload, authorized_root=authorized_root
        )
        ancestry.require_same()
        return identity


def _canonical_project_home_root(authorized_project_home_root: Path) -> Path:
    root = Path(os.path.abspath(authorized_project_home_root))
    try:
        resolved = root.resolve(strict=True)
    except OSError as error:
        raise ValueError("Canonical Project Home root is unavailable") from error
    if resolved != root or not root.is_dir():
        raise ValueError("Project Home evidence root must use the canonical path")
    return root


def _project_home_mirror_paths(
    submission_id: str, *, authorized_project_home_root: Path
) -> tuple[Path, Path]:
    root = _canonical_project_home_root(authorized_project_home_root)
    mirror_root = root / PROJECT_HOME_MIRROR_NAMESPACE / submission_id
    return (
        mirror_root / TERMINAL_RECEIPT_NAME,
        mirror_root / EXECUTION_RECEIPT_NAME,
    )


def _trampoline_completion_paths(
    submission_id: str,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> tuple[Path, Path]:
    relative = (
        TRAMPOLINE_COMPLETION_NAMESPACE
        / submission_id
        / TRAMPOLINE_COMPLETION_NAME
    )
    return authorized_pic_root / relative, authorized_project_home_root / relative


def _completion_binding(value: object, *, label: str) -> dict[str, object]:
    if (
        not isinstance(value, dict)
        or set(value) != {"path", "sha256", "byte_count", "filesystem_identity"}
        or not isinstance(value.get("path"), str)
        or not isinstance(value.get("byte_count"), int)
        or int(value["byte_count"]) < 0
        or not isinstance(value.get("filesystem_identity"), dict)
        or set(value["filesystem_identity"]) != {"device", "inode"}
        or any(
            type(value["filesystem_identity"][key]) is not int
            or int(value["filesystem_identity"][key]) < 0
            for key in ("device", "inode")
        )
    ):
        raise ValueError(f"{label} is malformed")
    return {
        "path": str(value["path"]),
        "sha256": _sha256(value["sha256"], label=f"{label}/sha256"),
        "byte_count": int(value["byte_count"]),
        "filesystem_identity": {
            "device": int(value["filesystem_identity"]["device"]),
            "inode": int(value["filesystem_identity"]["inode"]),
        },
    }


def _trusted_trampoline_completion(
    manifest: dict[str, object],
    event: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    submission_id = str(event.get("submission_id", ""))
    orion_path, project_home_path = _trampoline_completion_paths(
        submission_id,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    orion, orion_payload = _read_json(
        orion_path,
        root=authorized_pic_root,
        label="trusted Orion trampoline completion receipt",
    )
    project_home, project_home_payload = _read_json(
        project_home_path,
        root=authorized_project_home_root,
        label="trusted Project Home trampoline completion receipt",
    )
    if orion_payload != project_home_payload or orion != project_home:
        raise ValueError("Paired trampoline completion receipts differ")
    if (
        orion_path.parent.stat(follow_symlinks=False).st_mode & 0o222
        or project_home_path.parent.stat(follow_symlinks=False).st_mode & 0o222
    ):
        raise ValueError("Trampoline completion receipt directory is mutable")
    orion_identity = _read_exact_identity(
        orion_path, orion_payload, authorized_root=authorized_pic_root
    )
    project_home_identity = _read_exact_identity(
        project_home_path,
        project_home_payload,
        authorized_root=authorized_project_home_root,
    )
    if orion_identity == project_home_identity:
        raise ValueError("Paired trampoline completion receipts reuse one filesystem object")
    ledger_binding = event.get("trampoline_completion")
    if (
        not isinstance(ledger_binding, dict)
        or ledger_binding.get("receipt_sha256")
        != hashlib.sha256(orion_payload).hexdigest()
        or ledger_binding.get("receipt_byte_count") != len(orion_payload)
        or ledger_binding.get("paired_receipts")
        != {
            "orion": {
                "path": str(orion_path),
                "parent_identity": {
                    "device": orion_path.parent.stat(follow_symlinks=False).st_dev,
                    "inode": orion_path.parent.stat(follow_symlinks=False).st_ino,
                },
                "filesystem_identity": {
                    "device": orion_identity[0],
                    "inode": orion_identity[1],
                },
            },
            "project_home": {
                "path": str(project_home_path),
                "parent_identity": {
                    "device": project_home_path.parent.stat(
                        follow_symlinks=False
                    ).st_dev,
                    "inode": project_home_path.parent.stat(
                        follow_symlinks=False
                    ).st_ino,
                },
                "filesystem_identity": {
                    "device": project_home_identity[0],
                    "inode": project_home_identity[1],
                },
            },
        }
    ):
        raise ValueError(
            "Trampoline completion receipt differs from canonical mirrored-ledger anchor"
        )
    expected_keys = {
        "schema_version",
        "record_type",
        "receipt_role",
        "authority",
        "paired_paths",
        "receipt_parent_identities",
        "manifest_path",
        "manifest_sha256",
        "reservation_id",
        "submission_id",
        "slurm_job_id",
        "control_plane_version",
        "campaign",
        "test_id",
        "registered_science_authorization_id",
        "artifact_dir",
        "artifact_root_identity",
        "artifact_inventory",
        "artifact_records",
        "mandatory_stdout_stderr",
        "execution_binding",
    }
    authority = orion.get("authority")
    paired_paths = orion.get("paired_paths")
    receipt_parent_identities = orion.get("receipt_parent_identities")
    root_identity = orion.get("artifact_root_identity")
    inventory = orion.get("artifact_inventory")
    records = orion.get("artifact_records")
    mandatory = orion.get("mandatory_stdout_stderr")
    execution = orion.get("execution_binding")
    contract = validate_launch_contract(manifest.get("launch_contract"))
    mandatory_paths = sorted(
        {
            str(action[key])
            for action in contract["actions"]
            for key in ("stdout_artifact", "stderr_artifact")
        }
    )
    if (
        set(orion) != expected_keys
        or orion.get("schema_version") != 1
        or orion.get("record_type") != TRAMPOLINE_COMPLETION_RECORD_TYPE
        or orion.get("receipt_role")
        != "paired_immutable_pre_reconciliation_execution_anchor"
        or authority
        != {
            "launch_authorized": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
        }
        or paired_paths
        != {"orion": str(orion_path), "project_home": str(project_home_path)}
        or not isinstance(receipt_parent_identities, dict)
        or set(receipt_parent_identities) != {"orion", "project_home"}
        or any(
            not isinstance(receipt_parent_identities.get(root_name), dict)
            or set(receipt_parent_identities[root_name]) != {"device", "inode"}
            or any(
                type(receipt_parent_identities[root_name].get(key)) is not int
                for key in ("device", "inode")
            )
            for root_name in ("orion", "project_home")
        )
        or (
            orion_path.parent.stat(follow_symlinks=False).st_dev,
            orion_path.parent.stat(follow_symlinks=False).st_ino,
        )
        != (
            receipt_parent_identities["orion"]["device"],
            receipt_parent_identities["orion"]["inode"],
        )
        or (
            project_home_path.parent.stat(follow_symlinks=False).st_dev,
            project_home_path.parent.stat(follow_symlinks=False).st_ino,
        )
        != (
            receipt_parent_identities["project_home"]["device"],
            receipt_parent_identities["project_home"]["inode"],
        )
        or orion.get("manifest_path") != event.get("manifest_path")
        or orion.get("manifest_sha256") != event.get("manifest_sha256")
        or orion.get("reservation_id") != event.get("reservation_id")
        or orion.get("submission_id") != submission_id
        or orion.get("slurm_job_id") != event.get("job_id")
        or orion.get("control_plane_version") != event.get("control_plane_version")
        or orion.get("campaign") != event.get("campaign")
        or orion.get("test_id") != event.get("test_id")
        or orion.get("registered_science_authorization_id")
        != event.get("registered_science_authorization_id")
        or orion.get("artifact_dir") != event.get("artifact_dir")
        or not isinstance(root_identity, dict)
        or set(root_identity) != {"device", "inode"}
        or any(type(root_identity[key]) is not int for key in ("device", "inode"))
        or not isinstance(inventory, dict)
        or set(inventory)
        != {
            "path",
            "sha256",
            "byte_count",
            "filesystem_identity",
            "payload_base64",
        }
        or not isinstance(records, list)
        or not isinstance(mandatory, list)
        or not isinstance(execution, dict)
        or execution
        != {
            "launch_contract_sha256": launch_contract_sha256(contract),
            "job_script_sha256": record_for_role(manifest, "job-script")["sha256"],
            "executable_sha256": record_for_role(manifest, "executable")["sha256"],
            "input_deck_sha256": record_for_role(manifest, "input-deck")["sha256"],
        }
    ):
        raise ValueError("Trampoline completion receipt identity or schema drifted")
    normalized_records = [
        _completion_binding(record, label=f"trampoline completion artifact[{index}]")
        for index, record in enumerate(records)
    ]
    normalized_mandatory = [
        _completion_binding(record, label=f"trampoline completion stdio[{index}]")
        for index, record in enumerate(mandatory)
    ]
    indexed = {record["path"]: record for record in normalized_records}
    if (
        len(indexed) != len(normalized_records)
        or normalized_mandatory != [indexed[path] for path in mandatory_paths]
    ):
        raise ValueError("Trampoline completion mandatory stdout/stderr binding drifted")
    inventory_binding = _completion_binding(
        {
            key: inventory[key]
            for key in ("path", "sha256", "byte_count", "filesystem_identity")
        },
        label="trampoline completion artifact inventory",
    )
    try:
        inventory_payload = base64.b64decode(
            str(inventory["payload_base64"]), validate=True
        )
    except (ValueError, TypeError) as error:
        raise ValueError("Trampoline completion inventory payload is malformed") from error
    if (
        inventory_binding["path"]
        != str(Path(str(event["artifact_dir"])) / ARTIFACT_INVENTORY_NAME)
        or inventory_binding["sha256"]
        != hashlib.sha256(inventory_payload).hexdigest()
        or inventory_binding["byte_count"] != len(inventory_payload)
    ):
        raise ValueError("Trampoline completion inventory payload binding drifted")
    observed_ledger_binding = {
        "schema_version": 1,
        "record_type": TRAMPOLINE_COMPLETION_LEDGER_RECORD_TYPE,
        "authority": authority,
        "receipt_sha256": hashlib.sha256(orion_payload).hexdigest(),
        "receipt_byte_count": len(orion_payload),
        "paired_receipts": {
            "orion": {
                "path": str(orion_path),
                "parent_identity": receipt_parent_identities["orion"],
                "filesystem_identity": {
                    "device": orion_identity[0],
                    "inode": orion_identity[1],
                },
            },
            "project_home": {
                "path": str(project_home_path),
                "parent_identity": receipt_parent_identities["project_home"],
                "filesystem_identity": {
                    "device": project_home_identity[0],
                    "inode": project_home_identity[1],
                },
            },
        },
        "artifact_root_identity": root_identity,
        "artifact_inventory": {
            key: inventory_binding[key]
            for key in ("sha256", "byte_count", "filesystem_identity")
        },
        "artifact_records_sha256": _canonical_sha256(normalized_records),
        "mandatory_stdout_stderr_sha256": _canonical_sha256(normalized_mandatory),
    }
    if event.get("trampoline_completion") != observed_ledger_binding:
        raise ValueError(
            "Trampoline completion receipt differs from canonical mirrored-ledger anchor"
        )
    return {
        "orion_path": str(orion_path),
        "project_home_path": str(project_home_path),
        "sha256": hashlib.sha256(orion_payload).hexdigest(),
        "byte_count": len(orion_payload),
        "receipt": orion,
        "artifact_root_identity": root_identity,
        "artifact_inventory": inventory_binding,
        "artifact_inventory_payload": inventory_payload,
        "artifact_records": normalized_records,
        "mandatory_stdout_stderr": normalized_mandatory,
        "ledger_binding": observed_ledger_binding,
    }


def _read_open_regular(
    descriptor: int,
    parent_descriptor: int,
    name: str,
    *,
    relative: str,
    expected: dict[str, object] | None = None,
) -> bytes:
    os.lseek(descriptor, 0, os.SEEK_SET)
    before = os.fstat(descriptor)
    if (
        not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_mode & 0o222
    ):
        raise ValueError(f"Sealed launch artifact is not one read-only file: {relative}")
    payload = bytearray()
    while chunk := os.read(descriptor, 1024 * 1024):
        payload.extend(chunk)
    after = os.fstat(descriptor)
    current = os.stat(name, dir_fd=parent_descriptor, follow_symlinks=False)
    stable = ("st_dev", "st_ino", "st_mode", "st_nlink", "st_size", "st_mtime_ns", "st_ctime_ns")
    if (
        any(getattr(before, field) != getattr(after, field) for field in stable)
        or (after.st_dev, after.st_ino) != (current.st_dev, current.st_ino)
        or len(payload) != after.st_size
        or (
            expected is not None
            and (
                expected.get("path") != relative
                or expected.get("size") != len(payload)
                or expected.get("sha256")
                != hashlib.sha256(payload).hexdigest()
            )
        )
    ):
        raise ValueError(f"Sealed launch artifact changed while reading: {relative}")
    return bytes(payload)


def _open_artifact_member(
    artifact_dir_fd: int, relative: str
) -> tuple[int, int, str]:
    parts = PurePosixPath(_safe_relative(relative, label="artifact inventory path")).parts
    parent_fd = os.dup(artifact_dir_fd)
    try:
        for part in parts[:-1]:
            observed = os.stat(part, dir_fd=parent_fd, follow_symlinks=False)
            child_fd = os.open(part, _DIRECTORY_OPEN_FLAGS, dir_fd=parent_fd)
            child = os.fstat(child_fd)
            if (
                not stat.S_ISDIR(child.st_mode)
                or child.st_mode & 0o222
                or (observed.st_dev, observed.st_ino) != (child.st_dev, child.st_ino)
            ):
                os.close(child_fd)
                raise ValueError(f"Sealed launch artifact directory changed: {relative}")
            os.close(parent_fd)
            parent_fd = child_fd
        descriptor = os.open(
            parts[-1],
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        return parent_fd, descriptor, parts[-1]
    except BaseException:
        os.close(parent_fd)
        raise


def _scan_frozen_artifact_tree(
    directory_fd: int, *, prefix: tuple[str, ...] = ()
) -> set[str]:
    metadata = os.fstat(directory_fd)
    if not stat.S_ISDIR(metadata.st_mode) or metadata.st_mode & 0o222:
        raise ValueError("Launch artifact tree is not frozen")
    files: set[str] = set()
    for name in sorted(os.listdir(directory_fd)):
        if not prefix and name == "analysis":
            observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            if not stat.S_ISDIR(observed.st_mode) or stat.S_IMODE(observed.st_mode) != 0o700:
                raise ValueError("Launch artifact analysis directory is not private")
            continue
        if not prefix and name == ARTIFACT_INVENTORY_NAME:
            continue
        relative = "/".join((*prefix, name))
        observed = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if stat.S_ISREG(observed.st_mode):
            files.add(relative)
        elif stat.S_ISDIR(observed.st_mode):
            child_fd = os.open(name, _DIRECTORY_OPEN_FLAGS, dir_fd=directory_fd)
            try:
                child = os.fstat(child_fd)
                if (observed.st_dev, observed.st_ino) != (child.st_dev, child.st_ino):
                    raise ValueError(f"Launch artifact directory changed: {relative}")
                files.update(
                    _scan_frozen_artifact_tree(child_fd, prefix=(*prefix, name))
                )
            finally:
                os.close(child_fd)
        else:
            raise ValueError(f"Launch artifact tree has unsupported entry: {relative}")
    return files


def _sealed_artifact_inventory(
    artifact_dir: Path,
    *,
    authorized_pic_root: Path,
    trampoline_completion: dict[str, object],
) -> dict[str, object]:
    root = Path(os.path.abspath(authorized_pic_root))
    artifact_dir = require_canonical_path_below(artifact_dir, root / "runs")
    artifact_dir_fd = os.open(artifact_dir, _DIRECTORY_OPEN_FLAGS)
    inventory_fd: int | None = None
    retained: list[tuple[int, int, str, str, dict[str, object]]] = []
    captured_stdout: bytes | None = None
    try:
        require_same_directory(artifact_dir, artifact_dir_fd, root=root)
        root_metadata = os.fstat(artifact_dir_fd)
        if (
            root_metadata.st_mode & 0o222
            or (root_metadata.st_dev, root_metadata.st_ino)
            != (
                trampoline_completion["artifact_root_identity"]["device"],
                trampoline_completion["artifact_root_identity"]["inode"],
            )
        ):
            raise ValueError("Launch artifact root is not frozen")
        inventory_fd = os.open(
            ARTIFACT_INVENTORY_NAME,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=artifact_dir_fd,
        )
        inventory_payload = _read_open_regular(
            inventory_fd,
            artifact_dir_fd,
            ARTIFACT_INVENTORY_NAME,
            relative=ARTIFACT_INVENTORY_NAME,
        )
        inventory_metadata = os.fstat(inventory_fd)
        trusted_inventory = trampoline_completion["artifact_inventory"]
        if (
            inventory_payload != trampoline_completion["artifact_inventory_payload"]
            or hashlib.sha256(inventory_payload).hexdigest()
            != trusted_inventory["sha256"]
            or len(inventory_payload) != trusted_inventory["byte_count"]
            or (inventory_metadata.st_dev, inventory_metadata.st_ino)
            != (
                trusted_inventory["filesystem_identity"]["device"],
                trusted_inventory["filesystem_identity"]["inode"],
            )
        ):
            raise ValueError("Launch artifact inventory differs from trampoline completion")
        inventory = read_json_bytes(
            inventory_payload, label="sealed launch artifact inventory"
        )
        records = inventory.get("files")
        if (
            set(inventory) != {"schema_version", "files"}
            or inventory.get("schema_version") != 1
            or type(inventory.get("schema_version")) is not int
            or not isinstance(records, list)
        ):
            raise ValueError("Sealed launch artifact inventory schema is malformed")
        normalized: dict[str, dict[str, object]] = {}
        for raw in records:
            if not isinstance(raw, dict) or set(raw) != {"path", "sha256", "size"}:
                raise ValueError("Sealed launch artifact inventory record is malformed")
            relative = _safe_relative(raw["path"], label="artifact inventory path")
            if (
                relative in normalized
                or relative == ARTIFACT_INVENTORY_NAME
                or PurePosixPath(relative).parts[0] == "analysis"
                or type(raw["size"]) is not int
                or raw["size"] < 0
            ):
                raise ValueError("Sealed launch artifact inventory record is malformed")
            normalized[relative] = {
                "path": relative,
                "sha256": _sha256(
                    raw["sha256"], label=f"artifact inventory {relative} digest"
                ),
                "size": raw["size"],
            }
        if _scan_frozen_artifact_tree(artifact_dir_fd) != set(normalized):
            raise ValueError("Sealed launch artifact tree differs from its inventory")
        for relative, record in normalized.items():
            parent_fd, descriptor, name = _open_artifact_member(
                artifact_dir_fd, relative
            )
            retained.append((parent_fd, descriptor, name, relative, record))
            payload = _read_open_regular(
                descriptor,
                parent_fd,
                name,
                relative=relative,
                expected=record,
            )
            if relative == "athena_stdout.txt":
                captured_stdout = payload
        completion_records = {
            record["path"]: record for record in trampoline_completion["artifact_records"]
        }
        observed_completion_records = {
            relative: {
                "path": relative,
                "sha256": record["sha256"],
                "byte_count": record["size"],
                "filesystem_identity": {
                    "device": os.fstat(descriptor).st_dev,
                    "inode": os.fstat(descriptor).st_ino,
                },
            }
            for _, descriptor, _, relative, record in retained
        }
        if observed_completion_records != completion_records:
            raise ValueError("Launch artifact tree differs from trampoline completion")
        require_same_directory(artifact_dir, artifact_dir_fd, root=root)
        if _read_open_regular(
            inventory_fd,
            artifact_dir_fd,
            ARTIFACT_INVENTORY_NAME,
            relative=ARTIFACT_INVENTORY_NAME,
        ) != inventory_payload:
            raise ValueError("Sealed launch artifact inventory changed during verification")
        for parent_fd, descriptor, name, relative, record in retained:
            _read_open_regular(
                descriptor,
                parent_fd,
                name,
                relative=relative,
                expected=record,
            )
        require_same_directory(artifact_dir, artifact_dir_fd, root=root)
        if captured_stdout is None:
            raise ValueError("Sealed launch artifact tree lacks Q043 stdout")
        return {
            "path": str(artifact_dir / ARTIFACT_INVENTORY_NAME),
            "sha256": hashlib.sha256(inventory_payload).hexdigest(),
            "byte_count": len(inventory_payload),
            "records": normalized,
            "athena_stdout_payload": captured_stdout,
            "trampoline_completion": {
                key: trampoline_completion[key]
                for key in ("orion_path", "project_home_path", "sha256", "byte_count")
            },
        }
    finally:
        for parent_fd, descriptor, *_ in reversed(retained):
            os.close(descriptor)
            os.close(parent_fd)
        if inventory_fd is not None:
            os.close(inventory_fd)
        os.close(artifact_dir_fd)


def _raw_relative_path(case_id: str, *, field: str, cycle: int, rank: int, ranks: int) -> str:
    basename = case_id.replace("-", "_")
    filename = f"{basename}.{field}.{cycle:05d}.bin"
    if ranks > 1:
        return f"bin/rank_{rank:08d}/{filename}"
    return f"bin/{filename}"


def _trusted_wrapper_evidence(
    value: object, *, case_id: str, ranks: int, stdout_sha256: str
) -> dict[str, object]:
    if not isinstance(value, bytes):
        raise ValueError("Q043 retained stdout bytes are absent")
    try:
        lines = value.decode("utf-8").splitlines()
    except UnicodeDecodeError as error:
        raise ValueError("Q043 retained stdout is not UTF-8") from error
    expected_rank_line = (
        f"Q043_REGISTERED_EXECUTION case_id={case_id} "
        f"mpi_world_size={ranks} rank_ids={','.join(str(rank) for rank in range(ranks))}"
    )
    expected_exit_line = "Q043_REGISTERED_EXECUTION_EXIT exit_code=0 signal=0"
    if lines.count(expected_rank_line) != 1 or lines.count(expected_exit_line) != 1:
        raise ValueError("Q043 retained stdout lacks exact trusted rank/exit evidence")
    task_pattern = re.compile(
        r"^PIC trusted GPU launch: rank=([0-9]+) host=\S+ "
        r"ROCR_VISIBLE_DEVICES=[0-9]+ "
        r"linkage=libamdhip64,libmpi_amd,libmpi_gtl_hsa$"
    )
    observed_rank_ids = [
        int(match.group(1))
        for line in lines
        if (match := task_pattern.fullmatch(line)) is not None
    ]
    if sorted(observed_rank_ids) != list(range(ranks)) or len(observed_rank_ids) != ranks:
        raise ValueError("Q043 retained stdout task-rank evidence is incomplete or duplicated")
    finite = r"[+-]?(?:[0-9]+(?:[.][0-9]*)?|[.][0-9]+)(?:[eE][+-]?[0-9]+)?"
    time_values = [
        float(match.group(1))
        for line in lines
        if (match := re.fullmatch(rf"time=({finite}) cycle=1", line)) is not None
    ]
    tlim_values = [
        float(match.group(1))
        for line in lines
        if (match := re.fullmatch(rf"tlim=({finite}) nlim=1", line)) is not None
    ]
    if (
        lines.count("Terminating on cycle limit") != 1
        or len(time_values) != 1
        or not math.isfinite(time_values[0])
        or time_values[0] <= 0.0
        or len(tlim_values) != 1
        or not math.isfinite(tlim_values[0])
    ):
        raise ValueError("Q043 retained stdout lacks exact cycle-one termination evidence")
    return {
        "source": "installed_trampoline_retained_stdout_exact_bytes",
        "stdout_sha256": stdout_sha256,
        "required_exact_rank_line": expected_rank_line,
        "required_exact_exit_line": expected_exit_line,
        "observed_world_size": ranks,
        "observed_rank_ids": list(range(ranks)),
        "exit_code": 0,
        "signal": 0,
        "terminal_cycle": 1,
    }


def _structured_launch_evidence(
    manifest: dict[str, object],
    sealed_inventory: dict[str, object],
    *,
    case_id: str,
    artifact_dir: Path,
    producer: dict[str, object],
) -> dict[str, object]:
    contract = validate_launch_contract(manifest.get("launch_contract"))
    actions = contract["actions"]
    if len(actions) != 1:
        raise ValueError("Q043 requires one trusted trampoline Athena action")
    action = actions[0]
    arguments = action["arguments"]
    raw_directories = [
        str(arguments[index + 1]["artifact_directory"])
        for index, argument in enumerate(arguments[:-1])
        if argument == {"literal": "-d"}
        and isinstance(arguments[index + 1], dict)
        and set(arguments[index + 1]) == {"artifact_directory"}
    ]
    input_decks = sum(
        1
        for index, argument in enumerate(arguments[:-1])
        if argument == {"literal": "-i"}
        and arguments[index + 1] == {"snapshot_role": "input-deck"}
    )
    overrides = {
        str(argument["literal"]).partition("=")[0]: str(argument["literal"]).partition("=")[2]
        for argument in arguments
        if isinstance(argument, dict)
        and set(argument) == {"literal"}
        and "=" in str(argument["literal"])
    }
    if (
        len(raw_directories) != 1
        or input_decks != 1
        or overrides.get("time/nlim") != "1"
    ):
        raise ValueError(
            "Q043 launch contract lacks one raw directory, input deck, and cycle-one limit"
        )
    raw_relative = _safe_relative(
        raw_directories[0], label="Q043 raw artifact directory"
    )
    raw_root = artifact_dir.joinpath(*PurePosixPath(raw_relative).parts)
    ranks = action["resources"]["tasks"]
    if type(ranks) is not int or ranks <= 0:
        raise ValueError("Q043 launch contract has malformed MPI task count")
    records = sealed_inventory["records"]
    if not isinstance(records, dict):
        raise ValueError("Q043 sealed artifact records are malformed")
    stdout_record = records.get("athena_stdout.txt")
    if not isinstance(stdout_record, dict):
        raise ValueError("Q043 sealed artifact records lack trusted stdout")
    wrapper = _trusted_wrapper_evidence(
        sealed_inventory.get("athena_stdout_payload"),
        case_id=case_id,
        ranks=ranks,
        stdout_sha256=str(stdout_record["sha256"]),
    )
    expected = [
        _raw_relative_path(
            case_id, field=field, cycle=cycle, rank=rank, ranks=ranks
        )
        for cycle in REQUIRED_CYCLES
        for field in REQUIRED_FIELDS
        for rank in range(ranks)
    ]
    prefix = raw_relative + "/"
    observed = sorted(
        relative[len(prefix) :]
        for relative in records
        if relative.startswith(prefix)
    )
    if observed != sorted(expected):
        raise ValueError("Q043 sealed raw output differs from the exact expected inventory")
    raw_inventory = []
    for cycle in REQUIRED_CYCLES:
        for field in REQUIRED_FIELDS:
            for rank in range(ranks):
                relative = _raw_relative_path(
                    case_id, field=field, cycle=cycle, rank=rank, ranks=ranks
                )
                record = records[prefix + relative]
                raw_inventory.append(
                    {
                        "path": relative,
                        "sha256": record["sha256"],
                        "byte_count": record["size"],
                        "case_id": case_id,
                        "field": field,
                        "cycle": cycle,
                        "rank": rank,
                    }
                )
    command_evidence = {
        "source": "trusted_pre_submit_manifest_and_installed_trampoline",
        "executor": contract["executor"],
        "action": action,
        "launch_contract_sha256": launch_contract_sha256(contract),
        "launch_trampoline_entrypoint": TRAMPOLINE_ENTRYPOINT,
        "launch_trampoline_sha256": producer["launch_trampoline_sha256"],
        "trusted_wrapper_evidence": wrapper,
    }
    mpi_evidence = {
        "source": "trusted_pre_submit_manifest_and_installed_trampoline_stdout",
        **action["resources"],
        "observed_world_size": wrapper["observed_world_size"],
        "observed_rank_ids": wrapper["observed_rank_ids"],
    }
    return {
        "raw_root": raw_root,
        "raw_inventory": raw_inventory,
        "raw_inventory_sha256": _canonical_sha256(raw_inventory),
        "command_evidence": command_evidence,
        "mpi_evidence": mpi_evidence,
        "terminal_cycle": wrapper["terminal_cycle"],
    }


def _mirror_ack(
    event: dict[str, object],
    records: list[dict[str, object]],
    *,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    authorized_pic_root: Path,
) -> dict[str, object]:
    receipts = validate_receipts(
        receipts_jsonl,
        records,
        mirror_jsonl=mirror_jsonl,
        mirror_transport="filesystem_copy",
        root=authorized_pic_root,
    )
    matches = [
        receipt
        for receipt in receipts
        if receipt.get("mirrored_event_sha256") == event.get("event_sha256")
    ]
    if len(matches) != 1:
        raise ValueError("Q043 reconciliation lacks one canonical mirror acknowledgment")
    return matches[0]


def _q043_event(
    records: list[dict[str, object]],
    *,
    job_id: str,
    producer_control_plane_version: str,
) -> dict[str, object]:
    matches = [
        record
        for record in records
        if record.get("event_type") == "reconciliation"
        and record.get("job_id") == job_id
        and record.get("campaign") == REGISTERED_CAMPAIGN
    ]
    if len(matches) != 1:
        raise ValueError("Expected one Q043 reconciliation event for the scheduler job")
    event = matches[0]
    authorization = event.get("registered_science_authorization_id")
    if (
        event.get("submission_scope") != "registered_science"
        or not isinstance(authorization, str)
        or not authorization.startswith(AUTHORIZATION_PREFIX)
        or event.get("state") != "COMPLETED"
        or event.get("scheduler_exit_code") != "0:0"
        or event.get("reconciled") is not True
        or event.get("reconciled_by_control_plane_version")
        != producer_control_plane_version
        or not isinstance(event.get("event_sha256"), str)
        or _SHA256.fullmatch(str(event["event_sha256"])) is None
    ):
        raise ValueError("Q043 reconciliation event is not a successful registered execution")
    return event


def derive_q043_registered_execution_evidence(
    event: dict[str, object],
    mirror_ack: dict[str, object],
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> tuple[Path, bytes, Path, bytes, Path, Path]:
    """Derive exact Q043 receipts from controller-owned immutable evidence."""
    root = Path(os.path.abspath(authorized_pic_root))
    producer = _producer_binding(inventory)
    if (
        event.get("reconciled_by_control_plane_version")
        != producer["control_plane_version"]
        or mirror_ack.get("mirrored_event_sha256") != event.get("event_sha256")
        or not isinstance(mirror_ack.get("mirror_ack_sha256"), str)
        or _SHA256.fullmatch(str(mirror_ack["mirror_ack_sha256"])) is None
    ):
        raise ValueError("Q043 producer, reconciliation, and mirror acknowledgment differ")
    manifest_path = require_canonical_path_below(
        Path(str(event.get("manifest_path", ""))), root / "manifests"
    )
    manifest, manifest_payload = _read_json(
        manifest_path, root=root, label="Q043 pre-submit manifest"
    )
    if hashlib.sha256(manifest_payload).hexdigest() != event.get("manifest_sha256"):
        raise ValueError("Q043 pre-submit manifest digest differs from reconciliation")
    artifact_dir = _require_run_artifact_dir(manifest)
    if (
        manifest.get("submission_scope") != "registered_science"
        or manifest.get("campaign") != REGISTERED_CAMPAIGN
        or manifest.get("submission_id") != event.get("submission_id")
        or manifest.get("test_id") != event.get("test_id")
        or manifest.get("control_plane_version") != event.get("control_plane_version")
        or manifest.get("git_commit") != event.get("git_commit")
        or manifest.get("registered_science_authorization_id")
        != event.get("registered_science_authorization_id")
        or event.get("artifact_dir") != str(artifact_dir)
        or event.get("state") != "COMPLETED"
        or event.get("scheduler_exit_code") != "0:0"
    ):
        raise ValueError("Q043 pre-submit manifest differs from successful reconciliation")
    case_id = str(event["test_id"])
    submission_id = str(event["submission_id"])
    job_id = str(event["job_id"])
    if _UUID.fullmatch(submission_id) is None or _JOB_ID.fullmatch(job_id) is None:
        raise ValueError("Q043 reconciliation has malformed scheduler identity")
    trampoline_completion = _trusted_trampoline_completion(
        manifest,
        event,
        authorized_pic_root=root,
        authorized_project_home_root=authorized_project_home_root,
    )
    sealed_inventory = _sealed_artifact_inventory(
        artifact_dir,
        authorized_pic_root=root,
        trampoline_completion=trampoline_completion,
    )
    launch = _structured_launch_evidence(
        manifest,
        sealed_inventory,
        case_id=case_id,
        artifact_dir=artifact_dir,
        producer=producer,
    )
    terminal_mirror_path, receipt_mirror_path = _project_home_mirror_paths(
        submission_id,
        authorized_project_home_root=authorized_project_home_root,
    )
    terminal = {
        "schema_version": 1,
        "record_type": TERMINAL_RECEIPT_RECORD_TYPE,
        "campaign_id": CAMPAIGN_ID,
        "case_id": case_id,
        "submission_id": submission_id,
        "slurm_job_id": job_id,
        "slurm_terminal_state": event["state"],
        "slurm_exit_code": event["scheduler_exit_code"],
        "terminal_cycle": launch["terminal_cycle"],
        "registered_mpi_tasks": launch["mpi_evidence"]["tasks"],
        "artifact_inventory_sha256": sealed_inventory["sha256"],
        "trampoline_completion_receipt_sha256": trampoline_completion["sha256"],
        "raw_inventory_sha256": launch["raw_inventory_sha256"],
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "project_home_mirror_path": str(terminal_mirror_path),
        "producer": producer,
    }
    terminal_path = artifact_dir / "analysis" / TERMINAL_RECEIPT_NAME
    terminal_payload = _json_bytes(terminal)
    candidate_manifest_path = Path(str(manifest.get("clean_candidate_manifest_path", "")))
    candidate, candidate_payload = _read_json(
        candidate_manifest_path, root=root, label="Q043 clean-candidate manifest"
    )
    source = candidate.get("source")
    if (
        hashlib.sha256(candidate_payload).hexdigest()
        != manifest.get("clean_candidate_manifest_sha256")
        or not isinstance(source, dict)
        or source.get("git_commit") != event.get("git_commit")
    ):
        raise ValueError("Q043 clean candidate differs from pre-submit manifest")
    receipt = {
        "schema_version": 1,
        "record_type": EXECUTION_RECEIPT_RECORD_TYPE,
        "receipt_role": "immutable_reconciled_registered_execution",
        "registration_scope": "registered_science",
        "reconciled": True,
        "campaign_id": CAMPAIGN_ID,
        "case_id": case_id,
        "reservation_id": event["reservation_id"],
        "submission_id": submission_id,
        "reconciliation_event_sha256": event["event_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "control_plane_version": event["control_plane_version"],
        "project_home_mirrors": {
            "registered_execution_receipt_path": str(receipt_mirror_path),
            "terminal_receipt_path": str(terminal_mirror_path),
        },
        "producer": producer,
        "registered_science_authorization_id": event[
            "registered_science_authorization_id"
        ],
        "source_commit": event["git_commit"],
        "source_bundle_sha256": source["source_bundle_sha256"],
        "source_archive_sha256": source["archive_sha256"],
        "clean_candidate_manifest_sha256": manifest[
            "clean_candidate_manifest_sha256"
        ],
        "executable_sha256": record_for_role(manifest, "executable")["sha256"],
        "environment_sha256": record_for_role(manifest, "environment-profile")[
            "sha256"
        ],
        "deck_sha256": record_for_role(manifest, "input-deck")["sha256"],
        "command_evidence": launch["command_evidence"],
        "mpi_evidence": launch["mpi_evidence"],
        "slurm_job_id": job_id,
        "slurm_terminal_state": event["state"],
        "slurm_exit_code": event["scheduler_exit_code"],
        "terminal_cycle": launch["terminal_cycle"],
        "raw_output_root": str(launch["raw_root"]),
        "artifact_dir": str(artifact_dir),
        "artifact_inventory": {
            key: sealed_inventory[key] for key in ("path", "sha256", "byte_count")
        },
        "trampoline_completion_receipt": {
            key: trampoline_completion[key]
            for key in (
                "orion_path",
                "project_home_path",
                "sha256",
                "byte_count",
            )
        },
        "terminal_receipt_sha256": hashlib.sha256(terminal_payload).hexdigest(),
        "pre_submit_manifest_path": str(manifest_path),
        "pre_submit_manifest_sha256": event["manifest_sha256"],
        "raw_inventory": launch["raw_inventory"],
        "raw_inventory_sha256": launch["raw_inventory_sha256"],
    }
    receipt_path = artifact_dir / "analysis" / EXECUTION_RECEIPT_NAME
    return (
        terminal_path,
        terminal_payload,
        receipt_path,
        _json_bytes(receipt),
        terminal_mirror_path,
        receipt_mirror_path,
    )


def publish_q043_registered_execution_evidence(
    event: dict[str, object],
    mirror_ack: dict[str, object],
    inventory: dict[str, object],
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Publish or verify deterministic Q043 evidence under paired trusted roots."""
    (
        terminal_path,
        terminal_payload,
        receipt_path,
        receipt_payload,
        terminal_mirror_path,
        receipt_mirror_path,
    ) = derive_q043_registered_execution_evidence(
        event,
        mirror_ack,
        inventory,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    terminal_identity = _publish_exact(
        terminal_path, terminal_payload, authorized_root=authorized_pic_root
    )
    receipt_identity = _publish_exact(
        receipt_path, receipt_payload, authorized_root=authorized_pic_root
    )
    terminal_mirror_identity = _publish_exact(
        terminal_mirror_path,
        terminal_payload,
        authorized_root=authorized_project_home_root,
    )
    receipt_mirror_identity = _publish_exact(
        receipt_mirror_path,
        receipt_payload,
        authorized_root=authorized_project_home_root,
    )
    os.chmod(receipt_mirror_path.parent, 0o500, follow_symlinks=False)
    fsync_directory(receipt_mirror_path.parent)
    if receipt_mirror_path.parent.stat(follow_symlinks=False).st_mode & 0o222:
        raise ValueError("Q043 Project Home evidence mirror directory is mutable")
    if (
        terminal_identity == terminal_mirror_identity
        or receipt_identity == receipt_mirror_identity
    ):
        raise ValueError("Q043 Orion and Project Home evidence reuse one filesystem object")
    for path, payload, identity, root in (
        (terminal_path, terminal_payload, terminal_identity, authorized_pic_root),
        (receipt_path, receipt_payload, receipt_identity, authorized_pic_root),
        (
            terminal_mirror_path,
            terminal_payload,
            terminal_mirror_identity,
            authorized_project_home_root,
        ),
        (
            receipt_mirror_path,
            receipt_payload,
            receipt_mirror_identity,
            authorized_project_home_root,
        ),
    ):
        if _read_exact_identity(path, payload, authorized_root=root) != identity:
            raise ValueError("Q043 paired evidence identity changed during publication")
    receipt = read_json_bytes(
        receipt_payload, label="Q043 registered execution receipt"
    )
    return {
        "terminal_receipt": {
            "orion_path": str(terminal_path),
            "project_home_mirror_path": str(terminal_mirror_path),
            "sha256": hashlib.sha256(terminal_payload).hexdigest(),
        },
        "registered_execution_receipt": {
            "orion_path": str(receipt_path),
            "project_home_mirror_path": str(receipt_mirror_path),
            "sha256": hashlib.sha256(receipt_payload).hexdigest(),
        },
        "artifact_inventory_sha256": receipt["artifact_inventory"]["sha256"],
        "trampoline_completion_receipt_sha256": receipt[
            "trampoline_completion_receipt"
        ]["sha256"],
        "raw_inventory_sha256": receipt["raw_inventory_sha256"],
        "reconciliation_mirror_ack_sha256": mirror_ack["mirror_ack_sha256"],
        "producer": _producer_binding(inventory),
    }


def reconcile_q043(
    *,
    job_id: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    """Reconcile one Q043 job through the exact bootstrapped installed generation."""
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    _require_bootstrapped_installed_producer(
        control_plane_dir, inventory, authorized_pic_root=authorized_pic_root
    )
    paired = verify_installed_control_plane(
        Path(os.path.abspath(authorized_project_home_root))
        / "control_plane"
        / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    if paired != inventory:
        raise ValueError("Q043 producer installed control-plane pair differs")
    pre_reconciliation_records = validate_mirrored_state(
        ledger_jsonl,
        receipts_jsonl,
        mirror_jsonl,
        ledger_root=authorized_pic_root,
        receipts_root=authorized_pic_root,
        mirror_root=project_home_ledger_root(authorized_project_home_root),
    )
    require_explicit_genesis(pre_reconciliation_records)
    # Generic reconciliation remains available for legacy jobs, but a Q043 job
    # terminalized without this predecessor is intentionally stranded as
    # non-admissible evidence; it cannot acquire a completion anchor afterward.
    completion_events = [
        record
        for record in pre_reconciliation_records
        if record.get("event_type") == "trampoline_completion"
        and record.get("job_id") == job_id
        and record.get("campaign") == REGISTERED_CAMPAIGN
    ]
    if len(completion_events) != 1:
        raise ValueError(
            "Q043 reconciliation requires one canonical trampoline-completion anchor"
        )
    completion_event = completion_events[0]
    if (
        completion_event.get("control_plane_version") != inventory["version"]
        or completion_event.get("state") != "submitted"
        or completion_event.get("reconciled") is not False
    ):
        raise ValueError("Q043 trampoline-completion ledger anchor is invalid")
    _mirror_ack(
        completion_event,
        pre_reconciliation_records,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
    )
    completion_manifest_path = require_canonical_path_below(
        Path(str(completion_event.get("manifest_path", ""))),
        Path(os.path.abspath(authorized_pic_root)) / "manifests",
    )
    completion_manifest, completion_manifest_payload = _read_json(
        completion_manifest_path,
        root=Path(os.path.abspath(authorized_pic_root)),
        label="Q043 trampoline-completion pre-submit manifest",
    )
    if (
        hashlib.sha256(completion_manifest_payload).hexdigest()
        != completion_event.get("manifest_sha256")
    ):
        raise ValueError("Q043 trampoline-completion manifest digest differs")
    _trusted_trampoline_completion(
        completion_manifest,
        completion_event,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    event = reconcile(
        job_id=job_id,
        ledger_jsonl=ledger_jsonl,
        ledger_csv=ledger_csv,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    records = validate_mirrored_state(
        ledger_jsonl,
        receipts_jsonl,
        mirror_jsonl,
        ledger_root=authorized_pic_root,
        receipts_root=authorized_pic_root,
        mirror_root=project_home_ledger_root(authorized_project_home_root),
    )
    require_explicit_genesis(records)
    trusted_event = _q043_event(
        records,
        job_id=job_id,
        producer_control_plane_version=str(inventory["version"]),
    )
    if trusted_event != event:
        raise ValueError("Q043 reconciler result differs from canonical mirrored ledger")
    mirror_ack = _mirror_ack(
        event,
        records,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
    )
    evidence = publish_q043_registered_execution_evidence(
        event,
        mirror_ack,
        inventory,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    return {"reconciliation": event, "evidence": evidence}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    result = reconcile_q043(
        job_id=args.job_id,
        ledger_jsonl=args.ledger_jsonl,
        ledger_csv=args.ledger_csv,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )
    print(json.dumps(result["evidence"], sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
