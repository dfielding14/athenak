#!/usr/bin/env python3
"""Validate and freeze a PIC qualification manifest.

This gate is intentionally separate from the engineering-proxy publication
helpers.  It does not run a campaign or promote a claim.  It verifies that an
already prepared manifest is structurally complete, points at a clean source
candidate, references known claims, and binds existing artifacts by checksum.
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
import sys
from typing import Any
import uuid


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_PATH = READINESS_DIR / "schemas" / "validation_manifest.schema.json"
CLAIMS_PATH = READINESS_DIR / "claims_registry.json"
CONTROL_PLANE_DIR = Path(__file__).with_name("frontier_control_plane")
sys.path.insert(0, str(CONTROL_PLANE_DIR))

from control_plane_common import sha256_bytes  # noqa: E402
from control_plane_common import source_bundle_sha256  # noqa: E402
from control_plane_common import AUTHORIZED_PIC_ROOT  # noqa: E402
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT  # noqa: E402
from control_plane_common import BUILD_PROVENANCE_FILENAMES  # noqa: E402
from control_plane_common import PRODUCTION_RUNTIME_LOADED_MODULES  # noqa: E402
from control_plane_common import PRODUCTION_RUNTIME_MODULEFILES  # noqa: E402
from control_plane_common import PRODUCTION_RUNTIME_MODULEPATH  # noqa: E402
from control_plane_common import active_promotion_path  # noqa: E402
from control_plane_common import read_json_bytes  # noqa: E402
from control_plane_common import read_stable_regular_file  # noqa: E402
from control_plane_common import read_stable_regular_file_below  # noqa: E402
from control_plane_common import require_canonical_path_below  # noqa: E402
from control_plane_common import require_ledger_paths  # noqa: E402
from control_plane_common import require_storage_policy_unlock_snapshot  # noqa: E402
from control_plane_common import utc_datetime  # noqa: E402
from control_plane_common import validate_clean_candidate_bundle  # noqa: E402
from ledger import latest_reservations, require_explicit_genesis  # noqa: E402
from ledger import ledger_lock  # noqa: E402
from ledger import validate_mirrored_state  # noqa: E402

QUALIFYING_EVIDENCE_CLASSES = {
    "physics_validation",
    "sun_bai_2023_reproduction",
    "athenak_production_mode",
    "cross_code_comparison",
    "scoped_state_of_the_art",
}
RUNTIME_ENVIRONMENT_ALLOWLIST_KEYS = [
    "PIC_FRONTIER_PROFILE",
    "HSA_XNACK",
    "MPICH_ENV_DISPLAY",
    "MPICH_VERSION_DISPLAY",
    "MPICH_GPU_SUPPORT_ENABLED",
    "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED",
    "MPICH_OFI_NIC_POLICY",
    "MPICH_GPU_IPC_CACHE_MAX_SIZE",
    "MPICH_MPIIO_HINTS",
    "MPICH_OFI_NUM_CQ_ENTRIES",
    "FI_MR_CACHE_MONITOR",
    "FI_CXI_RX_MATCH_MODE",
    "OMP_NUM_THREADS",
    "SLURM_EXPORT_ENV",
    "ROCM_PATH",
    "LOADEDMODULES",
    "_LMFILES_",
    "MODULEPATH",
]


def _schema_matches(value: object, schema: dict[str, Any]) -> bool:
    try:
        validate_schema(value, schema)
    except ValueError:
        return False
    return True


def validate_schema(
    value: object, schema: dict[str, Any], path: str = "$"
) -> None:
    """Validate the JSON-schema keywords used by PIC readiness records."""
    for child_schema in schema.get("allOf", []):
        validate_schema(value, child_schema, path)
    conditional_schema = schema.get("if")
    if (
        isinstance(conditional_schema, dict)
        and _schema_matches(value, conditional_schema)
        and isinstance(schema.get("then"), dict)
    ):
        validate_schema(value, schema["then"], path)

    if "const" in schema and value != schema["const"]:
        raise ValueError(f"{path} does not match const")
    if "enum" in schema and value not in schema["enum"]:
        raise ValueError(f"{path} is not in enum")

    expected_type = schema.get("type")
    type_matches = {
        "object": isinstance(value, dict),
        "array": isinstance(value, list),
        "string": isinstance(value, str),
        "number": isinstance(value, (int, float)) and not isinstance(value, bool),
    }
    if expected_type is not None and not type_matches[str(expected_type)]:
        raise ValueError(f"{path} is not a {expected_type}")

    if isinstance(value, dict):
        required = schema.get("required", [])
        for key in required:
            if key not in value:
                raise ValueError(f"{path} is missing {key}")
        min_properties = schema.get("minProperties")
        if min_properties is not None and len(value) < min_properties:
            raise ValueError(f"{path} has too few properties")
        properties = schema.get("properties", {})
        if schema.get("additionalProperties") is False:
            extra = set(value) - set(properties)
            if extra:
                raise ValueError(f"{path} has unexpected properties: {extra}")
        for key, child_schema in properties.items():
            if key in value:
                validate_schema(value[key], child_schema, f"{path}.{key}")

    if isinstance(value, list):
        min_items = schema.get("minItems")
        if min_items is not None and len(value) < min_items:
            raise ValueError(f"{path} has too few items")
        item_schema = schema.get("items")
        if item_schema is not None:
            for index, item in enumerate(value):
                validate_schema(item, item_schema, f"{path}[{index}]")

    if isinstance(value, str):
        min_length = schema.get("minLength")
        if min_length is not None and len(value) < min_length:
            raise ValueError(f"{path} is too short")
        pattern = schema.get("pattern")
        if pattern is not None and re.search(str(pattern), value) is None:
            raise ValueError(f"{path} does not match pattern")
        if schema.get("format") == "date-time":
            utc_datetime(value, field=path)

    minimum = schema.get("minimum")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{path} is not finite")
    if minimum is not None and value < minimum:
        raise ValueError(f"{path} is below minimum")

    rejected_schema = schema.get("not")
    if rejected_schema is not None and _schema_matches(value, rejected_schema):
        raise ValueError(f"{path} matches rejected schema")


def _load_object(path: Path) -> dict[str, Any]:
    return read_json_bytes(read_stable_regular_file(path), label=str(path))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _validate_environment_allowlist(
    data: bytes, *, require_frontier_values: bool = False
) -> None:
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError("environment allowlist is not UTF-8") from error
    if "\r" in text:
        raise ValueError("environment allowlist must use LF line endings")
    if not text.endswith("\n"):
        raise ValueError("environment allowlist must end with a newline")
    records = text.splitlines()
    values = {}
    for record in records:
        key, separator, value = record.partition("=")
        if not separator or not key or not value:
            raise ValueError("environment allowlist contains a malformed record")
        if key in values:
            raise ValueError("environment allowlist contains a duplicate key")
        values[key] = value
    if list(values) != RUNTIME_ENVIRONMENT_ALLOWLIST_KEYS:
        raise ValueError("environment allowlist key set or order is not authorized")
    if not require_frontier_values:
        return
    profile = values["PIC_FRONTIER_PROFILE"]
    expected_by_profile = {
        "frontier_minimum_supported": {},
        "frontier_xnack1_experimental": {
            "HSA_XNACK": "1",
            "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED": "1",
        },
        "frontier_ofi_tuned_experimental": {
            "MPICH_OFI_NIC_POLICY": "GPU",
            "MPICH_GPU_IPC_CACHE_MAX_SIZE": "1000",
            "MPICH_MPIIO_HINTS": "*:romio_cb_write=disable",
            "MPICH_OFI_NUM_CQ_ENTRIES": "131072",
            "FI_MR_CACHE_MONITOR": "kdreg2",
            "FI_CXI_RX_MATCH_MODE": "software",
        },
    }
    if profile not in expected_by_profile:
        raise ValueError("environment allowlist selects an unsupported profile")
    expected = {
        "HSA_XNACK": "0",
        "MPICH_ENV_DISPLAY": "1",
        "MPICH_VERSION_DISPLAY": "1",
        "MPICH_GPU_SUPPORT_ENABLED": "1",
        "MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED": "<unset>",
        "MPICH_OFI_NIC_POLICY": "<unset>",
        "MPICH_GPU_IPC_CACHE_MAX_SIZE": "<unset>",
        "MPICH_MPIIO_HINTS": "<unset>",
        "MPICH_OFI_NUM_CQ_ENTRIES": "<unset>",
        "FI_MR_CACHE_MONITOR": "<unset>",
        "FI_CXI_RX_MATCH_MODE": "<unset>",
        "SLURM_EXPORT_ENV": "ALL",
    }
    expected.update(expected_by_profile[profile])
    for key, value in expected.items():
        if values[key] != value:
            raise ValueError(f"environment allowlist {key} is not authorized")
    if (
        values["OMP_NUM_THREADS"] != "<unset>"
        and not re.fullmatch(r"[1-9][0-9]*", values["OMP_NUM_THREADS"])
    ):
        raise ValueError("environment allowlist OMP_NUM_THREADS is malformed")
    rocm_path = PurePosixPath(values["ROCM_PATH"])
    if (
        not rocm_path.is_absolute()
        or rocm_path.as_posix() != values["ROCM_PATH"]
        or any(part in {"", ".", ".."} for part in rocm_path.parts)
    ):
        raise ValueError("environment allowlist ROCM_PATH is not canonical")
    if values["LOADEDMODULES"] != ":".join(PRODUCTION_RUNTIME_LOADED_MODULES):
        raise ValueError("environment allowlist LOADEDMODULES is not authorized")
    if values["_LMFILES_"] != ":".join(PRODUCTION_RUNTIME_MODULEFILES):
        raise ValueError("environment allowlist _LMFILES_ is not authorized")
    if values["MODULEPATH"] != PRODUCTION_RUNTIME_MODULEPATH:
        raise ValueError("environment allowlist MODULEPATH is not authorized")


def _source_bundle_sha256(
    source_archive_sha256: str,
    source_commit_sha256: str,
    submodules: list[dict[str, Any]],
) -> str:
    value = {
        "source_archive_sha256": source_archive_sha256,
        "source_commit_sha256": source_commit_sha256,
        "submodules": [
            {
                "path": record["path"],
                "archive_sha256": record["archive_sha256"],
                "commit_sha256": record["commit_sha256"],
                "git_commit": record["git_commit"],
                "git_tree": record["git_tree"],
            }
            for record in submodules
        ],
    }
    return source_bundle_sha256(
        str(value["source_archive_sha256"]),
        str(value["source_commit_sha256"]),
        value["submodules"],
    )


def _artifact_parts(artifact_root: Path, raw_path: str, label: str) -> tuple[str, ...]:
    if not isinstance(raw_path, str):
        raise ValueError(f"{label} path must be a string")
    pure = PurePosixPath(raw_path)
    if (
        not raw_path
        or raw_path != pure.as_posix()
        or any(part in {"", ".", ".."} for part in pure.parts)
    ):
        raise ValueError(f"{label} path is not canonical: {raw_path!r}")
    if pure.is_absolute():
        try:
            pure = PurePosixPath(Path(raw_path).relative_to(artifact_root).as_posix())
        except ValueError as error:
            raise ValueError(f"{label} escapes artifact root: {raw_path!r}") from error
    if not pure.parts:
        raise ValueError(f"{label} path is empty")
    return pure.parts


def _open_directory_component_at(directory_fd: int, name: str, *, label: str) -> int:
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(name, flags, dir_fd=directory_fd)
    except OSError as error:
        raise ValueError(f"{label} directory component cannot be opened: {name}") from error
    if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
        os.close(descriptor)
        raise ValueError(f"{label} component is not a directory: {name}")
    return descriptor


def _open_absolute_directory(path: Path, *, label: str) -> int:
    path = Path(os.path.abspath(path))
    if not path.is_absolute():
        raise ValueError(f"{label} must be absolute")
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open("/", flags)
    try:
        for part in path.parts[1:]:
            child_fd = _open_directory_component_at(descriptor, part, label=label)
            os.close(descriptor)
            descriptor = child_fd
        return descriptor
    except BaseException:
        os.close(descriptor)
        raise


def _artifact_bytes(
    artifact_root_fd: int, artifact_root: Path, raw_path: str, label: str
) -> bytes:
    parts = _artifact_parts(artifact_root, raw_path, label)
    directory_fd = os.dup(artifact_root_fd)
    try:
        for part in parts[:-1]:
            child_fd = _open_directory_component_at(directory_fd, part, label=label)
            os.close(directory_fd)
            directory_fd = child_fd
        flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(parts[-1], flags, dir_fd=directory_fd)
        except OSError as error:
            raise ValueError(f"{label} cannot be opened") from error
        try:
            if not stat.S_ISREG(os.fstat(descriptor).st_mode):
                raise ValueError(f"{label} is not a regular file")
            with os.fdopen(descriptor, "rb", closefd=False) as stream:
                return stream.read()
        finally:
            os.close(descriptor)
    finally:
        os.close(directory_fd)


def _bound_artifact_bytes(
    artifact_root_fd: int,
    artifact_root: Path,
    record: dict[str, Any],
    label: str,
) -> bytes:
    data = _artifact_bytes(artifact_root_fd, artifact_root, record["path"], label)
    if sha256_bytes(data) != record["sha256"]:
        raise ValueError(f"{label} checksum mismatch")
    return data


def _load_object_bytes(data: bytes, *, label: str) -> dict[str, Any]:
    return read_json_bytes(data, label=label)


def _project_submodules(records: list[dict[str, Any]]) -> list[dict[str, str]]:
    return [
        {
            "path": str(record["path"]),
            "archive_sha256": str(record["archive_sha256"]),
            "commit_sha256": str(record["commit_sha256"]),
            "git_commit": str(record["git_commit"]),
            "git_tree": str(record["git_tree"]),
        }
        for record in records
    ]


def _require_portable_candidate_layout(candidate: dict[str, Any]) -> None:
    """Require creator layout while allowing review-bundle artifact remapping."""
    source = candidate["source"]
    build = candidate["build"]
    archive = PurePosixPath(str(source["archive_path"]))
    candidate_dir = archive.parent
    fixed_paths = [
        source["archive_path"],
        source["commit_path"],
        build["profile_path"],
        build["profile_receipt_path"],
        build["executable_path"],
    ]
    fixed_paths.extend(
        path
        for record in source["submodules"]
        for path in (record["archive_path"], record["commit_path"])
    )
    if (
        any(str(path) != PurePosixPath(str(path)).as_posix() for path in fixed_paths)
        or any(not str(path).startswith("/") or str(path).startswith("//")
               for path in fixed_paths)
        or any(
            part in {".", ".."}
            for path in fixed_paths
            for part in PurePosixPath(str(path)).parts
        )
        or
        not archive.is_absolute()
        or candidate_dir.name != candidate["freeze_id"]
        or archive != candidate_dir / "source.tar"
        or PurePosixPath(str(source["commit_path"])) != candidate_dir / "source.commit"
        or PurePosixPath(str(build["profile_path"])) != candidate_dir / "build_profile.json"
        or PurePosixPath(str(build["profile_receipt_path"]))
        != candidate_dir / "profile_receipt.json"
        or PurePosixPath(str(build["executable_path"])) != candidate_dir / "athena"
    ):
        raise ValueError("clean-candidate manifest does not use the creator fixed layout")
    for index, record in enumerate(source["submodules"]):
        submodule_dir = candidate_dir / "submodules"
        if (
            PurePosixPath(str(record["archive_path"]))
            != submodule_dir / f"{index:04d}.tar"
            or PurePosixPath(str(record["commit_path"]))
            != submodule_dir / f"{index:04d}.commit"
        ):
            raise ValueError("clean-candidate submodule does not use the creator fixed layout")


def _live_candidate_authorization(candidate_sha256: str) -> dict[str, str]:
    """Bind review freezing to the exact clean candidate in the live mirrored policy."""
    promotion_path = active_promotion_path(AUTHORIZED_PIC_ROOT)
    promotion = read_json_bytes(
        read_stable_regular_file_below(
            promotion_path, AUTHORIZED_PIC_ROOT, require_read_only_mode=True
        ),
        label=str(promotion_path),
    )
    control_plane_version = str(promotion.get("control_plane_version", ""))
    policy, snapshot = require_storage_policy_unlock_snapshot(
        control_plane_version=control_plane_version,
    )
    science_freeze = policy["science_submission_freeze"]
    if science_freeze.get("status") != "authorized":
        raise ValueError("Live policy has not authorized a clean-candidate science freeze")
    manifest_path = require_canonical_path_below(
        Path(str(science_freeze["manifest_path"])),
        AUTHORIZED_PIC_ROOT / "clean_candidates",
    )
    manifest_bytes = read_stable_regular_file_below(
        manifest_path,
        AUTHORIZED_PIC_ROOT / "clean_candidates",
        require_read_only_mode=True,
    )
    if (
        science_freeze.get("manifest_sha256") != candidate_sha256
        or sha256_bytes(manifest_bytes) != candidate_sha256
    ):
        raise ValueError("Qualification candidate differs from the live authorized freeze")
    return {
        "control_plane_version": control_plane_version,
        "clean_candidate_manifest_path": str(manifest_path),
        "clean_candidate_manifest_sha256": candidate_sha256,
        "active_policy_sha256": snapshot["active_policy_sha256"],
        "active_promotion_sha256": snapshot["active_promotion_sha256"],
    }


def _require_frontier_ledger_binding(
    manifest: dict[str, Any],
    *,
    candidate_sha256: str,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    """Require a completed mirrored-ledger record for Frontier qualification."""
    resources = manifest["resources"]
    frontier_fields = {
        "submission_id",
        "reservation_id",
        "job_id",
        "pre_submit_manifest_path",
        "pre_submit_manifest_sha256",
        "run_artifact_dir",
        "node_hours",
    }
    if resources["platform"] != "Frontier":
        unexpected = frontier_fields & set(resources)
        if unexpected:
            raise ValueError(f"Host qualification resources claim Frontier fields: {unexpected}")
        return
    missing = frontier_fields - set(resources)
    if missing:
        raise ValueError(f"Frontier qualification resources are missing: {missing}")
    manifest_path = require_canonical_path_below(
        Path(str(resources["pre_submit_manifest_path"])),
        authorized_pic_root / "manifests",
    )
    manifest_bytes = read_stable_regular_file_below(
        manifest_path,
        authorized_pic_root / "manifests",
        require_read_only_mode=True,
    )
    manifest_sha256 = sha256_bytes(manifest_bytes)
    if manifest_sha256 != resources["pre_submit_manifest_sha256"]:
        raise ValueError("Frontier pre-submit manifest checksum mismatch")
    artifact_dir = require_canonical_path_below(
        Path(str(resources["run_artifact_dir"])),
        authorized_pic_root / "runs",
    )
    if Path(str(resources["artifact_root"])) != artifact_dir:
        raise ValueError("Frontier qualification artifact root differs from run artifact directory")
    ledger_jsonl = authorized_pic_root / "ledger" / "node_hours.jsonl"
    receipts_jsonl = authorized_pic_root / "ledger" / "mirror_receipts.jsonl"
    mirror_jsonl = authorized_project_home_root / "ledger" / "node_hours.jsonl"
    require_ledger_paths(
        ledger_jsonl,
        authorized_pic_root / "ledger" / "node_hours.csv",
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl):
        records = validate_mirrored_state(
            ledger_jsonl,
            receipts_jsonl,
            mirror_jsonl,
            ledger_root=authorized_pic_root / "ledger",
            receipts_root=authorized_pic_root / "ledger",
            mirror_root=authorized_project_home_root / "ledger",
        )
        require_explicit_genesis(records)
    reservation = latest_reservations(records).get(str(resources["reservation_id"]))
    if reservation is None:
        raise ValueError("Frontier qualification reservation is absent from the ledger")
    expected = {
        "event_type": "reconciliation",
        "reconciled": True,
        "state": "COMPLETED",
        "submission_scope": "registered_science",
        "submission_id": resources["submission_id"],
        "reservation_id": resources["reservation_id"],
        "job_id": resources["job_id"],
        "control_plane_version": control_plane_version,
        "clean_candidate_manifest_sha256": candidate_sha256,
        "manifest_path": str(manifest_path),
        "manifest_sha256": manifest_sha256,
        "artifact_dir": str(artifact_dir),
        "consumed_node_hours": resources["node_hours"],
    }
    if any(reservation.get(key) != value for key, value in expected.items()):
        raise ValueError("Frontier qualification does not match its reconciled ledger record")


def validate_qualification_manifest(
    manifest: dict[str, Any], *, verify_files: bool = True
) -> None:
    """Fail closed unless a manifest is ready for external scientific review."""
    validate_schema(manifest, _load_object(SCHEMA_PATH))
    resources = manifest["resources"]
    if resources["platform"] != "Frontier":
        _require_frontier_ledger_binding(
            manifest,
            candidate_sha256="",
            control_plane_version="",
        )

    evidence_class = manifest["evidence_class"]
    if evidence_class not in QUALIFYING_EVIDENCE_CLASSES:
        raise ValueError(
            f"evidence class is not qualifying evidence: {evidence_class}"
        )

    claim_registry = _load_object(CLAIMS_PATH)
    known_claim_ids = {
        claim["claim_id"] for claim in claim_registry["claims"]
    }
    unknown_claim_ids = set(manifest["claim_ids"]) - known_claim_ids
    if unknown_claim_ids:
        raise ValueError(f"manifest references unknown claims: {unknown_claim_ids}")

    if manifest["git"]["status"]:
        raise ValueError("qualification manifest requires a clean git status")
    submodules = manifest["git"]["submodules"]
    paths = [record["path"] for record in submodules]
    if paths != sorted(set(paths)):
        raise ValueError("qualification manifest submodules must be unique and sorted")
    expected_status = "clean_pinned_archived" if submodules else "absent"
    if manifest["git"]["submodule_status"] != expected_status:
        raise ValueError("qualification manifest submodule status does not match records")
    if manifest["git"]["source_bundle_sha256"] != _source_bundle_sha256(
        manifest["git"]["source_archive"]["sha256"],
        manifest["git"]["source_commit"]["sha256"],
        submodules,
    ):
        raise ValueError("qualification manifest source-bundle checksum mismatch")

    if not verify_files:
        return

    artifact_root = Path(manifest["resources"]["artifact_root"])
    if not artifact_root.is_absolute():
        raise ValueError("artifact root must be absolute")
    lexical_root = Path(os.path.abspath(artifact_root))
    if lexical_root != artifact_root or lexical_root.resolve() != lexical_root:
        raise ValueError("artifact root must use its canonical spelling")
    artifact_root = lexical_root
    if not artifact_root.is_dir():
        raise ValueError(f"artifact root is not a directory: {artifact_root}")

    artifact_root_fd = _open_absolute_directory(artifact_root, label="artifact root")
    try:
        source_archive = _bound_artifact_bytes(
            artifact_root_fd, artifact_root, manifest["git"]["source_archive"],
            "source archive",
        )
        source_commit = _bound_artifact_bytes(
            artifact_root_fd, artifact_root, manifest["git"]["source_commit"],
            "source commit object",
        )
        submodule_archives = []
        submodule_commits = []
        for index, submodule in enumerate(submodules):
            submodule_archives.append(
                _bound_artifact_bytes(
                    artifact_root_fd,
                    artifact_root,
                    {
                        "path": submodule["archive_path"],
                        "sha256": submodule["archive_sha256"],
                    },
                    f"submodule archive {index}",
                )
            )
            submodule_commits.append(
                _bound_artifact_bytes(
                    artifact_root_fd,
                    artifact_root,
                    {
                        "path": submodule["commit_path"],
                        "sha256": submodule["commit_sha256"],
                    },
                    f"submodule commit object {index}",
                )
            )
        candidate_bytes = _bound_artifact_bytes(
            artifact_root_fd, artifact_root,
            manifest["git"]["clean_candidate_manifest"], "clean-candidate manifest",
        )
        build_profile = _bound_artifact_bytes(
            artifact_root_fd, artifact_root,
            manifest["git"]["clean_candidate_build_profile"],
            "clean-candidate build profile",
        )
        build_profile_receipt = _bound_artifact_bytes(
            artifact_root_fd, artifact_root,
            manifest["git"]["clean_candidate_build_profile_receipt"],
            "clean-candidate build-profile receipt",
        )
        build_provenance = {
            label: _bound_artifact_bytes(
                artifact_root_fd,
                artifact_root,
                manifest["git"]["clean_candidate_build_provenance"][label],
                f"clean-candidate build provenance {label}",
            )
            for label in BUILD_PROVENANCE_FILENAMES
        }
        executable = _bound_artifact_bytes(
            artifact_root_fd, artifact_root, manifest["executable"], "executable"
        )
        _bound_artifact_bytes(
            artifact_root_fd,
            artifact_root,
            manifest["authorization"]["active_policy"],
            "active policy",
        )
        _bound_artifact_bytes(
            artifact_root_fd,
            artifact_root,
            manifest["authorization"]["active_promotion"],
            "active promotion",
        )
        for field in ("cmake_cache", "modules"):
            _bound_artifact_bytes(
                artifact_root_fd,
                artifact_root,
                manifest["executable"][field],
                f"executable provenance field {field}",
            )
        environment_allowlist = _bound_artifact_bytes(
            artifact_root_fd,
            artifact_root,
            manifest["executable"]["environment_allowlist"],
            "executable provenance field environment_allowlist",
        )
        for index, artifact in enumerate(manifest["artifacts"]):
            _bound_artifact_bytes(
                artifact_root_fd, artifact_root, artifact, f"artifact {index}"
            )
    finally:
        os.close(artifact_root_fd)

    executable_sha256 = sha256_bytes(executable)
    _validate_environment_allowlist(
        environment_allowlist,
        require_frontier_values=resources["platform"] == "Frontier",
    )
    candidate = _load_object_bytes(candidate_bytes, label="clean-candidate manifest")
    candidate_submodules = validate_clean_candidate_bundle(
        candidate,
        source_archive=source_archive,
        source_commit=source_commit,
        submodule_archives=submodule_archives,
        submodule_commits=submodule_commits,
        build_profile=build_profile,
        build_profile_receipt=build_profile_receipt,
        build_provenance=build_provenance,
        executable_sha256=executable_sha256,
        expected_control_plane_version=manifest["authorization"]["control_plane_version"],
    )
    _require_portable_candidate_layout(candidate)
    authorization = manifest["authorization"]
    expected_authorization = _live_candidate_authorization(
        manifest["git"]["clean_candidate_manifest"]["sha256"]
    )
    if {
        "control_plane_version": authorization["control_plane_version"],
        "clean_candidate_manifest_path": authorization["clean_candidate_manifest_path"],
        "clean_candidate_manifest_sha256": authorization["clean_candidate_manifest_sha256"],
        "active_policy_sha256": authorization["active_policy"]["sha256"],
        "active_promotion_sha256": authorization["active_promotion"]["sha256"],
    } != expected_authorization:
        raise ValueError("qualification manifest differs from the live authorization anchor")
    _require_frontier_ledger_binding(
        manifest,
        candidate_sha256=manifest["git"]["clean_candidate_manifest"]["sha256"],
        control_plane_version=authorization["control_plane_version"],
    )
    candidate_source = candidate["source"]
    if not isinstance(candidate_source, dict):
        raise ValueError("clean-candidate source attestation is missing")
    if (
        manifest["git"]["commit"] != candidate_source["git_commit"]
        or manifest["git"]["tree"] != candidate_source["git_tree"]
        or manifest["git"]["source_archive"]["sha256"]
        != candidate_source["archive_sha256"]
        or manifest["git"]["source_bundle_sha256"]
        != candidate_source["source_bundle_sha256"]
        or manifest["git"]["submodule_status"]
        != candidate_source["submodule_status"]
        or _project_submodules(submodules) != candidate_submodules
    ):
        raise ValueError("qualification manifest source projection differs from clean candidate")

def freeze_qualification_manifest(input_path: Path, output_path: Path) -> None:
    """Validate a prepared manifest and write a canonical immutable candidate."""
    manifest = _load_object(input_path)
    validate_qualification_manifest(manifest)
    data = (
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    output_path = Path(os.path.abspath(output_path))
    if not output_path.name or output_path.name in {".", ".."}:
        raise ValueError("Qualification output path must name one file")
    temporary_name = f".{output_path.name}.tmp-{uuid.uuid4()}"
    parent_fd = _open_absolute_directory(
        output_path.parent, label="qualification output parent"
    )
    descriptor = os.open(
        temporary_name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL,
        0o600,
        dir_fd=parent_fd,
    )
    temporary_exists = True
    published = False
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
        os.close(descriptor)
        descriptor = -1
        os.link(
            temporary_name,
            output_path.name,
            src_dir_fd=parent_fd,
            dst_dir_fd=parent_fd,
            follow_symlinks=False,
        )
        published = True
        try:
            os.fsync(parent_fd)
        except BaseException as error:
            try:
                os.unlink(output_path.name, dir_fd=parent_fd)
                published = False
                os.fsync(parent_fd)
            except BaseException as rollback_error:
                raise RuntimeError(
                    f"Failed to durably roll back qualification output: {output_path}"
                ) from rollback_error
            raise error
        os.unlink(temporary_name, dir_fd=parent_fd)
        temporary_exists = False
        try:
            os.fsync(parent_fd)
        except BaseException as error:
            try:
                os.unlink(output_path.name, dir_fd=parent_fd)
                published = False
                os.fsync(parent_fd)
            except BaseException as rollback_error:
                raise RuntimeError(
                    f"Failed to durably roll back qualification output: {output_path}"
                ) from rollback_error
            raise error
    finally:
        try:
            if descriptor >= 0:
                os.close(descriptor)
            if temporary_exists:
                try:
                    os.unlink(temporary_name, dir_fd=parent_fd)
                except FileNotFoundError:
                    pass
                try:
                    os.fsync(parent_fd)
                except BaseException as rollback_error:
                    raise RuntimeError(
                        f"Failed to durably remove qualification temporary: {output_path}"
                    ) from rollback_error
        finally:
            os.close(parent_fd)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--validate",
        type=Path,
        metavar="MANIFEST",
        help="validate an existing qualification manifest",
    )
    group.add_argument(
        "--input",
        type=Path,
        metavar="MANIFEST",
        help="validate and canonically freeze a prepared manifest",
    )
    parser.add_argument(
        "--output",
        type=Path,
        metavar="MANIFEST",
        help="exclusive output path required with --input",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    if args.validate is not None:
        if args.output is not None:
            raise ValueError("--output is only valid with --input")
        validate_qualification_manifest(_load_object(args.validate))
    else:
        if args.output is None:
            raise ValueError("--output is required with --input")
        freeze_qualification_manifest(args.input, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
