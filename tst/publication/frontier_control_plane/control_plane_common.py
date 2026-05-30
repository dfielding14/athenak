#!/opt/cray/pe/python/3.11.7/bin/python3
"""Shared helpers for immutable Frontier PIC submission snapshots."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tarfile
from typing import Iterable
import uuid


AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_PROJECT_HOME_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
AUTHORIZED_ACCOUNT = "AST207"
AUTHORIZED_PARTITION = "batch"
AUTHORIZED_NODE_HOUR_CAP = 10000.0
AUTHORIZED_LEDGER_MIRROR_TRANSPORT = "filesystem_copy"
AUTHORIZED_PROJECT_HOME_USAGE = [
    "small_append_only_ledger_and_control_plane_mirror",
]
AUTHORIZED_PROJECT_HOME_RETENTION_ROLE = "operational_ledger_mirror_only"
AUTHORIZED_ORION_BULK_EVIDENCE_USAGE = [
    "simulation_outputs",
    "immutable_bulk_signoff_bundles",
    "private_reference_artifact_staging",
    "restore_drill_evidence",
]
AUTHORIZED_ORION_RETENTION_ROLE = (
    "user_selected_sole_bulk_evidence_root_with_documented_durability_risk"
)
AUTHORIZED_LONG_TERM_STORAGE_STATUS = (
    "user_selected_orion_only_with_documented_durability_risk"
)
AUTHORIZED_LONG_TERM_STORAGE_RISK = (
    "Orion-only retention is user-directed and does not provide an institutional "
    "or approved off-site durable archive."
)
AUTHORIZED_LONG_TERM_STORAGE_BLOCKS = [
    "terminal_durable_retention_signoff_pending_external_review",
]
TRUSTED_GIT = "/usr/bin/git"
TRUSTED_PYTHON = "/opt/cray/pe/python/3.11.7/bin/python3"
TRUSTED_SQUEUE = "/usr/bin/squeue"
TRUSTED_SCONTROL = "/usr/bin/scontrol"
TRUSTED_SACCT = "/usr/bin/sacct"
TRUSTED_SBATCH = "/usr/bin/sbatch"
TRUSTED_SCANCEL = "/usr/bin/scancel"
REGISTERED_SCIENCE_SCOPE = "registered_science"
FRONTIER_ADMISSION_SMOKE_SCOPE = "frontier_admission_smoke"
SUBMISSION_SCOPES = {
    REGISTERED_SCIENCE_SCOPE,
    FRONTIER_ADMISSION_SMOKE_SCOPE,
}
PENDING_CLEAN_CANDIDATE_FREEZE = "pending_clean_candidate_freeze"
AUTHORIZED_CLEAN_CANDIDATE_FREEZE = "authorized"
AUTHORIZED_ADMISSION_SMOKE_STATUS = "authorized_f0_parser_contract_only"
PENDING_ADMISSION_SMOKE_STATUS = "pending_exact_executable_binding"
TRUSTED_LAUNCH_EXECUTOR = "trusted_trampoline_athena_argv_v1"
CANONICAL_POLICY_RELATIVE = Path("policy/storage_policy.json")
ACTIVE_PROMOTION_RELATIVE = Path("policy/active_promotion.json")
SITE_POLICY_MAX_AGE_SECONDS = 24 * 60 * 60
CONTROL_PLANE_FILES = [
    "clean_candidate.schema.json",
    "control_plane.schema.json",
    "control_plane_common.py",
    "create_clean_candidate_freeze.py",
    "create_pre_submit_manifest.py",
    "frontier_pic_environment.sh",
    "initialize_frontier_ledger.py",
    "launch_trampoline.py",
    "launch_with_frontier_profile.sh",
    "ledger.py",
    "promote_active_policy.py",
    "reconcile_frontier_job.py",
    "submit_frontier_job.sh",
    "validate_and_reserve_frontier_job.py",
    "verify_compute_node_snapshot.py",
]
PLACEHOLDER_PATTERN = re.compile(rb"REPLACE_[A-Z0-9_]+")
SENSITIVE_PATTERN = re.compile(
    rb"(TOKEN|PASSWORD|SECRET|PRIVATE[_-]?KEY|BEGIN [A-Z ]*PRIVATE KEY)",
    re.IGNORECASE,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _git_object_sha1(kind: str, data: bytes) -> bytes:
    header = f"{kind} {len(data)}\0".encode("ascii")
    return hashlib.sha1(header + data).digest()


def direct_submodule_gitlinks(
    records: Iterable[dict[str, object]], *, parent_path: str | None = None
) -> dict[str, str]:
    """Return the immediate Git links represented in one archived repository."""
    parsed: list[tuple[PurePosixPath, str]] = []
    seen: set[str] = set()
    for record in records:
        value = str(record["path"])
        path = PurePosixPath(value)
        commit = str(record["git_commit"])
        if (
            not value
            or path.is_absolute()
            or not path.parts
            or path.parts
            != tuple(part for part in path.parts if part not in {"", ".", ".."})
            or not re.fullmatch(r"[0-9a-f]{40}", commit)
        ):
            raise ValueError(f"Unsafe Git-link attestation: {value!r}")
        if value in seen:
            raise ValueError(f"Duplicate Git-link attestation: {value!r}")
        seen.add(value)
        parsed.append((path, commit))
    parent = PurePosixPath(parent_path) if parent_path is not None else None
    result: dict[str, str] = {}
    for path, commit in parsed:
        ancestors = [
            candidate
            for candidate, _ in parsed
            if len(candidate.parts) < len(path.parts)
            and path.parts[: len(candidate.parts)] == candidate.parts
        ]
        nearest = max(ancestors, key=lambda candidate: len(candidate.parts), default=None)
        if nearest != parent:
            continue
        relative = path.relative_to(parent).as_posix() if parent else path.as_posix()
        result[relative] = commit
    return result


def source_bundle_sha256(
    source_archive_sha256: str, submodules: Iterable[dict[str, object]]
) -> str:
    """Bind one parent archive and its ordered recursive submodule archives."""
    value = {
        "source_archive_sha256": source_archive_sha256,
        "submodules": list(submodules),
    }
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def git_tree_sha1_from_archive(
    archive: Path,
    *,
    gitlinks: dict[str, str] | None = None,
    reject_symlinks: bool = False,
) -> str:
    """Reconstruct the Git tree object ID represented by a git-archive tarball."""
    root: dict[str, object] = {}
    with tarfile.open(archive, "r:*") as stream:
        for member in stream.getmembers():
            path = PurePosixPath(member.name)
            if path.is_absolute() or not path.parts or any(
                part in {"", ".", ".."} for part in path.parts
            ):
                raise ValueError(f"Unsafe path in source archive: {member.name!r}")
            node = root
            for part in path.parts[:-1]:
                existing = node.setdefault(part, {})
                if not isinstance(existing, dict):
                    raise ValueError(f"Archive path collides with a file: {member.name!r}")
                node = existing
            name = path.parts[-1]
            if member.isdir():
                existing = node.setdefault(name, {})
                if not isinstance(existing, dict):
                    raise ValueError(f"Archive directory collides with a file: {member.name!r}")
                continue
            if name in node:
                raise ValueError(f"Duplicate path in source archive: {member.name!r}")
            if member.isfile():
                extracted = stream.extractfile(member)
                if extracted is None:
                    raise ValueError(f"Cannot read source-archive member: {member.name!r}")
                data = extracted.read()
                mode = "100755" if member.mode & 0o111 else "100644"
            elif member.issym():
                if reject_symlinks:
                    raise ValueError(f"Symlink is not allowed in source archive: {member.name!r}")
                data = member.linkname.encode("utf-8", "surrogateescape")
                mode = "120000"
            else:
                raise ValueError(f"Unsupported source-archive member: {member.name!r}")
            node[name] = (mode, _git_object_sha1("blob", data))

    for value, commit in (gitlinks or {}).items():
        path = PurePosixPath(value)
        if (
            not value
            or path.is_absolute()
            or not path.parts
            or path.parts
            != tuple(part for part in path.parts if part not in {"", ".", ".."})
            or not re.fullmatch(r"[0-9a-f]{40}", commit)
        ):
            raise ValueError(f"Unsafe Git-link attestation: {value!r}")
        node = root
        for part in path.parts[:-1]:
            existing = node.setdefault(part, {})
            if not isinstance(existing, dict):
                raise ValueError(f"Git-link path collides with a file: {value!r}")
            node = existing
        name = path.parts[-1]
        if name not in node or node[name] != {}:
            raise ValueError(f"Git-link path is not an empty archive directory: {value!r}")
        node[name] = ("160000", bytes.fromhex(commit))

    def tree_sha1(node: dict[str, object]) -> bytes:
        entries = []
        for name, value in node.items():
            encoded_name = name.encode("utf-8", "surrogateescape")
            if isinstance(value, dict):
                mode = "40000"
                digest = tree_sha1(value)
                sort_key = encoded_name + b"/"
            else:
                mode, digest = value
                sort_key = encoded_name
            entry = mode.encode("ascii") + b" " + encoded_name + b"\0" + digest
            entries.append((sort_key, entry))
        payload = b"".join(entry for _, entry in sorted(entries))
        return _git_object_sha1("tree", payload)

    return tree_sha1(root).hex()


def read_json(path: Path) -> dict[str, object]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected a JSON object in {path}")
    return value


def utc_datetime(value: object, *, field: str) -> datetime:
    text = str(value)
    if not text.endswith("Z"):
        raise ValueError(f"{field} must be an RFC-3339 UTC timestamp ending in Z")
    try:
        result = datetime.fromisoformat(text[:-1] + "+00:00")
    except ValueError as error:
        raise ValueError(f"Invalid {field}: {text}") from error
    if result.tzinfo != timezone.utc:
        raise ValueError(f"{field} must use UTC")
    return result


def write_json_exclusive(path: Path, value: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def atomic_write_bytes(path: Path, data: bytes, *, mode: int = 0o444) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{uuid.uuid4()}"
    try:
        with temporary.open("xb") as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        temporary.chmod(mode)
        temporary.replace(path)
    finally:
        if temporary.exists():
            temporary.unlink()


def atomic_write_json(path: Path, value: dict[str, object], *, mode: int = 0o444) -> None:
    data = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    atomic_write_bytes(path, data, mode=mode)


def make_tree_read_only(root: Path, *, executable_names: set[str] | None = None) -> None:
    executables = executable_names or set()
    for path in sorted(root.rglob("*"), reverse=True):
        if path.is_dir():
            path.chmod(0o555)
        else:
            path.chmod(0o555 if path.name in executables else 0o444)
    root.chmod(0o555)


def remove_tree(root: Path) -> None:
    """Remove a generated tree even if it was already staged read-only."""
    if not root.exists():
        return
    for path in root.rglob("*"):
        if path.is_dir():
            path.chmod(0o700)
    root.chmod(0o700)
    shutil.rmtree(root)


def require_read_only(path: Path) -> None:
    if path.stat().st_mode & 0o222:
        raise ValueError(f"Artifact is not read-only: {path}")


def require_not_symlink(path: Path) -> None:
    if path.is_symlink():
        raise ValueError(f"Path must not be a symlink alias: {path}")


def canonical_policy_path(authorized_pic_root: Path = AUTHORIZED_PIC_ROOT) -> Path:
    return authorized_pic_root.resolve() / CANONICAL_POLICY_RELATIVE


def active_promotion_path(root: Path) -> Path:
    return root.resolve() / ACTIVE_PROMOTION_RELATIVE


def require_below(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    try:
        resolved.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError(f"Path is outside authorized PIC root: {resolved}") from error
    return resolved


def require_no_symlink_components_below(path: Path, root: Path) -> Path:
    """Reject pre-existing symlink traversal below one trusted lexical root."""
    lexical_root = Path(os.path.abspath(root))
    lexical_path = Path(os.path.abspath(path))
    try:
        relative = lexical_path.relative_to(lexical_root)
    except ValueError as error:
        raise ValueError(f"Path is outside authorized lexical root: {lexical_path}") from error
    current = lexical_root
    for part in relative.parts:
        current /= part
        if current.is_symlink():
            raise ValueError(f"Path traverses a symlink below authorized root: {current}")
    return lexical_path


def require_canonical_path_below(path: Path, root: Path) -> Path:
    """Require one lexical path under root with no symlink aliases."""
    lexical = require_no_symlink_components_below(path, root)
    resolved = require_below(lexical, root)
    if lexical != resolved:
        raise ValueError(f"Path must use its canonical spelling: {lexical}; expected {resolved}")
    return lexical


def _relative_artifact_path(value: object, *, field: str) -> str:
    text = str(value)
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or not path.parts
        or path.parts != tuple(part for part in path.parts if part not in {"", ".", ".."})
    ):
        raise ValueError(f"{field} must be a non-empty relative artifact path")
    return text


def validate_launch_contract(value: object) -> dict[str, object]:
    """Validate the closed structured argv contract executed by the trampoline."""
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "executor",
        "pre_actions",
        "actions",
        "post_actions",
    }:
        raise ValueError(
            "Launch contract must contain only schema_version, executor, "
            "pre_actions, actions and post_actions"
        )
    if value.get("schema_version") != 1:
        raise ValueError("Unsupported launch-contract schema")
    if value.get("executor") != TRUSTED_LAUNCH_EXECUTOR:
        raise ValueError("Launch contract does not select the trusted Athena executor")
    actions = value.get("actions")
    if not isinstance(actions, list) or not 1 <= len(actions) <= 16:
        raise ValueError("Launch contract requires between one and sixteen Athena actions")
    identifiers = set()
    for phase in ["pre_actions", "post_actions"]:
        bounded_actions = value.get(phase)
        if not isinstance(bounded_actions, list) or len(bounded_actions) > 16:
            raise ValueError(f"Launch contract {phase} must contain at most sixteen actions")
        for action in bounded_actions:
            if not isinstance(action, dict):
                raise ValueError(f"Launch-contract {phase} action is malformed")
            identifier = str(action.get("action_id", ""))
            if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", identifier):
                raise ValueError(f"Launch-contract {phase} action ID is malformed")
            if identifier in identifiers:
                raise ValueError(f"Duplicate launch action ID: {identifier}")
            identifiers.add(identifier)
            kind = action.get("kind")
            if kind == "snapshot_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "snapshot_role",
                    "output_artifact",
                }:
                    raise ValueError(f"Launch-contract {phase} snapshot action is malformed")
                role = str(action.get("snapshot_role", ""))
                if role not in {
                    "job-script",
                    "executable",
                    "input-deck",
                    "environment-profile",
                    "timeout-margin",
                    "queue-snapshot",
                } and not re.fullmatch(r"analysis-script-[0-9]{3}", role):
                    raise ValueError(f"Launch-contract {phase} snapshot role is malformed")
                _relative_artifact_path(
                    action.get("output_artifact"), field="output_artifact"
                )
            elif kind == "artifact_sha256":
                if set(action) != {
                    "action_id",
                    "kind",
                    "artifact",
                    "output_artifact",
                }:
                    raise ValueError(f"Launch-contract {phase} artifact action is malformed")
                artifact = _relative_artifact_path(action.get("artifact"), field="artifact")
                output = _relative_artifact_path(
                    action.get("output_artifact"), field="output_artifact"
                )
                if artifact == output:
                    raise ValueError("Artifact checksum output must differ from its input")
            elif kind == "artifact_nonempty":
                if set(action) != {"action_id", "kind", "artifact"}:
                    raise ValueError(f"Launch-contract {phase} assertion action is malformed")
                _relative_artifact_path(action.get("artifact"), field="artifact")
            else:
                raise ValueError(f"Launch-contract {phase} accepts only bounded built-in actions")
    for action in actions:
        if not isinstance(action, dict) or set(action) != {
            "action_id",
            "kind",
            "resources",
            "arguments",
            "stdout_artifact",
            "stderr_artifact",
        }:
            raise ValueError("Launch action has unexpected or missing fields")
        identifier = str(action.get("action_id", ""))
        if not re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", identifier):
            raise ValueError("Launch action ID is malformed")
        if identifier in identifiers:
            raise ValueError(f"Duplicate launch action ID: {identifier}")
        identifiers.add(identifier)
        if action.get("kind") != "athena":
            raise ValueError("Trusted trampoline accepts only Athena actions")
        resources = action.get("resources")
        if not isinstance(resources, dict) or set(resources) != {
            "nodes",
            "tasks",
            "cpus_per_task",
            "gpus_per_task",
            "gpu_bind",
        }:
            raise ValueError("Launch-action resources are malformed")
        for field in ["nodes", "tasks", "cpus_per_task", "gpus_per_task"]:
            number = resources.get(field)
            if not isinstance(number, int) or isinstance(number, bool) or number <= 0:
                raise ValueError(f"Launch-action resources.{field} must be a positive integer")
        if resources.get("gpu_bind") != "closest":
            raise ValueError("Launch-action resources.gpu_bind must be closest")
        arguments = action.get("arguments")
        if not isinstance(arguments, list):
            raise ValueError("Launch-action arguments must be an array")
        for argument in arguments:
            if not isinstance(argument, dict) or len(argument) != 1:
                raise ValueError("Launch-action argument must be one structured token")
            if "literal" in argument:
                literal = argument["literal"]
                if (
                    not isinstance(literal, str)
                    or not literal
                    or "\0" in literal
                    or len(literal) > 4096
                ):
                    raise ValueError("Launch-action literal is malformed")
            elif argument.get("snapshot_role") == "input-deck":
                if set(argument) != {"snapshot_role"}:
                    raise ValueError("Launch-action snapshot token is malformed")
            elif "artifact_directory" in argument:
                if set(argument) != {"artifact_directory"}:
                    raise ValueError("Launch-action artifact token is malformed")
                _relative_artifact_path(
                    argument["artifact_directory"], field="artifact_directory"
                )
            else:
                raise ValueError("Launch-action argument type is not authorized")
        stdout = _relative_artifact_path(
            action.get("stdout_artifact"), field="stdout_artifact"
        )
        stderr = _relative_artifact_path(
            action.get("stderr_artifact"), field="stderr_artifact"
        )
        if stdout == stderr:
            raise ValueError("Launch-action stdout and stderr artifacts must differ")
    return value


def scrub_file(path: Path) -> None:
    data = path.read_bytes()
    if PLACEHOLDER_PATTERN.search(data):
        raise ValueError(f"Unresolved REPLACE_* placeholder in {path}")
    if SENSITIVE_PATTERN.search(data):
        raise ValueError(f"Potential sensitive string in submission artifact: {path}")


def snapshot_file(
    source: Path,
    destination: Path,
    *,
    role: str,
    destination_root: Path,
    scrub: bool = True,
) -> dict[str, str]:
    if not source.is_file():
        raise FileNotFoundError(f"Missing snapshot source for {role}: {source}")
    if scrub:
        scrub_file(source)
    try:
        destination.resolve().relative_to(destination_root.resolve())
    except ValueError as error:
        raise ValueError(
            f"Snapshot destination for {role} escapes staging root: {destination}"
        ) from error
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    if scrub:
        scrub_file(destination)
    return {
        "role": role,
        "path": str(destination),
        "sha256": sha256(destination),
        "source_path": str(source.resolve()),
        "source_sha256": sha256(source),
    }


def verify_snapshot_files(manifest: dict[str, object]) -> None:
    snapshot_files = manifest.get("snapshot_files")
    if not isinstance(snapshot_files, list) or not snapshot_files:
        raise ValueError("Manifest has no snapshot_files")
    for record in snapshot_files:
        if not isinstance(record, dict):
            raise ValueError("Malformed snapshot file record")
        path = Path(str(record["path"]))
        if not path.is_file():
            raise FileNotFoundError(f"Missing snapshotted dependency: {path}")
        if sha256(path) != record.get("sha256"):
            raise ValueError(f"Snapshot checksum mismatch: {path}")


def record_for_role(
    manifest: dict[str, object], role: str
) -> dict[str, object]:
    snapshot_files = manifest.get("snapshot_files", [])
    matches = [
        record for record in snapshot_files
        if isinstance(record, dict) and record.get("role") == role
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one snapshot role={role}")
    return matches[0]


def checksum_records(paths: Iterable[Path]) -> list[dict[str, str]]:
    return [{"path": str(path.resolve()), "sha256": sha256(path)} for path in paths]


def inventory_digest(records: list[dict[str, str]]) -> str:
    payload = json.dumps(records, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def verify_installed_control_plane(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> dict[str, object]:
    expected_parent = authorized_pic_root.resolve() / "control_plane"
    resolved = require_canonical_path_below(
        control_plane_dir, authorized_pic_root.resolve()
    )
    if resolved.parent != expected_parent:
        raise ValueError(
            f"Control plane is not installed under authorized root: {resolved}"
        )
    inventory_path = require_canonical_path_below(resolved / "inventory.json", resolved)
    if not inventory_path.is_file():
        raise ValueError(f"Missing installed control-plane inventory: {inventory_path}")
    inventory = read_json(inventory_path)
    if inventory.get("schema_version") != 1:
        raise ValueError("Unsupported control-plane inventory schema")
    raw_records = inventory.get("files")
    if not isinstance(raw_records, list):
        raise ValueError("Malformed control-plane inventory files")
    records: list[dict[str, str]] = []
    for raw in raw_records:
        if not isinstance(raw, dict):
            raise ValueError("Malformed control-plane inventory record")
        record = {"path": str(raw.get("path", "")), "sha256": str(raw.get("sha256", ""))}
        records.append(record)
    if [record["path"] for record in records] != CONTROL_PLANE_FILES:
        raise ValueError("Control-plane inventory file list differs from required list")
    version = inventory_digest(records)
    if inventory.get("version") != version or resolved.name != version:
        raise ValueError("Control-plane inventory digest mismatch")
    for record in records:
        path = require_canonical_path_below(resolved / record["path"], resolved)
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"Installed control-plane checksum mismatch: {path}")
        require_read_only(path)
    require_read_only(inventory_path)
    require_read_only(resolved)
    return inventory


def validate_storage_policy(
    policy: dict[str, object],
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
) -> dict[str, object]:
    frontier = policy.get("frontier")
    storage = policy.get("olcf_side_storage")
    long_term = policy.get("long_term_storage")
    science_freeze = policy.get("science_submission_freeze")
    admission_smoke = policy.get("frontier_admission_smoke")
    if (
        not isinstance(frontier, dict)
        or not isinstance(storage, dict)
        or not isinstance(long_term, dict)
        or not isinstance(science_freeze, dict)
        or not isinstance(admission_smoke, dict)
    ):
        raise ValueError(
            "Storage policy is missing Frontier, OLCF-side storage, science-freeze, "
            "or admission-smoke data"
        )
    expected_frontier = {
        "account": authorized_account,
        "partition": AUTHORIZED_PARTITION,
        "simulation_root": str(authorized_pic_root.resolve()),
        "maximum_node_hours": AUTHORIZED_NODE_HOUR_CAP,
        "serial_pic_submissions": True,
    }
    for key, expected in expected_frontier.items():
        if frontier.get(key) != expected:
            raise ValueError(f"Storage policy frontier.{key} must be {expected!r}")
    if Path(str(storage.get("project_home_mirror_root", ""))).resolve() != (
        authorized_project_home_root.resolve()
    ):
        raise ValueError("Storage policy Project Home mirror root is not authorized")
    if storage.get("project_home_usage") != AUTHORIZED_PROJECT_HOME_USAGE:
        raise ValueError("Storage policy Project Home usage is not authorized")
    if (
        storage.get("project_home_retention_role")
        != AUTHORIZED_PROJECT_HOME_RETENTION_ROLE
    ):
        raise ValueError("Storage policy Project Home retention role is not authorized")
    if ledger_mirror_transport != AUTHORIZED_LEDGER_MIRROR_TRANSPORT:
        raise ValueError("Only filesystem_copy Project Home ledger mirroring is authorized")
    if storage.get("project_home_ledger_mirror_transport") != ledger_mirror_transport:
        raise ValueError("Storage policy does not authorize the ledger mirror transport")
    if Path(str(storage.get("orion_bulk_evidence_root", ""))).resolve() != (
        authorized_pic_root.resolve()
    ):
        raise ValueError("Storage policy Orion bulk-evidence root is not authorized")
    if storage.get("orion_bulk_evidence_usage") != AUTHORIZED_ORION_BULK_EVIDENCE_USAGE:
        raise ValueError("Storage policy Orion bulk-evidence usage is not authorized")
    if storage.get("orion_retention_role") != AUTHORIZED_ORION_RETENTION_ROLE:
        raise ValueError("Storage policy Orion retention role is not authorized")
    if long_term.get("status") != AUTHORIZED_LONG_TERM_STORAGE_STATUS:
        raise ValueError("Storage policy must preserve the Orion-only durability risk")
    if Path(str(long_term.get("selected_destination", ""))).resolve() != (
        authorized_pic_root.resolve()
    ):
        raise ValueError("Storage policy long-term destination is not authorized")
    if long_term.get("risk") != AUTHORIZED_LONG_TERM_STORAGE_RISK:
        raise ValueError("Storage policy long-term risk statement is not authorized")
    if long_term.get("blocks") != AUTHORIZED_LONG_TERM_STORAGE_BLOCKS:
        raise ValueError("Storage policy long-term terminal block is not authorized")
    for key in [
        "orion_simulation_root_preflight",
        "project_home_preflight",
    ]:
        record = storage.get(key)
        if not isinstance(record, dict) or record.get("status") != "passed":
            raise ValueError(f"Storage policy {key} has not passed")
    if storage.get("ledger_genesis_allowed") is not True:
        raise ValueError("Storage policy does not allow Frontier PIC ledger genesis")
    if storage.get("installed_control_plane_version") != control_plane_version:
        raise ValueError("Storage policy does not authorize this control-plane version")
    freeze_status = science_freeze.get("status")
    if freeze_status == PENDING_CLEAN_CANDIDATE_FREEZE:
        if set(science_freeze) != {"status"}:
            raise ValueError("Pending science freeze must not retain candidate fields")
    elif freeze_status == AUTHORIZED_CLEAN_CANDIDATE_FREEZE:
        if set(science_freeze) != {"status", "manifest_path", "manifest_sha256"}:
            raise ValueError("Authorized science freeze must identify one exact manifest")
        manifest_path = Path(str(science_freeze["manifest_path"]))
        require_canonical_path_below(
            manifest_path, authorized_pic_root / "clean_candidates"
        )
        if manifest_path.name != "clean_candidate_manifest.json":
            raise ValueError("Authorized science-freeze manifest has an invalid path")
        if not re.fullmatch(r"[0-9a-f]{64}", str(science_freeze["manifest_sha256"])):
            raise ValueError("Authorized science-freeze manifest digest is malformed")
    else:
        raise ValueError("Storage policy does not declare a recognized science freeze state")
    if admission_smoke.get("status") == PENDING_ADMISSION_SMOKE_STATUS:
        if set(admission_smoke) != {"status"}:
            raise ValueError("Pending admission smoke must not retain executable fields")
        return policy
    if set(admission_smoke) != {
        "status",
        "campaign",
        "test_id",
        "evidence_class",
        "physical_mode",
        "selected_qos",
        "registered_short_nonproduction",
        "maximum_nodes",
        "maximum_walltime_seconds",
        "job_script_sha256",
        "input_deck_sha256",
        "environment_profile_sha256",
        "analysis_script_sha256",
        "executable_sha256",
    }:
        raise ValueError("Storage policy admission-smoke authorization is malformed")
    expected_admission_smoke = {
        "status": AUTHORIZED_ADMISSION_SMOKE_STATUS,
        "campaign": "f0_hipmpi_smoke",
        "test_id": "pic_parser_contract_guards",
        "evidence_class": "frontier_f0_admission_smoke_candidate",
        "physical_mode": "extended_mhd_pic_parser_contract",
        "selected_qos": "debug",
        "registered_short_nonproduction": True,
        "maximum_nodes": 1,
        "maximum_walltime_seconds": 15 * 60,
    }
    for key, expected in expected_admission_smoke.items():
        if admission_smoke.get(key) != expected:
            raise ValueError(f"Storage policy frontier_admission_smoke.{key} is invalid")
    for key in [
        "job_script_sha256",
        "input_deck_sha256",
        "environment_profile_sha256",
        "executable_sha256",
    ]:
        if not re.fullmatch(r"[0-9a-f]{64}", str(admission_smoke.get(key, ""))):
            raise ValueError(f"Storage policy frontier_admission_smoke.{key} is malformed")
    analysis_sha256 = admission_smoke.get("analysis_script_sha256")
    if (
        not isinstance(analysis_sha256, list)
        or len(analysis_sha256) != 1
        or not re.fullmatch(r"[0-9a-f]{64}", str(analysis_sha256[0]))
    ):
        raise ValueError(
            "Storage policy frontier_admission_smoke.analysis_script_sha256 is malformed"
        )
    return policy


def require_storage_policy_unlock(
    *,
    control_plane_version: str,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    ledger_mirror_transport: str = AUTHORIZED_LEDGER_MIRROR_TRANSPORT,
) -> dict[str, object]:
    policy_path = canonical_policy_path(authorized_pic_root)
    mirror_policy_path = canonical_policy_path(authorized_project_home_root)
    promotion_path = active_promotion_path(authorized_pic_root)
    mirror_promotion_path = active_promotion_path(authorized_project_home_root)
    for path, root in [
        (policy_path, authorized_pic_root),
        (mirror_policy_path, authorized_project_home_root),
        (promotion_path, authorized_pic_root),
        (mirror_promotion_path, authorized_project_home_root),
    ]:
        require_canonical_path_below(path, root)
        if not path.is_file():
            raise FileNotFoundError(f"Missing active policy artifact: {path}")
        require_read_only(path)
    promotion = read_json(promotion_path)
    mirror_promotion = read_json(mirror_promotion_path)
    if promotion != mirror_promotion:
        raise ValueError("Orion and Project Home active-policy promotion records differ")
    expected = {
        "schema_version": 1,
        "control_plane_version": control_plane_version,
        "policy_path": str(policy_path),
        "project_home_policy_path": str(mirror_policy_path),
        "policy_sha256": sha256(policy_path),
    }
    if promotion != expected:
        raise ValueError("Active-policy promotion record is not anchored to this control plane")
    if sha256(mirror_policy_path) != promotion["policy_sha256"]:
        raise ValueError("Project Home active-policy mirror checksum differs")
    return validate_storage_policy(
        read_json(policy_path),
        control_plane_version=control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        authorized_account=authorized_account,
        ledger_mirror_transport=ledger_mirror_transport,
    )


def require_ledger_paths(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    *,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    expected = [
        (
            ledger_jsonl,
            authorized_pic_root.resolve() / "ledger" / "node_hours.jsonl",
            authorized_pic_root.resolve(),
        ),
        (
            ledger_csv,
            authorized_pic_root.resolve() / "ledger" / "node_hours.csv",
            authorized_pic_root.resolve(),
        ),
        (
            receipts_jsonl,
            authorized_pic_root.resolve() / "ledger" / "mirror_receipts.jsonl",
            authorized_pic_root.resolve(),
        ),
        (
            mirror_jsonl,
            authorized_project_home_root.resolve() / "ledger" / "node_hours.jsonl",
            authorized_project_home_root.resolve(),
        ),
    ]
    for supplied, required, root in expected:
        actual = Path(os.path.abspath(supplied))
        if actual != required:
            raise ValueError(f"Unauthorized ledger path: {actual}; expected {required}")
        require_canonical_path_below(actual, root)
