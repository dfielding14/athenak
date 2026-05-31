#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Freeze clean source and Orion build provenance for Frontier PIC science."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, BUILD_PROVENANCE_FILENAMES, TRUSTED_GIT
from control_plane_common import PinnedStagingDirectory
from control_plane_common import direct_submodule_gitlinks, git_commit_tree_from_bytes
from control_plane_common import git_tree_sha1_from_archive
from control_plane_common import durable_mkdir_parents
from control_plane_common import make_tree_read_only
from control_plane_common import prepared_artifact_manifest_from_source_archive
from control_plane_common import read_json, read_json_bytes
from control_plane_common import read_stable_regular_file_below
from control_plane_common import require_production_build_provenance
from control_plane_common import require_below, require_canonical_path_below
from control_plane_common import require_no_symlink_components_below, sha256, sha256_bytes
from control_plane_common import source_bundle_sha256
from control_plane_common import trusted_git_command, trusted_git_environment
from control_plane_common import verify_installed_control_plane, write_json_exclusive


SCRIPT_DIR = Path(__file__).absolute().parent
AUTHORIZED_SOURCE_ROOT = Path("/ccs/home/dfielding/athenak-pic")
PROVENANCE_INPUT_FILENAMES = BUILD_PROVENANCE_FILENAMES
NONEMPTY_PROVENANCE_INPUTS = set(PROVENANCE_INPUT_FILENAMES) - {
    "git_status_preconfigure",
    "git_status",
    "submodule_status",
}


def _git(source_root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        trusted_git_command("-C", str(source_root), *arguments),
        text=True,
        env=trusted_git_environment(),
    ).rstrip("\n")


def _git_bytes(source_root: Path, *arguments: str) -> bytes:
    return subprocess.check_output(
        trusted_git_command("-C", str(source_root), *arguments),
        env=trusted_git_environment(),
    )


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def _archive_commit(archive: Path) -> str:
    with archive.open("rb") as stream:
        return subprocess.check_output(
            trusted_git_command("get-tar-commit-id"),
            stdin=stream,
            text=True,
            env=trusted_git_environment(),
        ).strip()


def _authorized_source_path(source_root: Path, authorized_source_root: Path) -> Path:
    source_root = Path(os.path.abspath(source_root))
    authorized_source_root = Path(os.path.abspath(authorized_source_root))
    if source_root != authorized_source_root:
        raise ValueError(
            f"Source root must use the authorized spelling: {authorized_source_root}"
        )
    if not source_root.is_dir():
        raise ValueError(f"Source root is not a directory: {source_root}")
    resolved_source_root = source_root.resolve(strict=True)
    if resolved_source_root != authorized_source_root.resolve(strict=True):
        raise ValueError(f"Source root differs from authorized source: {source_root}")
    return resolved_source_root


def _reviewed_utf8(data: bytes, *, label: str) -> str:
    try:
        text = data.decode("utf-8").strip()
    except UnicodeDecodeError as error:
        raise ValueError(f"{label} must be UTF-8 text") from error
    if not text:
        raise ValueError(f"{label} must not be blank")
    return text


def _require_exact_schema_version(value: dict[str, object], *, expected: int, label: str) -> None:
    if type(value.get("schema_version")) is not int or value["schema_version"] != expected:
        raise ValueError(f"{label} schema version is invalid")


def _documented_build_paths(
    *,
    authorized_pic_root: Path,
    git_commit: str,
    profile_id: str,
) -> dict[str, Path]:
    if re.fullmatch(r"[0-9a-f]{40}", git_commit) is None:
        raise ValueError("Build-profile Git commit must be a full lowercase hexadecimal commit")
    profile_id = profile_id.strip()
    if re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", profile_id) is None:
        raise ValueError("Build-profile ID must be a safe path component")
    authorized_pic_root = Path(os.path.abspath(authorized_pic_root))
    stem = f"{git_commit[:12]}.{profile_id}"
    artifact_dir = authorized_pic_root / "bin" / git_commit[:12] / profile_id
    log_dir = authorized_pic_root / "logs" / "build"
    return {
        "executable": artifact_dir / "athena",
        "build_profile": artifact_dir / "build_profile.json",
        "configure_log": log_dir / f"{stem}.configure.log",
        "build_log": log_dir / f"{stem}.build.log",
        "cmake_cache": artifact_dir / "CMakeCache.txt",
        "module_list": artifact_dir / "modules.txt",
        "toolchain": artifact_dir / "toolchain.txt",
        "build_invocations": artifact_dir / "build-invocations.json",
        "git_status_preconfigure": artifact_dir / "git_status.preconfigure.txt",
        "git_status": artifact_dir / "git_status.txt",
        "submodule_status": artifact_dir / "submodule_status.txt",
        "environment_allowlist": artifact_dir / "environment.allowlist.txt",
        "build_environment": artifact_dir / "build-environment.json",
    }


def _validate_documented_build_layout(
    *,
    authorized_pic_root: Path,
    git_commit: str,
    profile_id: str,
    executable: Path,
    build_profile: Path,
    provenance_inputs: dict[str, dict[str, str]],
) -> None:
    expected = _documented_build_paths(
        authorized_pic_root=authorized_pic_root,
        git_commit=git_commit,
        profile_id=profile_id,
    )
    actual = {
        "executable": executable,
        "build_profile": build_profile,
        **{
            label: Path(record["path"])
            for label, record in provenance_inputs.items()
        },
    }
    for label, expected_path in expected.items():
        actual_path = Path(os.path.abspath(actual[label]))
        if actual_path != expected_path:
            raise ValueError(
                f"Build {label} must use documented Orion layout: {expected_path}"
            )


def _snapshot_provenance_inputs(
    paths: dict[str, Path],
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, dict[str, str]], dict[str, bytes]]:
    if set(paths) != set(PROVENANCE_INPUT_FILENAMES):
        raise ValueError("Build provenance input set is incomplete")
    records: dict[str, dict[str, str]] = {}
    payloads: dict[str, bytes] = {}
    canonical_paths: set[Path] = set()
    for label in PROVENANCE_INPUT_FILENAMES:
        path = require_canonical_path_below(paths[label], authorized_pic_root)
        if path in canonical_paths:
            raise ValueError(f"Build provenance inputs must be distinct: {path}")
        canonical_paths.add(path)
        data = read_stable_regular_file_below(path, authorized_pic_root)
        if label in NONEMPTY_PROVENANCE_INPUTS and not data:
            raise ValueError(f"Build provenance input must not be empty: {label}")
        records[label] = {"path": str(path), "sha256": sha256_bytes(data)}
        payloads[label] = data
    return records, payloads


def _validate_recorded_provenance_inputs(
    value: object,
    *,
    authorized_pic_root: Path,
) -> tuple[dict[str, dict[str, str]], dict[str, bytes]]:
    if not isinstance(value, dict) or set(value) != set(PROVENANCE_INPUT_FILENAMES):
        raise ValueError("Build profile provenance input set is incomplete")
    paths: dict[str, Path] = {}
    expected_digests: dict[str, str] = {}
    for label in PROVENANCE_INPUT_FILENAMES:
        record = value[label]
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            raise ValueError(f"Malformed build provenance input: {label}")
        paths[label] = Path(str(record["path"]))
        expected_digests[label] = str(record["sha256"])
        if re.fullmatch(r"[0-9a-f]{64}", expected_digests[label]) is None:
            raise ValueError(f"Malformed build provenance checksum: {label}")
    records, payloads = _snapshot_provenance_inputs(
        paths, authorized_pic_root=authorized_pic_root
    )
    for label in PROVENANCE_INPUT_FILENAMES:
        if records[label]["sha256"] != expected_digests[label]:
            raise ValueError(f"Build provenance input checksum mismatch: {label}")
    return records, payloads


def _validate_source_status_inputs(source_root: Path, payloads: dict[str, bytes]) -> None:
    if payloads["git_status_preconfigure"]:
        raise ValueError("Recorded preconfigure Git status must be empty for a clean candidate")
    if payloads["git_status"]:
        raise ValueError("Recorded Git status must be empty for a clean candidate")
    try:
        recorded_submodules = payloads["submodule_status"].decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as error:
        raise ValueError("Recorded submodule status must be UTF-8 text") from error
    if recorded_submodules != _git(source_root, "submodule", "status", "--recursive"):
        raise ValueError("Recorded submodule status differs from the authorized source")


def _relative_submodule_path(value: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if (
        not value
        or path.is_absolute()
        or not path.parts
        or value != path.as_posix()
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError(f"Unsafe submodule path: {value!r}")
    return path


def _validated_submodules(source_root: Path) -> list[dict[str, str]]:
    source_root = Path(os.path.abspath(source_root))
    if source_root.resolve() != source_root or not source_root.is_dir():
        raise ValueError(f"Source root is not a canonical directory: {source_root}")
    output = _git(source_root, "submodule", "status", "--recursive")
    parsed: list[tuple[PurePosixPath, str]] = []
    seen: set[str] = set()
    for line in output.splitlines():
        match = re.fullmatch(r"([ +\-U])([0-9a-f]{40}) ([^ ]+)(?: .*)?", line)
        if match is None:
            raise ValueError(f"Cannot parse recursive submodule status: {line!r}")
        state, commit, value = match.groups()
        path = _relative_submodule_path(value)
        normalized = path.as_posix()
        if normalized in seen:
            raise ValueError(f"Duplicate recursive submodule path: {normalized}")
        if state != " ":
            raise ValueError(f"Submodule is not initialized at its pinned commit: {normalized}")
        seen.add(normalized)
        parsed.append((path, commit))

    records: list[dict[str, str]] = []
    for path, commit in sorted(parsed, key=lambda item: item[0].as_posix()):
        normalized = path.as_posix()
        module_root = source_root.joinpath(*path.parts)
        require_no_symlink_components_below(module_root, source_root)
        if module_root.resolve() != module_root or not module_root.is_dir():
            raise ValueError(f"Submodule path is not a canonical directory: {normalized}")
        if Path(_git(module_root, "rev-parse", "--show-toplevel")).resolve() != module_root:
            raise ValueError(f"Submodule path is not its own Git worktree: {normalized}")
        parent = max(
            (
                record
                for record in records
                if path.parts[: len(PurePosixPath(record["path"]).parts)]
                == PurePosixPath(record["path"]).parts
            ),
            key=lambda record: len(PurePosixPath(record["path"]).parts),
            default=None,
        )
        parent_path = PurePosixPath(parent["path"]) if parent else PurePosixPath()
        parent_root = source_root.joinpath(*parent_path.parts)
        relative_to_parent = path.relative_to(parent_path).as_posix()
        pinned_commit = _git(parent_root, "rev-parse", f"HEAD:{relative_to_parent}")
        actual_commit = _git(module_root, "rev-parse", "HEAD")
        if pinned_commit != commit or actual_commit != commit:
            raise ValueError(f"Submodule commit differs from its pinned commit: {normalized}")
        status = _git(
            module_root,
            "status",
            "--ignore-submodules=none",
            "--porcelain",
            "--untracked-files=all",
        )
        if status:
            raise ValueError(f"Submodule worktree is not clean: {normalized}")
        records.append(
            {
                "path": normalized,
                "git_commit": commit,
                "git_tree": _git(module_root, "rev-parse", "HEAD^{tree}"),
                "worktree_status": "clean",
            }
        )
    return records


def _source_identity(source_root: Path) -> tuple[str, str, list[dict[str, str]]]:
    status = _git(
        source_root,
        "status",
        "--ignore-submodules=none",
        "--porcelain",
        "--untracked-files=all",
    )
    if status:
        raise ValueError("Source worktree is not clean")
    return (
        _git(source_root, "rev-parse", "HEAD"),
        _git(source_root, "rev-parse", "HEAD^{tree}"),
        _validated_submodules(source_root),
    )


def _profile_submodules(records: list[dict[str, str]]) -> list[dict[str, str]]:
    return [
        {
            "path": record["path"],
            "archive_sha256": record["archive_sha256"],
            "commit_sha256": record["commit_sha256"],
            "git_commit": record["git_commit"],
            "git_tree": record["git_tree"],
        }
        for record in records
    ]


def create_freeze(
    *,
    source_root: Path,
    executable: Path,
    build_profile: Path,
    build_profile_id: str,
    prepared_artifact_inventory: str,
    freeze_id: str | None = None,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_source_root: Path = AUTHORIZED_SOURCE_ROOT,
) -> Path:
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    if (
        Path(os.path.abspath(authorized_pic_root))
        != Path(os.path.abspath(AUTHORIZED_PIC_ROOT))
        and Path(os.path.abspath(authorized_source_root))
        == Path(os.path.abspath(AUTHORIZED_SOURCE_ROOT))
    ):
        authorized_source_root = source_root
    source_root = _authorized_source_path(source_root, authorized_source_root)
    authorized_source_root = Path(os.path.abspath(authorized_source_root))
    if not build_profile_id.strip():
        raise ValueError("Build-profile ID must not be blank")
    executable = require_canonical_path_below(executable, authorized_pic_root)
    build_profile = require_canonical_path_below(build_profile, authorized_pic_root)
    profile_receipt = require_canonical_path_below(
        build_profile.with_name("profile_receipt.json"), authorized_pic_root
    )
    commit, tree, submodules = _source_identity(source_root)
    executable_payload = read_stable_regular_file_below(executable, authorized_pic_root)
    build_profile_payload = read_stable_regular_file_below(
        build_profile, authorized_pic_root
    )
    profile_receipt_payload = read_stable_regular_file_below(
        profile_receipt, authorized_pic_root
    )

    identifier = freeze_id or str(uuid.uuid4())
    uuid.UUID(identifier)
    freeze_dir = Path(os.path.abspath(authorized_pic_root)) / "clean_candidates" / identifier
    require_below(freeze_dir, authorized_pic_root / "clean_candidates")
    require_no_symlink_components_below(freeze_dir, authorized_pic_root)
    candidate_root = freeze_dir.parent
    durable_mkdir_parents(candidate_root, root=authorized_pic_root)
    if freeze_dir.exists():
        raise ValueError(f"Clean-candidate freeze already exists: {freeze_dir}")
    with PinnedStagingDirectory(
        candidate_root, prefix=f".tmp-{identifier}-", root=authorized_pic_root
    ) as staging:
        temporary = staging.path
        assert temporary is not None
        staged_archive = temporary / "source.tar"
        subprocess.run(
            trusted_git_command(
                "-C",
                str(source_root),
                "archive",
                "--format=tar",
                f"--output={staged_archive}",
                commit,
            ),
            check=True,
            env=trusted_git_environment(),
        )
        if _archive_commit(staged_archive) != commit:
            raise ValueError("Generated source archive does not identify HEAD")
        if (
            git_tree_sha1_from_archive(
                staged_archive,
                gitlinks=direct_submodule_gitlinks(submodules),
                reject_symlinks=True,
            )
            != tree
        ):
            raise ValueError("Generated source archive tree differs from HEAD")
        prepared_artifacts = prepared_artifact_manifest_from_source_archive(
            staged_archive.read_bytes(),
            inventory_path=prepared_artifact_inventory,
        )
        archive_sha256 = sha256(staged_archive)
        staged_commit = temporary / "source.commit"
        staged_commit.write_bytes(_git_bytes(source_root, "cat-file", "commit", commit))
        if git_commit_tree_from_bytes(staged_commit.read_bytes(), expected_commit=commit) != tree:
            raise ValueError("Generated source commit object tree differs from HEAD")
        commit_sha256 = sha256(staged_commit)
        staged_submodules: list[dict[str, str]] = []
        for index, record in enumerate(submodules):
            staged_submodule_archive = temporary / "submodules" / f"{index:04d}.tar"
            staged_submodule_archive.parent.mkdir(parents=True, exist_ok=True)
            module_root = source_root.joinpath(*PurePosixPath(record["path"]).parts)
            subprocess.run(
                trusted_git_command(
                    "-C",
                    str(module_root),
                    "archive",
                    "--format=tar",
                    f"--output={staged_submodule_archive}",
                    record["git_commit"],
                ),
                check=True,
                env=trusted_git_environment(),
            )
            if _archive_commit(staged_submodule_archive) != record["git_commit"]:
                raise ValueError(f"Generated submodule archive does not identify HEAD: {record['path']}")
            if (
                git_tree_sha1_from_archive(
                    staged_submodule_archive,
                    gitlinks=direct_submodule_gitlinks(
                        submodules, parent_path=record["path"]
                    ),
                    reject_symlinks=True,
                )
                != record["git_tree"]
            ):
                raise ValueError(f"Generated submodule archive tree differs from HEAD: {record['path']}")
            staged_submodule_commit = temporary / "submodules" / f"{index:04d}.commit"
            staged_submodule_commit.write_bytes(
                _git_bytes(module_root, "cat-file", "commit", record["git_commit"])
            )
            if (
                git_commit_tree_from_bytes(
                    staged_submodule_commit.read_bytes(),
                    expected_commit=record["git_commit"],
                )
                != record["git_tree"]
            ):
                raise ValueError(
                    f"Generated submodule commit object tree differs from HEAD: {record['path']}"
                )
            staged_submodules.append(
                {
                    **record,
                    "archive_path": str(
                        freeze_dir / "submodules" / staged_submodule_archive.name
                    ),
                    "archive_sha256": sha256(staged_submodule_archive),
                    "commit_path": str(
                        freeze_dir / "submodules" / staged_submodule_commit.name
                    ),
                    "commit_sha256": sha256(staged_submodule_commit),
                }
            )
        if _source_identity(source_root) != (commit, tree, submodules):
            raise ValueError("Source or submodule identity changed while creating freeze")
        profile_submodules = _profile_submodules(staged_submodules)
        bundle_sha256 = source_bundle_sha256(
            archive_sha256, commit_sha256, profile_submodules
        )
        executable_sha256 = sha256_bytes(executable_payload)
        profile = read_json_bytes(build_profile_payload, label=str(build_profile))
        _require_exact_schema_version(profile, expected=3, label="Build profile")
        provenance_inputs, provenance_payloads = _validate_recorded_provenance_inputs(
            profile.get("provenance_inputs"),
            authorized_pic_root=authorized_pic_root,
        )
        _validate_documented_build_layout(
            authorized_pic_root=authorized_pic_root,
            git_commit=commit,
            profile_id=build_profile_id,
            executable=executable,
            build_profile=build_profile,
            provenance_inputs=provenance_inputs,
        )
        _validate_source_status_inputs(source_root, provenance_payloads)
        build_invocations_sha256 = provenance_inputs["build_invocations"]["sha256"]
        expected_profile = {
            "schema_version": 3,
            "profile_id": build_profile_id.strip(),
            "authorized_source_root": str(authorized_source_root),
            "fresh_source_root": str(
                Path(os.path.abspath(authorized_pic_root))
                / "build"
                / commit[:12]
                / build_profile_id.strip()
                / "source"
            ),
            "git_commit": commit,
            "git_tree": tree,
            "source_archive_sha256": archive_sha256,
            "source_commit_sha256": commit_sha256,
            "source_bundle_sha256": bundle_sha256,
            "toolchain": _reviewed_utf8(
                provenance_payloads["toolchain"], label="Toolchain description"
            ),
            "build_invocations_sha256": build_invocations_sha256,
            "executable_sha256": executable_sha256,
            "provenance_inputs": provenance_inputs,
            "submodules": profile_submodules,
        }
        if not expected_profile["toolchain"]:
            raise ValueError("Build profile must declare toolchain")
        require_production_build_provenance(
            authorized_pic_root=authorized_pic_root,
            git_commit=commit,
            profile_id=build_profile_id.strip(),
            toolchain=str(expected_profile["toolchain"]),
            invocations=read_json_bytes(
                provenance_payloads["build_invocations"],
                label="Build invocations",
            ),
            module_list=provenance_payloads["module_list"],
            environment_allowlist=provenance_payloads["environment_allowlist"],
            build_environment=provenance_payloads["build_environment"],
        )
        if profile != expected_profile:
            raise ValueError("Build profile does not match clean source and executable")
        expected_receipt = {
            "schema_version": 1,
            "control_plane_version": inventory["version"],
            "profile_path": str(build_profile),
            "profile_sha256": sha256_bytes(build_profile_payload),
            "source_bundle_sha256": bundle_sha256,
            "fresh_source_root": expected_profile["fresh_source_root"],
            "build_invocations_sha256": build_invocations_sha256,
            "git_status_preconfigure_sha256": provenance_inputs[
                "git_status_preconfigure"
            ]["sha256"],
            "git_status_sha256": provenance_inputs["git_status"]["sha256"],
            "configure_log_sha256": provenance_inputs["configure_log"]["sha256"],
            "build_log_sha256": provenance_inputs["build_log"]["sha256"],
            "executable_path": str(executable),
            "executable_sha256": executable_sha256,
        }
        receipt = read_json_bytes(profile_receipt_payload, label=str(profile_receipt))
        _require_exact_schema_version(receipt, expected=1, label="Build-profile receipt")
        if receipt != expected_receipt:
            raise ValueError("Build-profile receipt does not match the trusted build")

        staged_profile = temporary / "build_profile.json"
        staged_receipt = temporary / "profile_receipt.json"
        staged_executable = temporary / "athena"
        staged_profile.write_bytes(build_profile_payload)
        staged_receipt.write_bytes(profile_receipt_payload)
        staged_executable.write_bytes(executable_payload)
        staged_provenance = temporary / "build_provenance"
        staged_provenance.mkdir()
        for label, filename in PROVENANCE_INPUT_FILENAMES.items():
            staged_input = staged_provenance / filename
            staged_input.write_bytes(provenance_payloads[label])
            if sha256(staged_input) != provenance_inputs[label]["sha256"]:
                raise ValueError(f"Frozen build provenance checksum mismatch: {label}")
        staged_profile_value = read_json(staged_profile)
        _require_exact_schema_version(staged_profile_value, expected=3, label="Staged build profile")
        if staged_profile_value != expected_profile:
            raise ValueError("Build profile changed while creating clean-candidate freeze")
        staged_receipt_value = read_json(staged_receipt)
        _require_exact_schema_version(
            staged_receipt_value, expected=1, label="Staged build-profile receipt"
        )
        if staged_receipt_value != expected_receipt:
            raise ValueError("Build-profile receipt changed while creating clean-candidate freeze")
        if sha256(staged_executable) != executable_sha256:
            raise ValueError("Executable changed while creating clean-candidate freeze")
        revalidated_inputs, revalidated_payloads = _validate_recorded_provenance_inputs(
            expected_profile["provenance_inputs"],
            authorized_pic_root=authorized_pic_root,
        )
        _validate_source_status_inputs(source_root, revalidated_payloads)
        if revalidated_inputs != provenance_inputs:
            raise ValueError("Build provenance changed while creating clean-candidate freeze")
        if _source_identity(source_root) != (commit, tree, submodules):
            raise ValueError("Source or submodule identity changed while creating freeze")
        archive = freeze_dir / staged_archive.name
        frozen_profile = freeze_dir / staged_profile.name
        frozen_executable = freeze_dir / staged_executable.name
        manifest = {
            "schema_version": 4,
            "freeze_id": identifier,
            "created_utc": _utc_now(),
            "prepared_artifacts": prepared_artifacts,
            "source": {
                "archive_path": str(archive),
                "archive_sha256": archive_sha256,
                "commit_path": str(freeze_dir / staged_commit.name),
                "commit_sha256": commit_sha256,
                "source_bundle_sha256": bundle_sha256,
                "git_commit": commit,
                "git_tree": tree,
                "worktree_status": "clean",
                "submodule_status": (
                    "clean_pinned_archived" if staged_submodules else "absent"
                ),
                "submodules": staged_submodules,
            },
            "build": {
                "profile_id": expected_profile["profile_id"],
                "profile_path": str(frozen_profile),
                "profile_sha256": sha256(staged_profile),
                "profile_receipt_path": str(freeze_dir / staged_receipt.name),
                "profile_receipt_sha256": sha256(staged_receipt),
                "source_archive_sha256": archive_sha256,
                "source_commit_sha256": commit_sha256,
                "source_bundle_sha256": bundle_sha256,
                "toolchain": expected_profile["toolchain"],
                "build_invocations_sha256": build_invocations_sha256,
                "executable_path": str(frozen_executable),
                "executable_sha256": executable_sha256,
            },
        }
        write_json_exclusive(temporary / "clean_candidate_manifest.json", manifest)
        make_tree_read_only(temporary, executable_names={"athena"})
        staging.publish_tree(freeze_dir)
    manifest_path = freeze_dir / "clean_candidate_manifest.json"
    print(f"{manifest_path} {sha256(manifest_path)}")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--build-profile", required=True, type=Path)
    parser.add_argument("--build-profile-id", required=True)
    parser.add_argument(
        "--prepared-artifact-inventory",
        required=True,
        help=(
            "committed source-relative JSON path; record paths are source-relative; format: "
            '{"schema_version":1,"paper_decks":[{"path":"inputs/tests/example.athinput",'
            '"sha256":"<64 lowercase hex>"}],"analyzers":[{"path":"tst/publication/'
            'analyze_example.py","sha256":"<64 lowercase hex>"}]}'
        ),
    )
    parser.add_argument("--freeze-id")
    args = parser.parse_args()
    create_freeze(
        source_root=args.source_root,
        executable=args.executable,
        build_profile=args.build_profile,
        build_profile_id=args.build_profile_id,
        prepared_artifact_inventory=args.prepared_artifact_inventory,
        freeze_id=args.freeze_id,
    )


if __name__ == "__main__":
    main()
