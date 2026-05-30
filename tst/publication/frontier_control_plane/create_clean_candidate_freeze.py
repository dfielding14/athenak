#!/opt/cray/pe/python/3.11.7/bin/python3
"""Freeze clean source and Orion build provenance for Frontier PIC science."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
import subprocess
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, TRUSTED_GIT
from control_plane_common import git_tree_sha1_from_archive
from control_plane_common import make_tree_read_only, read_json, remove_tree
from control_plane_common import require_below, require_canonical_path_below
from control_plane_common import require_no_symlink_components_below, sha256
from control_plane_common import verify_installed_control_plane, write_json_exclusive


SCRIPT_DIR = Path(__file__).absolute().parent


def _git(source_root: Path, *arguments: str) -> str:
    return subprocess.check_output(
        [TRUSTED_GIT, "-C", str(source_root), *arguments],
        text=True,
    ).strip()


def _utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace(
        "+00:00", "Z"
    )


def _archive_commit(archive: Path) -> str:
    with archive.open("rb") as stream:
        return subprocess.check_output(
            [TRUSTED_GIT, "get-tar-commit-id"], stdin=stream, text=True
        ).strip()


def create_freeze(
    *,
    source_root: Path,
    executable: Path,
    build_profile: Path,
    build_profile_id: str,
    freeze_id: str | None = None,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> Path:
    verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    source_root = source_root.resolve()
    if not source_root.is_dir():
        raise ValueError(f"Source root is not a directory: {source_root}")
    if not build_profile_id.strip():
        raise ValueError("Build-profile ID must not be blank")
    executable = require_canonical_path_below(executable, authorized_pic_root)
    build_profile = require_canonical_path_below(build_profile, authorized_pic_root)
    if not executable.is_file():
        raise FileNotFoundError(f"Missing Orion executable: {executable}")
    if not build_profile.is_file():
        raise FileNotFoundError(f"Missing Orion build profile: {build_profile}")

    status = _git(
        source_root,
        "status",
        "--ignore-submodules=none",
        "--porcelain",
        "--untracked-files=all",
    )
    if status:
        raise ValueError("Source worktree is not clean")
    submodules = _git(source_root, "submodule", "status", "--recursive")
    if submodules:
        raise ValueError("Clean-candidate freezes conservatively reject submodules")
    commit = _git(source_root, "rev-parse", "HEAD")
    tree = _git(source_root, "rev-parse", "HEAD^{tree}")

    identifier = freeze_id or str(uuid.uuid4())
    uuid.UUID(identifier)
    freeze_dir = authorized_pic_root.resolve() / "clean_candidates" / identifier
    require_below(freeze_dir, authorized_pic_root / "clean_candidates")
    require_no_symlink_components_below(freeze_dir, authorized_pic_root)
    candidate_root = freeze_dir.parent
    candidate_root.mkdir(parents=True, exist_ok=True)
    if freeze_dir.exists():
        raise ValueError(f"Clean-candidate freeze already exists: {freeze_dir}")
    temporary = candidate_root / f".tmp-{identifier}-{uuid.uuid4()}"
    temporary.mkdir()
    try:
        staged_archive = temporary / "source.tar"
        subprocess.run(
            [
                TRUSTED_GIT,
                "-C",
                str(source_root),
                "archive",
                "--format=tar",
                f"--output={staged_archive}",
                commit,
            ],
            check=True,
        )
        if _archive_commit(staged_archive) != commit:
            raise ValueError("Generated source archive does not identify HEAD")
        if git_tree_sha1_from_archive(staged_archive) != tree:
            raise ValueError("Generated source archive tree differs from HEAD")
        archive_sha256 = sha256(staged_archive)
        executable_sha256 = sha256(executable)
        profile = read_json(build_profile)
        expected_profile = {
            "schema_version": 1,
            "profile_id": build_profile_id.strip(),
            "source_archive_sha256": archive_sha256,
            "toolchain": str(profile.get("toolchain", "")).strip(),
            "build_command": str(profile.get("build_command", "")).strip(),
            "executable_sha256": executable_sha256,
        }
        if not expected_profile["toolchain"] or not expected_profile["build_command"]:
            raise ValueError("Build profile must declare toolchain and build_command")
        if profile != expected_profile:
            raise ValueError("Build profile does not match clean source and executable")

        staged_profile = temporary / "build_profile.json"
        staged_executable = temporary / "athena"
        shutil.copy2(build_profile, staged_profile)
        shutil.copy2(executable, staged_executable)
        if read_json(staged_profile) != expected_profile:
            raise ValueError("Build profile changed while creating clean-candidate freeze")
        if sha256(staged_executable) != executable_sha256:
            raise ValueError("Executable changed while creating clean-candidate freeze")
        archive = freeze_dir / staged_archive.name
        frozen_profile = freeze_dir / staged_profile.name
        frozen_executable = freeze_dir / staged_executable.name
        manifest = {
            "schema_version": 1,
            "freeze_id": identifier,
            "created_utc": _utc_now(),
            "source": {
                "archive_path": str(archive),
                "archive_sha256": archive_sha256,
                "git_commit": commit,
                "git_tree": tree,
                "worktree_status": "clean",
                "submodule_status": "absent",
            },
            "build": {
                "profile_id": expected_profile["profile_id"],
                "profile_path": str(frozen_profile),
                "profile_sha256": sha256(staged_profile),
                "source_archive_sha256": archive_sha256,
                "toolchain": expected_profile["toolchain"],
                "build_command": expected_profile["build_command"],
                "executable_path": str(frozen_executable),
                "executable_sha256": executable_sha256,
            },
        }
        write_json_exclusive(temporary / "clean_candidate_manifest.json", manifest)
        make_tree_read_only(temporary, executable_names={"athena"})
        os.replace(temporary, freeze_dir)
    finally:
        if temporary.exists():
            remove_tree(temporary)
    manifest_path = freeze_dir / "clean_candidate_manifest.json"
    print(f"{manifest_path} {sha256(manifest_path)}")
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--build-profile", required=True, type=Path)
    parser.add_argument("--build-profile-id", required=True)
    parser.add_argument("--freeze-id")
    args = parser.parse_args()
    create_freeze(
        source_root=args.source_root,
        executable=args.executable,
        build_profile=args.build_profile,
        build_profile_id=args.build_profile_id,
        freeze_id=args.freeze_id,
    )


if __name__ == "__main__":
    main()
