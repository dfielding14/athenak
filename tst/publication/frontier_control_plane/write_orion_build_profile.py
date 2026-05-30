#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Write one immutable Orion build profile for a clean Frontier PIC candidate."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess
import tempfile

from control_plane_common import AUTHORIZED_PIC_ROOT, TRUSTED_GIT
from control_plane_common import atomic_write_bytes, atomic_write_json
from control_plane_common import direct_submodule_gitlinks
from control_plane_common import durable_mkdir_parents
from control_plane_common import git_commit_tree_from_bytes, git_tree_sha1_from_archive
from control_plane_common import open_directory_below
from control_plane_common import PRODUCTION_BUILD_PROFILE
from control_plane_common import PRODUCTION_BUILD_ENVIRONMENT
from control_plane_common import PRODUCTION_TOOLCHAIN_DESCRIPTION
from control_plane_common import production_build_invocations
from control_plane_common import production_environment_allowlist_bytes
from control_plane_common import measured_production_module_list_bytes
from control_plane_common import production_module_list_bytes
from control_plane_common import read_stable_regular_file_below
from control_plane_common import require_production_build_environment
from control_plane_common import require_production_build_provenance
from control_plane_common import require_canonical_path_below
from control_plane_common import require_no_symlink_components_below
from control_plane_common import sha256, sha256_bytes, source_bundle_sha256
from control_plane_common import trusted_git_command, trusted_git_environment
from control_plane_common import verify_installed_control_plane
from create_clean_candidate_freeze import AUTHORIZED_SOURCE_ROOT
from create_clean_candidate_freeze import _archive_commit, _authorized_source_path
from create_clean_candidate_freeze import _documented_build_paths
from create_clean_candidate_freeze import _validate_documented_build_layout
from create_clean_candidate_freeze import _git_bytes, _reviewed_utf8
from create_clean_candidate_freeze import _snapshot_provenance_inputs
from create_clean_candidate_freeze import _source_identity, _validate_source_status_inputs


SCRIPT_DIR = Path(__file__).absolute().parent
def _archive_digest(
    *,
    repository: Path,
    commit: str,
    tree: str,
    archive: Path,
    gitlinks: dict[str, str],
) -> str:
    subprocess.run(
        trusted_git_command(
            "-C",
            str(repository),
            "archive",
            "--format=tar",
            f"--output={archive}",
            commit,
        ),
        check=True,
        env=trusted_git_environment(),
    )
    if _archive_commit(archive) != commit:
        raise ValueError(f"Generated source archive does not identify HEAD: {repository}")
    if (
        git_tree_sha1_from_archive(
            archive,
            gitlinks=gitlinks,
            reject_symlinks=True,
        )
        != tree
    ):
        raise ValueError(f"Generated source archive tree differs from HEAD: {repository}")
    return sha256(archive)


def _commit_digest(repository: Path, *, commit: str, tree: str) -> str:
    data = _git_bytes(repository, "cat-file", "commit", commit)
    if git_commit_tree_from_bytes(data, expected_commit=commit) != tree:
        raise ValueError(f"Generated source commit object tree differs from HEAD: {repository}")
    return sha256_bytes(data)


def write_profile(
    *,
    source_root: Path,
    fresh_source_root: Path,
    executable: Path,
    output: Path,
    profile_id: str,
    expected_git_commit: str,
    configure_log: Path,
    build_log: Path,
    cmake_cache: Path,
    module_list: Path,
    toolchain_file: Path,
    build_invocations_file: Path,
    git_status_preconfigure_file: Path,
    git_status_file: Path,
    submodule_status_file: Path,
    environment_allowlist_file: Path,
    build_environment_file: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_source_root: Path = AUTHORIZED_SOURCE_ROOT,
) -> Path:
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    authorized_pic_root = authorized_pic_root.resolve()
    source_root = _authorized_source_path(source_root, authorized_source_root)
    authorized_source_root = Path(os.path.abspath(authorized_source_root))
    if not profile_id.strip():
        raise ValueError("Build-profile ID must not be blank")
    if re.fullmatch(r"[0-9a-f]{40}", expected_git_commit) is None:
        raise ValueError("Expected Git commit must be a full lowercase hexadecimal commit")
    executable = require_canonical_path_below(executable, authorized_pic_root)
    output = require_canonical_path_below(
        Path(os.path.abspath(output)), authorized_pic_root
    )
    durable_mkdir_parents(output.parent, root=authorized_pic_root)
    require_no_symlink_components_below(output, authorized_pic_root)
    executable_sha256 = sha256_bytes(
        read_stable_regular_file_below(executable, authorized_pic_root)
    )
    provenance_inputs, provenance_payloads = _snapshot_provenance_inputs(
        {
            "configure_log": configure_log,
            "build_log": build_log,
            "cmake_cache": cmake_cache,
            "module_list": module_list,
            "toolchain": toolchain_file,
            "build_invocations": build_invocations_file,
            "git_status_preconfigure": git_status_preconfigure_file,
            "git_status": git_status_file,
            "submodule_status": submodule_status_file,
            "environment_allowlist": environment_allowlist_file,
            "build_environment": build_environment_file,
        },
        authorized_pic_root=authorized_pic_root,
    )
    _validate_documented_build_layout(
        authorized_pic_root=authorized_pic_root,
        git_commit=expected_git_commit,
        profile_id=profile_id,
        executable=executable,
        build_profile=output,
        provenance_inputs=provenance_inputs,
    )
    _validate_source_status_inputs(source_root, provenance_payloads)
    toolchain = _reviewed_utf8(
        provenance_payloads["toolchain"], label="Toolchain description"
    )
    invocations = json.loads(
        _reviewed_utf8(
            provenance_payloads["build_invocations"], label="Build invocations"
        )
    )
    if (
        not isinstance(invocations, dict)
        or set(invocations) != {"configure", "build"}
        or any(
            not isinstance(invocations[key], list)
            or not invocations[key]
            or any(not isinstance(token, str) or not token for token in invocations[key])
            for key in ["configure", "build"]
        )
    ):
        raise ValueError("Build invocations must contain exact configure and build argv")
    require_production_build_provenance(
        authorized_pic_root=authorized_pic_root,
        git_commit=expected_git_commit,
        profile_id=profile_id.strip(),
        toolchain=toolchain,
        invocations=invocations,
        module_list=provenance_payloads["module_list"],
        environment_allowlist=provenance_payloads["environment_allowlist"],
        build_environment=provenance_payloads["build_environment"],
    )
    commit, tree, submodules = _source_identity(source_root)
    if commit != expected_git_commit:
        raise ValueError("Authorized source HEAD differs from expected Git commit")
    fresh_source_root = Path(os.path.abspath(fresh_source_root))
    if _source_identity(fresh_source_root) != (commit, tree, submodules):
        raise ValueError("Fresh build source differs from the authorized source closure")
    with tempfile.TemporaryDirectory(prefix=".build-profile-") as temporary_text:
        temporary = Path(temporary_text)
        archive_sha256 = _archive_digest(
            repository=source_root,
            commit=commit,
            tree=tree,
            archive=temporary / "source.tar",
            gitlinks=direct_submodule_gitlinks(submodules),
        )
        commit_sha256 = _commit_digest(source_root, commit=commit, tree=tree)
        profile_submodules: list[dict[str, str]] = []
        for index, record in enumerate(submodules):
            module_root = source_root.joinpath(
                *PurePosixPath(record["path"]).parts
            )
            profile_submodules.append(
                {
                    "path": record["path"],
                    "archive_sha256": _archive_digest(
                        repository=module_root,
                        commit=record["git_commit"],
                        tree=record["git_tree"],
                        archive=temporary / f"submodule-{index:04d}.tar",
                        gitlinks=direct_submodule_gitlinks(
                            submodules, parent_path=record["path"]
                        ),
                    ),
                    "commit_sha256": _commit_digest(
                        module_root,
                        commit=record["git_commit"],
                        tree=record["git_tree"],
                    ),
                    "git_commit": record["git_commit"],
                    "git_tree": record["git_tree"],
                }
            )
        if _source_identity(source_root) != (commit, tree, submodules):
            raise ValueError("Source or submodule identity changed while writing profile")
    if sha256_bytes(
        read_stable_regular_file_below(executable, authorized_pic_root)
    ) != executable_sha256:
        raise ValueError("Executable changed while writing profile")
    revalidated_inputs, revalidated_payloads = _snapshot_provenance_inputs(
        {
            label: Path(record["path"]) for label, record in provenance_inputs.items()
        },
        authorized_pic_root=authorized_pic_root,
    )
    _validate_source_status_inputs(source_root, revalidated_payloads)
    if revalidated_inputs != provenance_inputs:
        raise ValueError("Build provenance changed while writing profile")
    profile = {
        "schema_version": 3,
        "profile_id": profile_id.strip(),
        "authorized_source_root": str(authorized_source_root),
        "fresh_source_root": str(fresh_source_root),
        "git_commit": commit,
        "git_tree": tree,
        "source_archive_sha256": archive_sha256,
        "source_commit_sha256": commit_sha256,
        "source_bundle_sha256": source_bundle_sha256(
            archive_sha256, commit_sha256, profile_submodules
        ),
        "toolchain": toolchain,
        "build_invocations_sha256": provenance_inputs["build_invocations"]["sha256"],
        "executable_sha256": executable_sha256,
        "provenance_inputs": provenance_inputs,
        "submodules": profile_submodules,
    }
    profile_sha256 = sha256_bytes(
        (json.dumps(profile, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
            "utf-8"
        )
    )
    atomic_write_json(output, profile, replace=False, root=authorized_pic_root)
    receipt = {
        "schema_version": 1,
        "control_plane_version": inventory["version"],
        "profile_path": str(output),
        "profile_sha256": profile_sha256,
        "source_bundle_sha256": profile["source_bundle_sha256"],
        "fresh_source_root": str(fresh_source_root),
        "build_invocations_sha256": provenance_inputs["build_invocations"]["sha256"],
        "git_status_preconfigure_sha256": provenance_inputs[
            "git_status_preconfigure"
        ]["sha256"],
        "git_status_sha256": provenance_inputs["git_status"]["sha256"],
        "configure_log_sha256": provenance_inputs["configure_log"]["sha256"],
        "build_log_sha256": provenance_inputs["build_log"]["sha256"],
        "executable_path": str(executable),
        "executable_sha256": executable_sha256,
    }
    atomic_write_json(
        output.with_name("profile_receipt.json"),
        receipt,
        replace=False,
        root=authorized_pic_root,
    )
    print(f"{output} {profile_sha256}")
    return output


def _write_exclusive(
    path: Path,
    data: bytes,
    *,
    authorized_pic_root: Path,
    mode: int = 0o444,
) -> None:
    atomic_write_bytes(path, data, mode=mode, replace=False, root=authorized_pic_root)


def _run_logged(
    command: list[str],
    *,
    log: Path,
    environment: dict[str, str],
    authorized_pic_root: Path,
) -> None:
    parent_descriptor = open_directory_below(log.parent, root=authorized_pic_root)
    descriptor = os.open(
        log.name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o444,
        dir_fd=parent_descriptor,
    )
    try:
        with os.fdopen(descriptor, "wb", closefd=False) as stream:
            _execute_logged_command(command, stream=stream, environment=environment)
            stream.flush()
            os.fsync(stream.fileno())
        os.fsync(parent_descriptor)
        require_canonical_path_below(log, authorized_pic_root)
    finally:
        os.close(descriptor)
        os.close(parent_descriptor)


def _execute_logged_command(
    command: list[str],
    *,
    stream: object,
    environment: dict[str, str],
) -> None:
    subprocess.run(
        command,
        stdout=stream,
        stderr=subprocess.STDOUT,
        check=True,
        env=environment,
    )


def _clone_fresh_source(
    source_root: Path,
    destination: Path,
    *,
    commit: str,
    submodules: list[dict[str, str]],
) -> None:
    subprocess.run(
        trusted_git_command(
            "clone", "--no-hardlinks", "--no-checkout", str(source_root), str(destination)
        ),
        check=True,
        env=trusted_git_environment(),
    )
    subprocess.run(
        trusted_git_command("-C", str(destination), "checkout", "--detach", commit),
        check=True,
        env=trusted_git_environment(),
    )
    cloned_paths: list[PurePosixPath] = []
    for record in sorted(
        submodules, key=lambda item: len(PurePosixPath(item["path"]).parts)
    ):
        relative = PurePosixPath(record["path"])
        parent = max(
            (
                path
                for path in cloned_paths
                if relative.parts[: len(path.parts)] == path.parts
            ),
            key=lambda path: len(path.parts),
            default=PurePosixPath(),
        )
        parent_root = destination.joinpath(*parent.parts)
        relative_to_parent = relative.relative_to(parent).as_posix()
        subprocess.run(
            trusted_git_command(
                "-C", str(parent_root), "submodule", "init", "--", relative_to_parent
            ),
            check=True,
            env=trusted_git_environment(),
        )
        original = source_root.joinpath(*relative.parts)
        cloned = destination.joinpath(*relative.parts)
        cloned.parent.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            trusted_git_command(
                "clone", "--no-hardlinks", "--no-checkout", str(original), str(cloned)
            ),
            check=True,
            env=trusted_git_environment(),
        )
        subprocess.run(
            trusted_git_command(
                "-C", str(cloned), "checkout", "--detach", record["git_commit"]
            ),
            check=True,
            env=trusted_git_environment(),
        )
        subprocess.run(
            trusted_git_command(
                "-C",
                str(parent_root),
                "submodule",
                "absorbgitdirs",
                relative_to_parent,
            ),
            check=True,
            env=trusted_git_environment(),
        )
        cloned_paths.append(relative)


def _production_build_environment() -> dict[str, str]:
    environment = dict(PRODUCTION_BUILD_ENVIRONMENT)
    require_production_build_environment(environment)
    return environment


def _git_status(source_root: Path) -> bytes:
    return subprocess.check_output(
        trusted_git_command(
            "-C",
            str(source_root),
            "status",
            "--ignore-submodules=none",
            "--porcelain",
            "--untracked-files=all",
        ),
        env=trusted_git_environment(),
    )


def _submodule_status(source_root: Path) -> bytes:
    return subprocess.check_output(
        trusted_git_command("-C", str(source_root), "submodule", "status", "--recursive"),
        env=trusted_git_environment(),
    )


def build_profile(
    *,
    source_root: Path,
    expected_git_commit: str,
    profile_id: str,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_source_root: Path = AUTHORIZED_SOURCE_ROOT,
) -> Path:
    verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    if (
        Path(os.path.abspath(authorized_pic_root)) == Path(os.path.abspath(AUTHORIZED_PIC_ROOT))
        and profile_id != PRODUCTION_BUILD_PROFILE
    ):
        raise ValueError(f"Unsupported installed Frontier build profile: {profile_id}")
    authorized_pic_root = Path(os.path.abspath(authorized_pic_root))
    production_root = authorized_pic_root == Path(os.path.abspath(AUTHORIZED_PIC_ROOT))
    measured_module_list = (
        measured_production_module_list_bytes()
        if production_root
        else production_module_list_bytes()
    )
    source_root = _authorized_source_path(source_root, authorized_source_root)
    commit, tree, submodules = _source_identity(source_root)
    if commit != expected_git_commit:
        raise ValueError("Authorized source HEAD differs from expected Git commit")
    with tempfile.TemporaryDirectory(prefix=".build-profile-preflight-") as temporary_text:
        temporary = Path(temporary_text)
        _archive_digest(
            repository=source_root,
            commit=commit,
            tree=tree,
            archive=temporary / "source.tar",
            gitlinks=direct_submodule_gitlinks(submodules),
        )
        for index, record in enumerate(submodules):
            module_root = source_root.joinpath(*PurePosixPath(record["path"]).parts)
            _archive_digest(
                repository=module_root,
                commit=record["git_commit"],
                tree=record["git_tree"],
                archive=temporary / f"submodule-{index:04d}.tar",
                gitlinks=direct_submodule_gitlinks(
                    submodules, parent_path=record["path"]
                ),
            )
    paths = _documented_build_paths(
        authorized_pic_root=authorized_pic_root,
        git_commit=commit,
        profile_id=profile_id,
    )
    build_dir = authorized_pic_root / "build" / commit[:12] / profile_id
    cmake_dir = build_dir / "cmake"
    fresh_source = build_dir / "source"
    for path in [build_dir, paths["executable"].parent, paths["configure_log"], paths["build_log"]]:
        require_no_symlink_components_below(path, authorized_pic_root)
        if path.exists():
            raise ValueError(f"Fresh Frontier build path already exists: {path}")
    durable_mkdir_parents(build_dir.parent, root=authorized_pic_root)
    durable_mkdir_parents(paths["executable"].parent.parent, root=authorized_pic_root)
    durable_mkdir_parents(paths["configure_log"].parent, root=authorized_pic_root)
    os.mkdir(build_dir, mode=0o700)
    os.mkdir(paths["executable"].parent)
    _clone_fresh_source(
        source_root,
        fresh_source,
        commit=commit,
        submodules=submodules,
    )
    if _source_identity(fresh_source) != (commit, tree, submodules):
        raise ValueError("Fresh detached source checkout differs from authorized source")
    invocations = production_build_invocations(
        authorized_pic_root=authorized_pic_root,
        git_commit=commit,
        profile_id=PRODUCTION_BUILD_PROFILE,
    )
    if profile_id != PRODUCTION_BUILD_PROFILE:
        invocations = {
            key: [
                token.replace(PRODUCTION_BUILD_PROFILE, profile_id)
                for token in value
            ]
            for key, value in invocations.items()
        }
    build_environment = _production_build_environment()
    invocation_bytes = (
        json.dumps(invocations, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    _write_exclusive(
        paths["build_invocations"], invocation_bytes, authorized_pic_root=authorized_pic_root
    )
    _write_exclusive(
        paths["build_environment"],
        (json.dumps(build_environment, indent=2, sort_keys=True, allow_nan=False) + "\n").encode(
            "utf-8"
        ),
        authorized_pic_root=authorized_pic_root,
    )
    preconfigure = _git_status(fresh_source)
    if preconfigure:
        raise ValueError("Fresh detached source is dirty before configuration")
    _write_exclusive(
        paths["git_status_preconfigure"], preconfigure, authorized_pic_root=authorized_pic_root
    )
    _run_logged(
        invocations["configure"],
        log=paths["configure_log"],
        environment=build_environment,
        authorized_pic_root=authorized_pic_root,
    )
    _run_logged(
        invocations["build"],
        log=paths["build_log"],
        environment=build_environment,
        authorized_pic_root=authorized_pic_root,
    )
    if production_root and measured_production_module_list_bytes() != measured_module_list:
        raise ValueError("Loaded Frontier module provenance changed during production build")
    post_build = _git_status(fresh_source)
    if post_build:
        raise ValueError("Fresh detached source is dirty after build")
    _write_exclusive(paths["git_status"], post_build, authorized_pic_root=authorized_pic_root)
    built_payload = None
    for candidate in [cmake_dir / "src" / "athena", cmake_dir / "athena"]:
        try:
            built_payload = read_stable_regular_file_below(candidate, authorized_pic_root)
            break
        except FileNotFoundError:
            pass
    if built_payload is None:
        raise ValueError("Frontier build did not produce the Athena executable")
    _write_exclusive(
        paths["executable"],
        built_payload,
        authorized_pic_root=authorized_pic_root,
        mode=0o555,
    )
    _write_exclusive(
        paths["cmake_cache"],
        read_stable_regular_file_below(cmake_dir / "CMakeCache.txt", authorized_pic_root),
        authorized_pic_root=authorized_pic_root,
    )
    _write_exclusive(
        paths["module_list"],
        measured_module_list,
        authorized_pic_root=authorized_pic_root,
    )
    _write_exclusive(
        paths["toolchain"],
        (PRODUCTION_TOOLCHAIN_DESCRIPTION + "\n").encode("utf-8"),
        authorized_pic_root=authorized_pic_root,
    )
    _write_exclusive(
        paths["submodule_status"],
        _submodule_status(fresh_source),
        authorized_pic_root=authorized_pic_root,
    )
    _write_exclusive(
        paths["environment_allowlist"],
        production_environment_allowlist_bytes(),
        authorized_pic_root=authorized_pic_root,
    )
    return write_profile(
        source_root=source_root,
        fresh_source_root=fresh_source,
        executable=paths["executable"],
        output=paths["build_profile"],
        profile_id=profile_id,
        expected_git_commit=expected_git_commit,
        configure_log=paths["configure_log"],
        build_log=paths["build_log"],
        cmake_cache=paths["cmake_cache"],
        module_list=paths["module_list"],
        toolchain_file=paths["toolchain"],
        build_invocations_file=paths["build_invocations"],
        git_status_preconfigure_file=paths["git_status_preconfigure"],
        git_status_file=paths["git_status"],
        submodule_status_file=paths["submodule_status"],
        environment_allowlist_file=paths["environment_allowlist"],
        build_environment_file=paths["build_environment"],
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_source_root=authorized_source_root,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--profile-id", required=True)
    parser.add_argument("--expected-git-commit", required=True)
    args = parser.parse_args()
    build_profile(
        source_root=args.source_root,
        profile_id=args.profile_id,
        expected_git_commit=args.expected_git_commit,
    )


if __name__ == "__main__":
    main()
