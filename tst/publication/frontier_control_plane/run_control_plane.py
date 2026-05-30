#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Verify an immutable installed control plane before importing sibling modules."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import runpy
import stat
import subprocess
import sys


TRUSTED_GIT = "/usr/bin/git"
TRUSTED_GIT_OPTIONS = [
    "-c",
    "core.fsmonitor=false",
    "-c",
    "core.hooksPath=/dev/null",
]
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
    "run_control_plane.py",
    "submit_frontier_job.sh",
    "validate_and_reserve_frontier_job.py",
    "verify_compute_node_snapshot.py",
    "write_orion_build_profile.py",
]
SOURCE_CONTROL_PLANE_FILES = [*CONTROL_PLANE_FILES, "install_control_plane.py"]


def _git(*arguments: str) -> list[str]:
    return [TRUSTED_GIT, *TRUSTED_GIT_OPTIONS, *arguments]


def _git_environment() -> dict[str, str]:
    return {
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "HOME": "/",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
    }


def _read_file_at(
    directory_descriptor: int, name: str, *, require_read_only: bool
) -> bytes:
    if not name or "/" in name or Path(name).name != name:
        raise ValueError(f"Invalid control-plane filename: {name!r}")
    descriptor = os.open(
        name,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_descriptor,
    )
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"Control-plane source is not a regular file: {name}")
        if require_read_only and metadata.st_mode & 0o222:
            raise ValueError(f"Installed control-plane file is not read-only: {name}")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def _inventory_digest(records: list[dict[str, str]]) -> str:
    payload = json.dumps(records, separators=(",", ":"), sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _verify_installed(script_dir: Path, directory_descriptor: int) -> None:
    if script_dir.parent.name != "control_plane":
        raise ValueError("Installed runner is not below a control_plane directory")
    if os.fstat(directory_descriptor).st_mode & 0o222:
        raise ValueError("Installed control-plane directory is not read-only")
    if set(os.listdir(directory_descriptor)) != {*CONTROL_PLANE_FILES, "inventory.json"}:
        raise ValueError("Installed control-plane entries differ from required list")
    inventory = json.loads(
        _read_file_at(
            directory_descriptor, "inventory.json", require_read_only=True
        ).decode("utf-8")
    )
    if not isinstance(inventory, dict) or inventory.get("schema_version") != 1:
        raise ValueError("Unsupported installed control-plane inventory")
    records = inventory.get("files")
    if not isinstance(records, list) or [
        record.get("path") if isinstance(record, dict) else None for record in records
    ] != CONTROL_PLANE_FILES:
        raise ValueError("Installed control-plane inventory file list differs from required list")
    if (
        inventory.get("version") != _inventory_digest(records)
        or script_dir.name != inventory["version"]
    ):
        raise ValueError("Installed control-plane inventory digest mismatch")
    for record in records:
        if set(record) != {"path", "sha256"}:
            raise ValueError("Malformed installed control-plane inventory record")
        data = _read_file_at(
            directory_descriptor, record["path"], require_read_only=True
        )
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(
                f"Installed control-plane checksum mismatch: {record['path']}"
            )


def _verify_source(script_dir: Path, directory_descriptor: int, target: str) -> None:
    if target != "install_control_plane.py":
        raise ValueError("Source runner may execute only install_control_plane.py")
    environment = _git_environment()
    repository = Path(
        subprocess.check_output(
            _git("-C", str(script_dir), "rev-parse", "--show-toplevel"),
            text=True,
            env=environment,
        ).strip()
    )
    paths = [
        str((script_dir / name).relative_to(repository))
        for name in SOURCE_CONTROL_PLANE_FILES
    ]
    subprocess.run(
        _git("-C", str(repository), "ls-files", "--error-unmatch", "--", *paths),
        check=True,
        stdout=subprocess.DEVNULL,
        env=environment,
    )
    head = subprocess.check_output(
        _git("-C", str(repository), "rev-parse", "HEAD"),
        text=True,
        env=environment,
    ).strip()
    status_command = _git(
        "-C",
        str(repository),
        "status",
        "--porcelain=v1",
        "--untracked-files=all",
        "--",
        *paths,
    )
    if subprocess.check_output(status_command, text=True, env=environment):
        raise ValueError("Production control-plane install requires clean tracked source files")
    for name, path in zip(SOURCE_CONTROL_PLANE_FILES, paths):
        tracked = subprocess.check_output(
            _git("-C", str(repository), "show", f"{head}:{path}"),
            env=environment,
        )
        if _read_file_at(directory_descriptor, name, require_read_only=False) != tracked:
            raise ValueError(f"Control-plane source differs from pinned HEAD blob: {name}")
    if (
        subprocess.check_output(
            _git("-C", str(repository), "rev-parse", "HEAD"),
            text=True,
            env=environment,
        ).strip()
        != head
        or subprocess.check_output(status_command, text=True, env=environment)
    ):
        raise ValueError("Control-plane source changed while verifying pinned HEAD blobs")


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit("usage: run_control_plane.py ENTRYPOINT [ARG ...]")
    script_dir = Path(os.path.abspath(__file__)).parent
    if Path(__file__).resolve(strict=True).parent != script_dir:
        raise ValueError("Refusing a symlink alias for the control-plane runner")
    target = sys.argv[1]
    if not target or "/" in target or not target.endswith(".py"):
        raise ValueError(f"Unsupported control-plane entrypoint: {target!r}")
    directory_descriptor = os.open(
        script_dir,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        if (script_dir / "inventory.json").exists():
            if target not in CONTROL_PLANE_FILES:
                raise ValueError(f"Unsupported installed control-plane entrypoint: {target!r}")
            _verify_installed(script_dir, directory_descriptor)
        else:
            _verify_source(script_dir, directory_descriptor, target)
    finally:
        os.close(directory_descriptor)
    sys.path.insert(0, str(script_dir))
    sys._pic_control_plane_bootstrapped = True
    sys.argv = [str(script_dir / target), *sys.argv[2:]]
    runpy.run_path(str(script_dir / target), run_name="__main__")


if __name__ == "__main__":
    main()
