#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Install an immutable, checksummed Frontier PIC control-plane version."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import json
import os
from pathlib import Path
import re
import stat
import subprocess

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import CONTROL_PLANE_FILES, TRUSTED_GIT, inventory_digest
from control_plane_common import PinnedStagingDirectory, durable_mkdir_parents
from control_plane_common import make_tree_read_only
from control_plane_common import read_json_bytes, read_stable_regular_file
from control_plane_common import require_read_only, sha256_bytes
from control_plane_common import require_no_symlink_components_below
from control_plane_common import trusted_git_command, trusted_git_environment
from control_plane_common import validate_control_plane_inventory
from control_plane_common import verify_installed_control_plane


SCRIPT_DIR = Path(__file__).resolve().parent


def _require_reviewed_source_for_production(
    pic_root: Path,
    *,
    expected_git_commit: str | None,
) -> list[tuple[str, bytes]] | None:
    if pic_root not in {
        Path(os.path.abspath(AUTHORIZED_PIC_ROOT)),
        Path(os.path.abspath(AUTHORIZED_PROJECT_HOME_ROOT)),
    }:
        return None
    if (
        not isinstance(expected_git_commit, str)
        or re.fullmatch(r"[0-9a-f]{40}", expected_git_commit) is None
    ):
        raise ValueError(
            "Production control-plane install requires an expected Git commit"
        )
    repository = Path(
        subprocess.check_output(
            trusted_git_command("-C", str(SCRIPT_DIR), "rev-parse", "--show-toplevel"),
            text=True,
            env=trusted_git_environment(),
        ).strip()
    )
    paths = [
        str((SCRIPT_DIR / name).relative_to(repository))
        for name in CONTROL_PLANE_FILES
    ]
    subprocess.run(
        trusted_git_command(
            "-C",
            str(repository),
            "ls-files",
            "--error-unmatch",
            "--",
            *paths,
        ),
        check=True,
        stdout=subprocess.DEVNULL,
        env=trusted_git_environment(),
    )
    status = subprocess.check_output(
        trusted_git_command(
            "-C",
            str(repository),
            "status",
            "--porcelain=v1",
            "--untracked-files=all",
            "--",
            *paths,
        ),
        text=True,
        env=trusted_git_environment(),
    )
    if status:
        raise ValueError(
            "Production control-plane install requires clean tracked source files"
        )
    head = subprocess.check_output(
        trusted_git_command("-C", str(repository), "rev-parse", "HEAD"),
        text=True,
        env=trusted_git_environment(),
    ).strip()
    if head != expected_git_commit:
        raise ValueError(
            "Production control-plane source HEAD differs from expected Git commit"
        )
    captured_sources = [
        (
            name,
            subprocess.check_output(
                trusted_git_command("-C", str(repository), "show", f"{head}:{path}"),
                env=trusted_git_environment(),
            ),
        )
        for name, path in zip(CONTROL_PLANE_FILES, paths)
    ]
    if (
        subprocess.check_output(
            trusted_git_command("-C", str(repository), "rev-parse", "HEAD"),
            text=True,
            env=trusted_git_environment(),
        ).strip()
        != head
        or subprocess.check_output(
            trusted_git_command(
                "-C",
                str(repository),
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
                "--",
                *paths,
            ),
            text=True,
            env=trusted_git_environment(),
        )
    ):
        raise ValueError("Production control-plane source changed while capturing pinned HEAD blobs")
    return captured_sources


def _snapshot_sources() -> list[tuple[str, bytes]]:
    def read_source(path: Path) -> bytes:
        descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
        try:
            before = os.fstat(descriptor)
            if not stat.S_ISREG(before.st_mode):
                raise ValueError(f"Control-plane source is not a regular file: {path}")
            with os.fdopen(descriptor, "rb", closefd=False) as stream:
                data = stream.read()
            after = os.fstat(descriptor)
            stable_fields = (
                "st_dev",
                "st_ino",
                "st_mode",
                "st_size",
                "st_mtime_ns",
                "st_ctime_ns",
            )
            if (
                any(
                    getattr(before, field) != getattr(after, field)
                    for field in stable_fields
                )
                or len(data) != after.st_size
            ):
                raise ValueError(
                    f"Control-plane source changed while snapshotting: {path}"
                )
            return data
        finally:
            os.close(descriptor)

    return [
        (name, read_source(SCRIPT_DIR / name))
        for name in CONTROL_PLANE_FILES
    ]


def _verify_staged_install(
    staging: Path,
    *,
    inventory: dict[str, object],
    records: list[dict[str, str]],
) -> None:
    expected_names = {*CONTROL_PLANE_FILES, "inventory.json"}
    actual_names = {path.name for path in staging.iterdir()}
    if actual_names != expected_names:
        raise ValueError("Staged control-plane file list differs from required list")
    staged_inventory = read_json_bytes(
        read_stable_regular_file(
            staging / "inventory.json",
            require_read_only_mode=True,
        ),
        label=str(staging / "inventory.json"),
    )
    if staged_inventory != inventory:
        raise ValueError(
            "Staged control-plane inventory differs from captured sources"
        )
    if validate_control_plane_inventory(staged_inventory) != records:
        raise ValueError("Staged control-plane inventory records differ from captured sources")
    for record in records:
        path = staging / record["path"]
        data = read_stable_regular_file(path, require_read_only_mode=True)
        if sha256_bytes(data) != record["sha256"]:
            raise ValueError(f"Staged control-plane checksum mismatch: {path}")
    require_read_only(staging)


def install(pic_root: Path, *, expected_git_commit: str | None = None) -> Path:
    pic_root = Path(os.path.abspath(pic_root))
    reviewed_sources = _require_reviewed_source_for_production(
        pic_root,
        expected_git_commit=expected_git_commit,
    )
    captured_sources = reviewed_sources if reviewed_sources is not None else _snapshot_sources()
    records = [
        {"path": name, "sha256": sha256_bytes(data)}
        for name, data in captured_sources
    ]
    digest = inventory_digest(records)
    destination = pic_root / "control_plane" / digest
    require_no_symlink_components_below(destination, pic_root)
    if destination.exists():
        raise ValueError(f"Control-plane version already exists: {destination}")
    parent = destination.parent
    durable_mkdir_parents(parent, root=pic_root)
    require_no_symlink_components_below(destination, pic_root)
    with PinnedStagingDirectory(
        parent, prefix=f".tmp-{digest}-", root=pic_root
    ) as staging:
        temporary = staging.path
        assert temporary is not None
        for name, data in captured_sources:
            (temporary / name).write_bytes(data)
        inventory = {"schema_version": 1, "version": digest, "files": records}
        (temporary / "inventory.json").write_text(
            json.dumps(inventory, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        make_tree_read_only(
            temporary,
            executable_names={
                path.name for path in temporary.iterdir()
                if path.suffix in {".py", ".sh"}
            },
        )
        _verify_staged_install(temporary, inventory=inventory, records=records)
        if (
            _require_reviewed_source_for_production(
                pic_root,
                expected_git_commit=expected_git_commit,
            )
            != reviewed_sources
        ):
            raise ValueError("Production control-plane pinned HEAD blobs changed before publication")
        staging.publish_tree(destination)
    published_inventory = verify_installed_control_plane(
        destination,
        authorized_pic_root=pic_root,
    )
    if published_inventory != inventory:
        raise ValueError(
            "Published control-plane inventory differs from captured sources"
        )
    if (
        _require_reviewed_source_for_production(
            pic_root,
            expected_git_commit=expected_git_commit,
        )
        != reviewed_sources
    ):
        raise ValueError("Production control-plane pinned HEAD blobs changed after publication")
    print(destination)
    return destination


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pic-root", required=True, type=Path)
    parser.add_argument("--expected-git-commit")
    args = parser.parse_args()
    install(args.pic_root, expected_git_commit=args.expected_git_commit)


if __name__ == "__main__":
    main()
