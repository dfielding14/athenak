#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Verify an immutable installed control plane before importing sibling modules."""

from __future__ import annotations

import hashlib
import importlib.abc
import importlib.util
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
from types import ModuleType


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
    "frontier_job.sh",
    "frontier_pic_environment.sh",
    "initialize_frontier_ledger.py",
    "launch_trampoline.py",
    "launch_with_frontier_profile.sh",
    "ledger.py",
    "operator_attestation.py",
    "promote_active_policy.py",
    "q011_pressure_review_packet_verifier.py",
    "reconcile_frontier_job.py",
    "reconcile_manual_frontier_allocations.py",
    "reconcile_q043_registered_execution.py",
    "revalidate_clean_candidate.py",
    "run_installed_control_plane_job.sh",
    "run_control_plane.py",
    "storage_preflight.schema.json",
    "submit_frontier_job.sh",
    "terminal_recovery_handoff.py",
    "validate_and_reserve_frontier_job.py",
    "verify_compute_node_snapshot.py",
    "write_orion_build_profile.py",
]
SOURCE_ONLY_ENTRYPOINTS = {
    "capture_storage_preflight_evidence.py",
    "install_control_plane.py",
}
SOURCE_CONTROL_PLANE_FILES = [
    *CONTROL_PLANE_FILES,
    "capture_storage_preflight_evidence.py",
    "install_control_plane.py",
]


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


def _reject_json_constant(value: str) -> None:
    raise ValueError(f"Non-finite JSON number is not allowed: {value}")


def _reject_duplicate_json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError(f"Duplicate JSON object key is not allowed: {key}")
        value[key] = item
    return value


def _inventory(data: bytes) -> dict[str, object]:
    try:
        inventory = json.loads(
            data.decode("utf-8"),
            parse_constant=_reject_json_constant,
            object_pairs_hook=_reject_duplicate_json_pairs,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("Installed control-plane inventory is not valid UTF-8 JSON") from error
    if (
        not isinstance(inventory, dict)
        or set(inventory) != {"schema_version", "version", "files"}
        or type(inventory.get("schema_version")) is not int
        or inventory.get("schema_version") != 1
        or not isinstance(inventory.get("version"), str)
        or re.fullmatch(r"[0-9a-f]{64}", str(inventory.get("version", ""))) is None
    ):
        raise ValueError("Unsupported installed control-plane inventory")
    records = inventory["files"]
    if not isinstance(records, list):
        raise ValueError("Malformed installed control-plane inventory files")
    for record in records:
        if (
            not isinstance(record, dict)
            or set(record) != {"path", "sha256"}
            or not isinstance(record["path"], str)
            or not isinstance(record["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", record["sha256"]) is None
        ):
            raise ValueError("Malformed installed control-plane inventory record")
    return inventory


def _verify_installed(
    script_dir: Path, directory_descriptor: int
) -> dict[str, bytes]:
    if script_dir.parent.name != "control_plane":
        raise ValueError("Installed runner is not below a control_plane directory")
    if os.fstat(directory_descriptor).st_mode & 0o222:
        raise ValueError("Installed control-plane directory is not read-only")
    if set(os.listdir(directory_descriptor)) != {*CONTROL_PLANE_FILES, "inventory.json"}:
        raise ValueError("Installed control-plane entries differ from required list")
    inventory = _inventory(
        _read_file_at(
            directory_descriptor, "inventory.json", require_read_only=True
        )
    )
    records = inventory["files"]
    if [record["path"] for record in records] != CONTROL_PLANE_FILES:
        raise ValueError("Installed control-plane inventory file list differs from required list")
    if (
        inventory.get("version") != _inventory_digest(records)
        or script_dir.name != inventory["version"]
    ):
        raise ValueError("Installed control-plane inventory digest mismatch")
    sources = {}
    for record in records:
        data = _read_file_at(
            directory_descriptor, record["path"], require_read_only=True
        )
        if hashlib.sha256(data).hexdigest() != record["sha256"]:
            raise ValueError(
                f"Installed control-plane checksum mismatch: {record['path']}"
            )
        sources[record["path"]] = data
    return sources


def _verify_source(
    script_dir: Path,
    directory_descriptor: int,
    target: str,
    *,
    expected_git_commit: str,
) -> dict[str, bytes]:
    if target not in SOURCE_ONLY_ENTRYPOINTS:
        raise ValueError(
            "Source runner may execute only authenticated source-only entrypoints"
        )
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
    if re.fullmatch(r"[0-9a-f]{40}", expected_git_commit) is None:
        raise ValueError("Expected source Git commit is malformed")
    head = subprocess.check_output(
        _git("-C", str(repository), "rev-parse", "HEAD"),
        text=True,
        env=environment,
    ).strip()
    if head != expected_git_commit:
        raise ValueError("Source control-plane HEAD differs from expected Git commit")
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
    sources = {}
    for name, path in zip(SOURCE_CONTROL_PLANE_FILES, paths):
        tracked = subprocess.check_output(
            _git("-C", str(repository), "show", f"{expected_git_commit}:{path}"),
            env=environment,
        )
        if _read_file_at(directory_descriptor, name, require_read_only=False) != tracked:
            raise ValueError(f"Control-plane source differs from pinned HEAD blob: {name}")
        sources[name] = tracked
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
    return sources


class _CapturedSourceLoader(importlib.abc.Loader):
    def __init__(self, name: str, source: bytes, filename: str) -> None:
        self.name = name
        self.source = source
        self.filename = filename

    def create_module(self, spec: object) -> ModuleType | None:
        return None

    def exec_module(self, module: ModuleType) -> None:
        module.__file__ = self.filename
        exec(compile(self.source, self.filename, "exec"), module.__dict__)


class _CapturedSourceFinder(importlib.abc.MetaPathFinder):
    def __init__(self, script_dir: Path, sources: dict[str, bytes]) -> None:
        self.script_dir = script_dir
        self.sources = {
            name[:-3]: source for name, source in sources.items() if name.endswith(".py")
        }

    def find_spec(
        self,
        fullname: str,
        path: object = None,
        target: ModuleType | None = None,
    ) -> object:
        del path, target
        source = self.sources.get(fullname) if "." not in fullname else None
        if source is None:
            return None
        filename = str(self.script_dir / f"{fullname}.py")
        return importlib.util.spec_from_loader(
            fullname,
            _CapturedSourceLoader(fullname, source, filename),
            origin=filename,
        )


def _execute_captured(script_dir: Path, target: str, sources: dict[str, bytes]) -> None:
    finder = _CapturedSourceFinder(script_dir, sources)
    filename = str(script_dir / target)
    globals_dict = {
        "__name__": "__main__",
        "__file__": filename,
        "__package__": None,
        "__cached__": None,
        "__builtins__": __builtins__,
    }
    sys.meta_path.insert(0, finder)
    try:
        exec(compile(sources[target], filename, "exec"), globals_dict)
    finally:
        sys.meta_path.remove(finder)


def main() -> None:
    if len(sys.argv) < 2:
        raise SystemExit(
            "usage: run_control_plane.py "
            "[--expected-git-commit FULL_GIT_COMMIT] ENTRYPOINT [ARG ...]"
        )
    script_dir = Path(os.path.abspath(__file__)).parent
    if Path(__file__).resolve(strict=True).parent != script_dir:
        raise ValueError("Refusing a symlink alias for the control-plane runner")
    arguments = sys.argv[1:]
    expected_git_commit: str | None = None
    if arguments[:1] == ["--expected-git-commit"]:
        if len(arguments) < 3:
            raise SystemExit(
                "usage: run_control_plane.py "
                "--expected-git-commit FULL_GIT_COMMIT ENTRYPOINT [ARG ...]"
            )
        expected_git_commit = arguments[1]
        arguments = arguments[2:]
    target = arguments[0]
    if not target or "/" in target or not target.endswith(".py"):
        raise ValueError(f"Unsupported control-plane entrypoint: {target!r}")
    directory_descriptor = os.open(
        script_dir,
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        entries = set(os.listdir(directory_descriptor))
        if "inventory.json" in entries:
            if expected_git_commit is not None:
                raise ValueError(
                    "Installed control-plane execution does not accept a source Git commit"
                )
            if target not in CONTROL_PLANE_FILES:
                raise ValueError(f"Unsupported installed control-plane entrypoint: {target!r}")
            sources = _verify_installed(script_dir, directory_descriptor)
        else:
            if expected_git_commit is None:
                raise ValueError(
                    "Source control-plane execution requires an expected Git commit"
                )
            sources = _verify_source(
                script_dir,
                directory_descriptor,
                target,
                expected_git_commit=expected_git_commit,
            )
        sys._pic_control_plane_bootstrapped = True
        sys.argv = [str(script_dir / target), *arguments[1:]]
        _execute_captured(script_dir, target, sources)
    finally:
        os.close(directory_descriptor)


if __name__ == "__main__":
    main()
