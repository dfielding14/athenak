#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Re-verify one reserved immutable snapshot before executing its job contract."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
import hashlib
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
from typing import Callable

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import durable_mkdir_parents, open_directory_below
from control_plane_common import record_for_role, require_ledger_paths
from control_plane_common import require_read_only, validate_launch_contract
from control_plane_common import verify_snapshot_files
from validate_and_reserve_frontier_job import _require_run_artifact_dir
from validate_and_reserve_frontier_job import reservation_bound_manifest


SRUN = "/usr/bin/srun"
_DIRECTORY_OPEN_FLAGS = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
_NEW_ARTIFACT_OPEN_FLAGS = (
    os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
)


def _artifact_path(artifact_dir: Path, relative: object) -> Path:
    text = str(relative)
    path = PurePosixPath(text)
    if (
        not text
        or path.is_absolute()
        or not path.parts
        or path.parts != tuple(part for part in path.parts if part not in {"", ".", ".."})
    ):
        raise ValueError("Artifact path must be a non-empty relative path")
    return artifact_dir.joinpath(*path.parts)


def _artifact_relative_parts(artifact_dir: Path, path: Path) -> tuple[str, ...]:
    try:
        parts = path.relative_to(artifact_dir).parts
    except ValueError as error:
        raise ValueError(f"Artifact path is outside launch artifact directory: {path}") from error
    if any(part in {"", ".", ".."} for part in parts):
        raise ValueError(f"Artifact path contains an unsafe component: {path}")
    return parts


def _open_artifact_directory(
    artifact_dir_fd: int,
    artifact_dir: Path,
    directory: Path,
    *,
    create: bool,
) -> int:
    descriptor = os.dup(artifact_dir_fd)
    try:
        for part in _artifact_relative_parts(artifact_dir, directory):
            created = False
            try:
                child_descriptor = os.open(part, _DIRECTORY_OPEN_FLAGS, dir_fd=descriptor)
            except FileNotFoundError:
                if not create:
                    raise
                try:
                    os.mkdir(part, mode=0o755, dir_fd=descriptor)
                    created = True
                except FileExistsError:
                    pass
                child_descriptor = os.open(part, _DIRECTORY_OPEN_FLAGS, dir_fd=descriptor)
            try:
                if created:
                    os.fsync(child_descriptor)
                    os.fsync(descriptor)
            except BaseException:
                os.close(child_descriptor)
                raise
            os.close(descriptor)
            descriptor = child_descriptor
    except BaseException:
        os.close(descriptor)
        raise
    return descriptor


def _mkdir_artifact_directory(
    artifact_dir_fd: int, artifact_dir: Path, directory: Path
) -> None:
    descriptor = _open_artifact_directory(
        artifact_dir_fd, artifact_dir, directory, create=True
    )
    os.close(descriptor)


def _open_artifact_file(
    artifact_dir_fd: int,
    artifact_dir: Path,
    path: Path,
    flags: int,
    *,
    create_parent: bool,
    mode: int = 0o600,
) -> int:
    parts = _artifact_relative_parts(artifact_dir, path)
    if not parts:
        raise ValueError("Artifact file path must name one file")
    parent_descriptor = _open_artifact_directory(
        artifact_dir_fd, artifact_dir, path.parent, create=create_parent
    )
    try:
        descriptor = os.open(parts[-1], flags, mode, dir_fd=parent_descriptor)
    finally:
        os.close(parent_descriptor)
    try:
        regular = stat.S_ISREG(os.fstat(descriptor).st_mode)
    except BaseException:
        os.close(descriptor)
        raise
    if not regular:
        os.close(descriptor)
        raise ValueError(f"Launch artifact is not a regular file: {path}")
    return descriptor


def _open_new_artifact(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> int:
    return _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        _NEW_ARTIFACT_OPEN_FLAGS,
        create_parent=True,
    )


def _read_artifact_bytes(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> bytes:
    descriptor = _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        create_parent=False,
    )
    try:
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def _artifact_size(artifact_dir_fd: int, artifact_dir: Path, path: Path) -> int:
    descriptor = _open_artifact_file(
        artifact_dir_fd,
        artifact_dir,
        path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        create_parent=False,
    )
    try:
        return os.fstat(descriptor).st_size
    finally:
        os.close(descriptor)


def _write_new_text_artifact(
    artifact_dir_fd: int, artifact_dir: Path, path: Path, text: str
) -> None:
    descriptor = _open_new_artifact(artifact_dir_fd, artifact_dir, path)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", closefd=False) as stream:
            stream.write(text)
    finally:
        os.close(descriptor)


def _profile_environment() -> dict[str, str]:
    return {
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
        "PIC_FRONTIER_PROFILE": "frontier_minimum_supported",
    }


def _create_artifact_directory(artifact_dir: Path, *, pic_root: Path) -> int:
    durable_mkdir_parents(artifact_dir.parent, root=pic_root)
    parent_fd = open_directory_below(artifact_dir.parent, root=pic_root)
    child_fd: int | None = None
    try:
        try:
            os.mkdir(artifact_dir.name, mode=0o755, dir_fd=parent_fd)
        except FileExistsError as error:
            raise ValueError(f"Launch artifact directory already exists: {artifact_dir}") from error
        child_fd = os.open(
            artifact_dir.name,
            os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        os.fsync(child_fd)
        os.fsync(parent_fd)
        return child_fd
    except BaseException:
        if child_fd is not None:
            os.close(child_fd)
        raise
    finally:
        os.close(parent_fd)


def _launch_actions(
    manifest: dict[str, object],
    *,
    executable_path: str,
    profile_launcher: str,
    slurm_job_id: str,
    runner: Callable[..., object],
) -> None:
    contract = validate_launch_contract(manifest.get("launch_contract"))
    pic_root = Path(str(manifest["pic_root"]))
    artifact_dir = _require_run_artifact_dir(manifest)
    artifact_dir_fd = _create_artifact_directory(artifact_dir, pic_root=pic_root)
    try:
        input_deck = Path(str(record_for_role(manifest, "input-deck")["path"])).resolve()
        _bounded_actions(contract["pre_actions"], manifest, artifact_dir, artifact_dir_fd)
        for action in contract["actions"]:
            resources = action["resources"]
            command = [
                profile_launcher,
                SRUN,
                f"--jobid={slurm_job_id}",
                f"-N{resources['nodes']}",
                f"-n{resources['tasks']}",
                f"-c{resources['cpus_per_task']}",
                f"--gpus-per-task={resources['gpus_per_task']}",
                f"--gpu-bind={resources['gpu_bind']}",
                executable_path,
            ]
            for argument in action["arguments"]:
                if "literal" in argument:
                    command.append(str(argument["literal"]))
                elif argument.get("snapshot_role") == "input-deck":
                    command.append(str(input_deck))
                else:
                    directory = _artifact_path(artifact_dir, argument["artifact_directory"])
                    _mkdir_artifact_directory(artifact_dir_fd, artifact_dir, directory)
                    command.append(str(directory))
            stdout_path = _artifact_path(artifact_dir, action["stdout_artifact"])
            stderr_path = _artifact_path(artifact_dir, action["stderr_artifact"])
            allowlist_path = _artifact_path(
                artifact_dir, f"{action['action_id']}.environment.allowlist.txt"
            )
            environment = _profile_environment()
            allowlist_fd: int | None = None
            try:
                allowlist_fd = _open_new_artifact(
                    artifact_dir_fd, artifact_dir, allowlist_path
                )
                environment["PIC_RUNTIME_ALLOWLIST_FD"] = str(allowlist_fd)
                environment["PIC_RUNTIME_ALLOWLIST_DIR_FD"] = str(artifact_dir_fd)
                with os.fdopen(
                    _open_new_artifact(artifact_dir_fd, artifact_dir, stdout_path), "wb"
                ) as stdout, os.fdopen(
                    _open_new_artifact(artifact_dir_fd, artifact_dir, stderr_path), "wb"
                ) as stderr:
                    runner(
                        command,
                        check=True,
                        stdout=stdout,
                        stderr=stderr,
                        env=environment,
                        pass_fds=(allowlist_fd, artifact_dir_fd),
                    )
            finally:
                if allowlist_fd is not None:
                    os.close(allowlist_fd)
        _bounded_actions(contract["post_actions"], manifest, artifact_dir, artifact_dir_fd)
    finally:
        os.close(artifact_dir_fd)


def _bounded_actions(
    actions: list[dict[str, object]],
    manifest: dict[str, object],
    artifact_dir: Path,
    artifact_dir_fd: int,
) -> None:
    for action in actions:
        kind = str(action["kind"])
        if kind == "snapshot_sha256":
            record = record_for_role(manifest, str(action["snapshot_role"]))
            source = Path(str(record["path"]))
            require_read_only(source)
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            if digest != record.get("sha256"):
                raise ValueError("Snapshot changed during bounded launch action")
            output = _artifact_path(artifact_dir, action["output_artifact"])
            _write_new_text_artifact(artifact_dir_fd, artifact_dir, output, digest + "\n")
        elif kind == "artifact_sha256":
            source = _artifact_path(artifact_dir, action["artifact"])
            try:
                source_bytes = _read_artifact_bytes(artifact_dir_fd, artifact_dir, source)
            except (FileNotFoundError, NotADirectoryError):
                raise ValueError(f"Missing artifact for checksum action: {source}")
            output = _artifact_path(artifact_dir, action["output_artifact"])
            _write_new_text_artifact(
                artifact_dir_fd,
                artifact_dir,
                output,
                hashlib.sha256(source_bytes).hexdigest() + "\n",
            )
        else:
            source = _artifact_path(artifact_dir, action["artifact"])
            try:
                size = _artifact_size(artifact_dir_fd, artifact_dir, source)
            except (FileNotFoundError, NotADirectoryError):
                size = 0
            if size <= 0:
                raise ValueError(f"Expected a non-empty launch artifact: {source}")


def launch(
    *,
    manifest_path: Path,
    manifest_sha256: str,
    job_script_sha256: str,
    executable_sha256: str,
    reservation_id: str,
    submission_id: str,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    runner: Callable[..., object] = subprocess.run,
    control_plane_dir: Path = Path(__file__).absolute().parent,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    slurm_job_id = os.environ.get("SLURM_JOB_ID", "")
    if not re.fullmatch(r"[0-9]+", slurm_job_id):
        raise ValueError("Trampoline requires a live numeric SLURM_JOB_ID")
    require_ledger_paths(
        ledger_jsonl,
        authorized_pic_root.resolve() / "ledger" / "node_hours.csv",
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    manifest, reservation = reservation_bound_manifest(
        manifest_path,
        reservation_id,
        ledger_jsonl=ledger_jsonl,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        executable_job_id=slurm_job_id,
    )
    if reservation.get("manifest_sha256") != manifest_sha256:
        raise ValueError("Scheduled manifest checksum differs from reservation ledger")
    if reservation.get("submission_id") != submission_id:
        raise ValueError("Scheduled submission ID differs from reservation ledger")
    verify_snapshot_files(manifest, root=authorized_pic_root)
    for record in manifest["snapshot_files"]:
        require_read_only(Path(str(record["path"])))

    job_script = record_for_role(manifest, "job-script")
    executable = record_for_role(manifest, "executable")
    if (
        job_script.get("sha256") != job_script_sha256
        or reservation.get("job_script_sha256") != job_script_sha256
    ):
        raise ValueError("Snapshotted job-script digest differs from scheduled binding")
    if (
        executable.get("sha256") != executable_sha256
        or reservation.get("executable_sha256") != executable_sha256
    ):
        raise ValueError("Snapshotted executable digest differs from scheduled binding")
    executable_path = str(Path(str(executable["path"])).resolve())
    if manifest.get("job_script_executable_env") != "PIC_EXECUTABLE":
        raise ValueError("Manifest does not declare the PIC_EXECUTABLE launch contract")
    declared = os.environ.get("PIC_EXECUTABLE")
    if declared and str(Path(declared).resolve()) != executable_path:
        raise ValueError("PIC_EXECUTABLE differs from verified executable snapshot")
    os.environ["PIC_EXECUTABLE"] = executable_path
    require_read_only(Path(executable_path))
    profile_launcher = control_plane_dir / "launch_with_frontier_profile.sh"
    require_read_only(profile_launcher)
    _launch_actions(
        manifest,
        executable_path=executable_path,
        profile_launcher=str(profile_launcher),
        slurm_job_id=slurm_job_id,
        runner=runner,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--job-script-sha256", required=True)
    parser.add_argument("--executable-sha256", required=True)
    parser.add_argument("--reservation-id", required=True)
    parser.add_argument("--submission-id", required=True)
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)
    args = parser.parse_args()
    launch(
        manifest_path=args.manifest,
        manifest_sha256=args.manifest_sha256,
        job_script_sha256=args.job_script_sha256,
        executable_sha256=args.executable_sha256,
        reservation_id=args.reservation_id,
        submission_id=args.submission_id,
        ledger_jsonl=args.ledger_jsonl,
        receipts_jsonl=args.receipts_jsonl,
        mirror_jsonl=args.mirror_jsonl,
    )


if __name__ == "__main__":
    main()
