#!/opt/cray/pe/python/3.11.7/bin/python3
"""Re-verify one reserved immutable snapshot before executing its job contract."""

from __future__ import annotations

import argparse
import hashlib
import os
from pathlib import Path
import subprocess
from typing import Callable

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import record_for_role, require_below
from control_plane_common import require_canonical_path_below, require_ledger_paths
from control_plane_common import require_read_only, validate_launch_contract
from control_plane_common import verify_snapshot_files
from validate_and_reserve_frontier_job import reservation_bound_manifest


SRUN = "/usr/bin/srun"


def _artifact_path(artifact_dir: Path, relative: object) -> Path:
    path = artifact_dir / str(relative)
    return require_canonical_path_below(path, artifact_dir)


def _write_new_text_artifact(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        stream.write(text)


def _launch_actions(
    manifest: dict[str, object],
    *,
    executable_path: str,
    profile_launcher: str,
    runner: Callable[..., object],
) -> None:
    contract = validate_launch_contract(manifest.get("launch_contract"))
    artifact_dir = require_canonical_path_below(
        Path(str(manifest["artifact_dir"])), Path(str(manifest["pic_root"]))
    )
    artifact_dir.mkdir(parents=True, exist_ok=True)
    input_deck = Path(str(record_for_role(manifest, "input-deck")["path"])).resolve()
    _bounded_actions(contract["pre_actions"], manifest, artifact_dir)
    for action in contract["actions"]:
        resources = action["resources"]
        command = [
            profile_launcher,
            SRUN,
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
                directory.mkdir(parents=True, exist_ok=True)
                command.append(str(directory))
        stdout_path = _artifact_path(artifact_dir, action["stdout_artifact"])
        stderr_path = _artifact_path(artifact_dir, action["stderr_artifact"])
        stdout_path.parent.mkdir(parents=True, exist_ok=True)
        stderr_path.parent.mkdir(parents=True, exist_ok=True)
        with stdout_path.open("xb") as stdout, stderr_path.open("xb") as stderr:
            runner(command, check=True, stdout=stdout, stderr=stderr)
    _bounded_actions(contract["post_actions"], manifest, artifact_dir)


def _bounded_actions(
    actions: list[dict[str, object]],
    manifest: dict[str, object],
    artifact_dir: Path,
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
            _write_new_text_artifact(output, digest + "\n")
        elif kind == "artifact_sha256":
            source = _artifact_path(artifact_dir, action["artifact"])
            if not source.is_file():
                raise ValueError(f"Missing artifact for checksum action: {source}")
            output = _artifact_path(artifact_dir, action["output_artifact"])
            _write_new_text_artifact(
                output,
                hashlib.sha256(source.read_bytes()).hexdigest() + "\n",
            )
        else:
            source = _artifact_path(artifact_dir, action["artifact"])
            if not source.is_file() or source.stat().st_size <= 0:
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
    if os.environ.get("PIC_MANIFEST_SHA256") != manifest_sha256:
        raise ValueError("Scheduled manifest digest does not match trampoline argument")
    if os.environ.get("PIC_RESERVATION_ID") != reservation_id:
        raise ValueError("Scheduled reservation ID does not match trampoline argument")
    if os.environ.get("PIC_SUBMISSION_ID") != submission_id:
        raise ValueError("Scheduled submission ID does not match trampoline argument")
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
    )
    if reservation.get("manifest_sha256") != manifest_sha256:
        raise ValueError("Scheduled manifest checksum differs from reservation ledger")
    if reservation.get("submission_id") != submission_id:
        raise ValueError("Scheduled submission ID differs from reservation ledger")
    verify_snapshot_files(manifest)
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
