#!/opt/cray/pe/python/3.11.7/bin/python3
"""Create an atomically promoted Frontier PIC pre-submit dependency snapshot."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, CONTROL_PLANE_FILES
from control_plane_common import REGISTERED_SCIENCE_SCOPE, SUBMISSION_SCOPES
from control_plane_common import make_tree_read_only, read_json, remove_tree, require_below
from control_plane_common import require_canonical_path_below
from control_plane_common import require_no_symlink_components_below
from control_plane_common import sha256, snapshot_file, verify_installed_control_plane
from control_plane_common import validate_launch_contract, write_json_exclusive


SCRIPT_DIR = Path(__file__).absolute().parent


def _required(config: dict[str, object], key: str) -> str:
    value = str(config.get(key, "")).strip()
    if not value:
        raise ValueError(f"Missing required pre-submit config key: {key}")
    return value


def _safe_filename_segment(config: dict[str, object], key: str) -> str:
    value = _required(config, key)
    if value in {".", ".."} or Path(value).name != value:
        raise ValueError(f"Pre-submit config key must be one filename segment: {key}")
    return value


def create_manifest(
    config_path: Path,
    *,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
) -> Path:
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    config = read_json(config_path)
    pic_root = Path(_required(config, "pic_root")).resolve()
    if pic_root != authorized_pic_root.resolve():
        raise ValueError(f"Manifest PIC root is not authorized: {pic_root}")
    campaign = _safe_filename_segment(config, "campaign")
    test_id = _safe_filename_segment(config, "test_id")
    submission_scope = _required(config, "submission_scope")
    if submission_scope not in SUBMISSION_SCOPES:
        raise ValueError(f"Unsupported submission scope: {submission_scope}")
    executable_env = _required(config, "job_script_executable_env")
    if executable_env != "PIC_EXECUTABLE":
        raise ValueError("Job scripts must declare job_script_executable_env=PIC_EXECUTABLE")
    launch_contract = validate_launch_contract(config.get("launch_contract"))
    submission_id = str(config.get("submission_id") or uuid.uuid4())
    uuid.UUID(submission_id)

    submission_dir = pic_root / "manifests" / campaign / submission_id
    require_below(submission_dir, pic_root)
    require_no_symlink_components_below(submission_dir, pic_root)
    if submission_dir.exists():
        raise ValueError(f"Submission directory already exists: {submission_dir}")
    campaign_dir = submission_dir.parent
    campaign_dir.mkdir(parents=True, exist_ok=True)
    temporary = campaign_dir / f".tmp-{submission_id}-{uuid.uuid4()}"
    temporary.mkdir()
    snapshot_dir = temporary / "snapshot"
    snapshot_dir.mkdir()
    try:
        snapshot_files = []
        snapshot_files.append(
            snapshot_file(
                config_path,
                snapshot_dir / "pre_submit_config.json",
                role="pre-submit-config",
                destination_root=snapshot_dir,
            )
        )
        clean_candidate_manifest = None
        if submission_scope == REGISTERED_SCIENCE_SCOPE:
            clean_candidate_manifest = require_canonical_path_below(
                Path(_required(config, "clean_candidate_manifest")),
                pic_root / "clean_candidates",
            )
            snapshot_files.append(
                snapshot_file(
                    clean_candidate_manifest,
                    snapshot_dir / "clean_candidate_manifest.json",
                    role="clean-candidate-manifest",
                    destination_root=snapshot_dir,
                    scrub=False,
                )
            )
        elif config.get("clean_candidate_manifest"):
            raise ValueError(
                "Admission-smoke configs must not claim a clean-candidate manifest"
            )
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "job_script")),
                snapshot_dir / "job.sh",
                role="job-script",
                destination_root=snapshot_dir,
            )
        )
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "executable")),
                snapshot_dir / "athena",
                role="executable",
                destination_root=snapshot_dir,
                scrub=False,
            )
        )
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "input_deck")),
                snapshot_dir / f"{test_id}.athinput",
                role="input-deck",
                destination_root=snapshot_dir,
            )
        )
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "environment_profile")),
                snapshot_dir / "frontier_pic_environment.sh",
                role="environment-profile",
                destination_root=snapshot_dir,
            )
        )
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "timeout_margin_artifact")),
                snapshot_dir / "timeout_margin.json",
                role="timeout-margin",
                destination_root=snapshot_dir,
            )
        )
        snapshot_files.append(
            snapshot_file(
                control_plane_dir / "verify_compute_node_snapshot.py",
                snapshot_dir / "verify_compute_node_snapshot.py",
                role="compute-node-verifier",
                destination_root=snapshot_dir,
            )
        )
        snapshot_files.append(
            snapshot_file(
                control_plane_dir / "control_plane_common.py",
                snapshot_dir / "control_plane_common.py",
                role="compute-node-verifier-library",
                destination_root=snapshot_dir,
                scrub=False,
            )
        )

        control_plane_inventory = []
        for name in CONTROL_PLANE_FILES:
            record = snapshot_file(
                control_plane_dir / name,
                snapshot_dir / "control_plane" / name,
                role=f"control-plane-{name}",
                destination_root=snapshot_dir,
                scrub=False,
            )
            snapshot_files.append(record)
            control_plane_inventory.append(
                {"path": name, "sha256": record["source_sha256"]}
            )
        control_plane_version = str(inventory["version"])
        if control_plane_inventory != inventory["files"]:
            raise ValueError("Installed control-plane inventory changed during snapshot")

        for index, raw_path in enumerate(config.get("analysis_scripts", [])):
            source = Path(str(raw_path))
            snapshot_files.append(
                snapshot_file(
                    source,
                    snapshot_dir / "analysis" / f"{index:03d}-{source.name}",
                    role=f"analysis-script-{index:03d}",
                    destination_root=snapshot_dir,
                )
            )

        queue_snapshot_source = Path(_required(config, "queue_snapshot"))
        queue_snapshot_record = snapshot_file(
            queue_snapshot_source,
            snapshot_dir / "queue_snapshot.txt",
            role="queue-snapshot",
            destination_root=snapshot_dir,
        )
        snapshot_files.append(queue_snapshot_record)
        timeout_margin = read_json(Path(_required(config, "timeout_margin_artifact")))
        for record in snapshot_files:
            staged_path = Path(record["path"])
            record["path"] = str(submission_dir / staged_path.relative_to(temporary))

        manifest: dict[str, object] = {
            "schema_version": 1,
            "control_plane_version": control_plane_version,
            "control_plane_inventory": control_plane_inventory,
            "submission_id": submission_id,
            "pic_root": str(pic_root),
            "campaign": campaign,
            "test_id": test_id,
            "submission_scope": submission_scope,
            "job_script_executable_env": executable_env,
            "launch_contract": launch_contract,
            "git_commit": _required(config, "git_commit"),
            "evidence_class": _required(config, "evidence_class"),
            "physical_mode": _required(config, "physical_mode"),
            "selected_qos": _required(config, "selected_qos"),
            "qos_selection_reason": _required(config, "qos_selection_reason"),
            "site_policy_checked_utc": _required(config, "site_policy_checked_utc"),
            "registered_short_nonproduction": bool(
                config.get("registered_short_nonproduction", False)
            ),
            "artifact_dir": _required(config, "artifact_dir"),
            "queue_snapshot_sha256": queue_snapshot_record["sha256"],
            "timeout_margin": timeout_margin,
            "snapshot_files": snapshot_files,
        }
        if clean_candidate_manifest is not None:
            manifest["clean_candidate_manifest_path"] = str(clean_candidate_manifest)
            manifest["clean_candidate_manifest_sha256"] = sha256(clean_candidate_manifest)
        write_json_exclusive(temporary / "pre_submit_manifest.json", manifest)
        make_tree_read_only(snapshot_dir, executable_names={"athena"})
        (temporary / "pre_submit_manifest.json").chmod(0o444)
        os.replace(temporary, submission_dir)
    finally:
        if temporary.exists():
            remove_tree(temporary)
    manifest_path = submission_dir / "pre_submit_manifest.json"
    print(manifest_path)
    return manifest_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    create_manifest(args.config)


if __name__ == "__main__":
    main()
