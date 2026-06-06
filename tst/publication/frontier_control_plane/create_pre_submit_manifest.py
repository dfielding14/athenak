#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Create an atomically promoted Frontier PIC pre-submit dependency snapshot."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from pathlib import Path
import uuid

from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import CONTROL_PLANE_FILES
from control_plane_common import PinnedStagingDirectory
from control_plane_common import REGISTERED_SCIENCE_SCOPE, SUBMISSION_SCOPES
from control_plane_common import durable_mkdir_parents
from control_plane_common import make_tree_read_only
from control_plane_common import project_home_ledger_root
from control_plane_common import read_json, require_below
from control_plane_common import require_canonical_path_below
from control_plane_common import require_no_symlink_components_below
from control_plane_common import sha256, snapshot_file, verify_installed_control_plane
from control_plane_common import launch_contract_sha256, validate_launch_contract
from control_plane_common import validate_planner_retention_binding, write_json_exclusive
from operator_attestation import validate_sealed_operator_attestation


SCRIPT_DIR = Path(__file__).absolute().parent
Q011_JOB_SCRIPT_SHA256 = (
    "3048493d376dfa7954595e586c7cebe6740460a7455110d0fc4020bedf697999"
)
Q011_INPUT_DECK_SHA256 = (
    "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
)
Q011_LAUNCH_CONTRACT_SHA256 = {
    "2413a91247d32fb6d93d4903bddab65ed518badc1be6c514ad29c03ca8eafb7b",
    "f9f615bfaa4cc18479dcd02ed4f688f3723fc3dd30a50e754e9f69e10616a8ea",
    "8ab4528b55bdf71e4a13047971aa5ccf8da35c28896ba0bf875a7dca368dfd2f",
    "d84b9217c8810f33ffda995f9611b44f7eeb519bb2e7667c4691ab5833895dd4",
}
Q011_PLANNER_RETENTION_REQUIRED = (
    "Registered Q011 science requires a planner-retention binding"
)


def _required(config: dict[str, object], key: str) -> str:
    value = config.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Missing required pre-submit config key: {key}")
    return value.strip()


def _safe_filename_segment(config: dict[str, object], key: str) -> str:
    value = _required(config, key)
    if value in {".", ".."} or Path(value).name != value:
        raise ValueError(f"Pre-submit config key must be one filename segment: {key}")
    return value


def _requires_q011_planner_retention(
    *,
    submission_scope: object,
    authorization_id: object,
    campaign: object,
    test_id: object,
    job_script_sha256: object,
    input_deck_sha256: object,
    launch_contract_sha256_value: object,
) -> bool:
    if submission_scope != REGISTERED_SCIENCE_SCOPE:
        return False
    return any(
        (
            isinstance(authorization_id, str)
            and authorization_id.startswith(("q011-", "q011_")),
            isinstance(campaign, str)
            and campaign.startswith(("q011-", "q011_")),
            isinstance(test_id, str)
            and test_id.startswith(("q011-", "q011_", "pic_parallel_shock_section54_")),
            job_script_sha256 == Q011_JOB_SCRIPT_SHA256,
            input_deck_sha256 == Q011_INPUT_DECK_SHA256,
            isinstance(launch_contract_sha256_value, str)
            and launch_contract_sha256_value in Q011_LAUNCH_CONTRACT_SHA256,
        )
    )


def create_manifest(
    config_path: Path,
    *,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
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
    physical_mode = _required(config, "physical_mode")
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
    durable_mkdir_parents(campaign_dir, root=pic_root)
    with PinnedStagingDirectory(
        campaign_dir, prefix=f".tmp-{submission_id}-", root=pic_root
    ) as staging:
        temporary = staging.path
        assert temporary is not None
        snapshot_dir = temporary / "snapshot"
        snapshot_dir.mkdir()
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
        pre_manifest_attestation = None
        authorization_id = None
        if submission_scope == REGISTERED_SCIENCE_SCOPE:
            authorization_id = _safe_filename_segment(
                config, "registered_science_authorization_id"
            )
            clean_candidate_manifest = require_canonical_path_below(
                Path(_required(config, "clean_candidate_manifest")),
                pic_root / "clean_candidates",
            )
            pre_manifest_attestation = validate_sealed_operator_attestation(
                Path(_required(config, "pre_manifest_attestation")),
                authorization_id=authorization_id,
                phase="pre_manifest",
                control_plane_version=str(inventory["version"]),
                authorized_pic_root=pic_root,
                authorized_project_home_root=project_home_ledger_root(
                    authorized_project_home_root
                ),
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
        job_script_record = snapshot_file(
            Path(_required(config, "job_script")),
            snapshot_dir / "job.sh",
            role="job-script",
            destination_root=snapshot_dir,
        )
        snapshot_files.append(job_script_record)
        snapshot_files.append(
            snapshot_file(
                Path(_required(config, "executable")),
                snapshot_dir / "athena",
                role="executable",
                destination_root=snapshot_dir,
                scrub=False,
            )
        )
        input_deck_record = snapshot_file(
            Path(_required(config, "input_deck")),
            snapshot_dir / f"{test_id}.athinput",
            role="input-deck",
            destination_root=snapshot_dir,
        )
        snapshot_files.append(input_deck_record)
        if (
            _requires_q011_planner_retention(
                submission_scope=submission_scope,
                authorization_id=authorization_id,
                campaign=campaign,
                test_id=test_id,
                job_script_sha256=job_script_record["sha256"],
                input_deck_sha256=input_deck_record["sha256"],
                launch_contract_sha256_value=launch_contract_sha256(launch_contract),
            )
            and config.get("planner_retention") is None
        ):
            raise ValueError(Q011_PLANNER_RETENTION_REQUIRED)
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

        analysis_destinations = set()
        for index, raw_path in enumerate(config.get("analysis_scripts", [])):
            source = Path(str(raw_path))
            snapshot_name = f"{index:03d}-{source.name}" if index == 0 else source.name
            if snapshot_name in analysis_destinations:
                raise ValueError("Analysis script snapshot basenames must be distinct")
            analysis_destinations.add(snapshot_name)
            snapshot_files.append(
                snapshot_file(
                    source,
                    snapshot_dir / "analysis" / snapshot_name,
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
        registered_short_nonproduction = config.get(
            "registered_short_nonproduction", False
        )
        if not isinstance(registered_short_nonproduction, bool):
            raise ValueError(
                "Pre-submit registered_short_nonproduction must be a boolean"
            )
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
            "physical_mode": physical_mode,
            "selected_qos": _required(config, "selected_qos"),
            "qos_selection_reason": _required(config, "qos_selection_reason"),
            "site_policy_checked_utc": _required(config, "site_policy_checked_utc"),
            "registered_short_nonproduction": registered_short_nonproduction,
            "artifact_dir": _required(config, "artifact_dir"),
            "queue_snapshot_sha256": queue_snapshot_record["sha256"],
            "timeout_margin": timeout_margin,
            "snapshot_files": snapshot_files,
        }
        if clean_candidate_manifest is not None:
            manifest["clean_candidate_manifest_path"] = str(clean_candidate_manifest)
            manifest["clean_candidate_manifest_sha256"] = sha256(clean_candidate_manifest)
            manifest["registered_science_authorization_id"] = authorization_id
            prior_case_closures = config.get("prior_case_closures", [])
            if not isinstance(prior_case_closures, list):
                raise ValueError("Registered-science prior-case closures must be a list")
            manifest["prior_case_closures"] = prior_case_closures
            assert pre_manifest_attestation is not None
            manifest["pre_manifest_attestation_path"] = pre_manifest_attestation["path"]
            manifest["pre_manifest_attestation_sha256"] = pre_manifest_attestation[
                "sha256"
            ]
            if config.get("planner_retention") is not None:
                manifest["planner_retention"] = validate_planner_retention_binding(
                    config["planner_retention"],
                    authorized_pic_root=pic_root,
                    expected_clean_candidate_manifest_sha256=sha256(
                        clean_candidate_manifest
                    ),
                )
        elif config.get("planner_retention") is not None:
            raise ValueError(
                "Admission-smoke configs must not claim a planner-retention binding"
            )
        write_json_exclusive(temporary / "pre_submit_manifest.json", manifest)
        make_tree_read_only(snapshot_dir, executable_names={"athena"})
        (temporary / "pre_submit_manifest.json").chmod(0o444)
        staging.publish_tree(submission_dir)
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
