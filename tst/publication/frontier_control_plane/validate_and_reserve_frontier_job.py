#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Validate and reserve a serialized Frontier PIC job before sbatch."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import io
import os
from pathlib import Path
import pwd
import re
import shlex
import subprocess
import stat
import tarfile
import types
import uuid

from control_plane_common import AUTHORIZED_ACCOUNT, AUTHORIZED_ADMISSION_SMOKE_STATUS
from control_plane_common import AUTHORIZED_CLEAN_CANDIDATE_FREEZE
from control_plane_common import AUTHORIZED_NODE_HOUR_CAP
from control_plane_common import AUTHORIZED_PIC_ROOT
from control_plane_common import AUTHORIZED_PROJECT_HOME_ROOT, SITE_POLICY_MAX_AGE_SECONDS
from control_plane_common import BUILD_PROVENANCE_FILENAMES
from control_plane_common import FRONTIER_ADMISSION_SMOKE_SCOPE, REGISTERED_SCIENCE_SCOPE
from control_plane_common import SUBMISSION_SCOPES
from control_plane_common import TRUSTED_GIT, TRUSTED_SCONTROL, TRUSTED_SQUEUE
from control_plane_common import atomic_write_bytes, atomic_write_json, fsync_directory
from control_plane_common import read_json
from control_plane_common import read_json_bytes
from control_plane_common import record_for_role, sha256_bytes
from control_plane_common import require_below, require_canonical_path_below
from control_plane_common import require_not_symlink, require_read_only, sha256
from control_plane_common import require_ledger_paths, require_storage_policy_unlock_snapshot
from control_plane_common import scheduler_account_matches_authorized
from control_plane_common import utc_datetime, validate_clean_candidate_bundle
from control_plane_common import verify_installed_control_plane
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import launch_contract_sha256, validate_launch_contract
from control_plane_common import trusted_slurm_environment
from control_plane_common import verify_snapshot_files
from ledger import accounting, append_primary_event_locked, ledger_lock
from ledger import latest_reservations, require_explicit_genesis, transition_payload
from ledger import repair_mirrored_state_locked, validate_mirrored_state
from ledger import slurm_walltime_seconds
from ledger import validated_read_only_mirrored_state_snapshot
from operator_attestation import validate_sealed_operator_attestation


DEBUG_MAX_SECONDS = 2 * 60 * 60
SCRIPT_DIR = Path(__file__).absolute().parent
ADMISSION_SMOKE_FIELDS = {
    "campaign": "f0_hipmpi_smoke",
    "test_id": "pic_parser_contract_guards",
    "evidence_class": "frontier_f0_admission_smoke_candidate",
    "physical_mode": "extended_mhd_pic_parser_contract",
    "selected_qos": "debug",
    "registered_short_nonproduction": True,
}
TIMEOUT_MARGIN_KEYS = {
    "athena_walltime_seconds",
    "scheduler_walltime_seconds",
    "environment_profile_sha256",
    "measured_utc",
    "expires_utc",
}
PRE_SUBMIT_MANIFEST_REQUIRED_KEYS = {
    "schema_version",
    "control_plane_version",
    "control_plane_inventory",
    "submission_id",
    "pic_root",
    "campaign",
    "test_id",
    "submission_scope",
    "job_script_executable_env",
    "launch_contract",
    "git_commit",
    "evidence_class",
    "physical_mode",
    "selected_qos",
    "qos_selection_reason",
    "site_policy_checked_utc",
    "registered_short_nonproduction",
    "artifact_dir",
    "queue_snapshot_sha256",
    "timeout_margin",
    "snapshot_files",
}
REGISTERED_SCIENCE_MANIFEST_KEYS = {
    "clean_candidate_manifest_path",
    "clean_candidate_manifest_sha256",
    "registered_science_authorization_id",
    "pre_manifest_attestation_path",
    "pre_manifest_attestation_sha256",
    "prior_case_closures",
}
Q011_PRESSURE_CASES = (
    (
        "ps_p0_1p00",
        "q011-section54-pressure-ps-p0-1p00-v2",
        "q011_section54_pressure_ps_p0_1p00",
        "pic_parallel_shock_section54_pressure_ps_p0_1p00",
    ),
    (
        "ps_p0_0p05",
        "q011-section54-pressure-ps-p0-0p05-v2",
        "q011_section54_pressure_ps_p0_0p05",
        "pic_parallel_shock_section54_pressure_ps_p0_0p05",
    ),
    (
        "ps_p0_0p10",
        "q011-section54-pressure-ps-p0-0p10-v2",
        "q011_section54_pressure_ps_p0_0p10",
        "pic_parallel_shock_section54_pressure_ps_p0_0p10",
    ),
    (
        "ps_p0_0p20",
        "q011-section54-pressure-ps-p0-0p20-v2",
        "q011_section54_pressure_ps_p0_0p20",
        "pic_parallel_shock_section54_pressure_ps_p0_0p20",
    ),
)
Q011_PRESSURE_BY_AUTHORIZATION = {
    authorization_id: (index, case_id, campaign, test_id)
    for index, (case_id, authorization_id, campaign, test_id) in enumerate(
        Q011_PRESSURE_CASES
    )
}
Q011_PRESSURE_BY_CAMPAIGN = {
    campaign: selected for selected in Q011_PRESSURE_CASES for campaign in [selected[2]]
}
Q011_PRESSURE_BY_TEST_ID = {
    test_id: selected for selected in Q011_PRESSURE_CASES for test_id in [selected[3]]
}
Q011_PRESSURE_BY_LAUNCH_CONTRACT_SHA256 = {
    "2413a91247d32fb6d93d4903bddab65ed518badc1be6c514ad29c03ca8eafb7b":
        Q011_PRESSURE_CASES[0],
    "f9f615bfaa4cc18479dcd02ed4f688f3723fc3dd30a50e754e9f69e10616a8ea":
        Q011_PRESSURE_CASES[1],
    "8ab4528b55bdf71e4a13047971aa5ccf8da35c28896ba0bf875a7dca368dfd2f":
        Q011_PRESSURE_CASES[2],
    "d84b9217c8810f33ffda995f9611b44f7eeb519bb2e7667c4691ab5833895dd4":
        Q011_PRESSURE_CASES[3],
}
Q011_JOB_SCRIPT_SHA256 = (
    "3048493d376dfa7954595e586c7cebe6740460a7455110d0fc4020bedf697999"
)
Q011_INPUT_DECK_SHA256 = (
    "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1"
)
Q011_REPAIRED_GENERATOR_PATH = "src/pgen/tests/pic_parallel_shock.cpp"
Q011_REPAIRED_GENERATOR_SHA256 = (
    "c0a01e4960f4ebb1a96bedc61fd918b4f7fc76addb59e0f9db3f9ab35eb8c9f9"
)
Q011_RAW_ANALYZER_SHA256 = (
    "8596d95c9b8952dcb760b10fbe00bf7ba8a2c895713cc3d36dd1efa1aac11ecc"
)
Q011_STRUCTURED_HELPER_SHA256 = (
    "cf090115bcdfd143f67b12339115102b57e74cf3205b1ebb1521c144c3415a5a"
)
Q011_FAILED_V1_GIT_COMMIT = "4972f998589e6f16d1ccff4f6b9facbcb127909d"
Q011_FAILED_V1_CLEAN_CANDIDATE_MANIFEST_SHA256 = (
    "e09b1ab8c7eb017a929a8c063fb53f176b43ca35050b13aa5a0db1f6432b7e97"
)
Q011_FAILED_V1_EXECUTABLE_SHA256 = (
    "c7c3a986fdbd1bb68dd2e8a0aa95df13934acacb33714b38bdb77ef73a1579c5"
)


def _q011_open_read_only_source(path: Path, *, label: str) -> tuple[int, bytes]:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        lexical = os.stat(path, follow_symlinks=False)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_mode & 0o222
            or (lexical.st_dev, lexical.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise ValueError(f"{label} is not a pinned read-only regular file")
        os.lseek(descriptor, 0, os.SEEK_SET)
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            payload = stream.read()
        os.lseek(descriptor, 0, os.SEEK_SET)
        after = os.fstat(descriptor)
        final_lexical = os.stat(path, follow_symlinks=False)
        stable = ("st_dev", "st_ino", "st_mode", "st_size", "st_mtime_ns", "st_ctime_ns")
        if (
            any(getattr(before, field) != getattr(after, field) for field in stable)
            or (final_lexical.st_dev, final_lexical.st_ino)
            != (after.st_dev, after.st_ino)
            or len(payload) != after.st_size
        ):
            raise ValueError(f"{label} changed while reading")
        return descriptor, payload
    except BaseException:
        os.close(descriptor)
        raise


def _strict_json_equal(left: object, right: object) -> bool:
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(
            _strict_json_equal(left[key], right[key]) for key in left
        )
    if isinstance(left, list):
        return len(left) == len(right) and all(
            _strict_json_equal(left_value, right_value)
            for left_value, right_value in zip(left, right)
        )
    return left == right


def _q011_selected_pressure(
    manifest: dict[str, object],
) -> tuple[int, str, str, str] | None:
    authorization = Q011_PRESSURE_BY_AUTHORIZATION.get(
        manifest.get("registered_science_authorization_id")
    )
    campaign = Q011_PRESSURE_BY_CAMPAIGN.get(manifest.get("campaign"))
    test_id = Q011_PRESSURE_BY_TEST_ID.get(manifest.get("test_id"))
    try:
        contract_sha256 = launch_contract_sha256(manifest.get("launch_contract"))
    except ValueError:
        contract_sha256 = None
    contract = Q011_PRESSURE_BY_LAUNCH_CONTRACT_SHA256.get(contract_sha256)
    try:
        job_script_sha256 = record_for_role(manifest, "job-script").get("sha256")
        input_deck_sha256 = record_for_role(manifest, "input-deck").get("sha256")
    except ValueError:
        job_script_sha256 = None
        input_deck_sha256 = None
    if not any(
        (
            authorization,
            campaign,
            test_id,
            contract,
            job_script_sha256 == Q011_JOB_SCRIPT_SHA256,
            input_deck_sha256 == Q011_INPUT_DECK_SHA256,
        )
    ):
        return None
    if (
        authorization is None
        or campaign != Q011_PRESSURE_CASES[authorization[0]]
        or test_id != Q011_PRESSURE_CASES[authorization[0]]
        or contract != Q011_PRESSURE_CASES[authorization[0]]
        or job_script_sha256 != Q011_JOB_SCRIPT_SHA256
        or input_deck_sha256 != Q011_INPUT_DECK_SHA256
    ):
        raise ValueError("Q011-equivalent registered-science aliases are forbidden")
    return authorization


def _verify_installed_control_plane_pair(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    verify_installed_control_plane(
        authorized_project_home_root / "control_plane" / str(inventory["version"]),
        authorized_pic_root=authorized_project_home_root,
    )
    return inventory


def _require_run_artifact_dir(manifest: dict[str, object]) -> Path:
    pic_root = Path(str(manifest["pic_root"])).resolve()
    runs_root = pic_root / "runs"
    artifact_dir = require_canonical_path_below(
        Path(str(manifest["artifact_dir"])), runs_root
    )
    relative = artifact_dir.relative_to(runs_root)
    if relative.parts != (
        str(manifest["campaign"]),
        str(manifest["submission_id"]),
    ):
        raise ValueError(
            "Manifest artifact directory must be runs/<campaign>/<submission-id>"
        )
    return artifact_dir


def _walltime_seconds(value: str) -> int:
    try:
        return slurm_walltime_seconds(value)
    except ValueError as error:
        raise ValueError(f"Unsupported Slurm walltime: {value}") from error


def _directives(job_script: Path) -> dict[str, str]:
    options: dict[str, str] = {}
    aliases = {
        "-p": "partition",
        "--partition": "partition",
        "-A": "account",
        "--account": "account",
        "-q": "qos",
        "--qos": "qos",
        "-N": "nodes",
        "--nodes": "nodes",
        "-t": "time",
        "--time": "time",
        "-o": "output",
        "--output": "output",
    }
    for line in job_script.read_text(encoding="utf-8").splitlines():
        if not line.startswith("#SBATCH"):
            continue
        tokens = shlex.split(line[len("#SBATCH"):].strip())
        if not tokens:
            continue
        option, separator, attached = tokens[0].partition("=")
        key = aliases.get(option)
        if key is None:
            continue
        value = attached if separator else (tokens[1] if len(tokens) == 2 else "")
        if not value:
            raise ValueError(f"Missing value for Slurm directive: {line}")
        if key in options:
            raise ValueError(f"Duplicate Slurm directive: {key}")
        options[key] = value
    missing = {"account", "partition", "qos", "nodes", "time", "output"} - options.keys()
    if missing:
        raise ValueError(f"Missing required Slurm directives: {sorted(missing)}")
    return options


def _queue_output() -> str:
    return subprocess.check_output(
        [
            TRUSTED_SQUEUE,
            "-u",
            pwd.getpwuid(os.getuid()).pw_name,
            "-h",
            "-o",
            "%i|%P|%q|%T|%j|%k",
        ],
        text=True,
        env=trusted_slurm_environment(),
    )


def _scheduler_job_output(job_id: str) -> str:
    return subprocess.check_output(
        [TRUSTED_SCONTROL, "show", "job", "--oneliner", job_id],
        text=True,
        env=trusted_slurm_environment(),
    )


def _verify_scheduler_job_binding(job_id: str, reservation_id: str) -> dict[str, str]:
    output = _scheduler_job_output(job_id).strip()
    fields: dict[str, str] = {}
    for item in output.split():
        key, separator, value = item.partition("=")
        if separator:
            fields[key] = value
    if fields.get("JobId") != job_id:
        raise ValueError("Slurm job identity does not match the requested attachment")
    if fields.get("Comment") != f"pic-reservation={reservation_id}":
        raise ValueError("Slurm job comment does not bind the PIC reservation")
    if not scheduler_account_matches_authorized(fields.get("Account")):
        raise ValueError("Slurm job account does not match the authorized PIC account")
    return fields


def _verify_scheduler_job(job_id: str, reservation_id: str) -> None:
    fields = _verify_scheduler_job_binding(job_id, reservation_id)
    if fields.get("JobState") not in {"PENDING", "CONFIGURING", "RUNNING", "COMPLETING"}:
        raise ValueError("Slurm job is not in an attachable scheduler state")


def _require_scheduler_output_path(value: str, authorized_pic_root: Path) -> Path:
    log_root = authorized_pic_root.resolve() / "logs" / "slurm"
    output = require_canonical_path_below(Path(value), log_root)
    expected = log_root / "%x.%j.log"
    if output != expected:
        raise ValueError(f"Slurm output must use the dedicated PIC log path: {expected}")
    return output


def _require_reservation_policy_snapshot(
    reservation: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> None:
    _, snapshot = require_storage_policy_unlock_snapshot(
        control_plane_version=str(reservation["control_plane_version"]),
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    for field in ["active_policy_sha256", "active_promotion_sha256"]:
        if reservation.get(field) != snapshot[field]:
            raise ValueError("Reservation belongs to another active-policy generation")


def _append_locked(
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    event: dict[str, object],
) -> dict[str, object]:
    return append_primary_event_locked(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        event,
        mirror_transport="filesystem_copy",
    )


def _pending_marker_path(authorized_pic_root: Path) -> Path:
    return authorized_pic_root.resolve() / "ledger" / "pending_submission.json"


def _write_pending_marker(path: Path, value: dict[str, object]) -> None:
    atomic_write_json(path, value)


def _matching_pending_marker(path: Path, reservation_id: str) -> dict[str, object] | None:
    if not path.is_file():
        return None
    marker = read_json(path)
    if marker.get("reservation_id") != reservation_id:
        raise ValueError("Pending marker belongs to another reservation")
    return marker


def _require_current_reservation_marker(
    marker: dict[str, object],
    reservation: dict[str, object],
    *,
    control_plane_version: str,
) -> None:
    expected = {
        "schema_version": 2,
        "state": marker.get("state"),
        "reservation_id": reservation["reservation_id"],
        "submission_id": reservation["submission_id"],
        "manifest_path": reservation["manifest_path"],
        "manifest_sha256": reservation["manifest_sha256"],
        "control_plane_version": control_plane_version,
    }
    if marker.get("state") in {"scheduler_job_id_received", "submitted_not_attached"}:
        expected["job_id"] = marker.get("job_id")
    if (
        reservation.get("control_plane_version") != control_plane_version
        or type(marker.get("schema_version")) is not int
        or marker != expected
    ):
        raise ValueError("Pending marker does not match the current reservation")


def _clear_matching_pending_marker(path: Path, reservation_id: str) -> None:
    if _matching_pending_marker(path, reservation_id) is not None:
        path.unlink()
        fsync_directory(path.parent)


def _write_reservation_attachments(
    manifest_path: Path,
    *,
    reservation_id: str,
    manifest_sha256: str,
    allow_existing: bool = False,
) -> None:
    attachments = {
        manifest_path.parent / "reservation_id.txt": reservation_id + "\n",
        manifest_path.parent / "manifest_sha256.txt": manifest_sha256 + "\n",
    }
    for path, text in attachments.items():
        require_not_symlink(path)
        if path.exists():
            if not allow_existing:
                raise ValueError(f"Reservation attachment already exists: {path}")
            if path.read_text(encoding="utf-8") != text:
                raise ValueError(f"Reservation attachment content mismatch: {path}")
            require_read_only(path)
        else:
            atomic_write_bytes(path, text.encode("utf-8"))


def _records_with_matching_mirror(
    ledger_jsonl: Path, receipts_jsonl: Path, mirror_jsonl: Path
) -> list[dict[str, object]]:
    records = validate_mirrored_state(ledger_jsonl, receipts_jsonl, mirror_jsonl)
    require_explicit_genesis(records)
    return records


def _verify_manifest(
    manifest_path: Path,
    *,
    control_plane_dir: Path,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> dict[str, object]:
    manifest_path = require_canonical_path_below(
        manifest_path, authorized_pic_root.resolve() / "manifests"
    )
    inventory = verify_installed_control_plane(
        control_plane_dir, authorized_pic_root=authorized_pic_root
    )
    manifest = read_json(manifest_path)
    if type(manifest.get("schema_version")) is not int or manifest.get("schema_version") != 1:
        raise ValueError("Unsupported pre-submit manifest schema")
    scope = manifest.get("submission_scope")
    if scope not in SUBMISSION_SCOPES:
        raise ValueError("Manifest does not declare a recognized submission scope")
    expected_manifest_keys = PRE_SUBMIT_MANIFEST_REQUIRED_KEYS | (
        REGISTERED_SCIENCE_MANIFEST_KEYS if scope == REGISTERED_SCIENCE_SCOPE else set()
    )
    if set(manifest) != expected_manifest_keys:
        raise ValueError("Pre-submit manifest root schema is malformed")
    for key in [
        "control_plane_version",
        "submission_id",
        "pic_root",
        "campaign",
        "test_id",
        "submission_scope",
        "job_script_executable_env",
        "git_commit",
        "evidence_class",
        "physical_mode",
        "selected_qos",
        "qos_selection_reason",
        "site_policy_checked_utc",
        "artifact_dir",
        "queue_snapshot_sha256",
    ]:
        if not isinstance(manifest.get(key), str) or not manifest[key]:
            raise ValueError(f"Pre-submit manifest {key} must be non-empty text")
    if type(manifest.get("registered_short_nonproduction")) is not bool:
        raise ValueError("Pre-submit manifest short-job flag must be a boolean")
    if re.fullmatch(r"[0-9a-f]{64}", manifest["queue_snapshot_sha256"]) is None:
        raise ValueError("Pre-submit manifest queue snapshot digest is malformed")
    pic_root = Path(manifest["pic_root"]).resolve()
    if pic_root != authorized_pic_root.resolve():
        raise ValueError(f"Manifest PIC root is not authorized: {pic_root}")
    if manifest.get("control_plane_version") != inventory["version"]:
        raise ValueError("Manifest control-plane version is not the active frozen version")
    if manifest.get("control_plane_inventory") != inventory["files"]:
        raise ValueError("Manifest control-plane inventory does not match installed files")
    require_canonical_path_below(manifest_path, pic_root / "manifests")
    require_read_only(manifest_path)
    _require_run_artifact_dir(manifest)
    verify_snapshot_files(manifest, root=authorized_pic_root)
    for record in manifest["snapshot_files"]:
        path = require_canonical_path_below(Path(str(record["path"])), pic_root)
        require_read_only(path)
    if manifest.get("job_script_executable_env") != "PIC_EXECUTABLE":
        raise ValueError("Manifest does not declare the PIC_EXECUTABLE launch contract")
    if record_for_role(manifest, "environment-profile").get("sha256") != sha256(
        control_plane_dir / "frontier_pic_environment.sh"
    ):
        raise ValueError("Manifest environment profile differs from installed trusted profile")
    validate_launch_contract(manifest.get("launch_contract"))
    clean_records = [
        record for record in manifest["snapshot_files"]
        if isinstance(record, dict) and record.get("role") == "clean-candidate-manifest"
    ]
    if scope == REGISTERED_SCIENCE_SCOPE:
        if len(clean_records) != 1:
            raise ValueError("Registered science requires one clean-candidate snapshot")
        if not manifest.get("clean_candidate_manifest_path"):
            raise ValueError("Registered science is missing the clean-candidate path")
        if not manifest.get("clean_candidate_manifest_sha256"):
            raise ValueError("Registered science is missing the clean-candidate digest")
        if not manifest.get("registered_science_authorization_id"):
            raise ValueError("Registered science is missing its authorization ID")
        if not isinstance(manifest.get("prior_case_closures"), list):
            raise ValueError("Registered science prior-case closures must be a list")
        pre_manifest_attestation = validate_sealed_operator_attestation(
            Path(str(manifest.get("pre_manifest_attestation_path", ""))),
            authorization_id=str(manifest["registered_science_authorization_id"]),
            phase="pre_manifest",
            control_plane_version=str(inventory["version"]),
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        if (
            pre_manifest_attestation["path"]
            != manifest.get("pre_manifest_attestation_path")
            or pre_manifest_attestation["sha256"]
            != manifest.get("pre_manifest_attestation_sha256")
        ):
            raise ValueError("Registered science pre-manifest attestation binding differs")
    else:
        if clean_records:
            raise ValueError("Admission smoke must not carry a clean-candidate snapshot")
        if (
            "clean_candidate_manifest_path" in manifest
            or "clean_candidate_manifest_sha256" in manifest
            or "registered_science_authorization_id" in manifest
            or "pre_manifest_attestation_path" in manifest
            or "pre_manifest_attestation_sha256" in manifest
        ):
            raise ValueError("Admission smoke must not claim a clean-candidate freeze")
    return manifest


def _mapping(record: dict[str, object], key: str) -> dict[str, object]:
    value = record.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Clean-candidate manifest has no {key} object")
    return value


def _text(record: dict[str, object], key: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Clean-candidate manifest has no {key}")
    return value.strip()


def _digest(record: dict[str, object], key: str) -> str:
    value = _text(record, key)
    if not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ValueError(f"Clean-candidate {key} is not a SHA-256 digest")
    return value


def _read_regular_file_at(directory_fd: int, name: str, *, label: str) -> bytes:
    if not name or "/" in name:
        raise ValueError(f"{label} has an invalid fixed-layout name")
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(name, flags, dir_fd=directory_fd)
    try:
        metadata = os.fstat(descriptor)
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError(f"{label} is not a regular file")
        if metadata.st_mode & 0o222:
            raise ValueError(f"{label} is not read-only")
        with os.fdopen(descriptor, "rb", closefd=False) as stream:
            return stream.read()
    finally:
        os.close(descriptor)


def _open_read_only_directory_at(directory_fd: int, name: str, *, label: str) -> int:
    if not name or "/" in name:
        raise ValueError(f"{label} has an invalid fixed-layout name")
    flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(name, flags, dir_fd=directory_fd)
    metadata = os.fstat(descriptor)
    if metadata.st_mode & 0o222:
        os.close(descriptor)
        raise ValueError(f"{label} is not read-only")
    return descriptor


def _require_exact_directory_entries(
    directory_fd: int, expected: set[str], *, label: str
) -> None:
    entries = set(os.listdir(directory_fd))
    if entries != expected:
        raise ValueError(f"{label} entries do not match the fixed layout")


def _require_exact_layout_path(record: dict[str, object], key: str, expected: Path) -> None:
    if _text(record, key) != str(expected):
        raise ValueError(f"Clean-candidate {key} does not match the fixed layout")


def _read_submodule_archives(
    source: dict[str, object], *, candidate_dir: Path, candidate_fd: int
) -> tuple[list[bytes], list[bytes]]:
    records = source.get("submodules")
    if not isinstance(records, list):
        raise ValueError("Clean-candidate submodules must be a list")
    if not records:
        return [], []
    submodules_fd = _open_read_only_directory_at(
        candidate_fd, "submodules", label="Clean-candidate submodule directory"
    )
    try:
        archives = []
        commits = []
        names = set()
        for index, record in enumerate(records):
            if not isinstance(record, dict):
                raise ValueError("Clean-candidate submodule attestation must be an object")
            name = f"{index:04d}.tar"
            commit_name = f"{index:04d}.commit"
            names.add(name)
            names.add(commit_name)
            _require_exact_layout_path(
                record, "archive_path", candidate_dir / "submodules" / name
            )
            archives.append(
                _read_regular_file_at(
                    submodules_fd, name, label=f"Clean-candidate submodule archive {index}"
                )
            )
            _require_exact_layout_path(
                record, "commit_path", candidate_dir / "submodules" / commit_name
            )
            commits.append(
                _read_regular_file_at(
                    submodules_fd,
                    commit_name,
                    label=f"Clean-candidate submodule commit object {index}",
                )
            )
        _require_exact_directory_entries(
            submodules_fd, names, label="Clean-candidate submodule directory"
        )
        return archives, commits
    finally:
        os.close(submodules_fd)


def _q011_generator_bytes(source_archive: bytes) -> bytes:
    try:
        with tarfile.open(fileobj=io.BytesIO(source_archive), mode="r:") as archive:
            members = [
                member
                for member in archive.getmembers()
                if member.name == Q011_REPAIRED_GENERATOR_PATH
            ]
            if len(members) != 1 or not members[0].isfile():
                raise ValueError("Q011 source archive lacks one regular repaired generator")
            stream = archive.extractfile(members[0])
            if stream is None:
                raise ValueError("Q011 repaired generator source cannot be read")
            return stream.read()
    except tarfile.TarError as error:
        raise ValueError("Q011 source archive is not a readable tar file") from error


def _verify_q011_repaired_clean_candidate(
    manifest: dict[str, object],
    candidate: dict[str, object],
    *,
    candidate_sha256: str,
    source_archive: bytes,
    executable_sha256: str,
) -> None:
    if _q011_selected_pressure(manifest) is None:
        return
    source = _mapping(candidate, "source")
    git_commit = _text(source, "git_commit")
    if (
        git_commit == Q011_FAILED_V1_GIT_COMMIT
        or candidate_sha256 == Q011_FAILED_V1_CLEAN_CANDIDATE_MANIFEST_SHA256
        or executable_sha256 == Q011_FAILED_V1_EXECUTABLE_SHA256
    ):
        raise ValueError("Q011 pressure retry rejects the failed v1 clean candidate")
    if sha256_bytes(_q011_generator_bytes(source_archive)) != Q011_REPAIRED_GENERATOR_SHA256:
        raise ValueError("Q011 pressure retry requires the frozen repaired generator source")


def _verify_clean_candidate(
    manifest: dict[str, object],
    policy: dict[str, object],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> str:
    freeze_policy = _mapping(policy, "science_submission_freeze")
    if freeze_policy.get("status") != AUTHORIZED_CLEAN_CANDIDATE_FREEZE:
        raise ValueError("Registered science is blocked pending an authorized clean freeze")
    candidate_root = authorized_pic_root.resolve() / "clean_candidates"
    candidate_path = require_canonical_path_below(
        Path(str(manifest["clean_candidate_manifest_path"])), candidate_root
    )
    if candidate_path.name != "clean_candidate_manifest.json":
        raise ValueError("Clean-candidate manifest has an invalid path")
    if candidate_path.parent.parent != candidate_root:
        raise ValueError("Clean-candidate manifest must be one level below clean_candidates")
    policy_candidate_path = require_canonical_path_below(
        Path(str(freeze_policy.get("manifest_path", ""))), candidate_root
    )
    if candidate_path != policy_candidate_path:
        raise ValueError("Science manifest does not reference the policy-authorized freeze")
    candidate_sha256 = manifest["clean_candidate_manifest_sha256"]
    if (
        not isinstance(candidate_sha256, str)
        or re.fullmatch(r"[0-9a-f]{64}", candidate_sha256) is None
    ):
        raise ValueError("Science manifest clean-candidate digest is malformed")
    if candidate_sha256 != freeze_policy.get("manifest_sha256"):
        raise ValueError("Science manifest digest is not authorized by storage policy")

    snapshot = record_for_role(manifest, "clean-candidate-manifest")
    snapshot_source = require_canonical_path_below(
        Path(str(snapshot.get("source_path", ""))), candidate_root
    )
    if snapshot_source != candidate_path:
        raise ValueError("Clean-candidate snapshot source path mismatch")
    if (
        snapshot.get("source_sha256") != candidate_sha256
        or snapshot.get("sha256") != candidate_sha256
    ):
        raise ValueError("Clean-candidate snapshot digest mismatch")

    root_flags = os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    candidate_root_fd = os.open(candidate_root, root_flags)
    try:
        candidate_fd = _open_read_only_directory_at(
            candidate_root_fd,
            candidate_path.parent.name,
            label="Clean-candidate directory",
        )
        try:
            candidate_bytes = _read_regular_file_at(
                candidate_fd,
                "clean_candidate_manifest.json",
                label="Clean-candidate manifest",
            )
            if sha256_bytes(candidate_bytes) != candidate_sha256:
                raise ValueError("Policy-authorized clean-candidate manifest checksum mismatch")
            candidate = read_json_bytes(candidate_bytes, label="Clean-candidate manifest")
            freeze_id = _text(candidate, "freeze_id")
            uuid.UUID(freeze_id)
            if candidate_path.parent.name != freeze_id:
                raise ValueError("Clean-candidate path does not match its freeze ID")
            utc_datetime(candidate.get("created_utc"), field="clean_candidate.created_utc")
            source = _mapping(candidate, "source")
            build = _mapping(candidate, "build")
            _require_exact_layout_path(source, "archive_path", candidate_path.parent / "source.tar")
            _require_exact_layout_path(source, "commit_path", candidate_path.parent / "source.commit")
            _require_exact_layout_path(
                build, "profile_path", candidate_path.parent / "build_profile.json"
            )
            _require_exact_layout_path(
                build,
                "profile_receipt_path",
                candidate_path.parent / "profile_receipt.json",
            )
            _require_exact_layout_path(
                build, "executable_path", candidate_path.parent / "athena"
            )
            source_archive = _read_regular_file_at(
                candidate_fd, "source.tar", label="Clean-candidate source archive"
            )
            source_commit = _read_regular_file_at(
                candidate_fd, "source.commit", label="Clean-candidate source commit object"
            )
            submodule_archives, submodule_commits = _read_submodule_archives(
                source, candidate_dir=candidate_path.parent, candidate_fd=candidate_fd
            )
            profile_bytes = _read_regular_file_at(
                candidate_fd, "build_profile.json", label="Clean-candidate build profile"
            )
            receipt_bytes = _read_regular_file_at(
                candidate_fd,
                "profile_receipt.json",
                label="Clean-candidate build-profile receipt",
            )
            executable_bytes = _read_regular_file_at(
                candidate_fd, "athena", label="Clean-candidate executable"
            )
            provenance_fd = _open_read_only_directory_at(
                candidate_fd, "build_provenance", label="Frozen build provenance directory"
            )
            try:
                _require_exact_directory_entries(
                    provenance_fd,
                    set(BUILD_PROVENANCE_FILENAMES.values()),
                    label="Frozen build provenance directory",
                )
                build_provenance = {
                    label: _read_regular_file_at(
                        provenance_fd,
                        filename,
                        label=f"Frozen build provenance {label}",
                    )
                    for label, filename in BUILD_PROVENANCE_FILENAMES.items()
                }
            finally:
                os.close(provenance_fd)
            entries = {
                "athena",
                "build_provenance",
                "build_profile.json",
                "clean_candidate_manifest.json",
                "profile_receipt.json",
                "source.commit",
                "source.tar",
            }
            if source.get("submodules"):
                entries.add("submodules")
            _require_exact_directory_entries(
                candidate_fd, entries, label="Clean-candidate directory"
            )
        finally:
            os.close(candidate_fd)
    finally:
        os.close(candidate_root_fd)
    executable_sha256 = sha256_bytes(executable_bytes)
    profile_submodules = validate_clean_candidate_bundle(
        candidate,
        source_archive=source_archive,
        source_commit=source_commit,
        submodule_archives=submodule_archives,
        submodule_commits=submodule_commits,
        build_profile=profile_bytes,
        build_profile_receipt=receipt_bytes,
        build_provenance=build_provenance,
        executable_sha256=executable_sha256,
    )
    _verify_q011_repaired_clean_candidate(
        manifest,
        candidate,
        candidate_sha256=candidate_sha256,
        source_archive=source_archive,
        executable_sha256=executable_sha256,
    )
    receipt = read_json_bytes(
        receipt_bytes, label="clean-candidate build-profile receipt"
    )
    receipt_control_plane_version = _text(receipt, "control_plane_version")
    if receipt_control_plane_version != freeze_policy.get(
        "build_profile_control_plane_version"
    ):
        raise ValueError(
            "Clean-candidate build-profile receipt belongs to an unauthorized "
            "control-plane version"
        )
    for root in [authorized_pic_root, authorized_project_home_root]:
        verify_historical_installed_control_plane(
            root / "control_plane" / receipt_control_plane_version,
            authorized_pic_root=root,
        )
    git_commit = _text(source, "git_commit")
    if manifest.get("git_commit") != git_commit:
        raise ValueError("Science manifest Git commit differs from clean candidate")
    executable = candidate_path.parent / "athena"
    executable_snapshot = record_for_role(manifest, "executable")
    if Path(str(executable_snapshot.get("source_path", ""))) != executable:
        raise ValueError("Science executable is not the clean-candidate executable")
    if (
        executable_snapshot.get("source_sha256") != executable_sha256
        or executable_snapshot.get("sha256") != executable_sha256
    ):
        raise ValueError("Science executable digest differs from clean candidate")
    return candidate_sha256


def _registered_science_authorization(
    manifest: dict[str, object],
    policy: dict[str, object],
    directives: dict[str, str],
    *,
    candidate_sha256: str,
) -> tuple[str, int]:
    identifier = manifest["registered_science_authorization_id"]
    if not isinstance(identifier, str):
        raise ValueError("Registered science authorization ID is malformed")
    records = policy.get("registered_science_slices")
    if not isinstance(records, list):
        raise ValueError("Storage policy does not carry registered-science slices")
    matches = [
        record for record in records
        if isinstance(record, dict) and record.get("authorization_id") == identifier
    ]
    if len(matches) != 1:
        raise ValueError("Registered science authorization ID is not active")
    authorization = matches[0]
    expected_fields = {
        "campaign": manifest.get("campaign"),
        "test_id": manifest.get("test_id"),
        "evidence_class": manifest.get("evidence_class"),
        "physical_mode": manifest.get("physical_mode"),
        "selected_qos": manifest.get("selected_qos"),
        "registered_short_nonproduction": manifest.get("registered_short_nonproduction"),
        "clean_candidate_manifest_sha256": candidate_sha256,
    }
    for key, expected in expected_fields.items():
        if authorization.get(key) != expected:
            raise ValueError(f"Registered science authorization {key} differs")
    if authorization.get("runtime_profile") != "frontier_minimum_supported":
        raise ValueError("Registered science runtime profile is not authorized")
    if int(directives["nodes"]) > int(authorization["maximum_nodes"]):
        raise ValueError("Registered science exceeds its authorized node ceiling")
    if _walltime_seconds(directives["time"]) > int(
        authorization["maximum_walltime_seconds"]
    ):
        raise ValueError("Registered science exceeds its authorized walltime ceiling")
    role_to_digest = {
        "job-script": "job_script_sha256",
        "input-deck": "input_deck_sha256",
        "environment-profile": "environment_profile_sha256",
        "executable": "executable_sha256",
    }
    for role, digest_key in role_to_digest.items():
        if record_for_role(manifest, role).get("sha256") != authorization.get(digest_key):
            raise ValueError(f"Registered-science {role} is not policy authorized")
    analysis_records = sorted(
        (
            record for record in manifest["snapshot_files"]
            if isinstance(record, dict)
            and str(record.get("role", "")).startswith("analysis-script-")
        ),
        key=lambda record: str(record["role"]),
    )
    if [record.get("sha256") for record in analysis_records] != authorization.get(
        "analysis_script_sha256"
    ):
        raise ValueError("Registered-science analysis scripts are not policy authorized")
    if launch_contract_sha256(manifest.get("launch_contract")) != authorization.get(
        "launch_contract_sha256"
    ):
        raise ValueError("Registered-science launch contract is not policy authorized")
    return identifier, int(authorization["maximum_attempts"])


def _q011_snapshot_analyzer(
    prior_manifest: dict[str, object], prior_manifest_path: Path
) -> object:
    if os.environ.get("PIC_F1_ANALYSIS_HELPER_FD") is not None:
        raise ValueError("Q011 descriptor verification rejects inherited helper overrides")
    analyzer = record_for_role(prior_manifest, "analysis-script-000")
    helper = record_for_role(prior_manifest, "analysis-script-001")
    snapshot_root = prior_manifest_path.parent / "snapshot"
    analyzer_path = require_canonical_path_below(
        Path(str(analyzer["path"])), snapshot_root
    )
    helper_path = require_canonical_path_below(Path(str(helper["path"])), snapshot_root)
    expected_analysis_entries = {
        "000-analyze_q011_section54_pressure_pilot_case.py",
        "frontier_f1_structured_artifacts.py",
    }
    if (
        analyzer_path
        != snapshot_root / "analysis" / "000-analyze_q011_section54_pressure_pilot_case.py"
        or helper_path
        != snapshot_root / "analysis" / "frontier_f1_structured_artifacts.py"
        or {entry.name for entry in analyzer_path.parent.iterdir()}
        != expected_analysis_entries
        or analyzer.get("sha256") != Q011_RAW_ANALYZER_SHA256
        or helper.get("sha256") != Q011_STRUCTURED_HELPER_SHA256
    ):
        raise ValueError("Q011 predecessor analyzer snapshot differs from the registered bytes")
    require_read_only(analyzer_path)
    require_read_only(helper_path)
    analyzer_descriptor, analyzer_bytes = _q011_open_read_only_source(
        analyzer_path, label="Q011 predecessor analyzer snapshot"
    )
    os.close(analyzer_descriptor)
    if sha256_bytes(analyzer_bytes) != Q011_RAW_ANALYZER_SHA256:
        raise ValueError("Q011 predecessor analyzer snapshot changed while opening source")
    helper_descriptor, helper_bytes = _q011_open_read_only_source(
        helper_path, label="Q011 predecessor helper snapshot"
    )
    try:
        if (
            sha256_bytes(helper_bytes) != Q011_STRUCTURED_HELPER_SHA256
            or os.environ.get("PIC_F1_ANALYSIS_HELPER_FD") is not None
        ):
            raise ValueError("Q011 predecessor helper snapshot cannot be pinned")
        module = types.ModuleType(
            f"_q011_snapshot_analyzer_{hashlib.sha256(str(analyzer_path).encode()).hexdigest()}"
        )
        module.__file__ = str(analyzer_path)
        os.environ["PIC_F1_ANALYSIS_HELPER_FD"] = str(helper_descriptor)
        code = compile(analyzer_bytes, str(analyzer_path), "exec", dont_inherit=True)
        exec(code, module.__dict__)
    finally:
        os.environ.pop("PIC_F1_ANALYSIS_HELPER_FD", None)
        os.close(helper_descriptor)
    return module


def _verify_q011_prior_case_closures(
    manifest: dict[str, object],
    records: list[dict[str, object]],
    *,
    authorized_pic_root: Path,
) -> None:
    selected = _q011_selected_pressure(manifest)
    if selected is None:
        return
    selected_index, _, campaign, test_id = selected
    if manifest.get("campaign") != campaign or manifest.get("test_id") != test_id:
        raise ValueError("Q011 pressure authorization does not match its exact campaign")
    closures = manifest.get("prior_case_closures")
    if not isinstance(closures, list) or len(closures) != selected_index:
        raise ValueError("Q011 pressure case requires exact ordered predecessor closures")
    current_executable_sha256 = record_for_role(manifest, "executable")["sha256"]
    manifest_root = authorized_pic_root.resolve() / "manifests"
    for closure, expected in zip(closures, Q011_PRESSURE_CASES[:selected_index]):
        case_id, expected_authorization, expected_campaign, expected_test_id = expected
        if not isinstance(closure, dict) or set(closure) != {
            "case_id",
            "submission_id",
            "artifact_dir",
            "descriptor_path",
            "descriptor_sha256",
            "reconciliation_event_sha256",
        }:
            raise ValueError("Q011 predecessor closure schema is malformed")
        if closure.get("case_id") != case_id:
            raise ValueError("Q011 predecessor closures are not in exact preregistered order")
        submission_id = str(closure.get("submission_id", ""))
        try:
            parsed_submission_id = uuid.UUID(submission_id)
        except ValueError as error:
            raise ValueError("Q011 predecessor submission ID is malformed") from error
        if str(parsed_submission_id) != submission_id:
            raise ValueError("Q011 predecessor submission ID is not canonical")
        artifact_dir = (
            authorized_pic_root.resolve() / "runs" / expected_campaign / submission_id
        )
        descriptor_path = artifact_dir / "analysis" / "analysis.json"
        if (
            closure.get("artifact_dir") != str(artifact_dir)
            or closure.get("descriptor_path") != str(descriptor_path)
            or re.fullmatch(r"[0-9a-f]{64}", str(closure.get("descriptor_sha256", "")))
            is None
            or re.fullmatch(
                r"[0-9a-f]{64}", str(closure.get("reconciliation_event_sha256", ""))
            )
            is None
        ):
            raise ValueError("Q011 predecessor closure path or digest binding is malformed")
        matches = [
            record
            for record in records
            if record.get("event_sha256") == closure["reconciliation_event_sha256"]
            and record.get("event_type") == "reconciliation"
            and record.get("submission_scope") == REGISTERED_SCIENCE_SCOPE
            and record.get("submission_id") == submission_id
            and record.get("campaign") == expected_campaign
            and record.get("test_id") == expected_test_id
            and record.get("registered_science_authorization_id") == expected_authorization
            and record.get("artifact_dir") == str(artifact_dir)
            and record.get("clean_candidate_manifest_sha256")
            == manifest["clean_candidate_manifest_sha256"]
            and record.get("git_commit") == manifest["git_commit"]
            and record.get("executable_sha256") == current_executable_sha256
            and record.get("reconciled") is True
            and record.get("state") == "COMPLETED"
        ]
        if len(matches) != 1:
            raise ValueError("Q011 predecessor lacks one matching completed reconciliation")
        prior_manifest_path = require_canonical_path_below(
            Path(str(matches[0].get("manifest_path", ""))), manifest_root
        )
        prior_manifest_descriptor, prior_manifest_bytes = _q011_open_read_only_source(
            prior_manifest_path, label="Q011 predecessor manifest"
        )
        os.close(prior_manifest_descriptor)
        if sha256_bytes(prior_manifest_bytes) != matches[0].get("manifest_sha256"):
            raise ValueError("Q011 predecessor manifest digest drifted")
        prior_manifest = read_json_bytes(
            prior_manifest_bytes, label="Q011 predecessor manifest"
        )
        if (
            prior_manifest.get("clean_candidate_manifest_sha256")
            != manifest["clean_candidate_manifest_sha256"]
            or prior_manifest.get("git_commit") != manifest["git_commit"]
            or record_for_role(prior_manifest, "executable").get("sha256")
            != current_executable_sha256
        ):
            raise ValueError("Q011 predecessor belongs to another clean candidate")
        module = _q011_snapshot_analyzer(prior_manifest, prior_manifest_path)
        with module.StructuredArtifactTree(artifact_dir) as tree:
            module.verify_published_case_descriptor(
                tree, case_id, str(closure["descriptor_sha256"])
            )


def _check_submission_scope(
    manifest: dict[str, object],
    policy: dict[str, object],
    directives: dict[str, str],
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
) -> tuple[str, str, int]:
    if manifest.get("submission_scope") == REGISTERED_SCIENCE_SCOPE:
        candidate_sha256 = _verify_clean_candidate(
            manifest,
            policy,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        identifier, maximum_attempts = _registered_science_authorization(
            manifest, policy, directives, candidate_sha256=candidate_sha256
        )
        return candidate_sha256, identifier, maximum_attempts
    for key, expected in ADMISSION_SMOKE_FIELDS.items():
        value = manifest.get(key)
        if type(value) is not type(expected) or value != expected:
            raise ValueError(f"Admission-smoke exemption requires {key}={expected!r}")
    if manifest.get("submission_scope") != FRONTIER_ADMISSION_SMOKE_SCOPE:
        raise ValueError("Submission scope is not authorized")
    authorization = _mapping(policy, "frontier_admission_smoke")
    if authorization.get("status") != AUTHORIZED_ADMISSION_SMOKE_STATUS:
        raise ValueError("Admission-smoke exemption is not authorized")
    for key, expected in ADMISSION_SMOKE_FIELDS.items():
        value = authorization.get(key)
        if type(value) is not type(expected) or value != expected:
            raise ValueError(f"Admission-smoke policy requires {key}={expected!r}")
    if int(directives["nodes"]) > int(authorization["maximum_nodes"]):
        raise ValueError("Admission smoke exceeds its authorized node ceiling")
    if _walltime_seconds(directives["time"]) > int(
        authorization["maximum_walltime_seconds"]
    ):
        raise ValueError("Admission smoke exceeds its authorized walltime ceiling")
    role_to_digest = {
        "job-script": "job_script_sha256",
        "input-deck": "input_deck_sha256",
        "environment-profile": "environment_profile_sha256",
    }
    for role, digest_key in role_to_digest.items():
        if record_for_role(manifest, role).get("sha256") != authorization.get(digest_key):
            raise ValueError(f"Admission-smoke {role} is not policy authorized")
    if record_for_role(manifest, "executable").get("sha256") != authorization.get(
        "executable_sha256"
    ):
        raise ValueError("Admission-smoke executable is not policy authorized")
    analysis_records = sorted(
        (
            record for record in manifest["snapshot_files"]
            if isinstance(record, dict)
            and str(record.get("role", "")).startswith("analysis-script-")
        ),
        key=lambda record: str(record["role"]),
    )
    if [record.get("sha256") for record in analysis_records] != authorization.get(
        "analysis_script_sha256"
    ):
        raise ValueError("Admission-smoke analysis scripts are not policy authorized")
    if launch_contract_sha256(manifest.get("launch_contract")) != authorization.get(
        "launch_contract_sha256"
    ):
        raise ValueError("Admission-smoke launch contract is not policy authorized")
    return "", "", 0


def _check_timeout_margin(
    manifest: dict[str, object],
    walltime_seconds: int,
    *,
    now: datetime,
) -> None:
    record = record_for_role(manifest, "timeout-margin")
    margin = read_json(Path(str(record["path"])))
    if (
        not isinstance(margin, dict)
        or set(margin) != TIMEOUT_MARGIN_KEYS
        or not _strict_json_equal(margin, manifest.get("timeout_margin"))
    ):
        raise ValueError("Timeout-margin snapshot differs from manifest record")
    scheduler = margin.get("scheduler_walltime_seconds")
    athena = margin.get("athena_walltime_seconds")
    if type(scheduler) is not int or type(athena) is not int:
        raise ValueError("Timeout-margin walltimes must be exact integers")
    if scheduler != walltime_seconds:
        raise ValueError("Timeout-margin artifact does not match Slurm walltime")
    if not 0 < athena < scheduler:
        raise ValueError("Athena timeout must be positive and below Slurm walltime")
    profile = record_for_role(manifest, "environment-profile")
    if margin.get("environment_profile_sha256") != profile.get("sha256"):
        raise ValueError("Timeout-margin artifact does not match environment profile")
    measured = utc_datetime(margin.get("measured_utc"), field="measured_utc")
    expires = utc_datetime(margin.get("expires_utc"), field="expires_utc")
    if not measured <= now < expires:
        raise ValueError("Timeout-margin artifact is stale or not yet valid")


def _check_site_policy(manifest: dict[str, object], *, now: datetime) -> None:
    checked = utc_datetime(
        manifest.get("site_policy_checked_utc"), field="site_policy_checked_utc"
    )
    age = (now - checked).total_seconds()
    if age < 0 or age > SITE_POLICY_MAX_AGE_SECONDS:
        raise ValueError("Recorded site-policy check is stale or in the future")


def _check_qos(
    manifest: dict[str, object],
    directives: dict[str, str],
    queue_output: str,
    walltime_seconds: int,
) -> None:
    qos = directives["qos"]
    if directives["partition"] != "batch":
        raise ValueError("Frontier PIC submissions must use partition=batch")
    if qos not in {"debug", "normal"}:
        raise ValueError("Frontier PIC submissions require qos=debug or qos=normal")
    if manifest.get("selected_qos") != qos:
        raise ValueError("Manifest QoS does not match frozen Slurm script")
    reason = str(manifest.get("qos_selection_reason", ""))
    debug_jobs = [
        line for line in queue_output.splitlines()
        if len(line.split("|")) >= 3 and line.split("|")[2] == "debug"
    ]
    short = manifest.get("registered_short_nonproduction", False)
    if not isinstance(short, bool):
        raise ValueError("Manifest short-job registration must be a boolean")
    if qos == "debug":
        if not short or walltime_seconds > DEBUG_MAX_SECONDS:
            raise ValueError("Debug QoS is restricted to registered short jobs")
        if debug_jobs:
            raise ValueError("A user debug job is already queued or running")
        if reason != "debug_available":
            raise ValueError("Debug QoS requires reason=debug_available")
    elif reason == "debug_slot_occupied":
        if not debug_jobs:
            raise ValueError("Normal fallback claims an occupied debug slot")
    elif reason == "debug_ineligible_walltime":
        if walltime_seconds <= DEBUG_MAX_SECONDS:
            raise ValueError("Normal fallback claims an ineligible walltime")
    elif reason == "debug_ineligible_production":
        if short:
            raise ValueError("Normal fallback claims a production registration")
    elif reason != "normal_required_by_registered_campaign":
        raise ValueError("Normal QoS requires an approved fallback reason")


def _check_launch_resources(
    manifest: dict[str, object], directives: dict[str, str]
) -> None:
    contract = validate_launch_contract(manifest.get("launch_contract"))
    for action in contract["actions"]:
        resources = action["resources"]
        if int(resources["nodes"]) > int(directives["nodes"]):
            raise ValueError("Launch action exceeds the reserved Slurm node count")


def reserve(
    *,
    manifest_path: Path,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    node_hour_cap: float,
    pending_marker: Path | None = None,
    pre_submit_wrapper_attestation: Path | None = None,
    reservation_id: str | None = None,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_account: str = AUTHORIZED_ACCOUNT,
    now: datetime | None = None,
) -> dict[str, object]:
    current_time = now or datetime.now(timezone.utc)
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    manifest_path = require_canonical_path_below(
        manifest_path, authorized_pic_root.resolve() / "manifests"
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        version = str(inventory["version"])
        policy, policy_snapshot = require_storage_policy_unlock_snapshot(
            control_plane_version=version,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
            authorized_account=authorized_account,
        )
        if not 0 < node_hour_cap <= AUTHORIZED_NODE_HOUR_CAP:
            raise ValueError("Requested PIC node-hour cap exceeds the authorized limit")
        manifest = _verify_manifest(
            manifest_path,
            control_plane_dir=control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        artifact_dir = _require_run_artifact_dir(manifest)
        if artifact_dir.exists():
            raise ValueError(f"Launch artifact directory already exists: {artifact_dir}")
        directives = _directives(
            Path(str(record_for_role(manifest, "job-script")["path"]))
        )
        (
            clean_candidate_manifest_sha256,
            registered_science_authorization_id,
            registered_science_maximum_attempts,
        ) = _check_submission_scope(
            manifest,
            policy,
            directives,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        pre_submit_wrapper_binding = None
        if registered_science_authorization_id:
            if pre_submit_wrapper_attestation is None:
                raise ValueError(
                    "Registered science requires a sealed pre-submit-wrapper attestation"
                )
            pre_submit_wrapper_binding = validate_sealed_operator_attestation(
                pre_submit_wrapper_attestation,
                authorization_id=registered_science_authorization_id,
                phase="pre_submit_wrapper",
                control_plane_version=version,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
                now=current_time,
            )
        elif pre_submit_wrapper_attestation is not None:
            raise ValueError(
                "Admission smoke must not claim a pre-submit-wrapper attestation"
            )
        if directives["account"] != authorized_account:
            raise ValueError(f"Frontier PIC submissions require account={authorized_account}")
        _require_scheduler_output_path(directives["output"], authorized_pic_root)
        queue_output = _queue_output()
        queue_sha256 = hashlib.sha256(queue_output.encode("utf-8")).hexdigest()
        if queue_sha256 != manifest.get("queue_snapshot_sha256"):
            raise ValueError("Fresh queue output differs from frozen queue snapshot")
        if "pic-reservation=" in queue_output:
            raise ValueError("A PIC-tagged Frontier job is already queued or running")
        expected_marker = _pending_marker_path(authorized_pic_root)
        if pending_marker is not None:
            if pending_marker.resolve() != expected_marker:
                raise ValueError(f"Unauthorized pending marker: {pending_marker}")
        if expected_marker.exists():
            raise ValueError(f"Pending PIC submission marker exists: {expected_marker}")

        walltime_seconds = _walltime_seconds(directives["time"])
        nodes = int(directives["nodes"])
        if nodes <= 0 or walltime_seconds <= 0:
            raise ValueError("Requested nodes and walltime must be positive")
        _check_timeout_margin(manifest, walltime_seconds, now=current_time)
        _check_site_policy(manifest, now=current_time)
        _check_qos(manifest, directives, queue_output, walltime_seconds)
        _check_launch_resources(manifest, directives)

        records = _records_with_matching_mirror(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        if registered_science_authorization_id:
            _verify_q011_prior_case_closures(
                manifest, records, authorized_pic_root=authorized_pic_root
            )
            prior_attempts = sum(
                record.get("event_type") == "reservation"
                and record.get("registered_science_authorization_id")
                == registered_science_authorization_id
                for record in records
            )
            if prior_attempts >= registered_science_maximum_attempts:
                raise ValueError("Registered-science authorization attempt ceiling is exhausted")
        totals = accounting(records)
        if totals["currently_reserved_node_hours"]:
            raise ValueError("An unreconciled PIC reservation already exists")
        reserved = nodes * walltime_seconds / 3600.0
        projected = (
            totals["cumulative_consumed_node_hours"]
            + totals["currently_reserved_node_hours"]
            + reserved
        )
        if projected > node_hour_cap:
            raise ValueError("PIC node-hour cap would be exceeded")

        reservation = reservation_id or str(uuid.uuid4())
        uuid.UUID(reservation)
        manifest_sha256 = sha256(manifest_path)
        for attachment in [
            manifest_path.parent / "reservation_id.txt",
            manifest_path.parent / "manifest_sha256.txt",
        ]:
            require_not_symlink(attachment)
            if attachment.exists():
                raise ValueError(f"Reservation attachment already exists: {attachment}")
        event = {
            "event_type": "reservation",
            "reservation_id": reservation,
            "submission_id": manifest["submission_id"],
            "control_plane_version": version,
            "active_policy_sha256": policy_snapshot["active_policy_sha256"],
            "active_promotion_sha256": policy_snapshot["active_promotion_sha256"],
            "git_commit": manifest["git_commit"],
            "campaign": manifest["campaign"],
            "test_id": manifest["test_id"],
            "manifest_path": str(manifest_path),
            "submission_scope": manifest["submission_scope"],
            "registered_science_authorization_id": registered_science_authorization_id,
            "clean_candidate_manifest_sha256": clean_candidate_manifest_sha256,
            "partition": directives["partition"],
            "qos": directives["qos"],
            "qos_selection_reason": manifest["qos_selection_reason"],
            "queue_snapshot_sha256": manifest["queue_snapshot_sha256"],
            "manifest_sha256": manifest_sha256,
            "job_script_sha256": record_for_role(manifest, "job-script")["sha256"],
            "executable_sha256": record_for_role(manifest, "executable")["sha256"],
            "site_policy_checked_utc": manifest["site_policy_checked_utc"],
            "requested_nodes": nodes,
            "requested_walltime": directives["time"],
            "reserved_node_hours": reserved,
            "artifact_dir": manifest["artifact_dir"],
            "state": "reserved",
            "reconciled": False,
        }
        if pre_submit_wrapper_binding is not None:
            event.update(
                {
                    "pre_manifest_attestation_path": manifest[
                        "pre_manifest_attestation_path"
                    ],
                    "pre_manifest_attestation_sha256": manifest[
                        "pre_manifest_attestation_sha256"
                    ],
                    "pre_submit_wrapper_attestation_path": pre_submit_wrapper_binding[
                        "path"
                    ],
                    "pre_submit_wrapper_attestation_sha256": pre_submit_wrapper_binding[
                        "sha256"
                    ],
                }
            )
        _write_pending_marker(
            expected_marker,
            {
                "schema_version": 2,
                "state": "reservation_intent",
                "reservation_id": reservation,
                "submission_id": manifest["submission_id"],
                "manifest_path": str(manifest_path),
                "manifest_sha256": manifest_sha256,
                "control_plane_version": version,
            },
        )
        result = _append_locked(
            ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl, event
        )
        _write_reservation_attachments(
            manifest_path,
            reservation_id=reservation,
            manifest_sha256=manifest_sha256,
        )
        _write_pending_marker(
            expected_marker,
            {
                "schema_version": 2,
                "state": "reserved_not_submitted",
                "reservation_id": reservation,
                "submission_id": manifest["submission_id"],
                "manifest_path": str(manifest_path),
                "manifest_sha256": manifest_sha256,
                "control_plane_version": version,
            },
        )
        return result


def transition(
    *,
    reservation_id: str,
    event_type: str,
    state: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    job_id: str = "",
    notes: str = "",
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, object]:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        records = _records_with_matching_mirror(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        latest = latest_reservations(records).get(reservation_id)
        if latest is None:
            raise ValueError(f"Unknown reservation: {reservation_id}")
        if latest.get("state") != "reserved":
            raise ValueError(f"Reservation is not attachable/cancellable: {reservation_id}")
        if latest.get("control_plane_version") != inventory["version"]:
            raise ValueError("Reservation belongs to another control-plane version")
        if event_type == "job_id_attached":
            _require_reservation_policy_snapshot(
                latest,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            marker = _matching_pending_marker(
                _pending_marker_path(authorized_pic_root), reservation_id
            )
            if (
                marker is None
                or marker.get("state") != "submitted_not_attached"
                or marker.get("job_id") != job_id
            ):
                raise ValueError("Scheduler attachment does not match the pending marker")
            _require_current_reservation_marker(
                marker,
                latest,
                control_plane_version=str(inventory["version"]),
            )
            _verify_scheduler_job(job_id, reservation_id)
        elif event_type == "reservation_cancelled":
            marker = _matching_pending_marker(
                _pending_marker_path(authorized_pic_root), reservation_id
            )
            if marker is None or marker.get("state") != "reserved_not_submitted":
                raise ValueError(
                    "Cannot cancel a reservation after scheduler submission"
                )
            _require_current_reservation_marker(
                marker,
                latest,
                control_plane_version=str(inventory["version"]),
            )
        event = transition_payload(latest)
        event.update({"event_type": event_type, "state": state})
        if job_id:
            event["job_id"] = job_id
        if notes:
            event["notes"] = notes
        result = _append_locked(
            ledger_jsonl, ledger_csv, receipts_jsonl, mirror_jsonl, event
        )
        if event_type in {"job_id_attached", "reservation_cancelled"}:
            _clear_matching_pending_marker(
                _pending_marker_path(authorized_pic_root), reservation_id
            )
        return result


def mark_dispatch_started(
    *,
    reservation_id: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        records = _records_with_matching_mirror(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        reservation = latest_reservations(records).get(reservation_id)
        if reservation is None or reservation.get("state") != "reserved":
            raise ValueError("Only a live reserved reservation may begin scheduler dispatch")
        if reservation.get("control_plane_version") != inventory["version"]:
            raise ValueError("Reservation belongs to another control-plane version")
        _require_reservation_policy_snapshot(
            reservation,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        marker_path = _pending_marker_path(authorized_pic_root)
        marker = _matching_pending_marker(marker_path, reservation_id)
        if marker is None or marker.get("state") != "reserved_not_submitted":
            raise ValueError("Expected a reserved-not-submitted recovery marker")
        _require_current_reservation_marker(
            marker,
            reservation,
            control_plane_version=str(inventory["version"]),
        )
        marker["state"] = "scheduler_dispatch_started"
        _write_pending_marker(marker_path, marker)


def mark_submitted(
    *,
    reservation_id: str,
    job_id: str,
    ledger_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> None:
    expected_ledger = authorized_pic_root.resolve() / "ledger" / "node_hours.jsonl"
    if require_canonical_path_below(ledger_jsonl, authorized_pic_root) != expected_ledger:
        raise ValueError(f"Unauthorized ledger path: {ledger_jsonl}")
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    mirror_jsonl = authorized_project_home_root / "ledger" / "node_hours.jsonl"
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        marker_path = _pending_marker_path(authorized_pic_root)
        records = _records_with_matching_mirror(
            ledger_jsonl,
            authorized_pic_root / "ledger" / "mirror_receipts.jsonl",
            mirror_jsonl,
        )
        reservation = latest_reservations(records).get(reservation_id)
        if reservation is None or reservation.get("state") != "reserved":
            raise ValueError("Expected a live reserved scheduler submission")
        marker = _matching_pending_marker(marker_path, reservation_id)
        if marker is None:
            raise ValueError("Expected a pending scheduler-submission recovery marker")
        _require_current_reservation_marker(
            marker,
            reservation,
            control_plane_version=str(inventory["version"]),
        )
        if marker.get("state") == "submitted_not_attached":
            if marker.get("job_id") != job_id:
                raise ValueError("Scheduler job differs from pending attachment marker")
            return
        if marker.get("state") == "scheduler_dispatch_started":
            marker["state"] = "scheduler_job_id_received"
            marker["job_id"] = job_id
            _write_pending_marker(marker_path, marker)
        elif (
            marker.get("state") != "scheduler_job_id_received"
            or marker.get("job_id") != job_id
        ):
            raise ValueError("Expected a scheduler-dispatch-started recovery marker")
        _verify_scheduler_job(job_id, reservation_id)
        marker["state"] = "submitted_not_attached"
        _write_pending_marker(marker_path, marker)


def repair_reservation_attachments(
    *,
    reservation_id: str,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> str:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        inventory = _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        repair_mirrored_state_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            mirror_transport="filesystem_copy",
        )
        records = _records_with_matching_mirror(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        latest = latest_reservations(records).get(reservation_id)
        marker_path = _pending_marker_path(authorized_pic_root)
        marker = _matching_pending_marker(marker_path, reservation_id)
        if marker is None:
            raise ValueError("Missing pending marker for attachment repair")
        if latest is not None:
            if latest.get("control_plane_version") != inventory["version"]:
                raise ValueError("Reservation belongs to another control-plane version")
            _require_current_reservation_marker(
                marker,
                latest,
                control_plane_version=str(inventory["version"]),
            )
            _require_reservation_policy_snapshot(
                latest,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
        if (
            latest is not None
            and latest.get("event_type") == "job_id_attached"
            and latest.get("state") == "submitted"
        ):
            if any(
                field in latest
                for field in [
                    "terminal_recovery_handoff_path",
                    "terminal_recovery_handoff_sha256",
                    "terminal_recovery_mode",
                ]
            ):
                raise ValueError(
                    "Terminal-recovery attachment must resume through its successor reconciler"
                )
            _clear_matching_pending_marker(marker_path, reservation_id)
            return "cleared_completed_attachment_pending_marker"
        if (
            latest is not None
            and latest.get("event_type") == "reservation_cancelled"
            and latest.get("state") == "cancelled"
        ):
            _clear_matching_pending_marker(marker_path, reservation_id)
            return "cleared_completed_cancellation_pending_marker"
        if (
            latest is not None
            and latest.get("event_type") == "reconciliation"
            and latest.get("reconciled") is True
        ):
            _clear_matching_pending_marker(marker_path, reservation_id)
            return "cleared_completed_reconciliation_pending_marker"
        if latest is None:
            if (
                type(marker.get("schema_version")) is not int
                or marker.get("schema_version") != 2
                or marker.get("control_plane_version") != inventory["version"]
                or set(marker) != {
                    "schema_version",
                    "state",
                    "reservation_id",
                    "submission_id",
                    "manifest_path",
                    "manifest_sha256",
                    "control_plane_version",
                }
            ):
                raise ValueError("Reservation intent belongs to another control-plane version")
            if marker.get("state") != "reservation_intent":
                raise ValueError("Only an unappended reservation intent may be cleared")
            for name in ["reservation_id.txt", "manifest_sha256.txt"]:
                attachment = Path(str(marker["manifest_path"])).parent / name
                require_not_symlink(attachment)
                if attachment.exists():
                    raise ValueError("Unappended reservation intent unexpectedly has attachments")
            marker_path.unlink()
            fsync_directory(marker_path.parent)
            return "cleared_unappended_reservation_intent"
        if latest.get("state") != "reserved":
            raise ValueError("Expected one recoverable reserved reservation")
        if marker.get("state") not in {"reservation_intent", "reserved_not_submitted"}:
            raise ValueError(
                "Scheduler dispatch may have started; reviewed scheduler reconciliation is required"
            )
        manifest_path = Path(str(latest["manifest_path"]))
        manifest_path = require_canonical_path_below(
            manifest_path, authorized_pic_root.resolve() / "manifests"
        )
        if sha256(manifest_path) != latest.get("manifest_sha256"):
            raise ValueError("Cannot repair attachments for a changed manifest")
        _write_reservation_attachments(
            manifest_path,
            reservation_id=reservation_id,
            manifest_sha256=str(latest["manifest_sha256"]),
            allow_existing=True,
        )
        marker["state"] = "reserved_not_submitted"
        _write_pending_marker(marker_path, marker)
        return "repaired_reserved_submission"


def repair_ledger_mirror(
    *,
    ledger_jsonl: Path,
    ledger_csv: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> dict[str, int]:
    require_ledger_paths(
        ledger_jsonl,
        ledger_csv,
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    with ledger_lock(ledger_jsonl, mirror_jsonl):
        _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        return repair_mirrored_state_locked(
            ledger_jsonl,
            ledger_csv,
            receipts_jsonl,
            mirror_jsonl,
            mirror_transport="filesystem_copy",
        )


def reservation_bound_manifest(
    manifest_path: Path,
    reservation_id: str,
    *,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    executable_job_id: str | None = None,
    require_reserved: bool = False,
    _compute_node_read_only_snapshot: bool = False,
) -> tuple[dict[str, object], dict[str, object]]:
    require_ledger_paths(
        ledger_jsonl,
        authorized_pic_root / "ledger" / "node_hours.csv",
        receipts_jsonl,
        mirror_jsonl,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    manifest_path = require_canonical_path_below(
        manifest_path, authorized_pic_root.resolve() / "manifests"
    )
    if _compute_node_read_only_snapshot and executable_job_id is None:
        raise ValueError("Compute-node stable snapshot requires a scheduler job ID")
    _verify_installed_control_plane_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
    )
    # Compute-node startup cannot rely on cross-root flock support. Pin and
    # recheck the complete authoritative mirrored state across verification
    # instead of entering the mutation hierarchy.
    state = (
        validated_read_only_mirrored_state_snapshot(
            ledger_jsonl,
            receipts_jsonl,
            mirror_jsonl,
            ledger_root=authorized_pic_root,
            receipts_root=authorized_pic_root,
            mirror_root=authorized_project_home_root,
        )
        if _compute_node_read_only_snapshot
        else ledger_lock(ledger_jsonl, mirror_jsonl)
    )
    with state as snapshot_records:
        _verify_installed_control_plane_pair(
            control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        records = (
            snapshot_records
            if _compute_node_read_only_snapshot
            else _records_with_matching_mirror(ledger_jsonl, receipts_jsonl, mirror_jsonl)
        )
        reservation = latest_reservations(records).get(reservation_id)
        if reservation is None:
            raise ValueError(f"Unknown reservation: {reservation_id}")
        if require_reserved and reservation.get("state") != "reserved":
            raise ValueError("Reservation is not available for pre-dispatch lookup")
        if executable_job_id is not None:
            _require_reservation_policy_snapshot(
                reservation,
                authorized_pic_root=authorized_pic_root,
                authorized_project_home_root=authorized_project_home_root,
            )
            _verify_scheduler_job(executable_job_id, reservation_id)
            if reservation.get("state") == "submitted":
                if reservation.get("job_id") != executable_job_id:
                    raise ValueError("Submitted reservation belongs to another Slurm job")
            else:
                raise ValueError("Reservation is not executable")
        manifest = _verify_manifest(
            manifest_path,
            control_plane_dir=control_plane_dir,
            authorized_pic_root=authorized_pic_root,
            authorized_project_home_root=authorized_project_home_root,
        )
        if Path(str(reservation.get("manifest_path", ""))) != manifest_path:
            raise ValueError("Manifest path differs from the reservation ledger")
        if reservation.get("submission_id") != manifest.get("submission_id"):
            raise ValueError("Submission ID differs from the reservation ledger")
        manifest_digest = str(reservation.get("manifest_sha256", ""))
        if sha256(manifest_path) != manifest_digest:
            raise ValueError("Pre-submit manifest checksum differs from the reservation ledger")
        for role, field in [
            ("job-script", "job_script_sha256"),
            ("executable", "executable_sha256"),
        ]:
            if record_for_role(manifest, role).get("sha256") != reservation.get(field):
                raise ValueError(f"{role} digest differs from the reservation ledger")
        attachment = manifest_path.parent / "reservation_id.txt"
        require_not_symlink(attachment)
        if not attachment.is_file():
            raise ValueError("Missing immutable reservation attachment")
        if attachment.read_bytes() != (reservation_id + "\n").encode("utf-8"):
            raise ValueError("Reservation ID does not match immutable attachment")
        require_read_only(attachment)
        manifest_digest_path = manifest_path.parent / "manifest_sha256.txt"
        require_not_symlink(manifest_digest_path)
        if not manifest_digest_path.is_file():
            raise ValueError("Missing immutable manifest checksum attachment")
        if manifest_digest_path.read_bytes() != (manifest_digest + "\n").encode("utf-8"):
            raise ValueError("Manifest checksum attachment differs from reservation ledger")
        require_read_only(manifest_digest_path)
        return manifest, reservation


def executable_reservation_bound_manifest(
    manifest_path: Path,
    reservation_id: str,
    *,
    ledger_jsonl: Path,
    receipts_jsonl: Path,
    mirror_jsonl: Path,
    executable_job_id: str,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
) -> tuple[dict[str, object], dict[str, object]]:
    if re.fullmatch(r"[0-9]+", executable_job_id) is None:
        raise ValueError("Executable snapshot lookup requires a numeric scheduler job ID")
    return reservation_bound_manifest(
        manifest_path,
        reservation_id,
        ledger_jsonl=ledger_jsonl,
        receipts_jsonl=receipts_jsonl,
        mirror_jsonl=mirror_jsonl,
        executable_job_id=executable_job_id,
        control_plane_dir=control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        _compute_node_read_only_snapshot=True,
    )


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--ledger-jsonl", required=True, type=Path)
    parser.add_argument("--ledger-csv", required=True, type=Path)
    parser.add_argument("--receipts-jsonl", required=True, type=Path)
    parser.add_argument("--mirror-jsonl", required=True, type=Path)


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("verify-control-plane")
    reserve_parser = subparsers.add_parser("reserve")
    _common(reserve_parser)
    reserve_parser.add_argument("--manifest", required=True, type=Path)
    reserve_parser.add_argument("--node-hour-cap", required=True, type=float)
    reserve_parser.add_argument("--pending-marker", type=Path)
    reserve_parser.add_argument("--pre-submit-wrapper-attestation", type=Path)

    attach_parser = subparsers.add_parser("attach-job-id")
    _common(attach_parser)
    attach_parser.add_argument("--reservation-id", required=True)
    attach_parser.add_argument("--job-id", required=True)

    cancel_parser = subparsers.add_parser("cancel-reservation")
    _common(cancel_parser)
    cancel_parser.add_argument("--reservation-id", required=True)
    cancel_parser.add_argument("--notes", required=True)

    submission_parser = subparsers.add_parser("submission-id")
    _common(submission_parser)
    submission_parser.add_argument("--manifest", required=True, type=Path)
    submission_parser.add_argument("--reservation-id", required=True)

    snapshot_parser = subparsers.add_parser("snapshot-path")
    _common(snapshot_parser)
    snapshot_parser.add_argument("--manifest", required=True, type=Path)
    snapshot_parser.add_argument("--reservation-id", required=True)
    snapshot_parser.add_argument("--role", required=True)
    manifest_sha_parser = subparsers.add_parser("manifest-sha256")
    _common(manifest_sha_parser)
    manifest_sha_parser.add_argument("--manifest", required=True, type=Path)
    manifest_sha_parser.add_argument("--reservation-id", required=True)
    snapshot_sha_parser = subparsers.add_parser("snapshot-sha256")
    _common(snapshot_sha_parser)
    snapshot_sha_parser.add_argument("--manifest", required=True, type=Path)
    snapshot_sha_parser.add_argument("--reservation-id", required=True)
    snapshot_sha_parser.add_argument("--role", required=True)
    directive_parser = subparsers.add_parser("directive")
    _common(directive_parser)
    directive_parser.add_argument("--manifest", required=True, type=Path)
    directive_parser.add_argument("--reservation-id", required=True)
    directive_parser.add_argument(
        "--key", required=True, choices=["account", "partition", "qos", "nodes", "time", "output"]
    )
    submitted_parser = subparsers.add_parser("mark-submitted")
    _common(submitted_parser)
    submitted_parser.add_argument("--reservation-id", required=True)
    submitted_parser.add_argument("--job-id", required=True)
    dispatch_parser = subparsers.add_parser("mark-dispatch-started")
    _common(dispatch_parser)
    dispatch_parser.add_argument("--reservation-id", required=True)
    repair_parser = subparsers.add_parser("repair-reservation-attachments")
    _common(repair_parser)
    repair_parser.add_argument("--reservation-id", required=True)
    mirror_repair_parser = subparsers.add_parser("repair-ledger-mirror")
    _common(mirror_repair_parser)

    args = parser.parse_args()
    if args.command == "verify-control-plane":
        inventory = verify_installed_control_plane(SCRIPT_DIR)
        print(inventory["version"])
        return
    if args.command in {
        "submission-id",
        "snapshot-path",
        "snapshot-sha256",
        "directive",
        "manifest-sha256",
    }:
        manifest, reservation = reservation_bound_manifest(
            args.manifest,
            args.reservation_id,
            ledger_jsonl=args.ledger_jsonl,
            receipts_jsonl=args.receipts_jsonl,
            mirror_jsonl=args.mirror_jsonl,
            require_reserved=True,
        )
        if args.command == "submission-id":
            print(reservation["submission_id"])
        elif args.command == "manifest-sha256":
            print(reservation["manifest_sha256"])
        elif args.command == "snapshot-path":
            print(record_for_role(manifest, args.role)["path"])
        elif args.command == "snapshot-sha256":
            fields = {
                "job-script": "job_script_sha256",
                "executable": "executable_sha256",
            }
            if args.role not in fields:
                raise ValueError("Only ledger-bound launch snapshot digests may be queried")
            print(reservation[fields[args.role]])
        else:
            script = Path(str(record_for_role(manifest, "job-script")["path"]))
            print(_directives(script)[args.key])
        return
    common = {
        "ledger_jsonl": args.ledger_jsonl,
        "ledger_csv": args.ledger_csv,
        "receipts_jsonl": args.receipts_jsonl,
        "mirror_jsonl": args.mirror_jsonl,
    }
    if args.command == "reserve":
        event = reserve(
            manifest_path=args.manifest,
            node_hour_cap=args.node_hour_cap,
            pending_marker=args.pending_marker,
            pre_submit_wrapper_attestation=args.pre_submit_wrapper_attestation,
            **common,
        )
    elif args.command == "mark-submitted":
        mark_submitted(
            reservation_id=args.reservation_id,
            job_id=args.job_id,
            ledger_jsonl=args.ledger_jsonl,
        )
        print(args.reservation_id)
        return
    elif args.command == "mark-dispatch-started":
        mark_dispatch_started(
            reservation_id=args.reservation_id,
            **common,
        )
        print(args.reservation_id)
        return
    elif args.command == "repair-reservation-attachments":
        result = repair_reservation_attachments(
            reservation_id=args.reservation_id,
            ledger_jsonl=args.ledger_jsonl,
            ledger_csv=args.ledger_csv,
            receipts_jsonl=args.receipts_jsonl,
            mirror_jsonl=args.mirror_jsonl,
        )
        print(result)
        return
    elif args.command == "repair-ledger-mirror":
        result = repair_ledger_mirror(**common)
        print(result)
        return
    elif args.command == "attach-job-id":
        event = transition(
            reservation_id=args.reservation_id,
            job_id=args.job_id,
            event_type="job_id_attached",
            state="submitted",
            **common,
        )
    else:
        event = transition(
            reservation_id=args.reservation_id,
            notes=args.notes,
            event_type="reservation_cancelled",
            state="cancelled",
            **common,
        )
    print(event["reservation_id"])


if __name__ == "__main__":
    main()
