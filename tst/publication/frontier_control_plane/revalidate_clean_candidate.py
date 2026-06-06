#!/opt/cray/pe/python/3.11.7/bin/python3 -I
"""Revalidate one immutable Orion clean candidate without changing control state."""

from __future__ import annotations

import sys as _sys
if __name__ == "__main__" and "/control_plane/" in __file__ and not getattr(
    _sys, "_pic_control_plane_bootstrapped", False
):
    raise SystemExit("Run installed control-plane tools through run_control_plane.py")

import argparse
from collections.abc import Callable
import hashlib
import json
import os
from pathlib import Path
import re
import sys

from control_plane_common import AUTHORIZED_CLEAN_CANDIDATE_SOURCE_ROOT
from control_plane_common import AUTHORIZED_PIC_ROOT, AUTHORIZED_PROJECT_HOME_ROOT
from control_plane_common import read_clean_candidate_tree, read_json_bytes
from control_plane_common import validate_clean_candidate_bundle
from control_plane_common import verify_historical_installed_control_plane
from control_plane_common import verify_installed_control_plane


SCRIPT_DIR = Path(__file__).absolute().parent
ENTRYPOINT_NAME = "revalidate_clean_candidate.py"
RECORD_TYPE = "frontier_pic_clean_candidate_read_only_revalidation"


InstalledVerifier = Callable[..., dict[str, object]]
CandidateReader = Callable[..., dict[str, object]]
BundleValidator = Callable[..., list[dict[str, str]]]


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def _sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _lowercase_sha256(value: object, *, label: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{64}", value) is None:
        raise ValueError(f"{label} is not a lowercase SHA-256 digest")
    return value


def _lowercase_git_commit(value: object, *, label: str) -> str:
    if not isinstance(value, str) or re.fullmatch(r"[0-9a-f]{40}", value) is None:
        raise ValueError(f"{label} is not a full lowercase Git commit")
    return value


def _expected_git_commit_argument(value: str) -> str:
    try:
        return _lowercase_git_commit(value, label="Expected source Git commit")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def _mapping(record: dict[str, object], key: str, *, label: str) -> dict[str, object]:
    value = record.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"{label} has no {key} object")
    return value


def _text(record: dict[str, object], key: str, *, label: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} has no {key}")
    return value.strip()


def _digest(record: dict[str, object], key: str, *, label: str) -> str:
    return _lowercase_sha256(
        _text(record, key, label=label),
        label=f"{label} {key}",
    )


def _verify_current_installed_pair(
    control_plane_dir: Path,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    verify_installed: InstalledVerifier,
) -> dict[str, object]:
    inventory = verify_installed(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
    )
    version = _digest(inventory, "version", label="Installed control-plane inventory")
    project_home_inventory = verify_installed(
        Path(os.path.abspath(authorized_project_home_root))
        / "control_plane"
        / version,
        authorized_pic_root=authorized_project_home_root,
    )
    if project_home_inventory != inventory:
        raise ValueError("Installed Orion and Project Home control-plane inventories differ")
    return inventory


def _verify_historical_installed_pair(
    control_plane_version: str,
    *,
    authorized_pic_root: Path,
    authorized_project_home_root: Path,
    verify_historical: InstalledVerifier,
) -> None:
    expected: dict[str, object] | None = None
    for root in [authorized_pic_root, authorized_project_home_root]:
        lexical_root = Path(os.path.abspath(root))
        inventory = verify_historical(
            lexical_root / "control_plane" / control_plane_version,
            authorized_pic_root=root,
        )
        if _digest(
            inventory,
            "version",
            label="Historical control-plane inventory",
        ) != control_plane_version:
            raise ValueError(
                "Historical control-plane inventory differs from build receipt"
            )
        if expected is not None and inventory != expected:
            raise ValueError(
                "Historical Orion and Project Home control-plane inventories differ"
            )
        expected = inventory


def revalidate_clean_candidate(
    candidate_manifest_path: Path,
    *,
    expected_manifest_sha256: str,
    expected_git_commit: str | None = None,
    expected_receipt_control_plane_version: str | None = None,
    control_plane_dir: Path = SCRIPT_DIR,
    authorized_pic_root: Path = AUTHORIZED_PIC_ROOT,
    authorized_project_home_root: Path = AUTHORIZED_PROJECT_HOME_ROOT,
    authorized_source_root: Path = AUTHORIZED_CLEAN_CANDIDATE_SOURCE_ROOT,
    verify_installed: InstalledVerifier = verify_installed_control_plane,
    verify_historical: InstalledVerifier = verify_historical_installed_control_plane,
    read_candidate_tree: CandidateReader = read_clean_candidate_tree,
    validate_bundle: BundleValidator = validate_clean_candidate_bundle,
) -> dict[str, object]:
    """Return one deterministic report after a read-only clean-candidate revalidation."""
    expected_candidate_sha256 = _lowercase_sha256(
        expected_manifest_sha256,
        label="Expected clean-candidate manifest SHA-256",
    )
    expected_source_git_commit = (
        _lowercase_git_commit(
            expected_git_commit,
            label="Expected source Git commit",
        )
        if expected_git_commit is not None
        else None
    )
    expected_receipt_version = (
        _lowercase_sha256(
            expected_receipt_control_plane_version,
            label="Expected build-receipt control-plane version",
        )
        if expected_receipt_control_plane_version is not None
        else None
    )
    current_inventory = _verify_current_installed_pair(
        control_plane_dir,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        verify_installed=verify_installed,
    )
    tree = read_candidate_tree(
        candidate_manifest_path,
        authorized_pic_root=authorized_pic_root,
    )
    candidate = tree.get("candidate")
    if not isinstance(candidate, dict):
        raise ValueError("Clean-candidate reader did not return one manifest object")
    candidate_bytes = tree.get("candidate_manifest_bytes")
    build_profile_receipt = tree.get("build_profile_receipt")
    executable = tree.get("executable")
    if (
        not isinstance(candidate_bytes, bytes)
        or not isinstance(build_profile_receipt, bytes)
        or not isinstance(executable, bytes)
    ):
        raise ValueError("Clean-candidate reader returned malformed captured bytes")
    actual_candidate_sha256 = _sha256(candidate_bytes)
    if actual_candidate_sha256 != expected_candidate_sha256:
        raise ValueError(
            "Captured clean-candidate manifest SHA-256 differs from expected binding"
        )
    source = _mapping(candidate, "source", label="Clean-candidate manifest")
    source_git_commit = _lowercase_git_commit(
        _text(source, "git_commit", label="Clean-candidate source"),
        label="Clean-candidate source git_commit",
    )
    if (
        expected_source_git_commit is not None
        and source_git_commit != expected_source_git_commit
    ):
        raise ValueError(
            "Captured clean-candidate source Git commit differs from expected binding"
        )
    receipt = read_json_bytes(
        build_profile_receipt,
        label="clean-candidate build-profile receipt",
    )
    receipt_control_plane_version = _digest(
        receipt,
        "control_plane_version",
        label="Clean-candidate build-profile receipt",
    )
    if (
        expected_receipt_version is not None
        and receipt_control_plane_version != expected_receipt_version
    ):
        raise ValueError(
            "Clean-candidate build-profile receipt belongs to a different "
            "control-plane version"
        )
    validated_submodules = validate_bundle(
        candidate,
        source_archive=tree["source_archive"],
        source_commit=tree["source_commit"],
        submodule_archives=tree["submodule_archives"],
        submodule_commits=tree["submodule_commits"],
        build_profile=tree["build_profile"],
        build_profile_receipt=build_profile_receipt,
        build_provenance=tree["build_provenance"],
        executable_sha256=_sha256(executable),
        authorized_pic_root=authorized_pic_root,
        authorized_source_root=authorized_source_root,
    )
    _verify_historical_installed_pair(
        receipt_control_plane_version,
        authorized_pic_root=authorized_pic_root,
        authorized_project_home_root=authorized_project_home_root,
        verify_historical=verify_historical,
    )
    build = _mapping(candidate, "build", label="Clean-candidate manifest")
    candidate_path = tree.get("candidate_manifest_path")
    if not isinstance(candidate_path, Path):
        raise ValueError("Clean-candidate reader did not return one manifest path")
    return {
        "build": {
            "executable_sha256": _sha256(executable),
            "profile_id": _text(build, "profile_id", label="Clean-candidate build"),
            "receipt_control_plane_version": receipt_control_plane_version,
        },
        "clean_candidate_manifest": {
            "expected_sha256": expected_candidate_sha256,
            "path": str(candidate_path),
            "sha256": actual_candidate_sha256,
        },
        "current_control_plane_version": _digest(
            current_inventory,
            "version",
            label="Installed control-plane inventory",
        ),
        "freeze_id": _text(candidate, "freeze_id", label="Clean-candidate manifest"),
        "record_type": RECORD_TYPE,
        "schema_version": 1,
        "source": {
            "git_commit": source_git_commit,
            "git_tree": _text(source, "git_tree", label="Clean-candidate source"),
            "source_bundle_sha256": _digest(
                source,
                "source_bundle_sha256",
                label="Clean-candidate source",
            ),
        },
        "status": "passed",
        "validated_submodules": validated_submodules,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Read-only revalidation for one fixed-root Orion clean candidate."
    )
    parser.add_argument(
        "--manifest",
        "--clean-candidate-manifest",
        dest="candidate_manifest_path",
        required=True,
        type=Path,
    )
    parser.add_argument("--expected-manifest-sha256", required=True)
    parser.add_argument("--expected-git-commit", type=_expected_git_commit_argument)
    parser.add_argument("--expected-receipt-control-plane-version")
    return parser


def main() -> None:
    if not getattr(sys, "_pic_control_plane_bootstrapped", False):
        raise SystemExit("Run installed control-plane tools through run_control_plane.py")
    args = _parser().parse_args()
    result = revalidate_clean_candidate(
        args.candidate_manifest_path,
        expected_manifest_sha256=args.expected_manifest_sha256,
        expected_git_commit=args.expected_git_commit,
        expected_receipt_control_plane_version=(
            args.expected_receipt_control_plane_version
        ),
    )
    print(_canonical_json_bytes(result).decode("utf-8"), end="")


if __name__ == "__main__":
    main()
