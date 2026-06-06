#!/usr/bin/env python3
"""Regression tests for source-controlled PIC readiness registries."""

from __future__ import annotations

import base64
from contextlib import contextmanager
import copy
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import tempfile
import unittest
from collections.abc import Iterator

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(REPO_ROOT / "tst" / "publication" / "frontier_control_plane"))

from control_plane_common import CONTROL_PLANE_FILES
from control_plane_common import PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
from control_plane_common import inventory_digest
from control_plane_common import launch_contract_sha256
from control_plane_common import validate_launch_contract
from control_plane_common import verify_historical_installed_control_plane
from ledger import record_sha256
from ledger import incomplete_manual_accounting_marker_paths
from ledger import validate_mirrored_state
import q011_pressure_review_packet_verifier as pressure_packet_verifier
from tst.publication import q011_section54_pressure_pilot_execution as pressure_execution
from tst.publication.pic_qualification_manifest import SCHEMA_PATH
from tst.publication.pic_qualification_manifest import validate_qualification_manifest
from tst.publication.pic_qualification_manifest import validate_schema


READINESS_DIR = REPO_ROOT / "tst" / "publication" / "readiness"
SCHEMA_DIR = READINESS_DIR / "schemas"
CONTROL_PLANE_DIR = REPO_ROOT / "tst" / "publication" / "frontier_control_plane"
VALIDATION_MANIFEST_SCHEMA = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
PAPER_TEX = (
    REPO_ROOT
    / "docs"
    / "reference_paper"
    / "arXiv-2304.10568v1"
    / "mnras_template.tex"
)

CLAIM_CLASSES = {
    "unit/regression",
    "engineering_proxy",
    "physics_validation",
    "sun_bai_2023_reproduction",
    "athenak_production_mode",
    "cross_code_comparison",
    "scoped_state_of_the_art",
    "unsupported",
}

REQUIRED_EXTENSION_CLAIMS = {
    "CLAIM-EXT-HALL-BELL-001",
    "CLAIM-EXT-CRSI-IN-DAMPING-001",
    "CLAIM-EXT-CRPAI-TRANSPORT-001",
    "CLAIM-STATEART-CRPAI-SCATTERING-001",
    "CLAIM-RELEASE-EXTENDED-MHD-PIC-001",
}

_REVIEWER_ID_PATTERN = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
)

_Q011_PACKET_GATE_SOURCE_TEST_CLOSURE_PATHS = (
    "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
    "tst/publication/frontier_control_plane/control_plane_common.py",
    "tst/publication/frontier_control_plane/run_control_plane.py",
    "tst/publication/q011_section54_pressure_selection.py",
    "tst/publication/q011_section54_historical_pressure_pilot_consumer.py",
    "tst/publication/q011_section54_qualifying_campaign_execution.py",
    "tst/publication/q011_section54_attempt_manifest_materializer.py",
    "tst/publication/analyze_q011_section54_campaign.py",
    "tst/publication/analyze_q011_section54_numerical_qualification.py",
    "tst/publication/analyze_q011_section54_pressure_pilot.py",
    "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
    "tst/publication/frontier_f1_structured_artifacts.py",
    "tst/publication/publish_q011_section54_campaign_attempt.py",
    "tst/publication/frontier_control_plane/test_q011_pressure_review_packet_verifier.py",
    "tst/publication/frontier_control_plane/test_control_plane.py",
    "tst/publication/test_q011_section54_pressure_selection.py",
    "tst/publication/test_q011_section54_historical_pressure_pilot_consumer.py",
    "tst/publication/test_q011_section54_qualifying_campaign_execution.py",
    "tst/publication/test_q011_section54_attempt_manifest_materializer.py",
    "tst/publication/test_analyze_q011_section54_campaign.py",
    "tst/publication/test_analyze_q011_section54_numerical_qualification.py",
    "tst/publication/test_analyze_q011_section54_pressure_pilot.py",
    "tst/publication/test_publish_q011_section54_pressure_pilot_bundle.py",
    "tst/publication/test_publish_q011_section54_campaign_attempt.py",
    "tst/publication/test_pic_readiness_registry.py",
    "tst/publication/frontier_q011_section54_pressure_gate_validation_job.sh",
    "tst/publication/readiness/q011_section54_pressure_pilot_postrun_aggregate_source_authorization_successor_v6_2026-06-05.json",
    "tst/publication/frontier_control_plane/prepared_pic_artifact_inventory.json",
)
_Q011_EXACT_PREDECESSOR_REPAIR_SOURCE_TEST_CLOSURE_PATHS = (
    "tst/publication/frontier_control_plane/README.md",
    "tst/publication/frontier_control_plane/control_plane_common.py",
    "tst/publication/frontier_control_plane/promote_active_policy.py",
    "tst/publication/frontier_control_plane/revalidate_clean_candidate.py",
    "tst/publication/frontier_control_plane/run_control_plane.py",
    "tst/publication/frontier_control_plane/install_control_plane.py",
    "tst/publication/frontier_control_plane/capture_storage_preflight_evidence.py",
    "tst/publication/frontier_control_plane/storage_preflight.schema.json",
    "tst/publication/frontier_control_plane/write_orion_build_profile.py",
    "tst/publication/frontier_control_plane/test_control_plane.py",
    "tst/publication/q011_section54_pressure_pilot_execution.py",
    "tst/publication/test_q011_section54_pressure_pilot_execution.py",
    "tst/publication/test_revalidate_clean_candidate.py",
    "tst/publication/test_capture_storage_preflight_evidence.py",
    "tst/publication/frontier_q011_clean_candidate_build_freeze_job.sh",
    "tst/publication/test_pic_readiness_registry.py",
)
_Q011_STAGE4_PRESSURE_SELECTION_CANDIDATE_CLOSURE_PATHS = (
    "tst/publication/publish_q011_section54_pressure_selection.py",
    "tst/publication/test_publish_q011_section54_pressure_selection.py",
    "tst/publication/frontier_control_plane/q011_pressure_review_packet_verifier.py",
    "tst/publication/frontier_control_plane/test_q011_pressure_review_packet_verifier.py",
    "tst/publication/frontier_q011_section54_pressure_gate_validation_job.sh",
    "tst/publication/frontier_control_plane/README.md",
    "tst/publication/test_pic_readiness_registry.py",
    "tst/publication/readiness/README.md",
    "tst/publication/PIC_SUN_BAI_RELEASE_HANDOFF.md",
    "tst/publication/PIC_PRODUCTION_READINESS_PLAN.md",
)


def _load(name: str) -> dict[str, object]:
    return json.loads((READINESS_DIR / name).read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _descriptor_bytes(descriptor: int) -> bytes:
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks = []
    while chunk := os.read(descriptor, 1024 * 1024):
        chunks.append(chunk)
    return b"".join(chunks)


def _require_regular_namespace_identity(path: Path, identity: tuple[int, int]) -> None:
    observed = os.stat(path, follow_symlinks=False)
    if not stat.S_ISREG(observed.st_mode) or (observed.st_dev, observed.st_ino) != identity:
        raise ValueError(f"pinned regular file namespace changed: {path}")


@contextmanager
def _pinned_regular_bytes(path: Path) -> Iterator[bytes]:
    descriptor = os.open(
        path,
        os.O_RDONLY | os.O_CLOEXEC | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise ValueError(f"pinned path is not a regular file: {path}")
        identity = (before.st_dev, before.st_ino)
        data = _descriptor_bytes(descriptor)
        after = os.fstat(descriptor)
        if (
            identity != (after.st_dev, after.st_ino)
            or before.st_size != after.st_size
            or len(data) != after.st_size
        ):
            raise ValueError(f"pinned regular file changed during read: {path}")
        _require_regular_namespace_identity(path, identity)
        yield data
        if _descriptor_bytes(descriptor) != data:
            raise ValueError(f"pinned regular file bytes changed during validation: {path}")
        _require_regular_namespace_identity(path, identity)
    finally:
        os.close(descriptor)


def _git_blob_sha256(commit: str, relative_path: str) -> str:
    contents = subprocess.check_output(
        ["git", "show", f"{commit}:{relative_path}"],
        cwd=REPO_ROOT,
    )
    return hashlib.sha256(contents).hexdigest()


def _git_archive_sha256(commit: str) -> str:
    contents = subprocess.check_output(
        ["git", "archive", "--format=tar", commit],
        cwd=REPO_ROOT,
    )
    return hashlib.sha256(contents).hexdigest()


def _canonical_decimal_integer(value: object) -> int:
    if (
        type(value) is not str
        or not value.isascii()
        or not value.isdecimal()
        or value != str(int(value))
    ):
        raise ValueError("identity is not a canonical decimal string")
    return int(value)


def _canonical_utc_second(value: object) -> datetime:
    if type(value) is not str or len(value) != 20:
        raise ValueError("timestamp is not a canonical UTC second")
    try:
        parsed = datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except ValueError as error:
        raise ValueError("timestamp is not a canonical UTC second") from error
    if parsed.strftime("%Y-%m-%dT%H:%M:%SZ") != value:
        raise ValueError("timestamp is not a canonical UTC second")
    return parsed


def _canonical_reviewer_id(value: object) -> str:
    if type(value) is not str or _REVIEWER_ID_PATTERN.fullmatch(value) is None:
        raise ValueError("reviewer identity is not a canonical lowercase UUID")
    return value


def _canonical_nonempty_fact_list(value: object) -> list[str]:
    if type(value) is not list or not value:
        raise ValueError("facts inspected must be a nonempty list")
    for fact in value:
        if (
            type(fact) is not str
            or not fact
            or not fact.isascii()
            or not fact.isprintable()
            or fact != " ".join(fact.split())
        ):
            raise ValueError("facts inspected contains a noncanonical fact")
    return value


def _exact_json_equal(value: object, expected: object) -> bool:
    if type(value) is not type(expected):
        return False
    if type(expected) is dict:
        if len(value) != len(expected):
            return False
        for expected_key, expected_value in expected.items():
            matching_keys = [
                key
                for key in value
                if type(key) is type(expected_key) and key == expected_key
            ]
            if len(matching_keys) != 1 or not _exact_json_equal(
                value[matching_keys[0]], expected_value
            ):
                return False
        return True
    if type(expected) is list:
        return len(value) == len(expected) and all(
            _exact_json_equal(observed, wanted)
            for observed, wanted in zip(value, expected)
        )
    return value == expected


def _regular_files_below(root: Path, pattern: str) -> list[Path]:
    def fail(error: OSError) -> None:
        raise ValueError(f"regular-file count failed below {root}: {error}")

    matches = []
    for current, _directories, files in os.walk(
        root,
        followlinks=False,
        onerror=fail,
    ):
        for name in files:
            path = Path(current) / name
            if Path(name).match(pattern) and stat.S_ISREG(path.lstat().st_mode):
                matches.append(path)
    return sorted(matches)


def _bounded_relative_matches(root: Path, pattern: str, max_depth: int) -> list[str]:
    if type(max_depth) is not int or max_depth < 0:
        raise ValueError("bounded absence search depth is not a nonnegative integer")

    def fail(error: OSError) -> None:
        raise ValueError(f"bounded absence search failed below {root}: {error}")

    matches = []
    for current, directories, files in os.walk(
        root,
        followlinks=False,
        onerror=fail,
    ):
        relative_current = Path(current).relative_to(root)
        depth = len(relative_current.parts)
        if depth >= max_depth:
            directories[:] = []
        for name in files:
            relative = relative_current / name
            root_level_double_star_match = (
                len(relative.parts) == 1
                and pattern.startswith("**/")
                and relative.match(pattern[3:])
            )
            if len(relative.parts) <= max_depth and (
                relative.match(pattern) or root_level_double_star_match
            ):
                matches.append(relative.as_posix())
    return sorted(matches)


def _installed_policy_unlock_snapshot(
    installed_control_plane_dir: Path,
    control_plane_version: str,
) -> dict[str, object]:
    if (
        type(control_plane_version) is not str
        or installed_control_plane_dir.name != control_plane_version
    ):
        raise ValueError("installed control-plane policy verifier binding is malformed")
    script = """
import json
import sys

sys.path.insert(0, sys.argv[1])
from control_plane_common import require_storage_policy_unlock_snapshot

policy, snapshot = require_storage_policy_unlock_snapshot(
    control_plane_version=sys.argv[2],
)
sys.stdout.write(json.dumps({"policy": policy, "snapshot": snapshot}, sort_keys=True))
"""
    payload = subprocess.check_output(
        [
            sys.executable,
            "-I",
            "-B",
            "-c",
            script,
            str(installed_control_plane_dir),
            control_plane_version,
        ],
        cwd="/",
        env={
            "LANG": "C",
            "LC_ALL": "C",
            "PATH": "/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
        },
        text=True,
    )
    value = json.loads(payload)
    if type(value) is not dict or set(value) != {"policy", "snapshot"}:
        raise ValueError("installed control-plane policy verifier returned malformed output")
    return value


def _validate_q011_post_publication_pressure_gate_status_successor(
    value: object,
) -> None:
    if type(value) is not dict:
        raise ValueError("post-publication pressure-gate status is not an object")
    if set(value) != {
        "schema_version",
        "record_type",
        "recorded_utc",
        "source_checkpoint_commit",
        "predecessor_record",
        "predecessor_sha256",
        "postrun_source_authorization_successor",
        "clean_snapshot_pressure_gate_validation_worker",
        "staged_packet_gate_control_plane",
        "live_operational_baseline",
        "source_test_closure",
        "published_pressure_evidence",
        "publication_acceptance_state",
        "scoped_absence_evidence",
        "advisory_pressure_options_memo",
        "pressure_selection",
        "packet_gate",
        "human_input_required_now",
        "first_required_human_input_after_packet_gate_repair",
        "frontier_launch_authorization",
        "status",
        "qualification_effect",
    }:
        raise ValueError("post-publication pressure-gate status shape drifted")
    expected_scalars = {
        "schema_version": 1,
        "record_type": "q011_section54_post_publication_pressure_gate_status_successor",
        "recorded_utc": "2026-06-05T09:23:34Z",
        "source_checkpoint_commit": "ce41d4b29bc646b4f0740468e1026f7308b29026",
        "predecessor_record": (
            "tst/publication/readiness/"
            "q011_section54_sixteenth_lustre_publication_rename_compatibility_"
            "transition_2026-06-05.json"
        ),
        "predecessor_sha256": (
            "8efd2b349da6c0c150448528b3656a84767cacc5bf08eae84f5b4848c9aeb422"
        ),
        "human_input_required_now": False,
        "first_required_human_input_after_packet_gate_repair": (
            "Review the ranked pressure options, select exactly one Section 5.4 "
            "problem/ps_p0 case, and seal the required reviewer attestation and "
            "schema-v3 receipt."
        ),
        "frontier_launch_authorization": (
            "none_no_bound_selection_receipt_and_packet_gate_blocked"
        ),
        "status": (
            "post_publication_pressure_evidence_bound_no_selection_receipt_bound_"
            "packet_gate_blocked"
        ),
        "qualification_effect": (
            "none_no_selection_no_execution_authorization_no_science_claim"
        ),
    }
    for key, expected in expected_scalars.items():
        if not _exact_json_equal(value[key], expected):
            raise ValueError(f"post-publication pressure-gate status {key} drifted")
    if (
        type(value["schema_version"]) is not int
        or type(value["human_input_required_now"]) is not bool
    ):
        raise ValueError("post-publication pressure-gate status scalar type drifted")
    if _canonical_utc_second(value["recorded_utc"]) <= _canonical_utc_second(
        "2026-06-05T05:06:20Z"
    ):
        raise ValueError("post-publication pressure-gate status chronology drifted")
    if not _exact_json_equal(value["postrun_source_authorization_successor"], {
        "path": (
            "tst/publication/readiness/"
            "q011_section54_pressure_pilot_postrun_aggregate_source_authorization_"
            "successor_v5_2026-06-05.json"
        ),
        "sha256": (
            "93c2b9174d1546b883ac4910afe3f7f021ed4324f2831f80b19aee9395e5bd3e"
        ),
    }):
        raise ValueError("postrun source authorization binding drifted")
    clean_snapshot_worker = value["clean_snapshot_pressure_gate_validation_worker"]
    if (
        type(clean_snapshot_worker) is not dict
        or any(
            type(clean_snapshot_worker.get(key)) is not int
            for key in (
                "expected_publication_python_files",
                "expected_publication_shell_files",
                "expected_publication_json_files",
                "expected_publication_test_modules",
            )
        )
        or not _exact_json_equal(clean_snapshot_worker, {
            "path": (
                "tst/publication/"
                "frontier_q011_section54_pressure_gate_validation_job.sh"
            ),
            "sha256": (
                "63ee6ca57229a96c0448134a7c105dcd0de820bbf46bd70182f49ff6514ed4d7"
            ),
            "expected_publication_python_files": 145,
            "expected_publication_shell_files": 17,
            "expected_publication_json_files": 290,
            "expected_publication_test_modules": 65,
            "status": (
                "no_retained_clean_snapshot_validation_evidence_bound_by_this_"
                "status_successor"
            ),
        })
    ):
        raise ValueError("clean-snapshot pressure-gate validation worker drifted")
    staged_packet_gate = value["staged_packet_gate_control_plane"]
    if (
        type(staged_packet_gate) is not dict
        or type(staged_packet_gate.get("inventoried_file_count")) is not int
        or not _exact_json_equal(staged_packet_gate, {
            "version": (
                "ccc9d8aef994bb64465f9b16de58236ff42ed2df028e8d9897ababb30f2cb7f1"
            ),
            "state": (
                "source_local_candidate_differs_from_active_live_generation"
            ),
            "inventoried_file_count": 24,
            "prepared_artifact_inventory": {
                "path": (
                    "tst/publication/frontier_control_plane/"
                    "prepared_pic_artifact_inventory.json"
                ),
                "sha256": (
                    "87f63ab91587b963347909b7d116ce0e1395f3488e5f61cd9249cb8a69dc02bc"
                ),
            },
            "pair_install_authorization_by_this_status": "none",
            "policy_promotion_authorization_by_this_status": "none",
        })
    ):
        raise ValueError("staged packet-gate control-plane candidate drifted")
    live = value["live_operational_baseline"]
    live_version = "821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721"
    inventory_sha256 = (
        "5e8764df23df793212a19f1f7365e077c5be3f6ab58c423560dcc7d557f7ae7b"
    )
    policy_sha256 = (
        "23a73b868146f63d4b363713f988d55e9dadffa15b26f2b2f1d07da66331f5c3"
    )
    promotion_sha256 = (
        "4824ea825e7b9e42becdca4b9a8b72c0454bd1a2e02d5e365b1878ed94c53243"
    )
    orion_root = "/lustre/orion/ast207/proj-shared/dfielding/PIC"
    project_root = "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
    if not _exact_json_equal(live, {
        "installed_control_plane_version": live_version,
        "inventoried_member_count": 23,
        "inventories": [
            {
                "path": f"{orion_root}/control_plane/{live_version}/inventory.json",
                "sha256": inventory_sha256,
            },
            {
                "path": f"{project_root}/control_plane/{live_version}/inventory.json",
                "sha256": inventory_sha256,
            },
        ],
        "inventories_byte_identical": True,
        "active_policies": [
            {
                "path": f"{orion_root}/policy/storage_policy.json",
                "sha256": policy_sha256,
            },
            {
                "path": f"{project_root}/policy/storage_policy.json",
                "sha256": policy_sha256,
            },
        ],
        "active_policies_byte_identical": True,
        "active_promotions": [
            {
                "path": f"{orion_root}/policy/active_promotion.json",
                "sha256": promotion_sha256,
            },
            {
                "path": f"{project_root}/policy/active_promotion.json",
                "sha256": promotion_sha256,
            },
        ],
        "active_promotions_byte_identical": True,
        "projected_policy_state": {
            "installed_control_plane_version": live_version,
            "staged_control_plane_candidate_version": live_version,
            "build_profile_control_plane_version": live_version,
            "registered_science_slices": [],
            "science_submission_freeze_status": "authorized",
        },
        "projected_promotion_state": {
            "control_plane_version": live_version,
            "policy_sha256": policy_sha256,
        },
    }):
        raise ValueError("live operational baseline drifted")
    closure = value["source_test_closure"]
    closure_files = closure.get("files") if type(closure) is dict else None
    if (
        type(closure) is not dict
        or set(closure)
        != {"status", "scope", "file_count", "closure_sha256", "files"}
        or closure["status"]
        != "source_local_change_and_validation_bytes_bound_by_this_status_successor"
        or closure["scope"]
        != (
            "pressure_selection_v3_sealed_reanalysis_reviewer_and_replay_change_"
            "validation_closure_not_transitive_runtime_closure"
        )
        or type(closure["file_count"]) is not int
        or closure["file_count"] != len(_Q011_PACKET_GATE_SOURCE_TEST_CLOSURE_PATHS)
        or type(closure["closure_sha256"]) is not str
        or re.fullmatch(r"[0-9a-f]{64}", closure["closure_sha256"]) is None
        or type(closure_files) is not list
        or any(type(binding) is not dict for binding in closure_files)
        or [binding.get("path") for binding in closure_files]
        != list(_Q011_PACKET_GATE_SOURCE_TEST_CLOSURE_PATHS)
    ):
        raise ValueError("packet-gate source/test closure drifted")
    for binding in closure_files:
        if (
            type(binding) is not dict
            or set(binding) != {"path", "sha256"}
            or type(binding["sha256"]) is not str
            or re.fullmatch(r"[0-9a-f]{64}", binding["sha256"]) is None
        ):
            raise ValueError("packet-gate source/test closure binding drifted")
    if closure["closure_sha256"] != inventory_digest(closure_files):
        raise ValueError("packet-gate source/test closure digest drifted")
    if not _exact_json_equal(value["published_pressure_evidence"], {
        "aggregate_manifest": {
            "path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
                "q011_section54_pressure_pilot_bundle/pressure_pilot_manifest.json"
            ),
            "sha256": (
                "7b3fb8de9e4dc8d6d8b2320a2c2aeaadf8bc6f076dab6b3e98ef8a415abdbf55"
            ),
        },
        "aggregate_analysis": {
            "path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
                "q011_section54_pressure_pilot_analysis.json"
            ),
            "sha256": (
                "d55b4c2020716df899c86dfe5a9018d48194d63590616ff60541067243daacb7"
            ),
        },
        "aggregate_receipt": {
            "path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
                "q011_section54_pressure_pilot_bundle_receipt.json"
            ),
            "sha256": (
                "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9"
            ),
        },
        "review_packet_receipt": {
            "path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
                "q011_section54_pressure_pilot_review_packet_receipt.json"
            ),
            "sha256": (
                "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d"
            ),
        },
        "review_packet_inventory": {
            "path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/publication/"
                "q011_section54_pressure_pilot_review_packet/packet_inventory.json"
            ),
            "sha256": (
                "ba38e86575baee720871e1f624e20edb5482f9a4a00afb4286ee136f821112c1"
            ),
        },
    }):
        raise ValueError("published pressure-evidence binding drifted")
    publication_root = f"{orion_root}/publication"
    acceptance_root = f"{orion_root}/publication_acceptance"
    aggregate_receipt_path = (
        f"{publication_root}/q011_section54_pressure_pilot_bundle_receipt.json"
    )
    packet_receipt_path = (
        f"{publication_root}/q011_section54_pressure_pilot_review_packet_receipt.json"
    )
    if not _exact_json_equal(value["publication_acceptance_state"], {
        "publication_root": {
            "path": publication_root,
            "identity": {"device": 135357496, "inode": 720587416766821486},
            "mode": "02755",
            "exact_entries": [
                "q011_section54_pressure_pilot_analysis.json",
                "q011_section54_pressure_pilot_bundle",
                "q011_section54_pressure_pilot_bundle_receipt.json",
                "q011_section54_pressure_pilot_review_packet",
                "q011_section54_pressure_pilot_review_packet_receipt.json",
            ],
        },
        "acceptance_root": {
            "path": acceptance_root,
            "identity": {"device": 135357496, "inode": 720587400627193972},
            "mode": "0700",
            "exact_entries": [
                ".q011_section54_pressure_pilot_bundle_receipt.json.publication-success",
                ".q011_section54_pressure_pilot_review_packet_receipt.json.publication-success",
            ],
        },
        "success_seals": [
            {
                "path": (
                    f"{acceptance_root}/"
                    ".q011_section54_pressure_pilot_bundle_receipt.json."
                    "publication-success"
                ),
                "sha256": (
                    "8eeb7eb24f9aa64652a27d619b155e959a302788eb99b8b4fa166052ae0158f8"
                ),
                "identity": {"device": 135357496, "inode": 720587443073450467},
                "mode": "0444",
                "receipt": {
                    "path": aggregate_receipt_path,
                    "sha256": (
                        "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9"
                    ),
                    "identity": {"device": 135357496, "inode": 720587443073450466},
                },
            },
            {
                "path": (
                    f"{acceptance_root}/"
                    ".q011_section54_pressure_pilot_review_packet_receipt.json."
                    "publication-success"
                ),
                "sha256": (
                    "109bc4522579feab92ef1778c0318a58e96aa8d059e481283c6faa83cd6c9e94"
                ),
                "identity": {"device": 135357496, "inode": 720587443090227212},
                "mode": "0444",
                "receipt": {
                    "path": packet_receipt_path,
                    "sha256": (
                        "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d"
                    ),
                    "identity": {"device": 135357496, "inode": 720587443090227211},
                },
            },
        ],
        "absent_publication_paths": [
            (
                f"{publication_root}/"
                ".q011_section54_pressure_pilot_bundle_receipt.json."
                "publication-invalid"
            ),
            (
                f"{publication_root}/"
                ".q011_section54_pressure_pilot_review_packet_receipt.json."
                "publication-invalid"
            ),
        ],
        "absent_namespace_globs": [
            {"root": publication_root, "glob": ".*.staging-*", "matches": []},
            {"root": acceptance_root, "glob": ".*.staging-*", "matches": []},
        ],
    }):
        raise ValueError("publication acceptance state drifted")
    if not _exact_json_equal(value["scoped_absence_evidence"], {
        "pressure_selection_receipt_searches": [
            {
                "root": "/ccs/home/dfielding/athenak-pic/tst/publication",
                "glob": "**/*pressure*selection*receipt*.json",
                "max_depth": 8,
                "matches": [],
            },
            {
                "root": orion_root,
                "glob": "**/*pressure*selection*receipt*.json",
                "max_depth": 6,
                "matches": [],
            },
            {
                "root": project_root,
                "glob": "**/*pressure*selection*receipt*.json",
                "max_depth": 6,
                "matches": [],
            },
        ],
        "clean_snapshot_validation_log_search": {
            "root": f"{orion_root}/logs/slurm",
            "glob": "pic-q011-pressure-gate-validate.*.log",
            "max_depth": 1,
            "matches": [],
        },
        "registered_science_slices": [],
    }):
        raise ValueError("scoped pressure-gate absence evidence drifted")
    if not _exact_json_equal(value["advisory_pressure_options_memo"], {
        "path": (
            "tst/publication/readiness/"
            "q011_section54_pressure_selection_options_2026-06-05.md"
        ),
        "sha256": (
            "330dba02daed0f3cf633d04e5e50da5d3208bbb54daba18f0828ae744e91ba86"
        ),
        "role": "human_review_memo_only_not_a_pressure_selection_receipt",
    }):
        raise ValueError("advisory pressure-options memo binding drifted")
    pressure_selection = value["pressure_selection"]
    if (
        type(pressure_selection) is not dict
        or type(pressure_selection.get("advisory_recommendation_is_selection"))
        is not bool
        or not _exact_json_equal(pressure_selection, {
            "status": "no_selection_receipt_bound_by_this_status_successor",
            "selection_receipt": None,
            "selected_case": None,
            "advisory_recommendation_is_selection": False,
        })
    ):
        raise ValueError("pressure-selection absence contract drifted")
    packet_gate = value["packet_gate"]
    if (
        type(packet_gate) is not dict
        or type(packet_gate.get("required_receipt_schema_version")) is not int
        or not _exact_json_equal(packet_gate, {
            "status": "blocked",
            "acceptance_authorization": "none",
            "required_receipt_schema_version": 3,
            "required_acceptance_and_replay_paths": [
                "source_local_pressure_selection_acceptance",
                "retained_qualifying_plan_replay",
                "installed_control_plane_pressure_selection_acceptance",
                "completed_attempt_replay",
            ],
            "blockers_bound_by_this_status_successor": [
                "pressure_selection_v3_gate_repair_commit_not_bound_by_"
                "this_status_successor",
                "pressure_selection_v3_gate_repair_independent_review_not_"
                "bound_by_this_status_successor",
                "pressure_selection_v3_gate_repair_clean_worker_validation_"
                "not_bound_by_this_status_successor",
                "pressure_selection_v3_gate_repair_pair_install_and_"
                "promotion_not_bound_by_this_status_successor",
                "sealed_authoritative_reanalysis_attestation_not_bound_by_this_"
                "status_successor",
                "sealed_human_reviewer_attestation_and_v3_selection_receipt_not_bound_"
                "by_this_status_successor",
            ],
        })
    ):
        raise ValueError("pressure-selection packet-gate contract drifted")


def _validate_q011_exact_predecessor_migration_repair_successor(
    value: object,
) -> None:
    if type(value) is not dict or set(value) != {
        "schema_version",
        "record_type",
        "recorded_utc",
        "source_checkpoint_commit",
        "predecessor_record",
        "predecessor_sha256",
        "retained_clean_worker",
        "retained_failed_migration_attempt",
        "retained_failed_build_freeze_attempt",
        "unchanged_live_operational_baseline",
        "source_local_exact_predecessor_repair",
        "next_clean_worker",
        "human_pressure_selection",
        "packet_gate",
        "frontier_launch_authorization",
        "qualification_effect",
        "status",
    }:
        raise ValueError("exact-predecessor migration repair successor shape drifted")
    if not _exact_json_equal(
        {
            key: value[key]
            for key in (
                "schema_version",
                "record_type",
                "recorded_utc",
                "source_checkpoint_commit",
                "predecessor_record",
                "predecessor_sha256",
                "frontier_launch_authorization",
                "qualification_effect",
                "status",
            )
        },
        {
            "schema_version": 1,
            "record_type": (
                "q011_section54_pressure_gate_exact_predecessor_migration_"
                "repair_successor"
            ),
            "recorded_utc": "2026-06-06T04:56:47Z",
            "source_checkpoint_commit": "f6471610a116ce5433550625c9dc61752315040f",
            "predecessor_record": (
                "tst/publication/readiness/"
                "q011_section54_post_publication_pressure_gate_status_successor_"
                "2026-06-05.json"
            ),
            "predecessor_sha256": (
                "0ba6ac19df478630830b13a68443e0560cfcf9dbd00d48f20ad25d59a39568be"
            ),
            "frontier_launch_authorization": "none_launch_prohibited",
            "qualification_effect": "none_no_science_claim",
            "status": (
                "accepted_build_freeze_timeout_bounded_successor_clone_repair_"
                "pending_commit_clean_worker_pair_install_and_promotion"
            ),
        },
    ):
        raise ValueError("exact-predecessor migration repair successor scalar drifted")
    if type(value["schema_version"]) is not int:
        raise ValueError("exact-predecessor migration repair schema type drifted")
    if not _exact_json_equal(
        value["retained_clean_worker"],
        {
            "job_id": "4769961",
            "job_name": "pic-q011-pressure-gate-validate",
            "state": "COMPLETED",
            "exit_code": "0:0",
            "submitted_utc": "2026-06-05T23:34:19",
            "started_utc": "2026-06-05T23:34:55",
            "ended_utc": "2026-06-05T23:46:20",
            "source_commit": "f6471610a116ce5433550625c9dc61752315040f",
            "worker": {
                "path": (
                    "tst/publication/"
                    "frontier_q011_section54_pressure_gate_validation_job.sh"
                ),
                "sha256": (
                    "3966f56d61bd73164873c88aac906da27cc002a11a88ea11d176329945a37d97"
                ),
                "expected_publication_python_files": 145,
                "expected_publication_shell_files": 17,
                "expected_publication_json_files": 291,
                "expected_publication_test_modules": 65,
            },
            "log": {
                "path": (
                    "/lustre/orion/ast207/proj-shared/dfielding/PIC/logs/slurm/"
                    "pic-q011-pressure-gate-validate.4769961.log"
                ),
                "sha256": (
                    "1bed9421b7e32d84c4b48ee0ada5d2189860e26363d640a632f51651772f3074"
                ),
            },
            "archive_focused_tests": 199,
            "publication_tests": 1508,
            "publication_test_skips": 2,
        },
    ):
        raise ValueError("retained clean-worker binding drifted")
    orion_root = "/lustre/orion/ast207/proj-shared/dfielding/PIC"
    project_root = "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
    failed_version = "ccc9d8aef994bb64465f9b16de58236ff42ed2df028e8d9897ababb30f2cb7f1"
    failed_inventory_sha256 = (
        "b843109601bda18cd5e9ebf0c7634480ec97d4d7feec07c1cc576964fb8f77eb"
    )
    if not _exact_json_equal(
        value["retained_failed_migration_attempt"],
        {
            "paired_control_plane": {
                "version": failed_version,
                "inventoried_file_count": 24,
                "orion_inventory": {
                    "path": f"{orion_root}/control_plane/{failed_version}/inventory.json",
                    "sha256": failed_inventory_sha256,
                },
                "project_home_inventory": {
                    "path": f"{project_root}/control_plane/{failed_version}/inventory.json",
                    "sha256": failed_inventory_sha256,
                },
                "inventories_byte_identical": True,
            },
            "fresh_storage_preflight": {
                "binding": {
                    "path": (
                        f"{orion_root}/policy/storage_preflight_bindings/"
                        "ba77f665-9bb1-4441-8475-f16ea6aa6ab6.json"
                    ),
                    "sha256": (
                        "25b4f01f92b52ef7e163e178f49fd2ab8519e5e2d035133752354b77aef24d36"
                    ),
                },
                "probe_id": "ba77f665-9bb1-4441-8475-f16ea6aa6ab6",
                "completed_utc": "2026-06-05T11:54:44.378099Z",
                "evidence_sha256": (
                    "6938aabe3da76c9dfffbfd6dce3b134cd65b0a9e052dc9bb01bedbcf152e7b48"
                ),
                "orion_evidence_path": (
                    f"{orion_root}/policy/storage_preflight_evidence/"
                    "ba77f665-9bb1-4441-8475-f16ea6aa6ab6.json"
                ),
                "project_home_evidence_path": (
                    f"{project_root}/policy/storage_preflight_evidence/"
                    "ba77f665-9bb1-4441-8475-f16ea6aa6ab6.json"
                ),
            },
            "reviewed_policy": {
                "path": (
                    f"{orion_root}/policy/"
                    "reviewed_launch_prohibited_pressure_gate_successor_"
                    "ccc9d8ae_ba77f665-9bb1-4441-8475-f16ea6aa6ab6.json"
                ),
                "sha256": (
                    "fcfc444d66b5c90986c3cc4e8ea6400858589e2eaea14e5b755703a89e040e0e"
                ),
            },
            "promotion_result": (
                "failed_closed_live_predecessor_source_authentication_not_authorized"
            ),
            "active_policy_changed": False,
            "active_promotion_changed": False,
            "authority": "none",
        },
    ):
        raise ValueError("retained failed migration-attempt binding drifted")
    if not _exact_json_equal(
        value["retained_failed_build_freeze_attempt"],
        {
            "job_id": "4769975",
            "job_name": "pic-q011-build-freeze",
            "top_level_scheduler_record": (
                "4769975|pic-q011-build-freeze|ast207|batch|debug|TIMEOUT|0:0|"
                "01:00:24|01:00:00|"
                "/autofs/nccs-svm1_home2/dfielding/athenak-pic"
            ),
            "batch_step_record": "4769975.batch|batch|CANCELLED|0:15|01:00:25",
            "worker_inputs": {
                "source_commit": "f6471610a116ce5433550625c9dc61752315040f",
                "control_plane_version": (
                    "b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48"
                ),
                "expected_active_policy_sha256": (
                    "48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1"
                ),
                "expected_active_promotion_sha256": (
                    "4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f"
                ),
                "expected_authorized_freeze_manifest": (
                    f"{orion_root}/clean_candidates/"
                    "98a372c9-2ea0-47e6-ad34-e66343e7eea1/"
                    "clean_candidate_manifest.json"
                ),
                "expected_authorized_freeze_manifest_sha256": (
                    "dee6be45657e99ec477eec513c45be4ba6ac43b4fd5deb5655750f51c668e42f"
                ),
                "expected_authorized_freeze_build_controller": (
                    "821d185856722bd0178acb9427f78ac82671a4b6670779ec8400fbac54c6d721"
                ),
            },
            "log": {
                "path": (
                    f"{orion_root}/logs/slurm/"
                    "pic-q011-build-freeze.4769975.log"
                ),
                "sha256": (
                    "d07a92eb89c5646bbb25977b9076b3e05e7a1b3935abdb1e1f4448741446cffa"
                ),
            },
            "failure_phase": "source_clone_before_checkout_configure_build_or_freeze",
            "published_outputs": {
                "clean_candidate_manifest": False,
                "candidate_only_policy": False,
                "science_launch": False,
            },
            "residue": {
                "build_root": {
                    "path": (
                        f"{orion_root}/build/f6471610a116/"
                        "hip-mpi-release-paper-pic"
                    ),
                    "top_level_entries": ["source"],
                    "file_count": 6912,
                    "directory_count": 209,
                    "symlink_count": 0,
                    "status": "preserved_in_place_non_authoritative_chronology",
                },
                "bin_root": {
                    "path": (
                        f"{orion_root}/bin/f6471610a116/"
                        "hip-mpi-release-paper-pic"
                    ),
                    "entry_count": 0,
                    "status": "preserved_in_place_non_authoritative_chronology",
                },
            },
            "authority": "none",
        },
    ):
        raise ValueError("retained failed build/freeze-attempt binding drifted")
    live_version = "b56d96b40f2c666b9fa5b421fac589d6a4d6500a716f20b479c354a3d239cb48"
    live_policy_sha256 = (
        "48e74f3151b51ba84b72e3214ef4a66c03441c745912b833395f98f6f60a0bd1"
    )
    live_promotion_sha256 = (
        "4ecb45bb399cee750ce69198286afc9fb903dde92cffc7507600b0a53f1c547f"
    )
    if not _exact_json_equal(
        value["unchanged_live_operational_baseline"],
        {
            "installed_control_plane_version": live_version,
            "active_policies": [
                {
                    "path": f"{orion_root}/policy/storage_policy.json",
                    "sha256": live_policy_sha256,
                },
                {
                    "path": f"{project_root}/policy/storage_policy.json",
                    "sha256": live_policy_sha256,
                },
            ],
            "active_policies_byte_identical": True,
            "active_promotions": [
                {
                    "path": f"{orion_root}/policy/active_promotion.json",
                    "sha256": live_promotion_sha256,
                },
                {
                    "path": f"{project_root}/policy/active_promotion.json",
                    "sha256": live_promotion_sha256,
                },
            ],
            "active_promotions_byte_identical": True,
            "registered_science_slices": [],
            "science_submission_freeze_status": "authorized",
        },
    ):
        raise ValueError("unchanged live operational baseline drifted")
    repair = value["source_local_exact_predecessor_repair"]
    closure = repair.get("source_test_closure") if type(repair) is dict else None
    closure_files = closure.get("files") if type(closure) is dict else None
    if (
        type(repair) is not dict
        or set(repair)
        != {
            "control_plane_version",
            "inventoried_file_count",
            "state",
            "final_binding_refresh",
            "exact_migration_contract",
            "source_test_closure",
        }
        or repair["control_plane_version"]
        != "930a04d1d39c873ea49abfcf500069011f6d5759240a8f5c5b3341a6d243b246"
        or type(repair["inventoried_file_count"]) is not int
        or repair["inventoried_file_count"] != 24
        or repair["state"]
        != "source_local_uncommitted_bounded_clone_and_schema_v2_predecessor_migration_repair"
        or not _exact_json_equal(
            repair["final_binding_refresh"],
            {
                "status": "completed_before_commit",
                "fields": ["control_plane_version", "source_test_closure"],
                "authority": "none",
            },
        )
        or repair["exact_migration_contract"]
        != [
            "exact_active_policy_and_promotion_hashes",
            "exact_active_controller_probe_evidence_and_source_authentication",
            "schema_v2_exact_predecessor_evidence_authorized",
            "new_control_plane_required",
            "empty_registered_science_allowlist_required",
            "authorized_science_freeze_preserved",
            "only_controller_and_different_strictly_newer_preflight_binding_change",
            "schema_v2_unique_promotion_identity_and_aba_resistant_compare_and_swap",
            "durable_mirrored_prepared_committed_four_anchor_transaction",
            "authorized_successor_revalidated_before_commit_after_both_committed_markers_on_normal_path_and_during_recovery",
            "preserved_authorized_freeze_revalidated_before_any_successor_anchor_publication",
            "exact_active_generation_verifier_authenticates_executing_installed_controller_and_holds_stable_serialization_anchor",
            "exact_successor_anchors_bracket_authorized_candidate_revalidation",
            "same_visible_transaction_generation_required_for_each_validation",
            "complete_successor_rollback_anchors_removed_before_marker_cleanup",
            "authorized_successor_revalidation_binds_current_and_freeze_build_receipt_controllers",
            "preserved_authorized_freeze_allows_paired_historical_build_receipt",
            "candidate_or_active_anchor_drift_at_complete_commit_boundary_retains_locked_recovery_evidence",
            "visible_ambiguous_postcommit_marker_state_retains_locked_recovery_evidence",
            "active_readers_fail_closed_while_transaction_marker_or_reserved_rollback_anchor_exists",
            "markerless_or_unexpected_rollback_anchor_evidence_requires_reviewed_manual_recovery",
            "complete_predecessor_absent_markers_leave_fail_closed_rollback_anchor_evidence",
            "complete_valid_successor_rolls_forward_regardless_of_prepared_or_mixed_marker_state",
            "complete_invalid_successor_retains_locked_recovery_evidence",
            "committed_recovery_validates_exact_successor_digests_semantics_and_controller_before_evidence_cleanup",
            "partial_successor_with_any_committed_marker_retains_all_transaction_evidence_for_reviewed_manual_recovery",
            "strict_partial_prepared_recovery_requires_reachable_publication_prefix_and_semantic_predecessor_before_mutation",
            "prepared_rollback_restores_reverse_publication_order",
            "prepared_rollback_revalidates_restored_predecessor_before_evidence_cleanup",
            "interrupted_transaction_recovery_requires_locked_promotion",
            "authorized_freeze_changes_require_exact_active_predecessor_compare_and_swap",
            "exact_freeze_replacement_revalidates_full_clean_candidate_bundle",
            "exact_freeze_replacement_requires_successor_controller_build_receipt",
            "active_predecessor_revalidated_immediately_before_transaction_setup",
            "exact_commit_required_for_source_runner_installer_and_q011_materializers",
            "preflight_capture_and_recovery_hold_stable_serialization_anchor",
            "preflight_capture_source_authentication_binds_executing_common_module",
            "preflight_recovery_requires_exact_pair_digest_source_and_existing_role",
            "preflight_evidence_publication_is_commit_forward_and_preserves_ambiguous_residue",
            "trusted_git_exact_clone_capabilities_required_before_build_path_creation",
            "bounded_top_level_exact_revision_and_full_independent_submodule_transport_required",
            "prebuild_exact_recursive_submodule_status_equivalence_required",
            "normal_unlock_and_promotion_remain_strict",
        ]
        or type(closure) is not dict
        or set(closure) != {"file_count", "closure_sha256", "files"}
        or type(closure["file_count"]) is not int
        or closure["file_count"]
        != len(_Q011_EXACT_PREDECESSOR_REPAIR_SOURCE_TEST_CLOSURE_PATHS)
        or type(closure_files) is not list
        or any(type(binding) is not dict for binding in closure_files)
        or [binding.get("path") for binding in closure_files]
        != list(_Q011_EXACT_PREDECESSOR_REPAIR_SOURCE_TEST_CLOSURE_PATHS)
    ):
        raise ValueError("source-local exact-predecessor repair binding drifted")
    historical_closure_commit = "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4"
    for binding in closure_files:
        if (
            type(binding) is not dict
            or set(binding) != {"path", "sha256"}
            or type(binding["path"]) is not str
            or type(binding["sha256"]) is not str
            or re.fullmatch(r"[0-9a-f]{64}", binding["sha256"]) is None
            or _git_blob_sha256(historical_closure_commit, binding["path"])
            != binding["sha256"]
        ):
            raise ValueError("exact-predecessor repair source/test closure drifted")
    if closure["closure_sha256"] != inventory_digest(closure_files):
        raise ValueError("exact-predecessor repair source/test closure digest drifted")
    next_worker = value["next_clean_worker"]
    if (
        type(next_worker) is not dict
        or any(
            type(next_worker.get(key)) is not int
            for key in (
                "expected_publication_python_files",
                "expected_publication_shell_files",
                "expected_publication_json_files",
                "expected_publication_test_modules",
            )
        )
        or not _exact_json_equal(
            next_worker,
            {
                "path": (
                    "tst/publication/"
                    "frontier_q011_section54_pressure_gate_validation_job.sh"
                ),
                "sha256": (
                    "3966f56d61bd73164873c88aac906da27cc002a11a88ea11d176329945a37d97"
                ),
                "expected_publication_python_files": 145,
                "expected_publication_shell_files": 17,
                "expected_publication_json_files": 291,
                "expected_publication_test_modules": 65,
                "status": "pending_clean_committed_worker_validation",
            },
        )
    ):
        raise ValueError("next clean-worker binding drifted")
    if not _exact_json_equal(
        value["human_pressure_selection"],
        {
            "status": "human_choice_recorded_not_yet_sealed_or_authoritative",
            "selected_case": {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0},
            "reviewer_id": "dfielding",
            "rationale": (
                "Selected p0=1.0 as the recommended baseline because it explicitly "
                "matches Bai et al. (2015), which uses P0=T0=1 and treats the choice "
                "as unimportant while thermal pressure is much smaller than ram "
                "pressure."
            ),
            "selection_receipt": None,
        },
    ):
        raise ValueError("human pressure-selection checkpoint drifted")
    if not _exact_json_equal(
        value["packet_gate"],
        {
            "status": "blocked",
            "required_receipt_schema_version": 3,
            "remaining_actions": [
                "commit_push_bounded_clone_and_schema_v2_predecessor_migration_repair",
                "clean_worker_validate_bounded_clone_repair_commit",
                "capture_fresh_preflight_pair_install_and_promote_repaired_controller",
                "build_freeze_revalidate_and_promote_fresh_clean_candidate",
                "seal_authoritative_reanalysis_reviewer_and_selection_receipt",
                "verify_all_pressure_selection_replay_boundaries",
            ],
            "authority": "none",
        },
    ):
        raise ValueError("exact-predecessor repair packet-gate state drifted")


def _validate_q011_stage4_pressure_selection_candidate_successor(
    value: object,
) -> None:
    if type(value) is not dict or set(value) != {
        "schema_version",
        "record_type",
        "recorded_utc",
        "source_checkpoint_commit",
        "historical_predecessor_final_binding_commit",
        "predecessor_record",
        "predecessor_sha256",
        "active_candidate_only_state",
        "human_pressure_selection",
        "stage4_publication_candidate",
        "source_local_validation",
        "next_clean_worker",
        "packet_gate",
        "frontier_launch_authorization",
        "qualification_effect",
        "status",
    }:
        raise ValueError("Stage-4 pressure-selection candidate successor shape drifted")
    if not _exact_json_equal(
        {
            key: value[key]
            for key in (
                "schema_version",
                "record_type",
                "recorded_utc",
                "source_checkpoint_commit",
                "historical_predecessor_final_binding_commit",
                "predecessor_record",
                "predecessor_sha256",
                "frontier_launch_authorization",
                "qualification_effect",
                "status",
            )
        },
        {
            "schema_version": 1,
            "record_type": (
                "q011_section54_pressure_selection_publication_candidate_successor"
            ),
            "recorded_utc": "2026-06-06T06:56:40Z",
            "source_checkpoint_commit": "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4",
            "historical_predecessor_final_binding_commit": (
                "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4"
            ),
            "predecessor_record": (
                "tst/publication/readiness/"
                "q011_section54_pressure_gate_exact_predecessor_migration_repair_"
                "successor_2026-06-05.json"
            ),
            "predecessor_sha256": (
                "71aa5f5e6c033db5e51e2a225d4bf233a9f52d4ab5d39f5148c19779ba16b661"
            ),
            "frontier_launch_authorization": "none_launch_prohibited",
            "qualification_effect": "none_no_science_claim",
            "status": (
                "stage4_pressure_selection_publication_candidate_pending_commit_"
                "push_clean_worker_rereview_and_no_science_publication"
            ),
        },
    ):
        raise ValueError("Stage-4 pressure-selection candidate successor scalar drifted")
    if type(value["schema_version"]) is not int:
        raise ValueError("Stage-4 pressure-selection candidate schema type drifted")
    orion_root = "/lustre/orion/ast207/proj-shared/dfielding/PIC"
    project_root = "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
    control_plane_version = (
        "930a04d1d39c873ea49abfcf500069011f6d5759240a8f5c5b3341a6d243b246"
    )
    if not _exact_json_equal(
        value["active_candidate_only_state"],
        {
            "control_plane_version": control_plane_version,
            "active_policy": {
                "orion_path": f"{orion_root}/policy/storage_policy.json",
                "project_home_path": f"{project_root}/policy/storage_policy.json",
                "sha256": (
                    "aeab7e4ef92c7cbbd5b84fcd139c046f2a96f21b4f280fa6dd89c80deca981d1"
                ),
            },
            "active_promotion": {
                "orion_path": f"{orion_root}/policy/active_promotion.json",
                "project_home_path": f"{project_root}/policy/active_promotion.json",
                "sha256": (
                    "ef11cb301ec4917cd32aaca8af56e4f2d753367682c6ba28fc048c004904613e"
                ),
            },
            "clean_candidate": {
                "manifest_path": (
                    f"{orion_root}/clean_candidates/"
                    "83dce7b7-0b03-4be2-b6da-17bb211d1fd4/"
                    "clean_candidate_manifest.json"
                ),
                "manifest_sha256": (
                    "ea5f295096b04d7e5f338873c2a29213f173677568fd33220e1e34ea239739e4"
                ),
                "git_commit": "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4",
                "source_archive_sha256": (
                    "7d4d84a11b5db5b1231c39a9e6a4fbe0f358c0187e0348dec4097f9a69de61cf"
                ),
            },
            "registered_science_slices": [],
            "frontier_admission_smoke": {"status": "closed_after_pass"},
            "pending_submission_marker": "absent",
            "manual_accounting_authorization_count": 4,
            "pending_manual_accounting_marker": "absent",
            "active_promotion_transaction": "absent",
        },
    ):
        raise ValueError("Stage-4 active candidate-only state drifted")
    rationale = (
        "Selected p0=1.0 as the recommended baseline because it explicitly matches "
        "Bai et al. (2015), which uses P0=T0=1 and treats the choice as unimportant "
        "while thermal pressure is much smaller than ram pressure."
    )
    if not _exact_json_equal(
        value["human_pressure_selection"],
        {
            "status": (
                "prior_p0_1p0_choice_recorded_post_reanalysis_human_decision_required"
            ),
            "selected_case": {"case_id": "ps_p0_1p00", "problem_ps_p0": 1.0},
            "reviewer_id": "dfielding",
            "rationale": rationale,
            "production_selection_receipt": None,
        },
    ):
        raise ValueError("Stage-4 human pressure selection drifted")
    candidate = value["stage4_publication_candidate"]
    closure = candidate.get("source_test_closure") if type(candidate) is dict else None
    closure_files = closure.get("files") if type(closure) is dict else None
    if (
        type(candidate) is not dict
        or set(candidate)
        != {
            "status",
            "publisher_mutation_execution_mode",
            "preparation_outputs",
            "publication_outputs",
            "no_science_contract",
            "source_test_closure",
        }
        or candidate["status"]
        != "source_local_uncommitted_candidate_pending_commit_push_clean_worker_and_publication"
        or candidate["publisher_mutation_execution_mode"]
        != (
            "authenticated_read_only_git_archive_expected_commit_and_active_"
            "candidate_archive_reanalysis"
        )
        or candidate["preparation_outputs"]
        != [
            "sealed_authoritative_reanalysis_attestation",
            "sealed_stage4_preparation_source_attestation",
            "empty_private_human_decision_root",
            "no_reviewer_attestation_or_candidate_receipt",
        ]
        or candidate["publication_outputs"]
        != [
            "explicit_post_reanalysis_human_decision_required",
            "sealed_post_reanalysis_human_reviewer_attestation",
            "read_only_non_authorizing_candidate_receipt",
            "read_only_candidate_publication_authorization",
            "canonical_schema_v3_pressure_selection_receipt",
            "sealed_launch_prohibited_controller_state_attestation",
            "inode_bound_publication_success_seal",
            "paired_live_active_state_reverification",
            "absent_publication_guard",
        ]
        or candidate["no_science_contract"]
        != [
            "no_registered_science_slice",
            "no_admission_smoke_authority",
            "no_scheduler_submission",
            "no_science_claim",
            "active_policy_and_promotion_unchanged",
        ]
        or type(closure) is not dict
        or set(closure) != {"file_count", "closure_sha256", "files"}
        or closure["file_count"]
        != len(_Q011_STAGE4_PRESSURE_SELECTION_CANDIDATE_CLOSURE_PATHS)
        or type(closure_files) is not list
        or [binding.get("path") for binding in closure_files if type(binding) is dict]
        != list(_Q011_STAGE4_PRESSURE_SELECTION_CANDIDATE_CLOSURE_PATHS)
    ):
        raise ValueError("Stage-4 publication candidate binding drifted")
    for binding in closure_files:
        if (
            type(binding) is not dict
            or set(binding) != {"path", "sha256"}
            or type(binding["path"]) is not str
            or type(binding["sha256"]) is not str
            or re.fullmatch(r"[0-9a-f]{64}", binding["sha256"]) is None
            or _sha256(REPO_ROOT / binding["path"]) != binding["sha256"]
        ):
            raise ValueError("Stage-4 publication candidate source/test closure drifted")
    if closure["closure_sha256"] != inventory_digest(closure_files):
        raise ValueError("Stage-4 publication candidate closure digest drifted")
    if not _exact_json_equal(
        value["source_local_validation"],
        {
            "focused_publisher_tests_passed": 39,
            "related_pressure_gate_tests_passed": 102,
            "readiness_registry_status": "pending_candidate_successor_finalization",
            "full_publication_suite_status": "pending_clean_committed_worker",
            "independent_rereview_status": "in_progress",
        },
    ):
        raise ValueError("Stage-4 source-local validation state drifted")
    worker = value["next_clean_worker"]
    worker_path = (
        "tst/publication/frontier_q011_section54_pressure_gate_validation_job.sh"
    )
    if (
        type(worker) is not dict
        or set(worker)
        != {
            "path",
            "sha256",
            "expected_publication_python_files",
            "expected_publication_shell_files",
            "expected_publication_json_files",
            "expected_publication_test_modules",
            "status",
        }
        or worker["path"] != worker_path
        or worker["sha256"] != _sha256(REPO_ROOT / worker_path)
        or worker["expected_publication_python_files"] != 147
        or worker["expected_publication_shell_files"] != 17
        or worker["expected_publication_json_files"] != 292
        or worker["expected_publication_test_modules"] != 66
        or worker["status"] != "pending_clean_committed_worker_validation"
    ):
        raise ValueError("Stage-4 next clean-worker binding drifted")
    if not _exact_json_equal(
        value["packet_gate"],
        {
            "status": "blocked_pending_stage4_no_science_publication",
            "required_receipt_schema_version": 3,
            "remaining_actions": [
                "commit_and_push_stage4_publication_candidate",
                "pass_exact_clean_committed_pressure_gate_worker",
                "close_independent_security_integration_and_science_rereviews",
                "prepare_machine_reanalysis_from_authenticated_committed_source_snapshot",
                "stop_for_explicit_post_reanalysis_human_pressure_selection",
                "seal_human_selection_and_publish_from_same_authenticated_source_snapshot",
                "independently_verify_receipt_seal_attestation_guard_and_active_state",
            ],
            "authority": "none",
        },
    ):
        raise ValueError("Stage-4 pressure-selection packet-gate state drifted")


def _validation_manifest_schema() -> dict[str, object]:
    return json.loads(
        (SCHEMA_DIR / "validation_manifest.schema.json").read_text(
            encoding="utf-8"
        )
    )


def _minimum_validation_manifest() -> dict[str, object]:
    return {
        "schema_version": 1,
        "manifest_id": "host-scaffold-001",
        "created_utc": "2026-05-30T12:00:00Z",
        "claim_ids": ["CLAIM-PAPER-GYRO-001"],
        "test_id": "pic_relativistic_gyro_paper",
        "evidence_class": "unit/regression",
        "physical_mode": "paper_mhd_pic",
        "git": {
            "commit": "0" * 40,
            "tree": "0" * 40,
            "status": [],
            "source_archive": {
                "path": "source.tar",
                "sha256": "0" * 64,
            },
            "source_commit": {
                "path": "source.commit",
                "sha256": "0" * 64,
            },
            "source_bundle_sha256": "0" * 64,
            "submodule_status": "absent",
            "submodules": [],
            "clean_candidate_manifest": {
                "path": "clean_candidate_manifest.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_profile": {
                "path": "build_profile.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_profile_receipt": {
                "path": "profile_receipt.json",
                "sha256": "0" * 64,
            },
            "clean_candidate_build_provenance": {
                label: {
                    "path": f"build_provenance/{filename}",
                    "sha256": "0" * 64,
                }
                for label, filename in {
                    "configure_log": "configure.log",
                    "build_log": "build.log",
                    "cmake_cache": "CMakeCache.txt",
                    "module_list": "modules.txt",
                    "toolchain": "toolchain.txt",
                    "build_invocations": "build-invocations.json",
                    "git_status_preconfigure": "git_status.preconfigure.txt",
                    "git_status": "git_status.txt",
                    "submodule_status": "submodule_status.txt",
                    "environment_allowlist": "environment.allowlist.txt",
                    "build_environment": "build-environment.json",
                }.items()
            },
        },
        "authorization": {
            "control_plane_version": "0" * 64,
            "build_profile_control_plane_version": "0" * 64,
            "clean_candidate_manifest_path": (
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/"
                "clean_candidates/00000000-0000-0000-0000-000000000000/"
                "clean_candidate_manifest.json"
            ),
            "clean_candidate_manifest_sha256": "0" * 64,
            "active_policy": {
                "path": "active_policy.json",
                "sha256": "0" * 64,
            },
            "active_promotion": {
                "path": "active_promotion.json",
                "sha256": "0" * 64,
            },
        },
        "executable": {
            "path": "/tmp/athena",
            "sha256": "0" * 64,
            "cmake_cache": {"path": "/tmp/CMakeCache.txt", "sha256": "0" * 64},
            "modules": {"path": "/tmp/modules.txt", "sha256": "0" * 64},
            "environment_allowlist": {
                "path": "/tmp/environment.txt",
                "sha256": "0" * 64,
            },
        },
        "parameters": {},
        "oracle": {
            "kind": "analytic",
            "reference": "bounded host scaffold",
            "tolerances": {"relative_error": 1.0e-6},
        },
        "metrics": [{"name": "relative_error", "value": 0.0}],
        "resources": {
            "platform": "host",
            "artifact_root": "/tmp/pic-readiness",
        },
        "artifacts": [{"path": "metrics.json", "sha256": "0" * 64}],
        "review": {
            "reviewer": "pending external review",
            "disposition": "pending external review",
        },
    }


def _replace_nested(value: dict[str, object], path: tuple[object, ...],
                    replacement: object) -> None:
    target = value
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = replacement


class PicReadinessRegistryTests(unittest.TestCase):
    def test_json_documents_parse(self) -> None:
        paths = sorted(READINESS_DIR.glob("*.json"))
        paths += sorted(SCHEMA_DIR.glob("*.json"))
        self.assertTrue(paths)
        for path in paths:
            with self.subTest(path=path):
                json.loads(path.read_text(encoding="utf-8"))

    def test_storage_policy_records_authorized_frontier_boundary(self) -> None:
        policy = _load("storage_policy.json")
        frontier = policy["frontier"]
        storage = policy["olcf_side_storage"]
        long_term = policy["long_term_storage"]
        self.assertEqual(frontier["maximum_node_hours"], 10000)
        self.assertEqual(frontier["partition"], "batch")
        self.assertTrue(frontier["serial_pic_submissions"])
        self.assertEqual(
            frontier["simulation_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        self.assertFalse(storage["ledger_genesis_allowed"])
        self.assertEqual(storage["ledger_genesis"]["status"], "initialized")
        self.assertEqual(
            storage["ledger_genesis"]["mirror_transport"],
            "filesystem_copy",
        )
        self.assertEqual(
            storage["orion_bulk_evidence_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC",
        )
        self.assertEqual(
            storage["project_home_retention_role"],
            "operational_ledger_mirror_only",
        )
        source_alias_candidate = _load(
            "q027_control_plane_source_alias_hardening_candidate_2026-05-30.json"
        )
        recovery_candidate = _load(
            "q027_frontier_f0_purged_submission_recovery_activation_2026-05-30.json"
        )
        candidate = _load(
            "q027_frontier_f0_compute_snapshot_activation_2026-05-30.json"
        )
        f1_candidate = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        manual_accounting_activation = _load(
            "q027_manual_frontier_accounting_activation_2026-05-30.json"
        )
        self.assertEqual(
            candidate["predecessor"]["control_plane_version"],
            recovery_candidate["active_successor"]["control_plane_version"],
        )
        self.assertEqual(
            f1_candidate["initial_live_predecessor"]["control_plane_version"],
            candidate["active_successor"]["control_plane_version"],
        )
        lifecycle = storage["installed_control_plane_lifecycle"]
        science_freeze = policy["science_submission_freeze"]
        if science_freeze == {"status": "pending_clean_candidate_freeze"}:
            phase0_successor = _load(
                "phase0_curated_candidate_successor_v14_2026-06-01.json"
            )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                phase0_successor["successor_source_control_plane_version"],
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                phase0_successor["successor_source_control_plane_version"],
            )
            self.assertEqual(policy["registered_science_slices"], [])
            return
        d720_promotion = _load(
            "phase0_clean_candidate_freeze_and_policy_promotion_"
            "successor_v2_2026-06-02.json"
        )
        c83e_promotion = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        strict_q011_promotion = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v6_2026-06-02.json"
        )
        q011_pressure_promotion = _load(
            "phase0_curated_candidate_successor_v20_2026-06-02.json"
        )
        if (
            science_freeze
            == q011_pressure_promotion["policy_promotion"][
                "science_submission_freeze"
            ]
        ):
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                q011_pressure_promotion["live_paired_control_plane_version"],
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                q011_pressure_promotion["live_paired_control_plane_version"],
            )
            self.assertEqual(
                [
                    record["authorization_id"]
                    for record in policy["registered_science_slices"]
                ],
                q011_pressure_promotion["policy_promotion"][
                    "registered_science_authorization_ids"
                ],
            )
            return
        if (
            science_freeze
            == d720_promotion["active_policy_promotion"]["science_submission_freeze"]
        ):
            expected_control_plane_version = d720_promotion["control_plane_version"]
            if (
                storage["installed_control_plane_version"]
                == c83e_promotion["control_plane_version"]
            ):
                expected_control_plane_version = c83e_promotion[
                    "control_plane_version"
                ]
                self.assertEqual(
                    science_freeze,
                    c83e_promotion["active_policy_promotion"][
                        "science_submission_freeze"
                    ],
                )
            if (
                storage["installed_control_plane_version"]
                == strict_q011_promotion["control_plane_version"]
            ):
                expected_control_plane_version = strict_q011_promotion[
                    "control_plane_version"
                ]
                self.assertEqual(
                    science_freeze,
                    strict_q011_promotion["active_policy_promotion"][
                        "science_submission_freeze"
                    ],
                )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                expected_control_plane_version,
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                expected_control_plane_version,
            )
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                d720_promotion["active_policy_promotion"][
                    "frontier_launch_authorization"
                ],
                "none_no_registered_science_slices",
            )
            return
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        if (
            storage["installed_control_plane_version"]
            == replay["control_plane_version"]
        ):
            accounting = _load(
                "phase0_scheduler_accounting_controller_successor_2026-06-01.json"
            )
            freeze = _load(
                "phase0_clean_candidate_freeze_and_policy_promotion_2026-06-01.json"
            )
            self.assertEqual(lifecycle, "paired_installed_reviewed_generation")
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            self.assertEqual(
                replay["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_scheduler_accounting_controller_successor_2026-06-01.json",
            )
            self.assertEqual(
                accounting["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_clean_candidate_freeze_and_policy_promotion_2026-06-01.json",
            )
            self.assertEqual(
                science_freeze["manifest_path"],
                freeze["clean_candidate"]["manifest_path"],
            )
            self.assertEqual(
                science_freeze["manifest_sha256"],
                freeze["clean_candidate"]["manifest_sha256"],
            )
            self.assertEqual(
                {
                    record["authorization_id"]: {
                        "campaign": record["campaign"],
                        "launch_contract_sha256": record["launch_contract_sha256"],
                    }
                    for record in policy["registered_science_slices"]
                },
                {
                    record["authorization_id"]: {
                        "campaign": record["campaign"],
                        "launch_contract_sha256": record["launch_contract_sha256"],
                    }
                    for record in replay["registered_science_slices"]
                },
            )
            terminal = replay["terminal_mirrored_ledger"]
            self.assertEqual(terminal["orion_ledger_records"], 66)
            self.assertEqual(
                terminal["orion_ledger_records"],
                terminal["project_home_ledger_records"],
            )
            self.assertEqual(
                terminal["orion_ledger_records"],
                terminal["orion_receipt_records"],
            )
            self.assertEqual(terminal["active_reservations"], 0)
            self.assertEqual(terminal["pending_submission_marker"], "absent")
            self.assertEqual(terminal["pending_manual_accounting_marker"], "absent")
            self.assertEqual(
                long_term["status"],
                "user_selected_orion_only_with_documented_durability_risk",
            )
            return
        if lifecycle == "live_active_generation_successor_staged_not_installed":
            prior_active = f1_candidate.get(
                "active_policy_transition",
                f1_candidate.get("prior_active_policy_transition"),
            )
            expected_active_version = (
                prior_active["control_plane_version"]
                if prior_active is not None
                else candidate["active_successor"]["control_plane_version"]
            )
            self.assertEqual(
                storage["installed_control_plane_version"],
                expected_active_version,
            )
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                f1_candidate["staged_successor"]["control_plane_version"],
            )
            self.assertNotEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            self.assertEqual(
                f1_candidate["paired_install_transition"]["status"],
                "pending_clean_commit_and_paired_install",
            )
            self.assertEqual(
                f1_candidate["staged_successor"]["paired_install_status"],
                f1_candidate["paired_install_transition"]["status"],
            )
        elif lifecycle == "paired_installed_reviewed_generation":
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            current_transition = manual_accounting_activation[
                "control_plane_transition"
            ]
            self.assertEqual(
                storage["installed_control_plane_version"],
                current_transition["control_plane_version"],
            )
            terminal_ledger = manual_accounting_activation["terminal_ledger"]
            self.assertEqual(terminal_ledger["orion_records"], 56)
            self.assertEqual(
                terminal_ledger["orion_records"],
                terminal_ledger["project_home_records"],
            )
            self.assertEqual(
                terminal_ledger["orion_records"],
                terminal_ledger["mirror_receipts"],
            )
            self.assertEqual(terminal_ledger["currently_reserved_node_hours"], 0.0)
            self.assertEqual(
                manual_accounting_activation["authorization"][
                    "scientific_evidence_eligible"
                ],
                False,
            )
            paired = f1_candidate["paired_install_transition"]
            self.assertEqual(
                storage["ledger_genesis"]["event_sha256"],
                paired["genesis_event_sha256"],
            )
            self.assertEqual(
                storage["ledger_genesis"]["mirror_ack_sha256"],
                paired["genesis_mirror_ack_sha256"],
            )
            self.assertEqual(
                paired["orion_ledger_records"],
                paired["project_home_ledger_records"],
            )
            self.assertEqual(paired["active_reservations"], 0)
            active_transition = f1_candidate.get("active_policy_transition")
            if active_transition is None:
                self.assertEqual(
                    f1_candidate["staged_successor"]["active_policy_promotion_status"],
                    "pending",
                )
            else:
                self.assertEqual(active_transition["status"], "pass")
            self.assertEqual(science_freeze["status"], "authorized")
            clean_candidate = candidate["inherited_clean_candidate_freeze"]
            clean_candidate_transition = source_alias_candidate[
                "clean_candidate_policy_transition"
            ]
            self.assertEqual(
                clean_candidate_transition["policy_sha256"],
                "e8909bf541d1c69d4d19c495ebfe983934291de37d6d1423ae4cbfca2e4bb155",
            )
            self.assertEqual(
                clean_candidate_transition["science_submission_freeze_status"],
                science_freeze["status"],
            )
            for key in ["manifest_path", "manifest_sha256"]:
                self.assertEqual(science_freeze[key], clean_candidate_transition[key])
                self.assertEqual(clean_candidate_transition[key], clean_candidate[key])
            self.assertEqual(clean_candidate_transition["orion_ledger_records"], 22)
            self.assertEqual(
                clean_candidate_transition["orion_ledger_records"],
                clean_candidate_transition["project_home_ledger_records"],
            )
            self.assertEqual(
                clean_candidate_transition["orion_ledger_records"],
                clean_candidate_transition["orion_receipt_records"],
            )

            self.assertEqual(clean_candidate_transition["active_reservations"], 0)
            admission_activation = _load(
                "q027_frontier_f0_clean_candidate_admission_policy_activation_2026-05-30.json"
            )
            self.assertEqual(
                admission_activation["control_plane_version"],
                recovery_candidate["predecessor"]["control_plane_version"],
            )
            self.assertEqual(
                admission_activation["policy_sha256"],
                recovery_candidate["predecessor"]["policy_sha256"],
            )
            self.assertEqual(
                admission_activation["science_submission_freeze"],
                {
                    key: science_freeze[key]
                    for key in ["status", "manifest_path", "manifest_sha256"]
                } | {
                    "executable_sha256": clean_candidate["executable_sha256"],
                },
            )
            self.assertEqual(
                admission_activation["frontier_admission_smoke"]["status"],
                "authorized_f0_parser_contract_only",
            )
            self.assertEqual(
                policy["frontier_admission_smoke"], {"status": "closed_after_pass"}
            )
            admission_ledger = admission_activation["ledger_validation"]
            self.assertEqual(admission_ledger["orion_ledger_records"], 22)
            self.assertEqual(
                admission_ledger["orion_ledger_records"],
                admission_ledger["project_home_ledger_records"],
            )
            self.assertEqual(
                admission_ledger["orion_ledger_records"],
                admission_ledger["orion_receipt_records"],
            )
            self.assertEqual(admission_ledger["active_reservations"], 0)
            self.assertEqual(admission_ledger["pending_submission_marker"], "absent")
            recovery_transition = candidate["active_policy_transition"]
            self.assertEqual(
                recovery_transition["control_plane_version"],
                candidate["active_successor"]["control_plane_version"],
            )
            self.assertEqual(
                recovery_transition["policy_sha256"],
                "cfe6610a4f7f54b153f2e20b30397d0417a2634b32997117f6d14a2eb4b771cc",
            )
            self.assertEqual(recovery_transition["orion_ledger_records"], 31)
            self.assertEqual(
                recovery_transition["orion_ledger_records"],
                recovery_transition["project_home_ledger_records"],
            )
            self.assertEqual(
                recovery_transition["orion_ledger_records"],
                recovery_transition["orion_receipt_records"],
            )
            self.assertEqual(recovery_transition["active_reservations"], 0)
            self.assertEqual(
                recovery_transition["pending_submission_marker"], "absent"
            )
        else:
            self.fail(f"Unknown installed-control-plane lifecycle: {lifecycle}")
        self.assertEqual(
            long_term["status"],
            "user_selected_orion_only_with_documented_durability_risk",
        )

    def test_registered_science_launch_contract_sidecars_match_policy(self) -> None:
        policy = _load("storage_policy.json")
        authorizations = {
            record["authorization_id"]: record
            for record in policy["registered_science_slices"]
        }
        sidecars = {
            "f1_gpu_relativistic_gyro": "frontier_f1_clean_gyro_launch_contract.json",
            "f1_gpu_paper_coupling": (
                "frontier_f1_clean_paper_coupling_launch_contract.json"
            ),
            "f2_multirank_runtime_metadata": (
                "frontier_f2_multirank_runtime_metadata_launch_contract.json"
            ),
        }
        if not authorizations:
            self.assertIn(
                policy["science_submission_freeze"]["status"],
                ("pending_clean_candidate_freeze", "authorized"),
            )
            self.assertEqual(authorizations, {})
            for campaign, filename in sidecars.items():
                with self.subTest(campaign=campaign):
                    validate_launch_contract(_load(filename))
            return
        authorizations_by_campaign = {
            record["campaign"]: record
            for record in policy["registered_science_slices"]
        }
        q011_sidecars = {
            case.campaign: case.launch_contract_path.name
            for case in pressure_execution.CASES
        }
        if set(authorizations_by_campaign) == set(q011_sidecars):
            for campaign, filename in q011_sidecars.items():
                authorization = authorizations_by_campaign[campaign]
                with self.subTest(campaign=campaign):
                    contract = _load(filename)
                    self.assertEqual(
                        pressure_execution._launch_contract_sha256(contract),
                        authorization["launch_contract_sha256"],
                    )
            return
        self.assertEqual(set(authorizations_by_campaign), set(sidecars))
        for campaign, filename in sidecars.items():
            authorization = authorizations_by_campaign[campaign]
            authorization_id = authorization["authorization_id"]
            with self.subTest(authorization_id=authorization_id):
                contract = _load(filename)
                validate_launch_contract(contract)
                self.assertEqual(
                    launch_contract_sha256(contract),
                    authorization["launch_contract_sha256"],
                )

    def test_q011_pressure_pilot_registered_execution_source_tranche_is_frozen(
        self,
    ) -> None:
        record = _load(
            "q011_section54_pressure_pilot_registered_execution_preregistration_2026-06-02.json"
        )
        boundary = record["execution_boundary"]
        self.assertFalse(boundary["frontier_execution_authorized_by_this_record"])
        self.assertFalse(boundary["scheduler_commands_authorized_by_this_record"])
        self.assertFalse(boundary["storage_policy_mutation_authorized_by_this_record"])
        self.assertEqual(
            record["source_bindings"]["materializer"]["sha256"],
            "57f6300b4d9a7cdf2327dc92137c9dcb040eeb815d5c79f33f01e646703b6ca4",
        )
        self.assertEqual(
            [entry["authorization_id"] for entry in record["launch_matrix"]],
            [
                f"q011-section54-pressure-{case.case_id.replace('_', '-')}-v1"
                for case in pressure_execution.CASES
            ],
        )

    def test_q011_pressure_pilot_registered_execution_retry_successor_is_frozen(
        self,
    ) -> None:
        record = _load(
            "q011_section54_pressure_pilot_registered_execution_retry_"
            "successor_v2_2026-06-02.json"
        )
        self.assertEqual(
            record["predecessor_sha256"],
            _sha256(REPO_ROOT / record["predecessor_record"]),
        )
        chronology = record["failed_attempt_chronology"]
        self.assertEqual(chronology["job_id"], "4754211")
        self.assertEqual(chronology["terminal_state"], "FAILED")
        self.assertEqual(
            chronology["terminal_reconciliation_event_sha256"],
            "fc0082bef800733c433395d48f551559abe083367e5549cc4bffbe8a0ab48bfa",
        )
        boundary = record["execution_boundary"]
        self.assertFalse(boundary["frontier_execution_authorized_by_this_record"])
        self.assertFalse(boundary["scheduler_commands_authorized_by_this_record"])
        self.assertFalse(boundary["storage_policy_mutation_authorized_by_this_record"])
        self.assertEqual(
            record["source_bindings"],
            pressure_execution._historical_v2_preregistration()["source_bindings"],
        )
        status = pressure_execution.historical_v2_source_tranche_status()
        self.assertEqual(status["state"], "historical_consumed_slice_non_authorizing")
        self.assertFalse(status["source_bindings_match_current_checkout"])
        self.assertFalse(status["consumed_slice_reauthorization_allowed"])
        with self.assertRaisesRegex(
            pressure_execution.ContractError,
            "historical v2 registered-execution tranche is consumed",
        ):
            pressure_execution.validate_source_tranche()
        self.assertEqual(len(record["launch_matrix"]), len(pressure_execution.CASES))
        binding_by_case = {
            binding["case_id"]: binding
            for binding in record["source_bindings"]["launch_contracts"]
        }
        for case in pressure_execution.CASES:
            with self.subTest(case_id=case.case_id):
                contract = _load(case.launch_contract_path.name)
                self.assertEqual(
                    contract, pressure_execution.expected_launch_contract(case)
                )
                validate_launch_contract(contract)
                binding = binding_by_case[case.case_id]
                self.assertEqual(binding["file_sha256"], _sha256(case.launch_contract_path))
                self.assertEqual(
                    binding["launch_contract_sha256"],
                    launch_contract_sha256(contract),
                )

    def test_f2_multirank_runtime_metadata_candidate_resolves_source_commit(self) -> None:
        candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        commit = candidate["implementation_source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        binding = candidate["registered_science_slice"]
        self.assertEqual(
            candidate["staged_policy_sha256"],
            candidate["accepted_v2_execution"]["active_policy_sha256"],
        )
        for digest_key, relative_path in {
            "job_script_sha256":
                "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh",
            "input_deck_sha256": "inputs/tests/pic_parser_contract_guards.athinput",
            "analysis_script_sha256":
                "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
            "analysis_support_sha256":
                "tst/publication/frontier_f1_structured_artifacts.py",
        }.items():
            self.assertEqual(binding[digest_key], _git_blob_sha256(commit, relative_path))
        contract = _load("frontier_f2_multirank_runtime_metadata_launch_contract.json")
        self.assertEqual(
            binding["launch_contract_sha256"],
            launch_contract_sha256(contract),
        )

    def _assert_current_registered_replay_bindings(
        self, policy: dict[str, object], staged_version: str
    ) -> None:
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        accounting = _load(
            "phase0_scheduler_accounting_controller_successor_2026-06-01.json"
        )
        closure = _load(
            "phase0_registered_prerequisite_replay_closure_2026-06-01.json"
        )
        storage = policy["olcf_side_storage"]
        self.assertEqual(staged_version, replay["control_plane_version"])
        self.assertEqual(
            storage["installed_control_plane_version"], staged_version
        )
        self.assertEqual(
            storage["staged_control_plane_candidate_version"], staged_version
        )
        promotion = replay["active_policy_promotion"]
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(READINESS_DIR / "storage_policy.json"),
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(Path(promotion["orion_policy_path"])),
        )
        self.assertEqual(
            promotion["project_home_policy_sha256"],
            _sha256(Path(promotion["project_home_policy_path"])),
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            _sha256(
                Path(promotion["orion_policy_path"]).with_name(
                    "active_promotion.json"
                )
            ),
        )
        self.assertEqual(
            promotion["project_home_promotion_sha256"],
            _sha256(
                Path(promotion["project_home_policy_path"]).with_name(
                    "active_promotion.json"
                )
            ),
        )
        clean_manifest_path = Path(
            policy["science_submission_freeze"]["manifest_path"]
        )
        self.assertEqual(
            _sha256(clean_manifest_path),
            replay["clean_candidate"]["manifest_sha256"],
        )
        clean_manifest = json.loads(clean_manifest_path.read_text(encoding="utf-8"))
        executable_path = Path(clean_manifest["build"]["executable_path"])
        self.assertEqual(
            _sha256(executable_path),
            replay["clean_candidate"]["executable_sha256"],
        )
        binding_paths = {
            "f1_gpu_relativistic_gyro": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_relativistic_gyro_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_relativistic_gyro_paper.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
            "f1_gpu_paper_coupling": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_paper_coupling_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_paper_coupling_conservation.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
            "f2_multirank_runtime_metadata": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_parser_contract_guards.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
        }
        expected_metadata = {
            "f1_gpu_relativistic_gyro": {
                "test_id": "pic_relativistic_gyro_paper",
                "evidence_class": "frontier_f1_clean_candidate_gpu_pusher_oracle",
                "physical_mode": "paper_test_particle",
            },
            "f1_gpu_paper_coupling": {
                "test_id": "pic_paper_coupling_conservation",
                "evidence_class": "frontier_f1_clean_candidate_gpu_paper_coupling_oracle",
                "physical_mode": "paper_mhd_pic",
            },
            "f2_multirank_runtime_metadata": {
                "test_id": "pic_parser_contract_guards",
                "evidence_class": "frontier_f2_clean_candidate_multirank_runtime_metadata",
                "physical_mode": "extended_mhd_pic_parser_contract",
            },
        }
        authorizations = {
            record["campaign"]: record
            for record in policy["registered_science_slices"]
        }
        replay_slices = {
            record["campaign"]: record for record in replay["registered_science_slices"]
        }
        self.assertEqual(set(authorizations), set(binding_paths))
        self.assertEqual(set(replay_slices), set(binding_paths))
        environment_path = CONTROL_PLANE_DIR / "frontier_pic_environment.sh"
        for campaign, authorization in authorizations.items():
            with self.subTest(campaign=campaign):
                paths = binding_paths[campaign]
                for key, expected in expected_metadata[campaign].items():
                    self.assertEqual(authorization[key], expected)
                self.assertEqual(authorization["runtime_profile"], "frontier_minimum_supported")
                self.assertEqual(authorization["selected_qos"], "debug")
                self.assertIs(authorization["registered_short_nonproduction"], True)
                self.assertEqual(authorization["maximum_nodes"], 1)
                self.assertEqual(authorization["maximum_walltime_seconds"], 900)
                self.assertEqual(authorization["maximum_attempts"], 1)
                self.assertEqual(
                    authorization["authorization_id"],
                    replay_slices[campaign]["authorization_id"],
                )
                self.assertEqual(
                    authorization["launch_contract_sha256"],
                    replay_slices[campaign]["launch_contract_sha256"],
                )
                self.assertEqual(
                    authorization["job_script_sha256"],
                    _sha256(paths["job_script_sha256"]),
                )
                self.assertEqual(
                    authorization["input_deck_sha256"],
                    _sha256(paths["input_deck_sha256"]),
                )
                self.assertEqual(
                    authorization["environment_profile_sha256"],
                    _sha256(environment_path),
                )
                self.assertEqual(
                    authorization["analysis_script_sha256"],
                    [_sha256(path) for path in paths["analysis_script_sha256"]],
                )
                self.assertEqual(
                    authorization["clean_candidate_manifest_sha256"],
                    _sha256(clean_manifest_path),
                )
                self.assertEqual(
                    authorization["executable_sha256"], _sha256(executable_path)
                )
        executions = {
            record["campaign"]: record for record in closure["registered_replays"]
        }
        self.assertEqual(set(executions), set(binding_paths))
        ledger_records = [
            json.loads(line)
            for line in Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/"
                "node_hours.jsonl"
            )
            .read_text(encoding="utf-8")
            .splitlines()
        ]
        terminal_events = {
            record["job_id"]: record
            for record in ledger_records
            if record.get("event_type") == "reconciliation"
            and record.get("job_id") in {
                execution["job_id"] for execution in executions.values()
            }
        }
        self.assertEqual(len(terminal_events), len(executions))
        for campaign, execution in executions.items():
            with self.subTest(campaign=campaign):
                self.assertEqual(
                    execution["authorization_id"],
                    authorizations[campaign]["authorization_id"],
                )
                manifest_path = Path(execution["manifest_path"])
                artifact_dir = Path(execution["artifact_dir"])
                self.assertEqual(_sha256(manifest_path), execution["manifest_sha256"])
                self.assertEqual(
                    _sha256(artifact_dir / "artifact_inventory.json"),
                    execution["artifact_inventory_sha256"],
                )
                result_path = artifact_dir / "analysis" / "analysis.json"
                self.assertEqual(
                    _sha256(result_path), execution["analysis_result_sha256"]
                )
                result = json.loads(result_path.read_text(encoding="utf-8"))
                self.assertEqual(result["status"], "pass")
                receipt_path = (
                    artifact_dir / "analysis" / "offline_analysis_receipt.json"
                )
                self.assertEqual(
                    _sha256(receipt_path),
                    execution["offline_analysis_receipt_sha256"],
                )
                receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
                self.assertEqual(
                    receipt["artifact_inventory"]["sha256"],
                    execution["artifact_inventory_sha256"],
                )
                self.assertEqual(
                    receipt["analysis_result"]["sha256"],
                    execution["analysis_result_sha256"],
                )
                for attestation_key in (
                    "pre_manifest_attestation",
                    "pre_submit_wrapper_attestation",
                ):
                    attestation = execution[attestation_key]
                    self.assertEqual(
                        _sha256(Path(attestation["path"])),
                        attestation["sha256"],
                    )
                terminal = terminal_events[execution["job_id"]]
                self.assertEqual(terminal["state"], "COMPLETED")
                self.assertEqual(
                    terminal["event_sha256"],
                    execution["terminal_ledger_event_sha256"],
                )
                self.assertEqual(
                    terminal["reservation_id"], execution["reservation_id"]
                )
                self.assertEqual(
                    terminal["submission_id"], execution["submission_id"]
                )
        ledger = closure["terminal_mirrored_ledger"]
        self.assertEqual(
            _sha256(Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl")),
            ledger["orion_node_hours_jsonl_sha256"],
        )
        self.assertEqual(
            _sha256(Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl")),
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        self.assertEqual(len(ledger_records), ledger["orion_ledger_records"])
        self.assertEqual(
            accounting["manual_accounting_activation"]["reviewed_job_count"], 10
        )

    def test_phase0_successor_v7_binds_completed_replay_boundary(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v7_2026-06-01.json")
        self.assertEqual(
            successor["status"],
            "canonical_clean_candidate_frozen_accounting_closed_"
            "registered_prerequisite_replays_authorized",
        )
        self.assertEqual(successor["remaining_phase0_actions"], [])
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v6_2026-06-01.json",
        )
        for key in (
            "clean_candidate_freeze_receipt",
            "scheduler_accounting_successor_receipt",
            "registered_prerequisite_replay_policy_promotion_receipt",
        ):
            receipt = successor[key]
            self.assertEqual(
                receipt["sha256"],
                _sha256(REPO_ROOT / receipt["path"]),
            )

    def test_phase0_successor_v8_binds_terminal_replay_closure(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v8_2026-06-01.json")
        self.assertEqual(
            successor["status"],
            "canonical_clean_candidate_frozen_registered_prerequisite_replays_pass",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v7_2026-06-01.json",
        )
        receipt = successor["registered_prerequisite_replay_closure_receipt"]
        self.assertEqual(receipt["sha256"], _sha256(REPO_ROOT / receipt["path"]))

    def test_phase0_successor_v15_binds_d720_policy_promotion(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v15_2026-06-02.json")
        self.assertEqual(
            successor["status"],
            "canonical_vl2_tsc_clean_candidate_authorized_policy_promoted_"
            "launch_prohibited_pending_registered_science_slices",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v14_2026-06-01.json",
        )
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        receipt_binding = successor[
            "clean_candidate_freeze_and_policy_promotion_receipt"
        ]
        self.assertEqual(
            receipt_binding["sha256"],
            _sha256(REPO_ROOT / receipt_binding["path"]),
        )
        receipt = json.loads(
            (REPO_ROOT / receipt_binding["path"]).read_text(encoding="utf-8")
        )
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        promotion = receipt["active_policy_promotion"]
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["frontier_launch_authorization"],
            "none_no_registered_science_slices",
        )
        self.assertEqual(
            promotion["repo_policy_sha256"],
            promotion["orion_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        attestation = receipt["pre_promotion_operator_attestation"]
        self.assertEqual(_sha256(Path(attestation["path"])), attestation["sha256"])
        ledger = receipt["mirrored_ledger_invariant"]
        self.assertEqual(
            ledger["orion_node_hours_jsonl_sha256"],
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        for key in (
            "orion_node_hours_jsonl_sha256",
            "project_home_node_hours_jsonl_sha256",
            "orion_node_hours_csv_sha256",
            "orion_mirror_receipts_jsonl_sha256",
        ):
            self.assertRegex(ledger[key], r"^[0-9a-f]{64}$")
        self.assertEqual(ledger["orion_ledger_records"], 108)
        self.assertEqual(ledger["active_reservations"], 0)
        self.assertEqual(receipt["external_review"]["reviewer"], "pending external review")

    def test_phase0_paired_successor_v4_binds_c83e_policy_promotion(self) -> None:
        receipt = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        self.assertEqual(
            receipt["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v3_2026-06-01.json",
        )
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        installed = receipt["paired_install"]
        for key in ("orion", "project_home"):
            self.assertEqual(
                installed[f"{key}_inventory_sha256"],
                _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
            )
        promotion = receipt["active_policy_promotion"]
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["frontier_launch_authorization"],
            "none_no_registered_science_slices",
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        attestation = receipt["pre_promotion_operator_attestation"]
        self.assertEqual(_sha256(Path(attestation["path"])), attestation["sha256"])
        ledger = receipt["terminal_mirrored_ledger"]
        self.assertEqual(
            ledger["orion_node_hours_jsonl_sha256"],
            ledger["project_home_node_hours_jsonl_sha256"],
        )
        for key in (
            "orion_node_hours_jsonl_sha256",
            "project_home_node_hours_jsonl_sha256",
            "orion_node_hours_csv_sha256",
            "orion_mirror_receipts_jsonl_sha256",
        ):
            self.assertRegex(ledger[key], r"^[0-9a-f]{64}$")
        self.assertEqual(ledger["ledger_records"], 108)
        self.assertEqual(ledger["active_reservations"], 0)

    def test_phase0_successor_v16_binds_c83e_transition(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v16_2026-06-02.json")
        self.assertEqual(
            successor["status"],
            "canonical_vl2_tsc_clean_candidate_authorized_successor_"
            "control_plane_promoted_launch_prohibited_pending_registered_"
            "science_slices",
        )
        self.assertEqual(
            successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v15_2026-06-02.json",
        )
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        for key in (
            "paired_control_plane_install_and_policy_promotion_receipt",
            "clean_candidate_freeze_and_policy_promotion_receipt",
        ):
            receipt = successor[key]
            self.assertEqual(
                receipt["sha256"],
                _sha256(REPO_ROOT / receipt["path"]),
            )
        baseline = successor["operational_baseline"]
        self.assertEqual(baseline["registered_science_slices"], [])
        self.assertEqual(
            successor["frontier_launch_authorization"],
            "none_no_registered_science_slices",
        )

    def test_registered_science_staged_bindings_recompute_from_exact_files(self) -> None:
        policy = _load("storage_policy.json")
        storage = policy["olcf_side_storage"]
        staged_version = inventory_digest(
            [
                {"path": name, "sha256": _sha256(CONTROL_PLANE_DIR / name)}
                for name in CONTROL_PLANE_FILES
            ]
        )
        current_repair = _load(
            "q011_section54_sixteenth_lustre_publication_rename_compatibility_"
            "transition_2026-06-05.json"
        )
        acceptance_repair = _load(
            "q011_section54_thirteenth_acceptance_root_setgid_repair_transition_"
            "2026-06-03.json"
        )
        repaired_staged = current_repair["repaired_staged_control_plane"]
        helper = acceptance_repair["acceptance_root_helper_repair"]
        self.assertEqual(
            current_repair["predecessor_sha256"],
            _sha256(REPO_ROOT / current_repair["predecessor_record"]),
        )
        self.assertGreater(
            _canonical_utc_second(current_repair["exact_patch_frozen_utc"]),
            _canonical_utc_second(current_repair["recorded_utc"]),
        )
        self.assertEqual(helper["sha256"], _sha256(REPO_ROOT / helper["path"]))
        runbook = (CONTROL_PLANE_DIR / "README.md").read_text(encoding="utf-8")
        recovery_receipt = json.loads(
            Path(
                current_repair["recovered_acceptance_root"]["recovery_receipt"]["path"]
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(
            helper["sha256"],
            _git_blob_sha256(recovery_receipt["validated_source_commit"], helper["path"]),
        )
        self.assertEqual(
            runbook.count(
                f"EXPECTED_ACCEPTANCE_HELPER_SHA256={helper['sha256']}"
            ),
            2,
        )
        self.assertEqual(
            runbook.count(
                "ACCEPTANCE_RECOVERY_VALIDATED_SOURCE_COMMIT="
                f"{recovery_receipt['validated_source_commit']}"
            ),
            1,
        )
        checkpoint = current_repair["recovered_acceptance_root"]
        acceptance_root = Path(checkpoint["path"])
        acceptance_status = acceptance_root.stat()
        self.assertTrue(stat.S_ISDIR(acceptance_status.st_mode))
        self.assertEqual(
            (acceptance_status.st_dev, acceptance_status.st_ino),
            (
                _canonical_decimal_integer(checkpoint["device"]),
                _canonical_decimal_integer(checkpoint["inode"]),
            ),
        )
        self.assertEqual(acceptance_status.st_uid, checkpoint["uid"])
        self.assertEqual(acceptance_status.st_gid, checkpoint["gid"])
        self.assertEqual(
            f"{stat.S_IMODE(acceptance_status.st_mode):05o}",
            checkpoint["mode"],
        )
        expected_acceptance_seals = {
            ".q011_section54_pressure_pilot_bundle_receipt.json.publication-success": {
                "seal_sha256": (
                    "8eeb7eb24f9aa64652a27d619b155e959a302788eb99b8b4fa166052ae0158f8"
                ),
                "receipt_name": "q011_section54_pressure_pilot_bundle_receipt.json",
                "receipt_sha256": (
                    "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9"
                ),
            },
            (
                ".q011_section54_pressure_pilot_review_packet_receipt.json."
                "publication-success"
            ): {
                "seal_sha256": (
                    "109bc4522579feab92ef1778c0318a58e96aa8d059e481283c6faa83cd6c9e94"
                ),
                "receipt_name": (
                    "q011_section54_pressure_pilot_review_packet_receipt.json"
                ),
                "receipt_sha256": (
                    "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d"
                ),
            },
        }
        self.assertEqual(
            {path.name for path in acceptance_root.iterdir()},
            set(expected_acceptance_seals),
        )
        for seal_name, expected in expected_acceptance_seals.items():
            with self.subTest(acceptance_seal=seal_name):
                seal_path = acceptance_root / seal_name
                seal_status = seal_path.lstat()
                self.assertTrue(stat.S_ISREG(seal_status.st_mode))
                self.assertEqual(stat.S_IMODE(seal_status.st_mode), 0o444)
                self.assertEqual(_sha256(seal_path), expected["seal_sha256"])
                seal = json.loads(seal_path.read_text(encoding="utf-8"))
                published_receipt = (
                    acceptance_root.parent
                    / "publication"
                    / expected["receipt_name"]
                )
                published_receipt_status = published_receipt.lstat()
                published_root_status = published_receipt.parent.lstat()
                self.assertEqual(
                    seal,
                    {
                        "schema_version": 1,
                        "record_type": (
                            "q011_receipt_inode_bound_publication_success_seal"
                        ),
                        "receipt_name": expected["receipt_name"],
                        "receipt_sha256": expected["receipt_sha256"],
                        "receipt_identity": {
                            "device": published_receipt_status.st_dev,
                            "inode": published_receipt_status.st_ino,
                        },
                        "publication_root_identity": {
                            "device": published_root_status.st_dev,
                            "inode": published_root_status.st_ino,
                        },
                    },
                )
        self.assertEqual(sorted(os.listxattr(acceptance_root)), checkpoint["xattrs"])
        publication_checkpoint = current_repair["publication_root_authority"]
        publication_root = Path(publication_checkpoint["path"])
        publication_status = publication_root.stat()
        self.assertTrue(stat.S_ISDIR(publication_status.st_mode))
        self.assertEqual(
            (publication_status.st_dev, publication_status.st_ino),
            (
                _canonical_decimal_integer(publication_checkpoint["device"]),
                _canonical_decimal_integer(publication_checkpoint["inode"]),
            ),
        )
        self.assertEqual(publication_status.st_uid, publication_checkpoint["uid"])
        self.assertEqual(publication_status.st_gid, publication_checkpoint["gid"])
        self.assertEqual(
            f"{stat.S_IMODE(publication_status.st_mode):05o}",
            publication_checkpoint["mode"],
        )
        expected_publication_entries = {
            "q011_section54_pressure_pilot_analysis.json": {
                "kind": "file",
                "mode": 0o444,
                "sha256": (
                    "d55b4c2020716df899c86dfe5a9018d48194d63590616ff60541067243daacb7"
                ),
            },
            "q011_section54_pressure_pilot_bundle": {
                "kind": "directory",
                "mode": 0o2500,
            },
            "q011_section54_pressure_pilot_bundle_receipt.json": {
                "kind": "file",
                "mode": 0o444,
                "sha256": (
                    "9117b3dbc7573187b2d080568e69bdbbee0642f2a965aa543273ab3ea3d67be9"
                ),
            },
            "q011_section54_pressure_pilot_review_packet": {
                "kind": "directory",
                "mode": 0o2500,
            },
            "q011_section54_pressure_pilot_review_packet_receipt.json": {
                "kind": "file",
                "mode": 0o444,
                "sha256": (
                    "3f20d3d26a479aa508439f9d038ec6510643bf407aa081fae22959a57571de5d"
                ),
            },
        }
        self.assertEqual(
            {path.name for path in publication_root.iterdir()},
            set(expected_publication_entries),
        )
        for entry_name, expected in expected_publication_entries.items():
            with self.subTest(publication_entry=entry_name):
                entry = publication_root / entry_name
                entry_status = entry.lstat()
                if expected["kind"] == "file":
                    self.assertTrue(stat.S_ISREG(entry_status.st_mode))
                    self.assertEqual(_sha256(entry), expected["sha256"])
                else:
                    self.assertTrue(stat.S_ISDIR(entry_status.st_mode))
                self.assertEqual(
                    stat.S_IMODE(entry_status.st_mode),
                    expected["mode"],
                )
        published_aggregate_receipt = json.loads(
            (
                publication_root
                / "q011_section54_pressure_pilot_bundle_receipt.json"
            ).read_text(encoding="utf-8")
        )
        aggregate_bundle = published_aggregate_receipt["aggregate_bundle"]
        self.assertEqual(
            _sha256(Path(aggregate_bundle["path"]) / "pressure_pilot_manifest.json"),
            aggregate_bundle["manifest_sha256"],
        )
        published_packet_receipt = json.loads(
            (
                publication_root
                / "q011_section54_pressure_pilot_review_packet_receipt.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(
            _sha256(
                Path(published_packet_receipt["packet_root"])
                / "packet_inventory.json"
            ),
            published_packet_receipt["inventory_sha256"],
        )
        self.assertEqual(
            sorted(os.listxattr(publication_root)),
            publication_checkpoint["xattrs"],
        )
        self.assertEqual(
            [
                name
                for name in sorted(os.listxattr(publication_root))
                if name in {"system.posix_acl_access", "system.posix_acl_default"}
            ],
            publication_checkpoint["acl_xattrs"],
        )
        self.assertTrue(publication_checkpoint["same_account_isolated_parent"])
        self.assertFalse(publication_checkpoint["group_or_other_write_bits"])
        anchor_checkpoint = current_repair[
            "stable_account_serialization_anchor_authority"
        ]
        anchor = Path(anchor_checkpoint["path"])
        anchor_status = anchor.stat()
        self.assertTrue(stat.S_ISDIR(anchor_status.st_mode))
        self.assertEqual(
            (anchor_status.st_dev, anchor_status.st_ino),
            (
                _canonical_decimal_integer(anchor_checkpoint["device"]),
                _canonical_decimal_integer(anchor_checkpoint["inode"]),
            ),
        )
        self.assertEqual(anchor_status.st_uid, anchor_checkpoint["uid"])
        self.assertEqual(anchor_status.st_gid, anchor_checkpoint["gid"])
        self.assertEqual(
            f"{stat.S_IMODE(anchor_status.st_mode):05o}",
            anchor_checkpoint["mode"],
        )
        self.assertEqual(sorted(os.listxattr(anchor)), anchor_checkpoint["xattrs"])
        self.assertEqual(
            [
                name
                for name in sorted(os.listxattr(anchor))
                if name in {"system.posix_acl_access", "system.posix_acl_default"}
            ],
            anchor_checkpoint["acl_xattrs"],
        )
        self.assertTrue(anchor_checkpoint["outside_replaceable_pic_root_name"])
        self.assertIn(anchor, acceptance_root.parent.parents)
        completed_validation = current_repair["completed_fifteenth_clean_worker_validation"]
        self.assertEqual(
            completed_validation["archived_log_sha256"],
            _sha256(Path(completed_validation["archived_log_path"])),
        )
        failed_publication = current_repair["failed_closed_aggregate_publication"]
        self.assertEqual(failed_publication["scheduler_job_token"], "4766456;frontier")
        self.assertEqual(failed_publication["terminal_state"], "FAILED")
        self.assertEqual(failed_publication["exit_code"], "1:0")
        self.assertEqual(
            failed_publication["archived_log_sha256"],
            _sha256(Path(failed_publication["archived_log_path"])),
        )
        self.assertEqual(
            failed_publication["source_archive_sha256"],
            _git_archive_sha256(failed_publication["validated_source_commit"]),
        )
        self.assertTrue(failed_publication["publication_root_empty_after_failure"])
        self.assertTrue(
            failed_publication["publication_acceptance_root_empty_after_failure"]
        )
        receipt_checkpoint = checkpoint["recovery_receipt"]
        receipt = Path(receipt_checkpoint["path"])
        receipt_status = receipt.stat()
        self.assertTrue(stat.S_ISREG(receipt_status.st_mode))
        self.assertEqual(
            (receipt_status.st_dev, receipt_status.st_ino),
            (
                _canonical_decimal_integer(receipt_checkpoint["device"]),
                _canonical_decimal_integer(receipt_checkpoint["inode"]),
            ),
        )
        self.assertEqual(receipt_status.st_uid, receipt_checkpoint["uid"])
        self.assertEqual(receipt_status.st_gid, receipt_checkpoint["gid"])
        self.assertEqual(
            f"{stat.S_IMODE(receipt_status.st_mode):05o}",
            receipt_checkpoint["mode"],
        )
        self.assertEqual(receipt_status.st_nlink, receipt_checkpoint["hard_link_count"])
        self.assertEqual(receipt_status.st_size, receipt_checkpoint["size_bytes"])
        self.assertEqual(_sha256(receipt), receipt_checkpoint["sha256"])
        self.assertEqual(sorted(os.listxattr(receipt)), receipt_checkpoint["xattrs"])
        self.assertEqual(
            sorted(receipt.parent.glob(".q011-acceptance-recovery.*")),
            [],
        )
        for label, path, identity, expected_mode, expected_xattrs in (
            (
                "pic",
                acceptance_root.parent,
                helper["reviewed_pic_root_identity"],
                helper["reviewed_pic_root_mode"],
                helper["reviewed_pic_root_xattrs"],
            ),
            (
                "policy",
                receipt.parent,
                helper["reviewed_policy_root_identity"],
                helper["reviewed_policy_root_mode"],
                helper["reviewed_policy_root_xattrs"],
            ),
        ):
            with self.subTest(authority=label):
                status = path.stat()
                self.assertEqual(
                    (status.st_dev, status.st_ino),
                    (identity["device"], identity["inode"]),
                )
                self.assertEqual(status.st_uid, checkpoint["uid"])
                self.assertEqual(status.st_gid, checkpoint["gid"])
                self.assertEqual(f"{stat.S_IMODE(status.st_mode):05o}", expected_mode)
                self.assertEqual(sorted(os.listxattr(path)), expected_xattrs)
        aggregate_repair = current_repair["lustre_publication_rename_compatibility_repair"]
        authorization = aggregate_repair["postrun_source_authorization_successor"]
        self.assertEqual(authorization["sha256"], _sha256(REPO_ROOT / authorization["path"]))
        historical_authorization = json.loads(
            (REPO_ROOT / authorization["path"]).read_text(encoding="utf-8")
        )
        historical_source_by_path = {
            record["path"]: record["sha256"]
            for record in historical_authorization["source_closure"]
        }
        compatibility = aggregate_repair["snapshot_time_compatibility_successor"]
        self.assertEqual(compatibility["sha256"], _sha256(REPO_ROOT / compatibility["path"]))
        for key in (
            "aggregate_analyzer",
            "aggregate_publisher",
            "review_packet_renderer",
            "repair_validation_worker",
        ):
            source_change = aggregate_repair[key]
            with self.subTest(source_change=key):
                self.assertEqual(
                    source_change["predecessor_sha256"],
                    _git_blob_sha256(
                        current_repair["source_checkpoint_commit"],
                        source_change["path"],
                    ),
                )
                self.assertEqual(
                    source_change["successor_sha256"],
                    historical_source_by_path.get(
                        source_change["path"],
                        _sha256(REPO_ROOT / source_change["path"]),
                    ),
                )
        current_authorization_path = (
            READINESS_DIR
            / "q011_section54_pressure_pilot_postrun_aggregate_source_authorization_"
            "successor_v6_2026-06-05.json"
        )
        current_authorization = json.loads(
            current_authorization_path.read_text(encoding="utf-8")
        )
        self.assertEqual(
            current_authorization["predecessor_record"],
            authorization["path"],
        )
        self.assertEqual(
            current_authorization["predecessor_sha256"],
            authorization["sha256"],
        )
        for record in current_authorization["source_closure"]:
            with self.subTest(current_postrun_source=record["role"]):
                self.assertEqual(
                    record["sha256"],
                    _sha256(REPO_ROOT / record["path"]),
                )
        self.assertTrue(
            aggregate_repair["repair_validation_worker"]["compute_node_lock_preflight"]
        )
        self.assertIn(
            "never automatically delete canonical artifacts or any still-present staging aliases",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "for lock-honoring reviewed workers, remove the retained guard inode",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "Preflight aggregate and renderer source authorization",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "exclusive advisory transaction lock",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "durably sync the publication parent",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "expose no public production verifier flags",
            aggregate_repair["accepted_contract"],
        )
        self.assertIn(
            "destructive_rollback_after_any_canonical_artifact_exposure",
            aggregate_repair["rejected_contracts"],
        )
        self.assertIn(
            "success_seal_is_committed_and_inode_bound_while_the_guard_remains_armed",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "guard_removal_is_the_final_publication_state_transition",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "reviewed_publication_workers_are_serialized_by_the_stable_account_anchor_and_exact_acceptance_root_descriptor_locks",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "clean_frontier_repair_validation_preflights_both_production_transaction_locks_on_a_compute_node",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "final_guard_removal_targets_the_retained_guard_inode_for_lock_honoring_reviewed_workers",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "staging_link_removal_targets_the_retained_source_inode_for_lock_honoring_reviewed_workers",
            aggregate_repair["compatibility_guarantees"],
        )
        self.assertIn(
            "public_production_receipt_verification_always_requires_guard_absence_and_the_inode_bound_success_seal",
            aggregate_repair["compatibility_guarantees"],
        )
        reviews = aggregate_repair["independent_read_only_forensic_reviews"]
        self.assertEqual(reviews["requested"], 2)
        self.assertEqual(
            reviews["requirement"],
            "publish_a_separate_exact_current_readiness_review_artifact_after_this_"
            "transition_is_frozen",
        )
        self.assertEqual(
            reviews["status"],
            "prior_findings_integrated_pending_separate_exact_current_review_artifact",
        )
        self.assertNotIn("completed", reviews)
        self.assertFalse(current_repair["human_input_required_now"])
        self.assertEqual(
            current_repair["first_required_human_input_after_automated_repairs"],
            "Select exactly one Section 5.4 problem/ps_p0 case from the verified "
            "immutable four-slice pressure-review packet.",
        )
        self.assertEqual(
            current_repair["release_blockers"],
            [
                "sixteenth_lustre_publication_compatibility_repair_clean_commit_and_push_pending",
                "replacement_full_clean_worker_validation_pending",
                "exact_latest_patch_independent_rereview_pending",
                "pressure_aggregate_and_review_packet_worker_publication_pending",
                "pressure_selection_review_packet_binding_pending_before_human_selection_acceptance",
                "qualifying_campaign_publishers_lustre_compatibility_repair_pending_before_any_qualifying_launch",
            ],
        )
        post_publication_status = _load(
            "q011_section54_post_publication_pressure_gate_status_successor_"
            "2026-06-05.json"
        )
        staged_packet_gate = post_publication_status[
            "staged_packet_gate_control_plane"
        ]
        exact_predecessor_repair_successor = _load(
            "q011_section54_pressure_gate_exact_predecessor_migration_repair_"
            "successor_2026-06-05.json"
        )
        _validate_q011_exact_predecessor_migration_repair_successor(
            exact_predecessor_repair_successor
        )
        exact_predecessor_repair = exact_predecessor_repair_successor[
            "source_local_exact_predecessor_repair"
        ]
        self.assertEqual(
            exact_predecessor_repair["final_binding_refresh"],
            {
                "status": "completed_before_commit",
                "fields": ["control_plane_version", "source_test_closure"],
                "authority": "none",
            },
        )
        exact_predecessor_repair_version = exact_predecessor_repair[
            "control_plane_version"
        ]
        if staged_version == exact_predecessor_repair_version:
            self.assertEqual(
                exact_predecessor_repair_successor["packet_gate"]["authority"],
                "none",
            )
            self.assertEqual(
                exact_predecessor_repair_successor["frontier_launch_authorization"],
                "none_launch_prohibited",
            )
            self.assertEqual(
                staged_packet_gate["version"],
                "ccc9d8aef994bb64465f9b16de58236ff42ed2df028e8d9897ababb30f2cb7f1",
            )
            self.assertNotEqual(staged_version, staged_packet_gate["version"])
            self.assertEqual(
                staged_packet_gate["pair_install_authorization_by_this_status"],
                "none",
            )
            self.assertEqual(
                staged_packet_gate["policy_promotion_authorization_by_this_status"],
                "none",
            )
            return
        if staged_version == staged_packet_gate["version"]:
            live_baseline = post_publication_status["live_operational_baseline"]
            self.assertNotEqual(
                staged_version,
                live_baseline["installed_control_plane_version"],
            )
            self.assertEqual(
                staged_packet_gate["pair_install_authorization_by_this_status"],
                "none",
            )
            self.assertEqual(
                staged_packet_gate["policy_promotion_authorization_by_this_status"],
                "none",
            )
            return
        if staged_version == repaired_staged["version"]:
            self.assertEqual(
                repaired_staged["state"],
                "installed_candidate_only_policy_promoted_acceptance_root_recovered_"
                "lustre_publication_compatibility_repair_staged_pending_clean_commit_push_"
                "worker_validation_exact_latest_patch_rereview_and_aggregate_retry",
            )
            prepared = repaired_staged["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
            self.assertEqual(prepared["inventory_sha256"], _sha256(prepared_path))
            self.assertEqual(
                prepared["paper_deck_count"], len(prepared_inventory["paper_decks"])
            )
            self.assertEqual(
                prepared["publication_analyzer_count"],
                len(prepared_inventory["analyzers"]),
            )
            self.assertNotEqual(
                storage["installed_control_plane_version"], staged_version
            )
            self.assertEqual(
                storage["installed_control_plane_version"],
                storage["staged_control_plane_candidate_version"],
            )
            return
        strict_q011 = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v6_2026-06-02.json"
        )
        if (
            staged_version == strict_q011["control_plane_version"]
            and storage["installed_control_plane_version"]
            == strict_q011["control_plane_version"]
        ):
            self.assertEqual(
                strict_q011["predecessor_sha256"],
                _sha256(REPO_ROOT / strict_q011["predecessor_record"]),
            )
            installed = strict_q011["paired_install"]
            for key in ("orion", "project_home"):
                self.assertEqual(
                    installed[f"{key}_inventory_sha256"],
                    _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
                )
            self.assertEqual(
                installed["orion_inventory_sha256"],
                installed["project_home_inventory_sha256"],
            )
            promotion = _load(
                "phase0_curated_candidate_successor_v20_2026-06-02.json"
            )["policy_promotion"]
            self.assertEqual(
                promotion["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    promotion[f"{key}_sha256"], _sha256(Path(promotion[f"{key}_path"]))
                )
            self.assertEqual(promotion["registered_science_slice_count"], 4)
            self.assertEqual(
                strict_q011["scheduler_isolation"],
                {
                    "status": "not_claimed_empty_allowlist_transition_only",
                    "reason": (
                        "unrelated_user_scheduler_job_4754394_active_"
                        "during_policy_transition"
                    ),
                },
            )
            terminal = strict_q011["terminal_mirrored_ledger"]
            self.assertEqual(
                terminal["orion_node_hours_jsonl_sha256"],
                "dfb6e24a431c4a44682f864be2dbd5520816488d58cada7f0672aeea39c204b0",
            )
            self.assertEqual(
                terminal["project_home_node_hours_jsonl_sha256"],
                terminal["orion_node_hours_jsonl_sha256"],
            )
            self.assertEqual(
                terminal["orion_node_hours_csv_sha256"],
                "92a833de7096c81954e9fed572e51c0706b61b54adc3d1704d7130652d47f946",
            )
            self.assertEqual(
                terminal["orion_mirror_receipts_jsonl_sha256"],
                "696e93105d386ea0cd77c14398552a0f17581a65e7714625e51d0b8d84b131d5",
            )
            self.assertEqual(terminal["ledger_records"], 111)
            self.assertEqual(terminal["latest_sequence_number"], 110)
            self.assertEqual(
                terminal["latest_event_sha256"],
                "fc0082bef800733c433395d48f551559abe083367e5549cc4bffbe8a0ab48bfa",
            )
            self.assertEqual(terminal["active_reservations"], 0)
            successor = _load("phase0_curated_candidate_successor_v18_2026-06-02.json")
            self.assertEqual(
                successor["predecessor_sha256"],
                _sha256(REPO_ROOT / successor["predecessor_record"]),
            )
            receipt = successor[
                "paired_control_plane_install_and_policy_promotion_receipt"
            ]
            self.assertEqual(receipt["sha256"], _sha256(REPO_ROOT / receipt["path"]))
            self.assertEqual(
                successor["live_paired_control_plane_version"], staged_version
            )
            self.assertEqual(
                successor["frontier_launch_authorization"],
                "none_no_registered_science_slices",
            )
            self.assertEqual(
                successor["operational_baseline"]["scheduler_isolation"],
                "not_claimed_unrelated_user_scheduler_job_4754394_active_"
                "empty_allowlist_transition_only",
            )
            return
        registered_pilot = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v5_2026-06-02.json"
        )
        if staged_version == registered_pilot["control_plane_version"]:
            self.assertEqual(
                registered_pilot["predecessor_sha256"],
                _sha256(REPO_ROOT / registered_pilot["predecessor_record"]),
            )
            installed = registered_pilot["paired_install"]
            for key in ("orion", "project_home"):
                self.assertEqual(
                    installed[f"{key}_inventory_sha256"],
                    _sha256(Path(installed[f"{key}_root"]) / "inventory.json"),
                )
            self.assertEqual(
                installed["orion_inventory_sha256"],
                installed["project_home_inventory_sha256"],
            )
            return
        q011_retry = _load("phase0_curated_candidate_successor_v17_2026-06-02.json")
        if staged_version == q011_retry["successor_source_control_plane_version"]:
            self.assertEqual(
                q011_retry["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_curated_candidate_successor_v16_2026-06-02.json",
            )
            self.assertEqual(
                q011_retry["predecessor_sha256"],
                _sha256(REPO_ROOT / q011_retry["predecessor_record"]),
            )
            self.assertEqual(
                q011_retry["live_paired_control_plane_version"],
                registered_pilot["control_plane_version"],
            )
            registration = q011_retry["q011_retry_registration"]
            self.assertEqual(
                registration["sha256"], _sha256(REPO_ROOT / registration["path"])
            )
            prepared = q011_retry["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
            self.assertEqual(prepared["inventory_sha256"], _sha256(prepared_path))
            self.assertEqual(prepared["paper_deck_count"], len(prepared_inventory["paper_decks"]))
            self.assertEqual(
                prepared["publication_analyzer_count"], len(prepared_inventory["analyzers"])
            )
            staged = q011_retry["staged_control_plane"]
            self.assertEqual(staged["version"], staged_version)
            self.assertEqual(staged["inventoried_file_count"], len(CONTROL_PLANE_FILES))
            self.assertEqual(q011_retry["qualification_effect"], "none")
            self.assertEqual(
                q011_retry["frontier_launch_authorization"],
                "none_live_8f0a9d7f_empty_registered_science_allowlist",
            )
            return
        replay = _load(
            "phase0_registered_prerequisite_replay_policy_promotion_2026-06-01.json"
        )
        if (
            storage["installed_control_plane_version"] == replay["control_plane_version"]
            and staged_version == replay["control_plane_version"]
        ):
            self._assert_current_registered_replay_bindings(policy, staged_version)
            return
        if storage["installed_control_plane_version"] == replay["control_plane_version"]:
            successor = _load("phase0_curated_candidate_successor_v11_2026-06-01.json")
            self.assertEqual(
                successor["predecessor_record"],
                "tst/publication/readiness/"
                "phase0_curated_candidate_successor_v10_2026-06-01.json",
            )
            self.assertEqual(
                successor["predecessor_sha256"],
                _sha256(READINESS_DIR / successor["predecessor_record"].split("/")[-1]),
            )
            self.assertEqual(
                successor["status"],
                "local_vl2_tsc_successor_staged_validation_rereview_"
                "install_build_and_freeze_pending",
            )
            self.assertEqual(successor["qualification_effect"], "none")
            self.assertEqual(
                successor["candidate_freeze_source_commit"],
                "pending_final_clean_receipt_commit",
            )
            self.assertEqual(
                successor["frontier_launch_authorization"],
                "none_pending_validation_rereview_paired_install_"
                "clean_build_and_freeze",
            )
            self.assertEqual(
                successor["live_paired_control_plane_version"],
                storage["installed_control_plane_version"],
            )
            self.assertEqual(
                successor["successor_source_control_plane_version"], staged_version
            )
            predecessor_commit = successor["curated_source_predecessor_commit"]
            self.assertRegex(predecessor_commit, r"^[0-9a-f]{40}$")
            subprocess.check_call(
                ["git", "cat-file", "-e", f"{predecessor_commit}^{{commit}}"],
                cwd=REPO_ROOT,
            )
            prepared = successor["prepared_artifacts"]
            prepared_path = REPO_ROOT / prepared["inventory_path"]
            self.assertEqual(_sha256(prepared_path), prepared["inventory_sha256"])
            prepared_inventory = json.loads(
                prepared_path.read_text(encoding="utf-8")
            )
            expected_decks = sorted(
                [
                    *(REPO_ROOT / "inputs" / "tests").glob("pic*.athinput"),
                    *(
                        REPO_ROOT / path
                        for path in PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
                    ),
                ]
            )
            expected_analyzers = sorted(
                (REPO_ROOT / "tst" / "publication").glob("analyze_*.py")
            )
            self.assertEqual(len(expected_decks), prepared["paper_deck_count"])
            self.assertEqual(
                len(expected_analyzers), prepared["publication_analyzer_count"]
            )
            for key, expected_paths in {
                "paper_decks": expected_decks,
                "analyzers": expected_analyzers,
            }.items():
                records = prepared_inventory[key]
                self.assertEqual(
                    [record["path"] for record in records],
                    [
                        path.relative_to(REPO_ROOT).as_posix()
                        for path in expected_paths
                    ],
                )
                for record in records:
                    self.assertEqual(
                        _sha256(REPO_ROOT / record["path"]), record["sha256"]
                    )
            baseline = successor["operational_baseline"]
            self.assertEqual(
                baseline["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"], _sha256(Path(baseline[f"{key}_path"]))
            )
            return
        c83e_successor = _load(
            "phase0_curated_candidate_successor_v16_2026-06-02.json"
        )
        c83e_receipt_binding = c83e_successor[
            "paired_control_plane_install_and_policy_promotion_receipt"
        ]
        c83e_receipt = json.loads(
            (REPO_ROOT / c83e_receipt_binding["path"]).read_text(encoding="utf-8")
        )
        if (
            storage["installed_control_plane_version"]
            == c83e_receipt["control_plane_version"]
            and policy["science_submission_freeze"]
            == c83e_receipt["active_policy_promotion"]["science_submission_freeze"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                c83e_receipt["control_plane_version"],
            )
            self.assertEqual(staged_version, c83e_receipt["control_plane_version"])
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                c83e_successor["predecessor_sha256"],
                _sha256(
                    READINESS_DIR
                    / c83e_successor["predecessor_record"].split("/")[-1]
                ),
            )
            self.assertEqual(
                c83e_receipt_binding["sha256"],
                _sha256(REPO_ROOT / c83e_receipt_binding["path"]),
            )
            baseline = c83e_successor["operational_baseline"]
            self.assertEqual(
                baseline["repo_policy_sha256"],
                _sha256(READINESS_DIR / "storage_policy.json"),
            )
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"],
                    _sha256(Path(baseline[f"{key}_path"])),
                )
            return
        d720_successor = _load(
            "phase0_curated_candidate_successor_v15_2026-06-02.json"
        )
        d720_receipt_binding = d720_successor[
            "clean_candidate_freeze_and_policy_promotion_receipt"
        ]
        d720_receipt = json.loads(
            (REPO_ROOT / d720_receipt_binding["path"]).read_text(encoding="utf-8")
        )
        if (
            storage["installed_control_plane_version"]
            == d720_receipt["control_plane_version"]
            and policy["science_submission_freeze"]
            == d720_receipt["active_policy_promotion"]["science_submission_freeze"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"],
                d720_receipt["control_plane_version"],
            )
            self.assertEqual(policy["registered_science_slices"], [])
            self.assertEqual(
                d720_successor["predecessor_sha256"],
                _sha256(
                    READINESS_DIR
                    / d720_successor["predecessor_record"].split("/")[-1]
                ),
            )
            self.assertEqual(
                d720_receipt_binding["sha256"],
                _sha256(REPO_ROOT / d720_receipt_binding["path"]),
            )
            baseline = d720_successor["operational_baseline"]
            for key in (
                "orion_policy",
                "project_home_policy",
                "orion_promotion",
                "project_home_promotion",
            ):
                self.assertEqual(
                    baseline[f"{key}_sha256"],
                    _sha256(Path(baseline[f"{key}_path"])),
                )
            return
        transition = _load("phase0_curated_candidate_successor_v14_2026-06-01.json")
        if (
            storage["installed_control_plane_version"]
            == transition["successor_source_control_plane_version"]
        ):
            self.assertEqual(
                storage["staged_control_plane_candidate_version"], staged_version
            )
            self.assertEqual(
                transition["successor_source_control_plane_version"], staged_version
            )
            self.assertEqual(
                transition["predecessor_sha256"],
                _sha256(READINESS_DIR / transition["predecessor_record"].split("/")[-1]),
            )
            self.assertEqual(
                policy["science_submission_freeze"],
                {"status": "pending_clean_candidate_freeze"},
            )
            self.assertEqual(policy["registered_science_slices"], [])
            clean_candidate = transition["clean_candidate_freeze_receipt"]
            self.assertEqual(
                clean_candidate["sha256"],
                _sha256(REPO_ROOT / clean_candidate["path"]),
            )
            prepared = transition["prepared_artifacts"]
            self.assertEqual(
                prepared["inventory_sha256"],
                _sha256(REPO_ROOT / prepared["inventory_path"]),
            )
            return
        phase0_successor = _load(
            "phase0_curated_candidate_successor_v6_2026-06-01.json"
        )
        self.assertEqual(
            phase0_successor["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_curated_candidate_successor_v4_2026-05-31.json",
        )
        predecessor_commit = phase0_successor["curated_source_predecessor_commit"]
        self.assertRegex(predecessor_commit, r"^[0-9a-f]{40}$")
        subprocess.check_call(
            ["git", "cat-file", "-e", f"{predecessor_commit}^{{commit}}"],
            cwd=REPO_ROOT,
        )
        self.assertEqual(
            staged_version, phase0_successor["successor_source_control_plane_version"]
        )
        prepared = phase0_successor["prepared_artifacts"]
        prepared_path = REPO_ROOT / prepared["inventory_path"]
        self.assertEqual(_sha256(prepared_path), prepared["inventory_sha256"])
        prepared_inventory = json.loads(prepared_path.read_text(encoding="utf-8"))
        expected_decks = sorted(
            [
                *(REPO_ROOT / "inputs" / "tests").glob("pic*.athinput"),
                *(
                    REPO_ROOT / path
                    for path in PREPARED_ARTIFACT_REQUIRED_PUBLICATION_DECK_PATHS
                ),
            ]
        )
        expected_analyzers = sorted(
            (REPO_ROOT / "tst" / "publication").glob("analyze_*.py")
        )
        self.assertEqual(len(expected_decks), prepared["paper_deck_count"])
        self.assertEqual(
            len(expected_analyzers), prepared["publication_analyzer_count"]
        )
        for key, expected_paths in {
            "paper_decks": expected_decks,
            "analyzers": expected_analyzers,
        }.items():
            records = prepared_inventory[key]
            self.assertEqual(
                [record["path"] for record in records],
                [path.relative_to(REPO_ROOT).as_posix() for path in expected_paths],
            )
            for record in records:
                self.assertEqual(
                    _sha256(REPO_ROOT / record["path"]),
                    record["sha256"],
                )
        paired_binding = phase0_successor[
            "paired_install_and_policy_promotion_receipt"
        ]
        paired_path = REPO_ROOT / paired_binding["path"]
        self.assertEqual(_sha256(paired_path), paired_binding["sha256"])
        paired = json.loads(paired_path.read_text(encoding="utf-8"))
        self.assertEqual(
            paired["predecessor_record"],
            "tst/publication/readiness/"
            "phase0_paired_control_plane_install_and_policy_promotion_successor_2026-05-31.json",
        )
        self.assertEqual(paired["control_plane_version"], staged_version)
        installed = paired["paired_install"]
        self.assertEqual(
            installed["orion_inventory_sha256"],
            installed["project_home_inventory_sha256"],
        )
        self.assertEqual(installed["inventoried_file_count"], len(CONTROL_PLANE_FILES))
        self.assertEqual(installed["generation_directory_mode"], "0555")
        self.assertTrue(installed["byte_identical_inventory"])
        self.assertEqual(installed["installed_pair_verify"], "pass")
        promotion = paired["active_policy_promotion"]
        self.assertEqual(promotion["maximum_node_hours"], 10000)
        self.assertEqual(promotion["registered_science_slices"], [])
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            _sha256(READINESS_DIR / "storage_policy.json"),
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        self.assertTrue(promotion["byte_identical_policy"])
        self.assertTrue(promotion["byte_identical_promotion"])
        self.assertEqual(
            storage["installed_control_plane_version"],
            storage["staged_control_plane_candidate_version"],
        )
        self.assertEqual(
            staged_version,
            storage["installed_control_plane_version"],
        )
        if policy["science_submission_freeze"] == {
            "status": "pending_clean_candidate_freeze"
        }:
            self.assertEqual(policy["registered_science_slices"], [])
            return
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        clean_manifest_path = Path(policy["science_submission_freeze"]["manifest_path"])
        clean_manifest = json.loads(clean_manifest_path.read_text(encoding="utf-8"))
        executable_path = Path(clean_manifest["build"]["executable_path"])
        binding_paths = {
            "f1-clean-gyro-mpich-stderr-v3": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_relativistic_gyro_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_relativistic_gyro_paper.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
            "f1-clean-paper-coupling-mpich-stderr-v2": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_gpu_paper_coupling_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_paper_coupling_conservation.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
                    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
            "f2-parser-multirank-runtime-metadata-v2": {
                "job_script_sha256": (
                    REPO_ROOT
                    / "tst/publication/frontier_f2_structured_multirank_runtime_metadata_job.sh"
                ),
                "input_deck_sha256": (
                    REPO_ROOT / "inputs/tests/pic_parser_contract_guards.athinput"
                ),
                "analysis_script_sha256": [
                    REPO_ROOT
                    / "tst/publication/frontier_f2_multirank_runtime_metadata_analysis.py",
                    REPO_ROOT
                    / "tst/publication/frontier_f1_structured_artifacts.py",
                ],
            },
        }
        environment_path = CONTROL_PLANE_DIR / "frontier_pic_environment.sh"
        f2_candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        successor_records = [
            *successor["registered_science_slices"],
            f2_candidate["registered_science_slice"],
        ]
        successors = {
            record["authorization_id"]: record
            for record in successor_records
        }
        self.assertEqual(len(successor_records), len(successors))
        self.assertEqual(set(successors), set(binding_paths))
        self.assertEqual(
            {
                record["authorization_id"]
                for record in policy["registered_science_slices"]
            },
            set(binding_paths),
        )
        for authorization in policy["registered_science_slices"]:
            authorization_id = authorization["authorization_id"]
            with self.subTest(authorization_id=authorization_id):
                paths = binding_paths[authorization_id]
                for key in (
                    "campaign",
                    "maximum_nodes",
                    "maximum_walltime_seconds",
                    "maximum_attempts",
                ):
                    self.assertEqual(
                        successors[authorization_id][key], authorization[key]
                    )
                self.assertEqual(
                    authorization["job_script_sha256"],
                    _sha256(paths["job_script_sha256"]),
                )
                self.assertEqual(
                    authorization["input_deck_sha256"],
                    _sha256(paths["input_deck_sha256"]),
                )
                self.assertEqual(
                    authorization["environment_profile_sha256"],
                    _sha256(environment_path),
                )
                analysis_sha256 = [
                    _sha256(path) for path in paths["analysis_script_sha256"]
                ]
                self.assertEqual(
                    authorization["analysis_script_sha256"], analysis_sha256
                )
                self.assertEqual(
                    successors[authorization_id]["analysis_script_sha256"],
                    analysis_sha256[0],
                )
                self.assertEqual(
                    successors[authorization_id]["analysis_support_sha256"],
                    analysis_sha256[1],
                )
                self.assertEqual(
                    authorization["clean_candidate_manifest_sha256"],
                    _sha256(clean_manifest_path),
                )
                self.assertEqual(
                    authorization["executable_sha256"], _sha256(executable_path)
                )

    def test_exact_current_identity_bindings_require_canonical_decimal_strings(
        self,
    ) -> None:
        self.assertEqual(_canonical_decimal_integer("0"), 0)
        self.assertEqual(
            _canonical_decimal_integer("720587400627193972"),
            720587400627193972,
        )
        for invalid in (
            1,
            True,
            None,
            1.0,
            "",
            "01",
            "+1",
            "-1",
            " 1",
            "1 ",
            "\u0661",
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    _canonical_decimal_integer(invalid)

    def test_exact_current_review_chronology_requires_canonical_utc_seconds(
        self,
    ) -> None:
        self.assertEqual(
            _canonical_utc_second("2026-06-05T05:17:39Z"),
            datetime(2026, 6, 5, 5, 17, 39, tzinfo=timezone.utc),
        )
        for invalid in (
            None,
            0,
            "2026-06-05T05:17:39",
            "2026-06-05T05:17:39+00:00",
            "2026-06-05T05:17:39.0Z",
            "2026-6-05T05:17:39Z",
            "2026-06-05t05:17:39Z",
            "2026-02-30T05:17:39Z",
        ):
            with self.subTest(invalid=invalid):
                with self.assertRaises(ValueError):
                    _canonical_utc_second(invalid)

    def test_exact_current_review_summary_requires_canonical_identity_and_facts(
        self,
    ) -> None:
        self.assertEqual(
            _canonical_reviewer_id("019e95d8-49db-7400-9aa8-b53aba417431"),
            "019e95d8-49db-7400-9aa8-b53aba417431",
        )
        self.assertEqual(
            _canonical_nonempty_fact_list(["frozen basis matched", "roots empty"]),
            ["frozen basis matched", "roots empty"],
        )
        for invalid in (None, False, "", " reviewer", "reviewer ", "review er"):
            with self.subTest(invalid_reviewer_id=invalid):
                with self.assertRaises(ValueError):
                    _canonical_reviewer_id(invalid)
        for invalid in (
            "019E95D8-49DB-7400-9AA8-B53ABA417431",
            "019e95d8-49db-7400-9aa8-b53aba417431\u200b",
            "019e95d8-49db-7400-9aa8-b53aba417431\0",
        ):
            with self.subTest(invalid_reviewer_id=invalid):
                with self.assertRaises(ValueError):
                    _canonical_reviewer_id(invalid)
        for invalid in (
            None,
            {},
            "fact",
            [],
            [""],
            [" fact"],
            ["fact "],
            ["two  spaces"],
            ["two\twords"],
            ["two\nlines"],
            ["fact\0"],
            ["fact\u200b"],
            ["valid", None],
        ):
            with self.subTest(invalid_facts=invalid):
                with self.assertRaises(ValueError):
                    _canonical_nonempty_fact_list(invalid)

    def test_exact_current_lustre_publication_rereview_artifact_recomputes(
        self,
    ) -> None:
        rereview = _load(
            "q011_section54_sixteenth_lustre_publication_exact_current_"
            "rereview_2026-06-05.json"
        )
        self.assertEqual(
            set(rereview),
            {
                "schema_version",
                "record_type",
                "recorded_utc",
                "source_checkpoint_commit",
                "transition",
                "reviewed_runtime",
                "required_scopes",
                "reviews",
                "status",
                "qualification_effect",
            },
        )
        self.assertEqual(rereview["schema_version"], 1)
        self.assertEqual(
            rereview["record_type"],
            "q011_section54_sixteenth_lustre_publication_exact_current_rereview",
        )
        self.assertEqual(
            rereview["status"],
            "exact_current_read_only_rereviews_complete_no_findings",
        )
        self.assertEqual(
            rereview["qualification_effect"],
            "none_no_execution_authorization_no_science_claim",
        )

        transition_binding = rereview["transition"]
        transition_path = (
            READINESS_DIR
            / "q011_section54_sixteenth_lustre_publication_rename_compatibility_"
            "transition_2026-06-05.json"
        )
        self.assertEqual(
            transition_binding,
            {
                "path": transition_path.relative_to(REPO_ROOT).as_posix(),
                "sha256": _sha256(transition_path),
            },
        )
        transition = json.loads(transition_path.read_text(encoding="utf-8"))
        self.assertEqual(
            rereview["source_checkpoint_commit"],
            transition["source_checkpoint_commit"],
        )
        self.assertGreater(
            _canonical_utc_second(rereview["recorded_utc"]),
            _canonical_utc_second(transition["exact_patch_frozen_utc"]),
        )

        reviewed_runtime = rereview["reviewed_runtime"]
        repair = transition["lustre_publication_rename_compatibility_repair"]
        expected_reviewed_runtime = {
            label: {
                "path": repair[label]["path"],
                "sha256": repair[label][hash_key],
            }
            for label, hash_key in (
                ("aggregate_publisher", "successor_sha256"),
                ("review_packet_renderer", "successor_sha256"),
                ("postrun_source_authorization_successor", "sha256"),
                ("repair_validation_worker", "successor_sha256"),
            )
        }
        self.assertEqual(
            reviewed_runtime,
            expected_reviewed_runtime,
        )
        historical_authorization = json.loads(
            (
                REPO_ROOT
                / reviewed_runtime["postrun_source_authorization_successor"]["path"]
            ).read_text(encoding="utf-8")
        )
        historical_source_by_path = {
            record["path"]: record["sha256"]
            for record in historical_authorization["source_closure"]
        }
        for label, binding in reviewed_runtime.items():
            with self.subTest(reviewed_runtime=label):
                self.assertEqual(
                    binding["sha256"],
                    historical_source_by_path.get(
                        binding["path"],
                        _sha256(REPO_ROOT / binding["path"]),
                    ),
                )

        required_scopes = {"filesystem_publication", "provenance_chronology"}
        self.assertEqual(set(rereview["required_scopes"]), required_scopes)
        self.assertEqual(len(rereview["required_scopes"]), len(required_scopes))
        reviews = rereview["reviews"]
        self.assertEqual(len(reviews), 2)
        self.assertEqual(
            {review["scope"] for review in reviews},
            required_scopes,
        )
        reviewer_ids = [
            _canonical_reviewer_id(review["reviewer_id"]) for review in reviews
        ]
        self.assertEqual(len(set(reviewer_ids)), 2)
        exact_hash_basis = {
            "source_checkpoint_commit": rereview["source_checkpoint_commit"],
            "transition_sha256": transition_binding["sha256"],
            **{
                f"{label}_sha256": binding["sha256"]
                for label, binding in reviewed_runtime.items()
            },
        }
        frozen_utc = _canonical_utc_second(transition["exact_patch_frozen_utc"])
        recorded_utc = _canonical_utc_second(rereview["recorded_utc"])
        for review in reviews:
            with self.subTest(reviewer=review["reviewer_id"]):
                self.assertEqual(
                    set(review),
                    {
                        "reviewer_id",
                        "reviewed_utc",
                        "scope",
                        "exact_hash_basis",
                        "facts_inspected",
                        "findings",
                        "disposition",
                        "read_only",
                        "no_files_edited",
                        "no_jobs_launched",
                    },
                )
                reviewed_utc = _canonical_utc_second(review["reviewed_utc"])
                self.assertGreater(reviewed_utc, frozen_utc)
                self.assertLessEqual(reviewed_utc, recorded_utc)
                self.assertEqual(review["exact_hash_basis"], exact_hash_basis)
                _canonical_nonempty_fact_list(review["facts_inspected"])
                self.assertEqual(review["findings"], [])
                self.assertEqual(review["disposition"], "no_findings")
                self.assertIs(review["read_only"], True)
                self.assertIs(review["no_files_edited"], True)
                self.assertIs(review["no_jobs_launched"], True)

    def test_q011_bounded_absence_search_matches_root_level_double_star(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            root_receipt = root / "human_pressure_selection_receipt.json"
            root_receipt.write_text("{}\n", encoding="utf-8")
            nested = root / "nested"
            nested.mkdir()
            nested_receipt = nested / "human_pressure_selection_receipt.json"
            nested_receipt.write_text("{}\n", encoding="utf-8")
            too_deep = nested / "deeper"
            too_deep.mkdir()
            (too_deep / "human_pressure_selection_receipt.json").write_text(
                "{}\n", encoding="utf-8"
            )
            self.assertEqual(
                _bounded_relative_matches(
                    root,
                    "**/*pressure*selection*receipt*.json",
                    2,
                ),
                [
                    "human_pressure_selection_receipt.json",
                    "nested/human_pressure_selection_receipt.json",
                ],
            )
            with self.assertRaisesRegex(ValueError, "nonnegative integer"):
                _bounded_relative_matches(root, "**/*.json", True)

    def test_q011_clean_snapshot_counts_exclude_symlinked_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            temporary_root = Path(temporary)
            root = temporary_root / "root"
            root.mkdir()
            regular = root / "test_regular.py"
            regular.write_text("pass\n", encoding="utf-8")
            (root / "test_symlink.py").symlink_to(regular)
            external = temporary_root / "external"
            external.mkdir()
            (external / "test_external.py").write_text("pass\n", encoding="utf-8")
            (root / "symlinked_directory").symlink_to(
                external,
                target_is_directory=True,
            )
            self.assertEqual(_regular_files_below(root, "test_*.py"), [regular])

    def test_q011_pressure_gate_validation_worker_isolates_python(self) -> None:
        wrapper = (
            REPO_ROOT
            / "tst/publication/frontier_q011_section54_pressure_gate_validation_job.sh"
        ).read_text(encoding="utf-8")
        self.assertIn(
            "TRUSTED_SYSTEM_SITE_PACKAGES="
            "/opt/cray/pe/python/3.11.7/lib/python3.11/site-packages",
            wrapper,
        )
        self.assertIn("PYTHON_ENV=(\n  /usr/bin/env -i", wrapper)
        self.assertIn('set_authenticated_python_roots "$SNAPSHOT_ROOT"', wrapper)
        self.assertIn('set_authenticated_python_roots "$FULL_SUITE_ROOT"', wrapper)
        self.assertIn(
            'clone --no-checkout --shared "$REPO_ROOT" "$FULL_SUITE_ROOT"',
            wrapper,
        )
        self.assertNotIn('set_authenticated_python_roots "$REPO_ROOT"', wrapper)
        self.assertNotIn("/usr/bin/git -c core.fsmonitor=false", wrapper)
        self.assertGreaterEqual(wrapper.count("--no-replace-objects"), 10)
        self.assertIn(
            "sys.path[:0] = [source_root, control_plane_root, site_packages]",
            wrapper,
        )
        self.assertEqual(
            wrapper.count(
                '"$AUTHENTICATED_SOURCE_ROOT" "$AUTHENTICATED_CONTROL_PLANE_ROOT" \\\n'
                '    "$TRUSTED_SYSTEM_SITE_PACKAGES"'
            ),
            2,
        )
        self.assertNotIn("export PYTHONPATH=", wrapper)
        python_launches = [
            line.strip() for line in wrapper.splitlines() if '"$PYTHON"' in line
        ]
        self.assertEqual(
            python_launches,
            ['"${PYTHON_ENV[@]}" "$PYTHON" -I -B -S -c \\'] * 2,
        )

    def test_q011_historical_status_binds_installed_controller_inventory(self) -> None:
        status = _load(
            "q011_section54_post_publication_pressure_gate_status_successor_"
            "2026-06-05.json"
        )
        live = status["live_operational_baseline"]
        inventory_payloads = []
        for binding in live["inventories"]:
            path = Path(binding["path"])
            payload = path.read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), binding["sha256"])
            inventory_payloads.append(payload)
            installed = verify_historical_installed_control_plane(
                path.parent,
                authorized_pic_root=path.parents[2],
            )
            self.assertEqual(installed["version"], live["installed_control_plane_version"])
        self.assertEqual(len(set(inventory_payloads)), 1)

    def test_q011_post_publication_pressure_gate_status_successor_recomputes(
        self,
    ) -> None:
        status_path = (
            READINESS_DIR
            / "q011_section54_post_publication_pressure_gate_status_successor_"
            "2026-06-05.json"
        )
        status_payload = status_path.read_text(encoding="utf-8")
        status = json.loads(status_payload)
        self.assertEqual(
            _sha256(status_path),
            "0ba6ac19df478630830b13a68443e0560cfcf9dbd00d48f20ad25d59a39568be",
        )
        _validate_q011_post_publication_pressure_gate_status_successor(status)
        self.assertEqual(status_payload, json.dumps(status, indent=2) + "\n")

        predecessor_path = REPO_ROOT / status["predecessor_record"]
        self.assertEqual(_sha256(predecessor_path), status["predecessor_sha256"])
        predecessor = json.loads(predecessor_path.read_text(encoding="utf-8"))
        self.assertGreater(
            _canonical_utc_second(status["recorded_utc"]),
            _canonical_utc_second(predecessor["exact_patch_frozen_utc"]),
        )

        authorization = status["postrun_source_authorization_successor"]
        self.assertEqual(
            _sha256(REPO_ROOT / authorization["path"]),
            authorization["sha256"],
        )
        # This successor is an immutable point-in-time status. Its source
        # closure, validation worker, staged controller, active anchors, and
        # scoped absence observations are expected to be superseded later.
        # Recompute only bindings that the record declares immutable.
        evidence = status["published_pressure_evidence"]
        for label, binding in evidence.items():
            with self.subTest(published_pressure_evidence=label):
                path = Path(binding["path"])
                self.assertTrue(stat.S_ISREG(path.lstat().st_mode))
                self.assertEqual(_sha256(path), binding["sha256"])
                self.assertEqual(stat.S_IMODE(path.stat().st_mode) & 0o222, 0)

        aggregate_receipt = json.loads(
            Path(evidence["aggregate_receipt"]["path"]).read_text(encoding="utf-8")
        )
        packet_receipt = json.loads(
            Path(evidence["review_packet_receipt"]["path"]).read_text(
                encoding="utf-8"
            )
        )
        self.assertEqual(
            packet_receipt["aggregate_receipt"],
            evidence["aggregate_receipt"],
        )
        self.assertEqual(
            packet_receipt["inventory_sha256"],
            evidence["review_packet_inventory"]["sha256"],
        )
        self.assertEqual(
            aggregate_receipt["aggregate_analysis"],
            evidence["aggregate_analysis"],
        )
        self.assertEqual(
            aggregate_receipt["aggregate_bundle"]["manifest_sha256"],
            evidence["aggregate_manifest"]["sha256"],
        )
        for label, receipt in (
            ("aggregate", aggregate_receipt),
            ("review_packet", packet_receipt),
        ):
            with self.subTest(published_receipt=label):
                self.assertEqual(
                    receipt["source_bindings"]["postrun_aggregate_source_authorization"],
                    authorization,
                )
                self.assertEqual(
                    receipt["source_bindings"]["runtime_source_archive"]["git_commit"],
                    status["source_checkpoint_commit"],
                )
        consumed = (
            pressure_packet_verifier.consume_published_pressure_pilot_review_packet(
                evidence["review_packet_receipt"]["path"],
                aggregate_receipt_binding=evidence["aggregate_receipt"],
                authorized_pic_root=Path(evidence["aggregate_receipt"]["path"]).parents[1],
            )
        )
        self.assertEqual(
            consumed["receipt_binding"],
            evidence["review_packet_receipt"],
        )
        self.assertEqual(
            consumed["aggregate_receipt_binding"],
            evidence["aggregate_receipt"],
        )
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", status["source_checkpoint_commit"]],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )

        memo = status["advisory_pressure_options_memo"]
        memo_path = REPO_ROOT / memo["path"]
        self.assertEqual(_sha256(memo_path), memo["sha256"])
        memo_text = memo_path.read_text(encoding="utf-8")
        self.assertIn("no pressure selected; no execution", memo_text)
        self.assertIn(
            "recommendation is not a human pressure-selection receipt",
            memo_text,
        )

    def test_q011_post_publication_pressure_gate_status_successor_rejects_drift(
        self,
    ) -> None:
        status = _load(
            "q011_section54_post_publication_pressure_gate_status_successor_"
            "2026-06-05.json"
        )
        rejection_cases = [
            ("non-successor chronology", ("recorded_utc",), "2026-06-05T05:06:20Z"),
            ("later unbound timestamp", ("recorded_utc",), "2026-06-05T07:55:57Z"),
            ("schema boolean", ("schema_version",), True),
            ("source checkpoint", ("source_checkpoint_commit",), "0" * 40),
            ("predecessor", ("predecessor_sha256",), "0" * 64),
            (
                "v5 authorization",
                ("postrun_source_authorization_successor", "sha256"),
                "0" * 64,
            ),
            (
                "clean-snapshot worker",
                ("clean_snapshot_pressure_gate_validation_worker", "sha256"),
                "0" * 64,
            ),
            (
                "clean-snapshot Python count",
                (
                    "clean_snapshot_pressure_gate_validation_worker",
                    "expected_publication_python_files",
                ),
                142,
            ),
            (
                "clean-snapshot shell count type",
                (
                    "clean_snapshot_pressure_gate_validation_worker",
                    "expected_publication_shell_files",
                ),
                True,
            ),
            (
                "clean-snapshot worker launch status",
                ("clean_snapshot_pressure_gate_validation_worker", "status"),
                "passed",
            ),
            (
                "unscoped clean-snapshot worker status",
                ("clean_snapshot_pressure_gate_validation_worker", "status"),
                "no_retained_clean_snapshot_validation_evidence",
            ),
            (
                "staged packet-gate control-plane version",
                ("staged_packet_gate_control_plane", "version"),
                "0" * 64,
            ),
            (
                "staged packet-gate file count type",
                ("staged_packet_gate_control_plane", "inventoried_file_count"),
                True,
            ),
            (
                "staged packet-gate prepared inventory",
                (
                    "staged_packet_gate_control_plane",
                    "prepared_artifact_inventory",
                    "sha256",
                ),
                "0" * 64,
            ),
            (
                "staged packet-gate pair install authorization",
                (
                    "staged_packet_gate_control_plane",
                    "pair_install_authorization_by_this_status",
                ),
                "authorized",
            ),
            (
                "live inventory",
                ("live_operational_baseline", "inventories", 0, "sha256"),
                "0" * 64,
            ),
            (
                "live inventoried count numeric alias",
                ("live_operational_baseline", "inventoried_member_count"),
                23.0,
            ),
            (
                "live byte-identical boolean alias",
                ("live_operational_baseline", "inventories_byte_identical"),
                1,
            ),
            (
                "live policy projection",
                (
                    "live_operational_baseline",
                    "projected_policy_state",
                    "registered_science_slices",
                ),
                ["q011"],
            ),
            (
                "live promotion",
                ("live_operational_baseline", "active_promotions", 0, "sha256"),
                "0" * 64,
            ),
            (
                "source/test closure path",
                ("source_test_closure", "files", 0, "path"),
                "tst/publication/other.py",
            ),
            (
                "unscoped source/test closure status",
                ("source_test_closure", "status"),
                "source_local_change_and_validation_bytes_bound_uncommitted",
            ),
            (
                "source/test closure malformed hash",
                ("source_test_closure", "files", 0, "sha256"),
                "0",
            ),
            (
                "source/test closure valid-looking member hash drift",
                ("source_test_closure", "files", 0, "sha256"),
                "0" * 64,
            ),
            (
                "source/test closure valid-looking digest drift",
                ("source_test_closure", "closure_sha256"),
                "0" * 64,
            ),
            (
                "aggregate receipt",
                ("published_pressure_evidence", "aggregate_receipt", "sha256"),
                "0" * 64,
            ),
            (
                "aggregate manifest",
                ("published_pressure_evidence", "aggregate_manifest", "sha256"),
                "0" * 64,
            ),
            (
                "packet receipt",
                ("published_pressure_evidence", "review_packet_receipt", "sha256"),
                "0" * 64,
            ),
            (
                "packet inventory",
                ("published_pressure_evidence", "review_packet_inventory", "sha256"),
                "0" * 64,
            ),
            (
                "publication root identity",
                (
                    "publication_acceptance_state",
                    "publication_root",
                    "identity",
                    "inode",
                ),
                0,
            ),
            (
                "publication root identity numeric alias",
                (
                    "publication_acceptance_state",
                    "publication_root",
                    "identity",
                    "device",
                ),
                135357496.0,
            ),
            (
                "success seal",
                ("publication_acceptance_state", "success_seals", 0, "sha256"),
                "0" * 64,
            ),
            (
                "publication guard absence",
                ("publication_acceptance_state", "absent_publication_paths"),
                [],
            ),
            (
                "selection-receipt absence",
                (
                    "scoped_absence_evidence",
                    "pressure_selection_receipt_searches",
                    0,
                    "matches",
                ),
                ["human_pressure_selection_receipt.json"],
            ),
            (
                "selection-receipt search depth numeric alias",
                (
                    "scoped_absence_evidence",
                    "pressure_selection_receipt_searches",
                    0,
                    "max_depth",
                ),
                8.0,
            ),
            (
                "validation-log absence",
                (
                    "scoped_absence_evidence",
                    "clean_snapshot_validation_log_search",
                    "matches",
                ),
                ["pic-q011-pressure-gate-validate.1.log"],
            ),
            (
                "advisory memo",
                ("advisory_pressure_options_memo", "sha256"),
                "0" * 64,
            ),
            (
                "implied selection status",
                ("pressure_selection", "status"),
                "human_selection_complete",
            ),
            (
                "selection receipt",
                ("pressure_selection", "selection_receipt"),
                {},
            ),
            (
                "selected case",
                ("pressure_selection", "selected_case"),
                {"problem_ps_p0": 1.0},
            ),
            (
                "advisory recommendation promoted",
                ("pressure_selection", "advisory_recommendation_is_selection"),
                True,
            ),
            ("packet gate unblocked", ("packet_gate", "status"), "ready"),
            (
                "packet gate acceptance authorization",
                ("packet_gate", "acceptance_authorization"),
                "human_review_only",
            ),
            (
                "packet gate schema downgrade",
                ("packet_gate", "required_receipt_schema_version"),
                1,
            ),
            (
                "packet gate schema type drift",
                ("packet_gate", "required_receipt_schema_version"),
                True,
            ),
            (
                "acceptance path removed",
                ("packet_gate", "required_acceptance_and_replay_paths"),
                status["packet_gate"]["required_acceptance_and_replay_paths"][:-1],
            ),
            (
                "blocker removed",
                ("packet_gate", "blockers_bound_by_this_status_successor"),
                status["packet_gate"]["blockers_bound_by_this_status_successor"][:-1],
            ),
            (
                "unscoped blocker claim",
                ("packet_gate", "blockers_bound_by_this_status_successor"),
                ["packet_bound_pressure_selection_gate_repair_not_committed"],
            ),
            ("human input enabled", ("human_input_required_now",), True),
            (
                "human input type drift",
                ("human_input_required_now",),
                0,
            ),
            (
                "launch authorization granted",
                ("frontier_launch_authorization",),
                "authorized",
            ),
        ]
        for label, path, replacement in rejection_cases:
            with self.subTest(label=label):
                candidate = copy.deepcopy(status)
                _replace_nested(candidate, path, replacement)
                with self.assertRaises(ValueError):
                    _validate_q011_post_publication_pressure_gate_status_successor(
                        candidate
                    )
        status_with_extra_key = copy.deepcopy(status)
        status_with_extra_key["unexpected"] = True
        with self.assertRaises(ValueError):
            _validate_q011_post_publication_pressure_gate_status_successor(
                status_with_extra_key
            )

    def test_q011_exact_predecessor_repair_source_test_closure_contract(
        self,
    ) -> None:
        newly_critical_paths = {
            "tst/publication/frontier_control_plane/README.md",
            "tst/publication/frontier_control_plane/run_control_plane.py",
            "tst/publication/frontier_control_plane/install_control_plane.py",
            "tst/publication/frontier_control_plane/capture_storage_preflight_evidence.py",
            "tst/publication/frontier_control_plane/storage_preflight.schema.json",
            "tst/publication/frontier_control_plane/write_orion_build_profile.py",
            "tst/publication/test_capture_storage_preflight_evidence.py",
        }
        closure_paths = _Q011_EXACT_PREDECESSOR_REPAIR_SOURCE_TEST_CLOSURE_PATHS
        self.assertEqual(len(closure_paths), 16)
        self.assertEqual(len(set(closure_paths)), len(closure_paths))
        self.assertTrue(newly_critical_paths.issubset(closure_paths))
        for relative in closure_paths:
            with self.subTest(closure_path=relative):
                self.assertTrue((REPO_ROOT / relative).is_file())

    def test_q011_exact_predecessor_migration_repair_successor_recomputes(
        self,
    ) -> None:
        successor_path = (
            READINESS_DIR
            / "q011_section54_pressure_gate_exact_predecessor_migration_repair_"
            "successor_2026-06-05.json"
        )
        successor_payload = successor_path.read_text(encoding="utf-8")
        successor = json.loads(successor_payload)
        _validate_q011_exact_predecessor_migration_repair_successor(successor)
        self.assertEqual(successor_payload, json.dumps(successor, indent=2) + "\n")

        predecessor_path = REPO_ROOT / successor["predecessor_record"]
        self.assertEqual(_sha256(predecessor_path), successor["predecessor_sha256"])
        predecessor = json.loads(predecessor_path.read_text(encoding="utf-8"))
        self.assertGreater(
            _canonical_utc_second(successor["recorded_utc"]),
            _canonical_utc_second(predecessor["recorded_utc"]),
        )

        retained_worker = successor["retained_clean_worker"]
        self.assertEqual(
            _sha256(Path(retained_worker["log"]["path"])),
            retained_worker["log"]["sha256"],
        )

        failed_attempt = successor["retained_failed_migration_attempt"]
        paired = failed_attempt["paired_control_plane"]
        installed_inventory_payloads = []
        for key in ("orion_inventory", "project_home_inventory"):
            binding = paired[key]
            payload = Path(binding["path"]).read_bytes()
            self.assertEqual(hashlib.sha256(payload).hexdigest(), binding["sha256"])
            installed_inventory_payloads.append(payload)
        self.assertEqual(len(set(installed_inventory_payloads)), 1)
        self.assertTrue(paired["inventories_byte_identical"])

        preflight = failed_attempt["fresh_storage_preflight"]
        self.assertEqual(
            _sha256(Path(preflight["binding"]["path"])),
            preflight["binding"]["sha256"],
        )
        for key in ("orion_evidence_path", "project_home_evidence_path"):
            self.assertEqual(
                _sha256(Path(preflight[key])),
                preflight["evidence_sha256"],
            )
        reviewed_policy = failed_attempt["reviewed_policy"]
        self.assertEqual(
            _sha256(Path(reviewed_policy["path"])),
            reviewed_policy["sha256"],
        )

        failed_build = successor["retained_failed_build_freeze_attempt"]
        self.assertEqual(
            _sha256(Path(failed_build["log"]["path"])),
            failed_build["log"]["sha256"],
        )
        build_residue = failed_build["residue"]["build_root"]
        build_root = Path(build_residue["path"])
        self.assertEqual(
            sorted(path.name for path in build_root.iterdir()),
            build_residue["top_level_entries"],
        )
        self.assertEqual(
            sum(1 for path in build_root.rglob("*") if path.is_file()),
            build_residue["file_count"],
        )
        self.assertEqual(
            1 + sum(1 for path in build_root.rglob("*") if path.is_dir()),
            build_residue["directory_count"],
        )
        self.assertEqual(
            sum(1 for path in build_root.rglob("*") if path.is_symlink()),
            build_residue["symlink_count"],
        )
        bin_residue = failed_build["residue"]["bin_root"]
        self.assertEqual(
            sum(1 for _path in Path(bin_residue["path"]).iterdir()),
            bin_residue["entry_count"],
        )
        self.assertNotIn(
            "clean_candidate_manifest=",
            Path(failed_build["log"]["path"]).read_text(encoding="utf-8"),
        )

        live = successor["unchanged_live_operational_baseline"]
        for bindings_key, identical_key in (
            ("active_policies", "active_policies_byte_identical"),
            ("active_promotions", "active_promotions_byte_identical"),
        ):
            self.assertEqual(
                len({binding["sha256"] for binding in live[bindings_key]}),
                1,
            )
            self.assertTrue(live[identical_key])
        stage4_successor = _load(
            "q011_section54_pressure_selection_publication_candidate_"
            "successor_2026-06-05.json"
        )
        self.assertEqual(
            stage4_successor["predecessor_record"],
            (
                "tst/publication/readiness/"
                "q011_section54_pressure_gate_exact_predecessor_migration_repair_"
                "successor_2026-06-05.json"
            ),
        )
        self.assertNotEqual(
            live["active_policies"][0]["sha256"],
            stage4_successor["active_candidate_only_state"]["active_policy"]["sha256"],
        )
        active_predecessor_controller = failed_build["worker_inputs"][
            "control_plane_version"
        ]
        preserved_freeze_controller = failed_build["worker_inputs"][
            "expected_authorized_freeze_build_controller"
        ]
        self.assertNotEqual(active_predecessor_controller, preserved_freeze_controller)

        historical_closure_commit = "67a418c432e2d424aa9e6cf5ed16316ea40fc0a4"
        repair = successor["source_local_exact_predecessor_repair"]
        checkpoint_control_plane_bindings = [
            {
                "path": name,
                "sha256": _git_blob_sha256(
                    historical_closure_commit,
                    f"tst/publication/frontier_control_plane/{name}",
                ),
            }
            for name in CONTROL_PLANE_FILES
        ]
        checkpoint_control_plane_version = inventory_digest(
            checkpoint_control_plane_bindings
        )
        self.assertRegex(checkpoint_control_plane_version, r"^[0-9a-f]{64}$")
        self.assertEqual(
            repair["control_plane_version"],
            checkpoint_control_plane_version,
        )
        self.assertEqual(
            repair["inventoried_file_count"],
            len(checkpoint_control_plane_bindings),
        )
        checkpoint_closure = [
            {
                "path": relative,
                "sha256": _git_blob_sha256(historical_closure_commit, relative),
            }
            for relative in _Q011_EXACT_PREDECESSOR_REPAIR_SOURCE_TEST_CLOSURE_PATHS
        ]
        closure = repair["source_test_closure"]
        checkpoint_closure_sha256 = inventory_digest(checkpoint_closure)
        self.assertRegex(checkpoint_closure_sha256, r"^[0-9a-f]{64}$")
        self.assertEqual(
            closure["files"],
            checkpoint_closure,
        )
        self.assertEqual(
            [binding["path"] for binding in closure["files"]],
            [binding["path"] for binding in checkpoint_closure],
        )
        self.assertEqual(
                repair["final_binding_refresh"],
                {
                    "status": "completed_before_commit",
                    "fields": ["control_plane_version", "source_test_closure"],
                    "authority": "none",
                },
        )

        next_worker = successor["next_clean_worker"]
        self.assertEqual(
            _git_blob_sha256(historical_closure_commit, next_worker["path"]),
            next_worker["sha256"],
        )
        publication_paths = subprocess.check_output(
            [
                "git",
                "ls-tree",
                "-r",
                "--name-only",
                historical_closure_commit,
                "--",
                "tst/publication",
            ],
            cwd=REPO_ROOT,
            text=True,
        ).splitlines()
        expected_counts = {
            "expected_publication_python_files": sum(
                path.endswith(".py") for path in publication_paths
            ),
            "expected_publication_shell_files": sum(
                path.endswith(".sh") for path in publication_paths
            ),
            "expected_publication_json_files": sum(
                path.endswith(".json") for path in publication_paths
            ),
            "expected_publication_test_modules": sum(
                Path(path).name.startswith("test_") and path.endswith(".py")
                for path in publication_paths
            ),
        }
        for key, observed in expected_counts.items():
            with self.subTest(next_clean_worker_count=key):
                self.assertEqual(next_worker[key], observed)

    def test_q011_exact_predecessor_migration_repair_successor_rejects_drift(
        self,
    ) -> None:
        successor = _load(
            "q011_section54_pressure_gate_exact_predecessor_migration_repair_"
            "successor_2026-06-05.json"
        )
        rejection_cases = [
            ("schema numeric alias", ("schema_version",), 1.0),
            ("predecessor", ("predecessor_sha256",), "0" * 64),
            (
                "retained worker log",
                ("retained_clean_worker", "log", "sha256"),
                "0" * 64,
            ),
            (
                "failed-attempt authority",
                ("retained_failed_migration_attempt", "authority"),
                "pair_install",
            ),
            (
                "failed build/freeze authority",
                ("retained_failed_build_freeze_attempt", "authority"),
                "build_evidence",
            ),
            (
                "live active policy",
                (
                    "unchanged_live_operational_baseline",
                    "active_policies",
                    0,
                    "sha256",
                ),
                "0" * 64,
            ),
            (
                "repair controller",
                ("source_local_exact_predecessor_repair", "control_plane_version"),
                "0" * 64,
            ),
            (
                "final-refresh status",
                (
                    "source_local_exact_predecessor_repair",
                    "final_binding_refresh",
                    "status",
                ),
                "complete",
            ),
            (
                "source/test closure path",
                (
                    "source_local_exact_predecessor_repair",
                    "source_test_closure",
                    "files",
                    0,
                    "path",
                ),
                "tst/publication/other.py",
            ),
            (
                "source/test closure member hash",
                (
                    "source_local_exact_predecessor_repair",
                    "source_test_closure",
                    "files",
                    0,
                    "sha256",
                ),
                "0" * 64,
            ),
            (
                "source/test closure digest",
                (
                    "source_local_exact_predecessor_repair",
                    "source_test_closure",
                    "closure_sha256",
                ),
                "0" * 64,
            ),
            ("next worker", ("next_clean_worker", "sha256"), "0" * 64),
            (
                "human pressure numeric alias",
                (
                    "human_pressure_selection",
                    "selected_case",
                    "problem_ps_p0",
                ),
                1,
            ),
            (
                "selection receipt invented",
                ("human_pressure_selection", "selection_receipt"),
                {},
            ),
            ("packet gate unblocked", ("packet_gate", "status"), "ready"),
            (
                "launch authorization granted",
                ("frontier_launch_authorization",),
                "authorized",
            ),
        ]
        for label, path, replacement in rejection_cases:
            with self.subTest(label=label):
                candidate = copy.deepcopy(successor)
                _replace_nested(candidate, path, replacement)
                with self.assertRaises(ValueError):
                    _validate_q011_exact_predecessor_migration_repair_successor(
                        candidate
                    )

        coordinated_closure_drift = copy.deepcopy(successor)
        coordinated_closure = coordinated_closure_drift[
            "source_local_exact_predecessor_repair"
        ]["source_test_closure"]
        coordinated_closure["files"][0]["sha256"] = "0" * 64
        coordinated_closure["closure_sha256"] = inventory_digest(
            coordinated_closure["files"]
        )
        with self.assertRaises(ValueError):
            _validate_q011_exact_predecessor_migration_repair_successor(
                coordinated_closure_drift
            )

        successor_with_extra_key = copy.deepcopy(successor)
        successor_with_extra_key["unexpected"] = True
        with self.assertRaises(ValueError):
            _validate_q011_exact_predecessor_migration_repair_successor(
                successor_with_extra_key
            )

    def test_q011_stage4_pressure_selection_candidate_successor_recomputes(
        self,
    ) -> None:
        successor_path = (
            READINESS_DIR
            / "q011_section54_pressure_selection_publication_candidate_"
            "successor_2026-06-05.json"
        )
        successor_payload = successor_path.read_text(encoding="utf-8")
        successor = json.loads(successor_payload)
        _validate_q011_stage4_pressure_selection_candidate_successor(successor)
        self.assertEqual(successor_payload, json.dumps(successor, indent=2) + "\n")

        predecessor_path = REPO_ROOT / successor["predecessor_record"]
        self.assertEqual(_sha256(predecessor_path), successor["predecessor_sha256"])
        predecessor = json.loads(predecessor_path.read_text(encoding="utf-8"))
        self.assertGreater(
            _canonical_utc_second(successor["recorded_utc"]),
            _canonical_utc_second(predecessor["recorded_utc"]),
        )

        active = successor["active_candidate_only_state"]
        for key in ("active_policy", "active_promotion"):
            binding = active[key]
            payloads = [
                Path(binding[path_key]).read_bytes()
                for path_key in ("orion_path", "project_home_path")
            ]
            self.assertEqual(len(set(payloads)), 1)
            self.assertEqual(hashlib.sha256(payloads[0]).hexdigest(), binding["sha256"])
        manifest_path = Path(active["clean_candidate"]["manifest_path"])
        self.assertEqual(
            _sha256(manifest_path),
            active["clean_candidate"]["manifest_sha256"],
        )
        source_archive = manifest_path.parent / "source.tar"
        self.assertEqual(
            _sha256(source_archive),
            active["clean_candidate"]["source_archive_sha256"],
        )
        self.assertFalse(
            (
                Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/publication")
                / "q011_section54_pressure_selection_receipt.json"
            ).exists()
        )
        for marker in (
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/pending_submission.json"),
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/pending_manual_accounting.json"),
            Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC/ledger/pending_manual_accounting.json"),
        ):
            self.assertFalse(marker.exists())

        closure = successor["stage4_publication_candidate"]["source_test_closure"]
        current_closure = [
            {"path": relative, "sha256": _sha256(REPO_ROOT / relative)}
            for relative in _Q011_STAGE4_PRESSURE_SELECTION_CANDIDATE_CLOSURE_PATHS
        ]
        self.assertEqual(closure["files"], current_closure)
        self.assertEqual(closure["closure_sha256"], inventory_digest(current_closure))
        worker = successor["next_clean_worker"]
        self.assertEqual(_sha256(REPO_ROOT / worker["path"]), worker["sha256"])

    def test_q011_stage4_pressure_selection_candidate_successor_rejects_drift(
        self,
    ) -> None:
        successor = _load(
            "q011_section54_pressure_selection_publication_candidate_"
            "successor_2026-06-05.json"
        )
        rejection_cases = [
            ("schema numeric alias", ("schema_version",), 1.0),
            ("predecessor", ("predecessor_sha256",), "0" * 64),
            (
                "active policy",
                ("active_candidate_only_state", "active_policy", "sha256"),
                "0" * 64,
            ),
            (
                "active slices",
                ("active_candidate_only_state", "registered_science_slices"),
                [{}],
            ),
            (
                "admission smoke",
                ("active_candidate_only_state", "frontier_admission_smoke", "status"),
                "authorized",
            ),
            (
                "selected pressure",
                (
                    "human_pressure_selection",
                    "selected_case",
                    "problem_ps_p0",
                ),
                0.1,
            ),
            (
                "production receipt invented",
                ("human_pressure_selection", "production_selection_receipt"),
                {},
            ),
            (
                "publisher execution mode",
                (
                    "stage4_publication_candidate",
                    "publisher_mutation_execution_mode",
                ),
                "trusted_checkout",
            ),
            (
                "source closure member",
                (
                    "stage4_publication_candidate",
                    "source_test_closure",
                    "files",
                    0,
                    "sha256",
                ),
                "0" * 64,
            ),
            ("worker", ("next_clean_worker", "sha256"), "0" * 64),
            ("packet gate", ("packet_gate", "status"), "ready"),
            ("launch authority", ("frontier_launch_authorization",), "authorized"),
        ]
        for label, path, replacement in rejection_cases:
            with self.subTest(label=label):
                candidate = copy.deepcopy(successor)
                _replace_nested(candidate, path, replacement)
                with self.assertRaises(ValueError):
                    _validate_q011_stage4_pressure_selection_candidate_successor(
                        candidate
                    )

        successor_with_extra_key = copy.deepcopy(successor)
        successor_with_extra_key["unexpected"] = True
        with self.assertRaises(ValueError):
            _validate_q011_stage4_pressure_selection_candidate_successor(
                successor_with_extra_key
            )

    def test_reviewed_mpich_stderr_fixture_matches_failed_attempt_provenance(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        superseded = successor["superseded_registered_execution"]
        provenance_path = (
            REPO_ROOT
            / superseded["reviewed_stderr_transcript_provenance"]
        )
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        failed = provenance["failed_attempt"]
        for key in (
            "job_id",
            "submission_id",
            "reservation_id",
            "registered_science_authorization_id",
            "run_artifact_dir",
        ):
            self.assertEqual(failed[key], superseded[key])
        fixture = provenance["local_review_fixture"]
        fixture_path = REPO_ROOT / fixture["path"]
        decoded = base64.b64decode(
            fixture_path.read_bytes().replace(b"\n", b""),
            validate=True,
        )
        stderr_entry = failed["stderr_inventory_entry"]
        inventory_path = REPO_ROOT / failed["artifact_inventory_fixture_path"]
        self.assertEqual(_sha256(inventory_path), failed["artifact_inventory_sha256"])
        inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
        inventory_records = {
            record["path"]: record for record in inventory["files"]
        }
        self.assertEqual(_sha256(fixture_path), fixture["encoded_file_sha256"])
        self.assertEqual(hashlib.sha256(decoded).hexdigest(), fixture["decoded_sha256"])
        self.assertEqual(len(decoded), fixture["decoded_size"])
        self.assertEqual(fixture["decoded_sha256"], stderr_entry["sha256"])
        self.assertEqual(fixture["decoded_size"], stderr_entry["size"])
        self.assertEqual(inventory_records[stderr_entry["path"]], stderr_entry)
        self.assertTrue(fixture["verified_byte_identical_to_live_immutable_stderr"])
        binding = provenance["reviewed_validator_binding"]
        commit = binding["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        for key, relative_key in {
            "support_module_sha256": "support_module",
            "gyro_analyzer_sha256": "gyro_analyzer",
            "paper_coupling_analyzer_sha256": "paper_coupling_analyzer",
        }.items():
            self.assertEqual(binding[key], _git_blob_sha256(commit, binding[relative_key]))
        self.assertEqual(
            provenance["disposition"],
            "pass_historical_transcript_bound_to_local_fixture_retry_requires_separate_v2_policy_activation",
        )

    def test_rejected_pre_reservation_manifest_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["rejected_pre_reservation_manifest_chronology"]
        fixture = _load(
            "q027_frontier_f1_rejected_pre_reservation_manifest_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["rejected_pre_reservation_manifest_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_rejected_pre_reservation_manifest_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        self.assertEqual(
            set(chronology),
            {
                "status",
                "submission_id",
                "registered_science_authorization_id",
                "control_plane_version",
                "manifest_path",
                "manifest_sha256",
                "manifest_analysis_support_sha256",
                "active_policy_analysis_support_sha256",
                "rejection",
                "pre_manifest_attestation_sha256",
                "pre_submit_wrapper_attestation_sha256",
                "reservation_attachments",
                "pending_submission_marker",
                "ledger_uuid_hits",
                "orion_ledger_records",
                "orion_receipt_records",
                "project_home_ledger_records",
                "active_reservations",
            },
        )
        self.assertEqual(
            chronology["status"],
            "fail_closed_before_reservation_intent_and_scheduler_submission",
        )
        self.assertEqual(
            chronology["active_policy_analysis_support_sha256"],
            "c5c77b6a952ed319498f08c91b9adc40101c090f91dc57d798dafdc455a38c01",
        )
        self.assertNotEqual(
            chronology["manifest_analysis_support_sha256"],
            chronology["active_policy_analysis_support_sha256"],
        )
        binding = fixture["manifest_binding"]
        self.assertEqual(binding["mode"], "0444")
        for key in (
            "submission_id",
            "control_plane_version",
            "registered_science_authorization_id",
            "manifest_sha256",
        ):
            self.assertEqual(binding[key], chronology[key])
        self.assertEqual(binding["analysis_support_role"], "analysis-script-001")
        self.assertEqual(
            binding["analysis_support_sha256"],
            chronology["manifest_analysis_support_sha256"],
        )
        self.assertEqual(chronology["reservation_attachments"], "absent")
        self.assertEqual(chronology["pending_submission_marker"], "absent")
        self.assertEqual(chronology["ledger_uuid_hits"], 0)
        self.assertEqual(chronology["orion_ledger_records"], 34)
        self.assertEqual(chronology["orion_receipt_records"], 34)
        self.assertEqual(chronology["project_home_ledger_records"], 34)
        self.assertEqual(chronology["active_reservations"], 0)
        template = _load(
            "q027_frontier_registered_science_same_account_isolation_attestation_template_2026-05-30.json"
        )
        attestations = fixture["attestations"]
        self.assertEqual(
            [record["contents"]["phase"] for record in attestations],
            ["pre_manifest", "pre_submit_wrapper"],
        )
        self.assertLess(
            attestations[0]["contents"]["recorded_utc"],
            attestations[1]["contents"]["recorded_utc"],
        )
        digest_keys = (
            "pre_manifest_attestation_sha256",
            "pre_submit_wrapper_attestation_sha256",
        )
        for record, digest_key in zip(attestations, digest_keys):
            contents = record["contents"]
            rendered = (json.dumps(contents, indent=2, sort_keys=True) + "\n").encode()
            self.assertEqual(
                hashlib.sha256(rendered).hexdigest(),
                chronology[digest_key],
            )
            self.assertEqual(record["attestation_sha256"], chronology[digest_key])
            self.assertEqual(contents["control_plane_version"], chronology["control_plane_version"])
            self.assertEqual(
                contents["registered_science_authorization_id"],
                chronology["registered_science_authorization_id"],
            )
            self.assertEqual(contents["operator_statement"], template["operator_statement"])
            self.assertEqual(contents["pending_submission_marker"]["value"], "absent")
            self.assertEqual(
                sorted(contents["mirrored_ledger_line_counts"]["counts"].values()),
                [34, 34, 34],
            )

    def test_rejected_pre_reservation_manifest_live_preflight_contract_is_explicit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["rejected_pre_reservation_manifest_chronology"]
        preflight = successor["required_live_preflight"]
        self.assertEqual(
            preflight["status"],
            "required_immediately_before_policy_promotion_and_each_registered_science_submission_boundary",
        )
        queue_chronology = successor[
            "rejected_pre_reservation_operator_queue_format_chronology"
        ]
        self.assertEqual(
            preflight["historical_submission_ids_must_remain_absent_from_live_ledgers"],
            [
                chronology["submission_id"],
                queue_chronology["submission_id"],
                successor["coupling_v2_submission_artifact_scrub_transition"][
                    "failed_manifest_creation_submission_id"
                ],
            ],
        )
        self.assertEqual(
            preflight["checks"],
            [
                "current Orion ledger, Orion mirror-receipt and Project Home mirror chains are coherent",
                "current Orion pending_submission.json marker is absent",
                "current mirrored ledger state has zero active reservations",
                "each historical rejected pre-reservation submission UUID has zero hits in current Orion JSONL, Orion CSV, Orion receipts and Project Home mirror streams",
                "same-account process and scheduler snapshots are reviewed for the current boundary",
            ],
        )
        self.assertIn("Legitimate later reservations", preflight["evidence_rule"])

    @unittest.skipUnless(
        os.environ.get("PIC_RUN_LIVE_PREFLIGHT") == "1",
        "set PIC_RUN_LIVE_PREFLIGHT=1 for Orion live-state checks",
    )
    def test_rejected_pre_reservation_manifest_live_preflight(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        preflight = successor["required_live_preflight"]
        ledger_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"
        )
        live_surfaces = (
            ledger_path,
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv"),
            Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl"),
            Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"),
        )
        line_counts = [len(path.read_text().splitlines()) for path in live_surfaces]
        self.assertEqual(line_counts[0], line_counts[2])
        self.assertEqual(line_counts[0], line_counts[3])
        self.assertEqual(line_counts[0] + 1, line_counts[1])
        self.assertFalse(
            Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/pending_submission.json"
            ).exists()
        )
        for marker in incomplete_manual_accounting_marker_paths(
            ledger_path, live_surfaces[3]
        ):
            self.assertFalse(marker.exists())
        ledger_events = validate_mirrored_state(
            ledger_path,
            live_surfaces[2],
            live_surfaces[3],
        )
        active_reservations = set()
        for event in ledger_events:
            if event["event_type"] == "reservation":
                active_reservations.add(event["reservation_id"])
            elif event["event_type"] in {"reconciliation", "reservation_cancelled"}:
                active_reservations.discard(event["reservation_id"])
        self.assertEqual(active_reservations, set())
        for surface in live_surfaces:
            contents = surface.read_text()
            for submission_id in preflight[
                "historical_submission_ids_must_remain_absent_from_live_ledgers"
            ]:
                self.assertNotIn(submission_id, contents)
        queue_fixture = _load(
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json"
        )
        rejected_manifest = Path(queue_fixture["chronology"]["manifest_path"])
        self.assertEqual(
            sorted(path.name for path in rejected_manifest.parent.iterdir()),
            ["pre_submit_manifest.json", "snapshot"],
        )
        gyro_fixture = _load(
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json"
        )
        self.assertIn(gyro_fixture["terminal_reconciliation_event"], ledger_events)

    def test_accepted_registered_f1_source_local_closure_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        closure_path = successor["source_local_accepted_closure_fixture"]
        closure = json.loads((REPO_ROOT / closure_path).read_text(encoding="utf-8"))
        executions = closure["registered_executions"]
        self.assertEqual(
            set(executions),
            {
                "accepted_gyro_v3_registered_execution",
                "accepted_paper_coupling_v2_registered_execution",
            },
        )
        for key, fixture in executions.items():
            execution = successor[key]
            documents = {}
            for name, record in fixture.items():
                path = REPO_ROOT / record["path"]
                self.assertEqual(_sha256(path), record["sha256"])
                documents[name] = json.loads(path.read_text(encoding="utf-8"))

            self.assertEqual(
                fixture["pre_submit_manifest"]["sha256"],
                execution["pre_submit_manifest_sha256"],
            )
            self.assertEqual(
                fixture["artifact_inventory"]["sha256"],
                execution["artifact_inventory_sha256"],
            )
            self.assertEqual(
                fixture["analysis_result"]["sha256"],
                execution["analysis_result_sha256"],
            )
            self.assertEqual(
                fixture["offline_analysis_receipt"]["sha256"],
                execution["offline_analysis_receipt_sha256"],
            )
            self.assertEqual(
                fixture["terminal_qualification_manifest"]["sha256"],
                execution["qualification_manifest_sha256"],
            )
            pre_submit_manifest = documents["pre_submit_manifest"]
            reconciliation = documents["terminal_reconciliation_event"]
            receipt = documents["terminal_reconciliation_mirror_receipt"]
            offline_receipt = documents["offline_analysis_receipt"]
            qualification = documents["terminal_qualification_manifest"]
            result = documents["analysis_result"]
            for field, expected in {
                "submission_id": execution["submission_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "artifact_dir": execution["run_artifact_dir"],
            }.items():
                self.assertEqual(pre_submit_manifest[field], expected)
            for field, expected in {
                "submission_id": execution["submission_id"],
                "reservation_id": execution["reservation_id"],
                "job_id": execution["job_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "manifest_sha256": execution["pre_submit_manifest_sha256"],
                "artifact_dir": execution["run_artifact_dir"],
                "consumed_node_hours": execution["consumed_node_hours"],
                "cumulative_consumed_node_hours":
                    execution["cumulative_consumed_node_hours"],
            }.items():
                self.assertEqual(reconciliation[field], expected)
            self.assertEqual(
                reconciliation["event_sha256"],
                record_sha256(reconciliation, "event_sha256"),
            )
            self.assertEqual(
                receipt["mirror_ack_sha256"],
                record_sha256(receipt, "mirror_ack_sha256"),
            )
            self.assertEqual(
                receipt["mirrored_event_sha256"],
                reconciliation["event_sha256"],
            )
            self.assertEqual(
                offline_receipt["artifact_inventory"]["sha256"],
                execution["artifact_inventory_sha256"],
            )
            self.assertEqual(
                offline_receipt["analysis_result"]["sha256"],
                execution["analysis_result_sha256"],
            )
            resources = qualification["resources"]
            for field, expected in {
                "submission_id": execution["submission_id"],
                "reservation_id": execution["reservation_id"],
                "job_id": execution["job_id"],
                "registered_science_authorization_id":
                    execution["registered_science_authorization_id"],
                "pre_submit_manifest_sha256": execution["pre_submit_manifest_sha256"],
                "artifact_inventory_sha256": execution["artifact_inventory_sha256"],
                "analysis_result_sha256": execution["analysis_result_sha256"],
                "offline_analysis_receipt_sha256":
                    execution["offline_analysis_receipt_sha256"],
                "node_hours": execution["consumed_node_hours"],
            }.items():
                self.assertEqual(resources[field], expected)
            self.assertEqual(
                qualification["authorization"]["active_policy"]["sha256"],
                execution["qualification_active_policy_sha256"],
            )
            self.assertEqual(
                qualification["authorization"]["active_promotion"]["sha256"],
                execution["qualification_active_promotion_sha256"],
            )
            if key == "accepted_gyro_v3_registered_execution":
                for field in [
                    "status",
                    "cycle",
                    "particle_count",
                    "max_abs_velocity_error",
                    "velocity_tolerance",
                ]:
                    self.assertEqual(result[field], execution["analysis"][field])
                self.assertEqual(
                    fixture["initial_qualification_manifest"]["sha256"],
                    execution["initial_qualification_manifest_sha256"],
                )
            else:
                coeff0 = result["cases"]["coeff0"]
                coeff7 = result["cases"]["coeff7"]
                self.assertEqual(
                    coeff0["max_abs_momentum_conservation_error"],
                    execution["analysis"]["coeff0_max_abs_momentum_conservation_error"],
                )
                self.assertEqual(
                    coeff0["abs_energy_conservation_error"],
                    execution["analysis"]["coeff0_abs_energy_conservation_error"],
                )
                self.assertEqual(
                    coeff7["max_abs_momentum_conservation_error"],
                    execution["analysis"]["coeff7_max_abs_momentum_conservation_error"],
                )
                self.assertEqual(
                    coeff7["abs_energy_conservation_error"],
                    execution["analysis"]["coeff7_abs_energy_conservation_error"],
                )
                self.assertEqual(
                    result["coefficient_invariance"]["max_abs_particle_momentum_error"],
                    execution["analysis"][
                        "max_abs_particle_momentum_coefficient_invariance_error"
                    ],
                )

    def test_accepted_registered_f2_source_local_closure_is_bound(self) -> None:
        candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        execution = candidate["accepted_v2_execution"]
        closure_path = candidate["source_local_accepted_closure_fixture"]
        closure = json.loads((REPO_ROOT / closure_path).read_text(encoding="utf-8"))
        fixture = closure["accepted_v2_registered_execution"]
        self.assertEqual(fixture["submission_id"], execution["submission_id"])
        self.assertEqual(fixture["reservation_id"], execution["reservation_id"])
        self.assertEqual(fixture["job_id"], execution["job_id"])
        self.assertEqual(
            fixture["registered_science_authorization_id"],
            execution["authorization_id"],
        )
        self.assertEqual(
            set(fixture["documents"]),
            {
                "pre_submit_manifest",
                "artifact_inventory",
                "analysis_result",
                "offline_analysis_receipt",
                "terminal_qualification_manifest",
                "terminal_reconciliation_event",
                "terminal_reconciliation_mirror_receipt",
                "pre_policy_promotion_attestation",
                "pre_manifest_attestation",
                "pre_submit_wrapper_attestation",
            },
        )
        documents = {}
        for name, record in fixture["documents"].items():
            path = REPO_ROOT / record["path"]
            self.assertEqual(_sha256(path), record["sha256"])
            documents[name] = json.loads(path.read_text(encoding="utf-8"))

        for name, digest_key in {
            "pre_submit_manifest": "pre_submit_manifest_sha256",
            "artifact_inventory": "artifact_inventory_sha256",
            "analysis_result": "analysis_result_sha256",
            "offline_analysis_receipt": "offline_analysis_receipt_sha256",
            "terminal_qualification_manifest": "qualification_manifest_sha256",
            "pre_policy_promotion_attestation":
                "pre_policy_promotion_attestation_sha256",
            "pre_manifest_attestation": "pre_manifest_attestation_sha256",
            "pre_submit_wrapper_attestation":
                "pre_submit_wrapper_attestation_sha256",
        }.items():
            self.assertEqual(
                fixture["documents"][name]["sha256"],
                execution[digest_key],
            )

        self.assertEqual(
            fixture["contained_output_path"], execution["contained_output_path"]
        )
        self.assertEqual(
            fixture["contained_output_sha256"], execution["contained_output_sha256"]
        )
        manifest = documents["pre_submit_manifest"]
        self.assertEqual(manifest["submission_id"], execution["submission_id"])
        self.assertEqual(
            manifest["registered_science_authorization_id"],
            execution["authorization_id"],
        )
        self.assertEqual(manifest["artifact_dir"], fixture["run_artifact_dir"])
        inventory = {
            record["path"]: record
            for record in documents["artifact_inventory"]["files"]
        }
        contained_output = inventory[fixture["contained_output_path"]]
        self.assertEqual(
            contained_output["sha256"], fixture["contained_output_sha256"]
        )
        analysis = documents["analysis_result"]
        self.assertEqual(analysis["status"], "pass")
        self.assertEqual(analysis["parallel_ranks"], 8)
        self.assertEqual(len(analysis["hosts"]), 1)
        self.assertEqual(len(analysis["rank_gpu_bindings"]), 8)
        self.assertEqual(
            len({
                binding["rocr_visible_device"]
                for binding in analysis["rank_gpu_bindings"]
            }),
            8,
        )
        self.assertEqual(
            analysis["runtime_artifacts"][fixture["contained_output_path"]],
            fixture["contained_output_sha256"],
        )
        receipt = documents["offline_analysis_receipt"]
        self.assertEqual(
            receipt["artifact_inventory"]["sha256"],
            execution["artifact_inventory_sha256"],
        )
        self.assertEqual(
            receipt["analysis_result"]["sha256"],
            execution["analysis_result_sha256"],
        )
        qualification = documents["terminal_qualification_manifest"]
        resources = qualification["resources"]
        for field, expected in {
            "submission_id": execution["submission_id"],
            "reservation_id": execution["reservation_id"],
            "job_id": execution["job_id"],
            "registered_science_authorization_id": execution["authorization_id"],
            "pre_submit_manifest_sha256": execution["pre_submit_manifest_sha256"],
            "artifact_inventory_sha256": execution["artifact_inventory_sha256"],
            "analysis_result_sha256": execution["analysis_result_sha256"],
            "offline_analysis_receipt_sha256":
                execution["offline_analysis_receipt_sha256"],
            "node_hours": execution["consumed_node_hours"],
        }.items():
            self.assertEqual(resources[field], expected)
        self.assertEqual(
            qualification["authorization"]["active_policy"]["sha256"],
            execution["active_policy_sha256"],
        )
        self.assertEqual(
            qualification["authorization"]["active_promotion"]["sha256"],
            execution["active_promotion_sha256"],
        )
        reconciliation = documents["terminal_reconciliation_event"]
        for field, expected in {
            "event_type": "reconciliation",
            "state": execution["state"],
            "submission_id": execution["submission_id"],
            "reservation_id": execution["reservation_id"],
            "job_id": execution["job_id"],
            "registered_science_authorization_id": execution["authorization_id"],
            "manifest_sha256": execution["pre_submit_manifest_sha256"],
            "artifact_dir": fixture["run_artifact_dir"],
            "active_policy_sha256": execution["active_policy_sha256"],
            "active_promotion_sha256": execution["active_promotion_sha256"],
            "consumed_node_hours": execution["consumed_node_hours"],
            "cumulative_consumed_node_hours":
                execution["cumulative_consumed_node_hours"],
        }.items():
            self.assertEqual(reconciliation[field], expected)
        self.assertEqual(
            reconciliation["event_sha256"],
            record_sha256(reconciliation, "event_sha256"),
        )
        self.assertEqual(
            reconciliation["event_sha256"],
            execution["terminal_reconciliation_event_sha256"],
        )
        mirror_receipt = documents["terminal_reconciliation_mirror_receipt"]
        self.assertEqual(
            mirror_receipt["mirror_ack_sha256"],
            record_sha256(mirror_receipt, "mirror_ack_sha256"),
        )
        self.assertEqual(
            mirror_receipt["mirror_ack_sha256"],
            execution["terminal_reconciliation_mirror_ack_sha256"],
        )
        self.assertEqual(
            mirror_receipt["mirrored_event_sha256"],
            reconciliation["event_sha256"],
        )

        chronology = candidate["rejected_v2_pre_reservation_chronology"]
        rejected = closure["rejected_pre_reservation_chronology"]
        rejected_manifest = rejected["queue_format_manifest"]
        self.assertEqual(rejected_manifest["submission_id"], chronology["submission_id"])
        self.assertEqual(
            rejected_manifest["pre_submit_manifest"]["sha256"],
            chronology["pre_submit_manifest_sha256"],
        )
        self.assertEqual(
            rejected_manifest["attestations"][0]["sha256"],
            chronology["pre_manifest_attestation_sha256"],
        )
        self.assertEqual(
            rejected_manifest["attestations"][1]["sha256"],
            chronology["pre_submit_wrapper_attestation_sha256"],
        )
        queue_drift_attestation = rejected["transient_queue_drift"]["attestation"]
        self.assertEqual(
            queue_drift_attestation["sha256"],
            chronology["transient_queue_drift_retry_attestation_sha256"],
        )
        for record in [
            rejected_manifest["pre_submit_manifest"],
            *rejected_manifest["attestations"],
            queue_drift_attestation,
        ]:
            self.assertEqual(_sha256(REPO_ROOT / record["path"]), record["sha256"])
        self.assertEqual(chronology["reservation_attachments"], "absent")
        self.assertEqual(chronology["ledger_intent"], "absent")
        self.assertEqual(chronology["scheduler_submission"], "absent")
        terminal = closure["active_terminal_ledger"]
        self.assertEqual(terminal["records"], execution["terminal_ledger_records"])
        self.assertEqual(
            terminal["cumulative_consumed_node_hours"],
            execution["cumulative_consumed_node_hours"],
        )
        self.assertEqual(terminal["active_reservations"], 0)
        self.assertEqual(terminal["pending_submission_marker"], "absent")

    @unittest.skipUnless(
        os.environ.get("PIC_RUN_LIVE_PREFLIGHT") == "1",
        "set PIC_RUN_LIVE_PREFLIGHT=1 for Orion live-state checks",
    )
    def test_registered_terminal_review_manifests_replay_live(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        f2_candidate = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )
        manual_accounting_activation = _load(
            "q027_manual_frontier_accounting_activation_2026-05-30.json"
        )
        f2 = f2_candidate["accepted_v2_execution"]
        active_policy_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/policy/storage_policy.json"
        )
        active_promotion_path = Path(
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/policy/active_promotion.json"
        )
        current_policy_sha256 = _sha256(active_policy_path)
        current_promotion_sha256 = _sha256(active_promotion_path)
        current_transition = manual_accounting_activation["control_plane_transition"]
        paired_transition = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v3_2026-06-01.json"
        )
        paired_promotion = paired_transition["active_policy_promotion"]
        c83e_transition = _load(
            "phase0_paired_control_plane_install_and_policy_promotion_"
            "successor_v4_2026-06-02.json"
        )
        c83e_promotion = c83e_transition["active_policy_promotion"]
        d720_transition = _load(
            "phase0_clean_candidate_freeze_and_policy_promotion_"
            "successor_v2_2026-06-02.json"
        )
        d720_promotion = d720_transition["active_policy_promotion"]
        if current_policy_sha256 == c83e_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, c83e_promotion["orion_promotion_sha256"]
            )
            terminal = c83e_transition["terminal_mirrored_ledger"]
        elif current_policy_sha256 == d720_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, d720_promotion["orion_promotion_sha256"]
            )
            terminal = d720_transition["mirrored_ledger_invariant"]
        elif current_policy_sha256 == paired_promotion["orion_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, paired_promotion["orion_promotion_sha256"]
            )
            terminal = paired_transition["terminal_mirrored_ledger"]
        elif current_policy_sha256 == current_transition["active_policy_sha256"]:
            self.assertEqual(
                current_promotion_sha256, current_transition["active_promotion_sha256"]
            )
            terminal = manual_accounting_activation["terminal_ledger"]
        else:
            active_policy = json.loads(active_policy_path.read_text(encoding="utf-8"))
            active_promotion = json.loads(
                active_promotion_path.read_text(encoding="utf-8")
            )
            self.assertEqual(
                active_promotion["policy_sha256"],
                current_policy_sha256,
            )
            project_home_policy_path = Path(active_promotion["project_home_policy_path"])
            self.assertEqual(_sha256(project_home_policy_path), current_policy_sha256)
            self.assertEqual(
                json.loads(project_home_policy_path.read_text(encoding="utf-8")),
                active_policy,
            )
            validate_mirrored_state(
                Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"),
                Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl"),
                Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"),
            )
            terminal = {
                "orion_node_hours_jsonl_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl")
                ),
                "orion_node_hours_csv_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv")
                ),
                "orion_mirror_receipts_jsonl_sha256": _sha256(
                    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/mirror_receipts.jsonl")
                ),
                "project_home_node_hours_jsonl_sha256": _sha256(
                    Path("/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl")
                ),
            }
        for digest_key, path in {
            "orion_node_hours_jsonl_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.jsonl"
            ),
            "orion_node_hours_csv_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/node_hours.csv"
            ),
            "orion_mirror_receipts_jsonl_sha256": Path(
                "/lustre/orion/ast207/proj-shared/dfielding/PIC/ledger/"
                "mirror_receipts.jsonl"
            ),
            "project_home_node_hours_jsonl_sha256": Path(
                "/ccs/proj/ast207/proj-shared/PIC/ledger/node_hours.jsonl"
            ),
        }.items():
            expected_key = digest_key
            if terminal is manual_accounting_activation["terminal_ledger"]:
                expected_key = {
                    "orion_node_hours_jsonl_sha256": "orion_jsonl_sha256",
                    "orion_node_hours_csv_sha256": "orion_csv_sha256",
                    "orion_mirror_receipts_jsonl_sha256":
                        "orion_mirror_receipts_sha256",
                    "project_home_node_hours_jsonl_sha256":
                        "project_home_jsonl_sha256",
                }[digest_key]
            self.assertEqual(_sha256(path), terminal[expected_key])
        for key in [
            "accepted_gyro_v3_registered_execution",
            "accepted_paper_coupling_v2_registered_execution",
        ]:
            with self.subTest(execution=key):
                execution = successor[key]
                manifest_path = Path(execution["qualification_manifest_path"])
                with _pinned_regular_bytes(manifest_path) as manifest_bytes:
                    self.assertEqual(
                        hashlib.sha256(manifest_bytes).hexdigest(),
                        execution["qualification_manifest_sha256"],
                    )
                    manifest = json.loads(manifest_bytes)
                for path_key, sha256_key in {
                    "pre_submit_manifest_path": "pre_submit_manifest_sha256",
                    "artifact_inventory_path": "artifact_inventory_sha256",
                    "analysis_result_path": "analysis_result_sha256",
                    "offline_analysis_receipt_path":
                        "offline_analysis_receipt_sha256",
                }.items():
                    with _pinned_regular_bytes(
                        Path(manifest["resources"][path_key])
                    ) as resource_bytes:
                        self.assertEqual(
                            hashlib.sha256(resource_bytes).hexdigest(),
                            manifest["resources"][sha256_key],
                        )
                manifest_policy_sha256 = manifest["authorization"][
                    "active_policy"
                ]["sha256"]
                manifest_promotion_sha256 = manifest["authorization"][
                    "active_promotion"
                ]["sha256"]
                self.assertEqual(
                    manifest_policy_sha256,
                    execution["qualification_active_policy_sha256"],
                )
                self.assertEqual(
                    manifest_promotion_sha256,
                    execution["qualification_active_promotion_sha256"],
                )
                validate_schema(manifest, VALIDATION_MANIFEST_SCHEMA)
        with _pinned_regular_bytes(
            Path(f2["qualification_manifest_path"])
        ) as manifest_bytes:
            self.assertEqual(
                hashlib.sha256(manifest_bytes).hexdigest(),
                f2["qualification_manifest_sha256"],
            )
            validate_schema(json.loads(manifest_bytes), VALIDATION_MANIFEST_SCHEMA)

    def test_registered_parser_policy_transition_resolves_source_commit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        transition = successor["gyro_v3_parser_policy_transition"]
        commit = transition["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        slices = {
            record["authorization_id"]: record
            for record in transition["registered_science_slices"]
        }
        for authorization_id, analyzer in {
            "f1-clean-gyro-mpich-stderr-v3":
                "tst/publication/frontier_f1_gpu_relativistic_gyro_analysis.py",
            "f1-clean-paper-coupling-mpich-stderr-v2":
                "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
        }.items():
            self.assertEqual(
                slices[authorization_id]["analysis_script_sha256"],
                _git_blob_sha256(commit, analyzer),
            )
            self.assertEqual(
                slices[authorization_id]["analysis_support_sha256"],
                _git_blob_sha256(
                    commit,
                    "tst/publication/frontier_f1_structured_artifacts.py",
                ),
            )

    def test_rejected_operator_queue_format_manifest_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor[
            "rejected_pre_reservation_operator_queue_format_chronology"
        ]
        fixture = _load(
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["rejected_pre_reservation_operator_queue_format_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_rejected_operator_queue_format_manifest_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        binding = fixture["manifest_binding"]
        manifest_path = REPO_ROOT / binding["local_fixture_path"]
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(_sha256(manifest_path), chronology["manifest_sha256"])
        self.assertEqual(manifest["submission_id"], chronology["submission_id"])
        self.assertEqual(
            manifest["queue_snapshot_sha256"], binding["queue_snapshot_sha256"]
        )
        self.assertEqual(
            binding["queue_snapshot_format"], "%i|%a|%P|%q|%T|%j|%k"
        )
        self.assertEqual(
            binding["validator_queue_snapshot_format"], "%i|%P|%q|%T|%j|%k"
        )
        self.assertNotEqual(
            binding["queue_snapshot_sha256"],
            binding["validator_queue_snapshot_sha256"],
        )
        self.assertEqual(chronology["reservation_attachments"], "absent")
        for attestation, digest_key in zip(
            fixture["attestations"],
            (
                "pre_manifest_attestation_sha256",
                "pre_submit_wrapper_attestation_sha256",
            ),
        ):
            path = REPO_ROOT / attestation["local_fixture_path"]
            self.assertEqual(attestation["sha256"], chronology[digest_key])
            self.assertEqual(_sha256(path), chronology[digest_key])
            contents = json.loads(path.read_text(encoding="utf-8"))
            self.assertEqual(contents["phase"], attestation["phase"])
            self.assertEqual(contents["pending_submission_marker"]["value"], "absent")
            self.assertEqual(
                sorted(contents["mirrored_ledger_line_counts"]["counts"].values()),
                [
                    chronology["orion_ledger_records"],
                    chronology["orion_receipt_records"],
                    chronology["project_home_ledger_records"],
                ],
            )
            self.assertEqual(
                contents["queue_snapshot"]["sha256"],
                binding["queue_snapshot_sha256"],
            )

    def test_coupling_scrub_transition_resolves_source_commit(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        transition = successor["coupling_v2_submission_artifact_scrub_transition"]
        commit = transition["source_commit"]
        self.assertEqual(
            subprocess.check_output(
                ["git", "cat-file", "-t", commit],
                cwd=REPO_ROOT,
                text=True,
            ).strip(),
            "commit",
        )
        self.assertEqual(
            transition["scrub_safe_analysis_script_sha256"],
            _git_blob_sha256(
                commit,
                "tst/publication/frontier_f1_gpu_paper_coupling_analysis.py",
            ),
        )
        self.assertEqual(
            transition["active_policy_sha256"],
            _git_blob_sha256(
                commit,
                "tst/publication/readiness/storage_policy.json",
            ),
        )
        gyro = successor["accepted_gyro_v3_registered_execution"]
        coupling = successor["accepted_paper_coupling_v2_registered_execution"]
        f2 = _load(
            "q027_frontier_f2_multirank_runtime_metadata_candidate_2026-05-30.json"
        )["accepted_v2_execution"]
        self.assertNotEqual(
            gyro["initial_qualification_manifest_sha256"],
            gyro["qualification_manifest_sha256"],
        )
        self.assertEqual(
            coupling["active_policy_sha256"],
            transition["active_policy_sha256"],
        )
        self.assertEqual(
            coupling["active_promotion_sha256"],
            transition["active_promotion_sha256"],
        )
        for execution in (gyro, coupling):
            self.assertEqual(
                execution["qualification_active_policy_sha256"],
                f2["active_policy_sha256"],
            )
            self.assertEqual(
                execution["qualification_active_promotion_sha256"],
                f2["active_promotion_sha256"],
            )

    def test_reconciled_gyro_v2_analysis_rejection_chronology_is_bound(self) -> None:
        successor = _load(
            "q027_frontier_f1_registered_science_successor_candidate_2026-05-30.json"
        )
        chronology = successor["superseded_gyro_v2_registered_execution"]
        fixture = _load(
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json"
        )
        self.assertEqual(
            successor["superseded_gyro_v2_registered_execution_fixture"],
            "tst/publication/readiness/"
            "q027_frontier_f1_gyro_v2_analysis_rejection_fixture_2026-05-30.json",
        )
        self.assertEqual(chronology, fixture["chronology"])
        inventory_path = REPO_ROOT / fixture["artifact_inventory_fixture_path"]
        self.assertEqual(_sha256(inventory_path), chronology["artifact_inventory_sha256"])
        self.assertEqual(
            json.loads(inventory_path.read_text(encoding="utf-8")),
            fixture["artifact_inventory"],
        )
        output_paths = [
            record["path"]
            for record in fixture["artifact_inventory"]["files"]
            if record["path"].startswith("output/")
        ]
        self.assertEqual(len(output_paths), 21)
        self.assertEqual(
            output_paths,
            sorted(output_paths),
        )
        event = fixture["terminal_reconciliation_event"]
        self.assertEqual(event["submission_id"], chronology["submission_id"])
        self.assertEqual(event["reservation_id"], chronology["reservation_id"])
        self.assertEqual(event["job_id"], chronology["job_id"])
        self.assertEqual(event["manifest_sha256"], chronology["pre_submit_manifest_sha256"])
        self.assertEqual(event["consumed_node_hours"], chronology["consumed_node_hours"])
        self.assertEqual(
            event["event_sha256"],
            record_sha256(event, "event_sha256"),
        )
        receipt_fixture = fixture["terminal_reconciliation_mirror_receipt_fixture"]
        receipt_path = REPO_ROOT / receipt_fixture["path"]
        self.assertEqual(_sha256(receipt_path), receipt_fixture["sha256"])
        receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
        self.assertEqual(
            receipt["mirror_ack_sha256"],
            record_sha256(receipt, "mirror_ack_sha256"),
        )
        self.assertEqual(receipt["mirrored_event_sha256"], event["event_sha256"])
        self.assertEqual(chronology["analysis_result"], "absent")
        self.assertEqual(chronology["offline_analysis_receipt"], "absent")

    def test_registered_science_same_account_isolation_attestation_template(self) -> None:
        template = _load(
            "q027_frontier_registered_science_same_account_isolation_attestation_template_2026-05-30.json"
        )
        self.assertEqual(
            template["archive_root"],
            "/lustre/orion/ast207/proj-shared/dfielding/PIC/operator_attestations",
        )
        self.assertEqual(
            template["status"],
            "required_before_each_registered_science_submission",
        )
        self.assertIn("same_account_process_snapshot", template["required_fields"])
        self.assertIn("operator_statement", template["required_fields"])
        self.assertIn(
            "test -e \"${PIC_ROOT}/ledger/pending_submission.json\" "
            "&& printf 'present\\n' || printf 'absent\\n'",
            template["required_snapshot_commands"],
        )

    def test_pic_vl2_tsc_phase1_stage_ordering_review_receipt_is_bounded(self) -> None:
        receipt = _load(
            "pic_vl2_tsc_phase1_stage_ordering_review_receipt_successor_"
            "2026-06-02.json"
        )
        predecessor = receipt["predecessor"]
        predecessor_path = REPO_ROOT / predecessor["path"]
        self.assertEqual(predecessor["sha256"], _sha256(predecessor_path))
        registration = json.loads(predecessor_path.read_text(encoding="utf-8"))
        self.assertEqual(
            registration["identity_split"]["successor_mode"],
            receipt["physical_mode"],
        )
        self.assertEqual(receipt["phase"], "phase1")
        self.assertEqual(receipt["finding_id"], "PIC-P1-002")
        self.assertEqual(receipt["verification_gate"], "Q-004")
        self.assertIn("Bounded source-local", receipt["scope"])

        trace = receipt["reused_stage_trace"]
        self.assertEqual(trace["sha256"], _sha256(REPO_ROOT / trace["path"]))
        self.assertEqual(
            trace["sha256"],
            registration["active_successor_files"][trace["path"]],
        )
        self.assertIn(
            trace["listed_pass_result"],
            registration["local_validation"]["serial_oracles_passed"],
        )
        self.assertTrue(trace["reuse_only_no_new_dynamic_execution_or_artifact"])

        for binding_group in ("source_bindings", "script_bindings", "fixture_bindings"):
            for relative, expected_sha256 in receipt[binding_group].items():
                with self.subTest(binding_group=binding_group, relative=relative):
                    self.assertEqual(expected_sha256, _sha256(REPO_ROOT / relative))
                    predecessor_sha256 = registration["active_successor_files"].get(relative)
                    if predecessor_sha256 is not None:
                        self.assertEqual(expected_sha256, predecessor_sha256)

        review = receipt["stage_ordering_review"]
        expected_sequences = {
            "push": [
                "stage_1_inserted_push_is_predictor_no_op_at_x_ini",
                "stage_1_post_deposit_half_step_drift_reaches_x_mid",
                "stage_2_inserted_push_applies_midpoint_boris_kick_at_x_mid",
                "stage_2_post_deposit_half_step_drift_reaches_x_end",
            ],
            "deposition": [
                "Particles::Push",
                "Particles::SaveOldPositions",
                "Particles::ZeroMoments",
                "Particles::InitRecvMoments",
                "Particles::DepositMoments",
            ],
            "feedback_placement": [
                "MHD::RKUpdate",
                "paper_vl2_staged_particle_wrapper_chain",
                "MHD::MHDSrcTerms",
            ],
            "boundary_synchronization": [
                "Particles::DepositMoments",
                "Particles::RestrictMoments",
                "Particles::SendMoments",
                "Particles::RecvMoments",
                "Particles::ClearRecvMoments",
                "Particles::ClearSendMoments",
                "Particles::ApplyMomentPhysicalBCs",
                "Particles::ProlongateMoments",
            ],
            "migration_communication_ordering": [
                "Particles::DriftPaperCosmicRaysHalfStep",
                "after_stagen",
                "Particles::NewGID",
                "Particles::SendCnt",
                "Particles::InitRecv",
                "Particles::SendP",
                "Particles::RecvP",
                "Particles::ClearRecv",
                "Particles::ClearSend",
            ],
            "ct_ordering": [
                "MHD::CornerE",
                "MHD::EFieldSrc",
                "MHD::SendE",
                "MHD::RecvE",
                "MHD::CT",
            ],
        }
        self.assertEqual(
            set(review),
            set(expected_sequences),
        )
        for topic_name, topic in review.items():
            self.assertEqual(topic["status"], "covered_bounded_source_local_review")
            self.assertEqual(topic["sequence"], expected_sequences[topic_name])
        self.assertIn(
            "excluded from AddsCRCurrentToCT",
            review["ct_ordering"]["interpretation"],
        )

        boundary = receipt["qualification_boundary"]
        self.assertTrue(boundary["bounded_source_local_scope"])
        for key in (
            "adds_dynamic_evidence",
            "frontier_execution_authorized",
            "frontier_qualified",
            "mpi_qualified",
            "hip_qualified",
            "claim_closure",
        ):
            self.assertFalse(boundary[key])

    def test_claim_classes_and_required_extensions(self) -> None:
        registry = _load("claims_registry.json")
        self.assertEqual(registry["default_reviewer"], "pending external review")
        claims = registry["claims"]
        ids = {claim["claim_id"] for claim in claims}
        self.assertEqual(len(ids), len(claims))
        self.assertTrue(REQUIRED_EXTENSION_CLAIMS.issubset(ids))
        for claim in claims:
            self.assertIn(claim["claim_class"], CLAIM_CLASSES)
            self.assertEqual(claim["disposition"], "open")
            self.assertTrue(claim["required_gates"])
            if claim["claim_id"] in REQUIRED_EXTENSION_CLAIMS:
                self.assertTrue(claim["authorized_extension_required"])

    def test_initial_findings_are_unique_and_have_valid_status(self) -> None:
        registry = _load("findings_registry.json")
        findings = registry["findings"]
        ids = {finding["finding_id"] for finding in findings}
        self.assertEqual(len(ids), len(findings))
        self.assertEqual(ids, {f"PIC-P0-00{i}" for i in range(1, 7)} | {
            f"PIC-P1-{i:03d}" for i in range(1, 11)
        })
        expected_verifying = {
            "PIC-P0-001", "PIC-P0-002", "PIC-P0-003", "PIC-P0-004",
            "PIC-P0-005", "PIC-P0-006", "PIC-P1-002", "PIC-P1-004",
            "PIC-P1-008", "PIC-P1-009", "PIC-P1-010",
        }
        for finding in findings:
            expected = "verifying" if finding["finding_id"] in expected_verifying else "open"
            self.assertEqual(finding["status"], expected)
            self.assertRegex(finding["verification_gate"], r"^Q-\d{3}$")

    def test_paper_source_checksum_is_frozen(self) -> None:
        inventory = _load("external_artifacts.json")
        artifacts = {
            artifact["artifact_id"]: artifact for artifact in inventory["artifacts"]
        }
        tex = artifacts["SUN_BAI_2023_ARXIV_V1_TEX"]
        self.assertEqual(tex["sha256"], _sha256(PAPER_TEX))
        entity = artifacts["ENTITY_TOOLKIT_REPLACEMENT_CANDIDATE"]
        self.assertEqual(
            entity["git_commit"],
            "512998c471bf3fdec292cb4a64150c4f0aeea539",
        )
        self.assertEqual(entity["historical_unavailable_commit"], "a59065fc")
        snapshot = _load("entity_snapshot.json")
        self.assertEqual(entity["snapshot_manifest"],
                         "tst/publication/readiness/entity_snapshot.json")
        self.assertEqual(snapshot["git_commit"], entity["git_commit"])
        self.assertEqual(snapshot["git_tree"], entity["git_tree"])
        self.assertEqual(snapshot["git_archive_sha256"],
                         entity["git_archive_sha256"])
        self.assertEqual(snapshot["worktree_status"], "clean")
        self.assertTrue(snapshot["files"])
        for source_file in snapshot["files"]:
            self.assertRegex(source_file["sha256"], r"^[0-9a-f]{64}$")

    def test_validation_manifest_schema_has_release_minimum(self) -> None:
        schema = _validation_manifest_schema()
        required = set(schema["required"])
        self.assertTrue(
            {
                "claim_ids",
                "evidence_class",
                "physical_mode",
                "git",
                "executable",
                "oracle",
                "metrics",
                "resources",
                "artifacts",
                "review",
            }.issubset(required)
        )
        classes = set(schema["properties"]["evidence_class"]["enum"])
        self.assertEqual(classes, CLAIM_CLASSES)

    def test_q011_repaired_clean_candidate_freeze_receipt_is_frozen(self) -> None:
        receipt = _load("phase0_clean_candidate_freeze_successor_v3_2026-06-02.json")
        self.assertEqual(
            receipt["predecessor_sha256"],
            _sha256(REPO_ROOT / receipt["predecessor_record"]),
        )
        candidate = receipt["clean_candidate"]
        manifest_path = Path(candidate["manifest_path"])
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        self.assertEqual(candidate["manifest_sha256"], _sha256(manifest_path))
        self.assertEqual(candidate["manifest_sha256"], "ef527ed467995bd60fda07b5a3b09b56ea871595ace12fd64a948e246720dbe3")
        self.assertEqual(stat.S_IMODE(manifest_path.parent.stat().st_mode) & 0o222, 0)
        descendants = list(manifest_path.parent.rglob("*"))
        self.assertEqual(
            receipt["local_validation"]["expected_tree_entries"],
            len(descendants),
        )
        for descendant in descendants:
            self.assertFalse(descendant.is_symlink())
            self.assertEqual(stat.S_IMODE(descendant.stat().st_mode) & 0o222, 0)
        self.assertEqual(candidate["source_git_commit"], manifest["source"]["git_commit"])
        self.assertEqual(candidate["source_git_tree"], manifest["source"]["git_tree"])
        self.assertEqual(candidate["source_archive_sha256"], manifest["source"]["archive_sha256"])
        self.assertEqual(candidate["source_commit_sha256"], manifest["source"]["commit_sha256"])
        self.assertEqual(candidate["source_bundle_sha256"], manifest["source"]["source_bundle_sha256"])
        self.assertEqual(candidate["prepared_artifact_inventory_sha256"], manifest["prepared_artifacts"]["inventory_sha256"])
        self.assertEqual(candidate["paper_deck_count"], len(manifest["prepared_artifacts"]["paper_decks"]))
        self.assertEqual(candidate["publication_analyzer_count"], len(manifest["prepared_artifacts"]["analyzers"]))
        self.assertEqual(candidate["build_profile_sha256"], _sha256(Path(manifest["build"]["profile_path"])))
        self.assertEqual(candidate["profile_receipt_sha256"], _sha256(Path(manifest["build"]["profile_receipt_path"])))
        self.assertEqual(candidate["executable_sha256"], _sha256(Path(manifest["build"]["executable_path"])))
        fragment = receipt["post_freeze_policy_fragment"]
        self.assertEqual(fragment["sha256"], _sha256(Path(fragment["path"])))
        self.assertEqual(fragment["status"], "non_authorizing_review_fragment_only")
        fragment_payload = json.loads(Path(fragment["path"]).read_text(encoding="utf-8"))
        self.assertEqual(fragment["registered_science_slice_count"], len(fragment_payload["registered_science_slices"]))
        self.assertEqual(fragment["registered_science_slice_count"], 4)
        policy = receipt["policy_promotion"]
        self.assertEqual(policy["registered_science_slices"], [])
        self.assertEqual(
            policy["live_orion_policy_sha256"],
            "1a331137c7d83717a890fa64046890074101f0fec3de7b8b22ca41b8bb644d28",
        )
        self.assertEqual(
            policy["live_project_home_policy_sha256"],
            policy["live_orion_policy_sha256"],
        )
        self.assertEqual(
            policy["live_orion_promotion_sha256"],
            "e64db1ef1ca755e1fcf4ff86e078856e381ff31f6aab926f01658b7aa69cd73e",
        )
        self.assertEqual(
            policy["live_project_home_promotion_sha256"],
            policy["live_orion_promotion_sha256"],
        )
        successor = _load("phase0_curated_candidate_successor_v19_2026-06-02.json")
        self.assertEqual(successor["predecessor_sha256"], _sha256(REPO_ROOT / successor["predecessor_record"]))
        successor_receipt = successor["clean_candidate_freeze_receipt"]
        self.assertEqual(
            successor_receipt["sha256"],
            _sha256(REPO_ROOT / successor_receipt["path"]),
        )
        self.assertEqual(successor["frontier_launch_authorization"], "none_no_registered_science_slices")
        self.assertEqual(successor["operational_baseline"]["registered_science_slices"], [])

    def test_q011_pressure_pilot_four_slice_policy_promotion_is_historical(self) -> None:
        successor = _load("phase0_curated_candidate_successor_v20_2026-06-02.json")
        self.assertEqual(
            successor["predecessor_sha256"],
            _sha256(REPO_ROOT / successor["predecessor_record"]),
        )
        promotion = successor["policy_promotion"]
        self.assertEqual(
            promotion["repo_policy_sha256"],
            _sha256(REPO_ROOT / promotion["repo_policy_path"]),
        )
        self.assertEqual(
            promotion["orion_policy_sha256"],
            promotion["project_home_policy_sha256"],
        )
        self.assertEqual(
            promotion["orion_promotion_sha256"],
            promotion["project_home_promotion_sha256"],
        )
        policy = _load("storage_policy.json")
        self.assertEqual(
            [
                record["authorization_id"]
                for record in policy["registered_science_slices"]
            ],
            promotion["registered_science_authorization_ids"],
        )
        self.assertEqual(promotion["registered_science_slice_count"], 4)
        current = _load(
            "q011_section54_sixteenth_lustre_publication_rename_compatibility_"
            "transition_2026-06-05.json"
        )["candidate_only_policy_promotion"]
        self.assertEqual(
            current["reviewed_policy_sha256"],
            _sha256(Path(current["reviewed_policy_path"])),
        )
        self.assertEqual(
            current["orion_active_policy_sha256"],
            current["reviewed_policy_sha256"],
        )
        self.assertEqual(
            current["project_home_active_policy_sha256"],
            current["reviewed_policy_sha256"],
        )
        self.assertEqual(
            current["orion_active_promotion_sha256"],
            current["project_home_active_promotion_sha256"],
        )
        current_policy = json.loads(
            Path(current["reviewed_policy_path"]).read_text(encoding="utf-8")
        )
        self.assertEqual(current_policy["registered_science_slices"], [])
        attestation = successor["pre_policy_promotion_operator_attestation"]
        self.assertEqual(attestation["sha256"], _sha256(Path(attestation["path"])))
        self.assertEqual(
            attestation["queue_snapshot_sha256"],
            hashlib.sha256(b"").hexdigest(),
        )
        self.assertEqual(attestation["active_reservation_count"], 0)
        self.assertEqual(attestation["pending_submission_marker"], "absent")
        self.assertEqual(attestation["pending_manual_accounting_marker"], "absent")

    def test_validation_manifest_schema_accepts_reviewable_scaffolds(self) -> None:
        schema = _validation_manifest_schema()
        pending = _minimum_validation_manifest()
        validate_schema(pending, schema)

        qualified = copy.deepcopy(pending)
        qualified["review"] = {
            "reviewer": "Named External Reviewer",
            "disposition": "qualified",
        }
        validate_schema(qualified, schema)

    def test_validation_manifest_schema_rejects_incomplete_evidence(self) -> None:
        schema = _validation_manifest_schema()
        cases = [
            ("empty metrics", ("metrics",), []),
            ("empty metric record", ("metrics",), [{}]),
            ("empty artifacts", ("artifacts",), []),
            ("empty tolerances", ("oracle", "tolerances"), {}),
            ("blank manifest ID", ("manifest_id",), " "),
            ("blank claim ID", ("claim_ids", 0), " "),
            ("blank test ID", ("test_id",), " "),
            ("blank physical mode", ("physical_mode",), " "),
            ("blank executable path", ("executable", "path"), " "),
            ("blank CMake cache path", ("executable", "cmake_cache", "path"), " "),
            ("blank modules path", ("executable", "modules", "path"), " "),
            (
                "blank environment allowlist path",
                ("executable", "environment_allowlist", "path"),
                " ",
            ),
            ("blank oracle kind", ("oracle", "kind"), " "),
            ("blank oracle reference", ("oracle", "reference"), " "),
            ("blank platform", ("resources", "platform"), " "),
            ("blank artifact root", ("resources", "artifact_root"), " "),
            ("blank artifact path", ("artifacts", 0, "path"), " "),
            ("blank reviewer", ("review", "reviewer"), " "),
            (
                "pending reviewer qualified disposition",
                ("review", "disposition"),
                "qualified",
            ),
        ]
        for label, path, replacement in cases:
            with self.subTest(label=label):
                manifest = _minimum_validation_manifest()
                _replace_nested(manifest, path, replacement)
                with self.assertRaises(ValueError):
                    validate_schema(manifest, schema)


if __name__ == "__main__":
    unittest.main()
