#!/opt/cray/pe/python/3.11.7/bin/python3
"""Promote one reviewed Q019 carrier-calibration retirement policy."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import time

from tst.publication import (
    q019_q023_carrier_calibration_retirement_v1 as retirement,
)


PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
PROJECT_HOME_ROOT = Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC")
PYTHON = Path("/opt/cray/pe/python/3.11.7/bin/python3")


class RetirementDriverError(RuntimeError):
    """Reject drifted bindings or non-empty queues at policy mutation."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RetirementDriverError(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _run(argv: list[str | Path], *, cwd: Path | None = None) -> str:
    command = [str(item) for item in argv]
    result = subprocess.run(
        command,
        cwd=str(cwd) if cwd is not None else None,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        env={
            **os.environ,
            "HOME": "/",
            "LANG": "C",
            "LC_ALL": "C",
            "SLURM_CLUSTERS": "frontier",
        },
    )
    if result.returncode != 0:
        raise RetirementDriverError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout.strip()


def _read_bound_json(
    path: Path, *, expected_sha256: str, label: str
) -> dict[str, object]:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and before.st_nlink == 1
            and not before.st_mode & 0o222,
            f"{label} is not one read-only regular file",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
        current = path.stat(follow_symlinks=False)
        _require(
            (
                before.st_dev,
                before.st_ino,
                before.st_mode,
                before.st_nlink,
                before.st_size,
                before.st_mtime_ns,
                before.st_ctime_ns,
            )
            == (
                after.st_dev,
                after.st_ino,
                after.st_mode,
                after.st_nlink,
                after.st_size,
                after.st_mtime_ns,
                after.st_ctime_ns,
            )
            and (after.st_dev, after.st_ino) == (current.st_dev, current.st_ino)
            and len(payload) == after.st_size
            and hashlib.sha256(payload).hexdigest() == expected_sha256,
            f"{label} changed or its digest drifted",
        )
    finally:
        os.close(descriptor)
    value = json.loads(payload)
    _require(type(value) is dict, f"{label} must be one JSON object")
    return value


def _active_pair() -> tuple[dict[str, object], str, str]:
    policy = PIC_ROOT / "policy/storage_policy.json"
    promotion = PIC_ROOT / "policy/active_promotion.json"
    project_policy = PROJECT_HOME_ROOT / "policy/storage_policy.json"
    project_promotion = PROJECT_HOME_ROOT / "policy/active_promotion.json"
    _require(
        policy.read_bytes() == project_policy.read_bytes()
        and promotion.read_bytes() == project_promotion.read_bytes(),
        "active policy or promotion mirrors differ",
    )
    return (
        json.loads(policy.read_text(encoding="utf-8")),
        _sha256(policy),
        _sha256(promotion),
    )


def _queue() -> str:
    return _run(
        [
            "/usr/bin/squeue",
            "--clusters=frontier",
            "-u",
            os.environ.get("USER", "dfielding"),
            "-h",
            "-o",
            "%i|%a|%P|%q|%T|%j|%k",
        ]
    )


def _wait_for_empty_queue() -> None:
    while True:
        queue = _queue()
        if not queue:
            return
        print(f"waiting for empty same-user queue:\n{queue}", flush=True)
        time.sleep(30)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--control-plane-version", required=True)
    parser.add_argument("--final-bindings", required=True, type=Path)
    parser.add_argument("--final-bindings-sha256", required=True)
    parser.add_argument("--engineering-qualification", required=True, type=Path)
    parser.add_argument("--engineering-qualification-sha256", required=True)
    parser.add_argument("--storage-preflight-binding", required=True, type=Path)
    parser.add_argument("--reviewed-retired-policy", required=True, type=Path)
    parser.add_argument("--reviewed-retired-policy-sha256", required=True)
    parser.add_argument("--expected-active-policy-sha256", required=True)
    parser.add_argument("--expected-active-promotion-sha256", required=True)
    return parser


def main() -> int:
    arguments = _parser().parse_args()
    source_root = arguments.source_root.resolve(strict=True)
    _require(
        _run(["/usr/bin/git", "-C", source_root, "rev-parse", "HEAD"])
        == arguments.source_commit
        and not _run(
            [
                "/usr/bin/git",
                "-C",
                source_root,
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
            ]
        ),
        "carrier retirement source is not the exact clean selected commit",
    )
    final = _read_bound_json(
        arguments.final_bindings.resolve(strict=True),
        expected_sha256=arguments.final_bindings_sha256,
        label="carrier final bindings",
    )
    qualified = _read_bound_json(
        arguments.engineering_qualification.resolve(strict=True),
        expected_sha256=arguments.engineering_qualification_sha256,
        label="carrier engineering qualification",
    )
    reviewed = _read_bound_json(
        arguments.reviewed_retired_policy.resolve(strict=True),
        expected_sha256=arguments.reviewed_retired_policy_sha256,
        label="reviewed carrier retired policy",
    )
    active, active_sha256, promotion_sha256 = _active_pair()
    _require(
        active_sha256 == arguments.expected_active_policy_sha256
        and promotion_sha256 == arguments.expected_active_promotion_sha256,
        "active carrier predecessor digest drifted",
    )
    expected = retirement.materialize_retired_policy(
        active_policy=active,
        engineering_qualification=qualified,
        final_bindings=final,
        successor_control_plane_version=arguments.control_plane_version,
        storage_preflight_binding=arguments.storage_preflight_binding.resolve(
            strict=True
        ),
    )
    _require(
        expected == reviewed,
        "reviewed retired carrier policy differs from exact materialization",
    )
    _wait_for_empty_queue()
    active, active_sha256, promotion_sha256 = _active_pair()
    _require(
        active_sha256 == arguments.expected_active_policy_sha256
        and promotion_sha256 == arguments.expected_active_promotion_sha256,
        "active carrier predecessor changed while waiting for retirement",
    )
    control_plane_root = PIC_ROOT / "control_plane" / arguments.control_plane_version
    _run(
        [
            PYTHON,
            "-I",
            "-B",
            control_plane_root / "run_control_plane.py",
            "promote_active_policy.py",
            "--reviewed-policy",
            arguments.reviewed_retired_policy,
            "--retire-completed-q019-carrier-calibration-slices",
            "--q019-carrier-qualification",
            arguments.engineering_qualification,
            "--q019-carrier-qualification-sha256",
            arguments.engineering_qualification_sha256,
            "--expected-active-policy-sha256",
            active_sha256,
            "--expected-active-promotion-sha256",
            promotion_sha256,
        ]
    )
    observed, observed_sha256, observed_promotion_sha256 = _active_pair()
    _require(
        observed == reviewed
        and observed_sha256 == arguments.reviewed_retired_policy_sha256
        and observed["registered_science_slices"] == [],
        "live carrier retirement produced the wrong active successor",
    )
    print(
        json.dumps(
            {
                "retired_slice_count": 6,
                "successor_policy_sha256": observed_sha256,
                "successor_promotion_sha256": observed_promotion_sha256,
                "registered_science_slices": 0,
                "production_resource_freeze_authorized": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
