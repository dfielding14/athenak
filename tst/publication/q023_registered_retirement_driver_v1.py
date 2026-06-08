#!/opt/cray/pe/python/3.11.7/bin/python3
"""Retire one completed registered Q023 allowlist into a new controller."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import sys
import uuid


PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
PROJECT_HOME_ROOT = Path("/autofs/nccs-svm1_proj/ast207/proj-shared/PIC")
PYTHON = Path("/opt/cray/pe/python/3.11.7/bin/python3")
MATRIX_PATH = (
    PIC_ROOT
    / "analysis/q023_paper_bell_linear_joverc_registered_successor_v1"
    / "q023_registered_matrix_qualification.json"
)
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
COMMIT_PATTERN = re.compile(r"[0-9a-f]{40}")


class RetirementError(RuntimeError):
    """Reject drifted source, evidence, controller, or active policy state."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RetirementError(message)


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
        env={**os.environ, "HOME": "/", "LANG": "C", "LC_ALL": "C"},
    )
    if result.returncode != 0:
        raise RetirementError(
            f"command failed ({result.returncode}): {' '.join(command)}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result.stdout.strip()


def _stable_read_only_json(path: Path) -> dict[str, object]:
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        before = os.fstat(descriptor)
        _require(
            stat.S_ISREG(before.st_mode)
            and not before.st_mode & 0o222
            and before.st_nlink == 1,
            f"expected one read-only regular file: {path}",
        )
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        payload = b"".join(chunks)
        after = os.fstat(descriptor)
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
            and len(payload) == after.st_size,
            f"file changed while reading: {path}",
        )
    finally:
        os.close(descriptor)
    try:
        value = json.loads(payload)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RetirementError(f"invalid JSON object: {path}") from error
    _require(type(value) is dict, f"expected one JSON object: {path}")
    return value


def _atomic_json(path: Path, value: object) -> None:
    payload = (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.parent / f".{path.name}.tmp-{uuid.uuid4()}"
    descriptor = os.open(
        temporary,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o444,
    )
    try:
        view = memoryview(payload)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError(f"short write: {temporary}")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)
    try:
        os.link(temporary, path, follow_symlinks=False)
    finally:
        os.unlink(temporary)
    directory = os.open(path.parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(directory)
    finally:
        os.close(directory)


def _active_paths(root: Path) -> tuple[Path, Path]:
    return root / "policy/storage_policy.json", root / "policy/active_promotion.json"


def validate_retired_active_pair(
    expected_policy: dict[str, object],
    *,
    control_plane_version: str,
    pic_root: Path = PIC_ROOT,
    project_home_root: Path = PROJECT_HOME_ROOT,
) -> tuple[str, str]:
    orion_policy, orion_promotion = _active_paths(pic_root)
    project_policy, project_promotion = _active_paths(project_home_root)
    _require(
        orion_policy.read_bytes() == project_policy.read_bytes()
        and orion_promotion.read_bytes() == project_promotion.read_bytes(),
        "active policy or promotion mirrors differ",
    )
    observed = _stable_read_only_json(orion_policy)
    _require(
        observed == expected_policy,
        "active retired Q023 policy differs from reviewed successor",
    )
    storage = observed.get("olcf_side_storage")
    _require(
        observed.get("registered_science_slices") == []
        and type(storage) is dict
        and storage.get("installed_control_plane_version")
        == control_plane_version,
        "active retired Q023 generation has wrong lifecycle state",
    )
    return _sha256(orion_policy), _sha256(orion_promotion)


class Driver:
    def __init__(self, arguments: argparse.Namespace) -> None:
        self.source_root = arguments.source_root.resolve(strict=True)
        self.source_commit = arguments.source_commit
        self.control_plane_version = arguments.control_plane_version
        self.storage_preflight_binding = arguments.storage_preflight_binding.resolve(
            strict=True
        )
        self.expected_predecessor_policy_sha256 = (
            arguments.expected_predecessor_policy_sha256
        )
        self.expected_predecessor_promotion_sha256 = (
            arguments.expected_predecessor_promotion_sha256
        )
        self.matrix_path = arguments.q023_registered_matrix.resolve(strict=True)
        self.control_plane_root = PIC_ROOT / "control_plane" / self.control_plane_version
        self.run_control_plane = self.control_plane_root / "run_control_plane.py"
        self._validate_identities()
        sys.path.insert(0, str(self.source_root))
        from tst.publication import (  # noqa: PLC0415
            q023_registered_execution_linear_qualification_successor_v1
            as qualification,
        )
        from tst.publication import (  # noqa: PLC0415
            q023_registered_launch_policy_preparation_successor_v1
            as preparation,
        )

        self.qualification = qualification
        self.preparation = preparation

    def _validate_identities(self) -> None:
        for label, value, pattern in (
            ("source commit", self.source_commit, COMMIT_PATTERN),
            ("control-plane version", self.control_plane_version, SHA256_PATTERN),
            (
                "predecessor policy SHA-256",
                self.expected_predecessor_policy_sha256,
                SHA256_PATTERN,
            ),
            (
                "predecessor promotion SHA-256",
                self.expected_predecessor_promotion_sha256,
                SHA256_PATTERN,
            ),
        ):
            _require(pattern.fullmatch(value) is not None, f"malformed Q023 {label}")
        _require(
            self.matrix_path == MATRIX_PATH,
            "Q023 retirement matrix path is not canonical",
        )
        _require(
            _run(["/usr/bin/git", "-C", self.source_root, "rev-parse", "HEAD"])
            == self.source_commit
            and not _run(
                [
                    "/usr/bin/git",
                    "-C",
                    self.source_root,
                    "status",
                    "--porcelain=v1",
                    "--untracked-files=all",
                ]
            ),
            "Q023 retirement source is not the exact clean selected commit",
        )
        _require(
            self.run_control_plane.is_file(),
            "Q023 retirement installed control-plane runner is absent",
        )

    def run(self) -> int:
        matrix = self.qualification.validate_downstream_q019_prerequisite(
            _stable_read_only_json(self.matrix_path),
            q043_artifact_root=PIC_ROOT,
        )
        matrix_sha256 = _sha256(self.matrix_path)
        active_policy_path, active_promotion_path = _active_paths(PIC_ROOT)
        active_policy = _stable_read_only_json(active_policy_path)
        active_policy_sha256 = _sha256(active_policy_path)
        active_promotion_sha256 = _sha256(active_promotion_path)
        _require(
            active_policy_sha256 == self.expected_predecessor_policy_sha256
            and active_promotion_sha256
            == self.expected_predecessor_promotion_sha256,
            "active Q023 predecessor policy or promotion digest drifted",
        )
        successor = self.preparation.materialize_q023_retired_policy(
            active_policy=active_policy,
            registered_matrix=matrix,
            successor_control_plane_version=self.control_plane_version,
            storage_preflight_binding=self.storage_preflight_binding,
            q043_artifact_root=PIC_ROOT,
        )
        reviewed_path = (
            PIC_ROOT
            / "policy"
            / (
                "reviewed_q023_retired_successor_"
                f"{self.source_commit[:8]}_{self.control_plane_version[:8]}_"
                f"{matrix_sha256[:8]}.json"
            )
        )
        if reviewed_path.exists():
            _require(
                _stable_read_only_json(reviewed_path) == successor,
                "existing reviewed Q023 retirement policy differs",
            )
        else:
            _atomic_json(reviewed_path, successor)
        _run(
            [
                PYTHON,
                "-I",
                "-B",
                self.run_control_plane,
                "promote_active_policy.py",
                "--reviewed-policy",
                reviewed_path,
                "--retire-completed-q023-registered-slices",
                "--q023-registered-matrix",
                self.matrix_path,
                "--q023-registered-matrix-sha256",
                matrix_sha256,
                "--expected-active-policy-sha256",
                active_policy_sha256,
                "--expected-active-promotion-sha256",
                active_promotion_sha256,
            ]
        )
        successor_policy_sha256, successor_promotion_sha256 = (
            validate_retired_active_pair(
                successor, control_plane_version=self.control_plane_version
            )
        )
        print(
            json.dumps(
                {
                    "matrix_path": str(self.matrix_path),
                    "matrix_sha256": matrix_sha256,
                    "reviewed_policy_path": str(reviewed_path),
                    "reviewed_policy_sha256": _sha256(reviewed_path),
                    "successor_policy_sha256": successor_policy_sha256,
                    "successor_promotion_sha256": successor_promotion_sha256,
                    "installed_control_plane_version": self.control_plane_version,
                    "registered_science_slices": 0,
                },
                sort_keys=True,
            )
        )
        return 0


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--source-commit", required=True)
    parser.add_argument("--control-plane-version", required=True)
    parser.add_argument("--storage-preflight-binding", required=True, type=Path)
    parser.add_argument("--expected-predecessor-policy-sha256", required=True)
    parser.add_argument("--expected-predecessor-promotion-sha256", required=True)
    parser.add_argument(
        "--q023-registered-matrix", type=Path, default=MATRIX_PATH
    )
    return parser


def main() -> int:
    return Driver(_parser().parse_args()).run()


if __name__ == "__main__":
    raise SystemExit(main())
