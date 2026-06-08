#!/usr/bin/env python3
"""Materialize Q-011 pressure-pilot registered-execution review artifacts.

This source-local tool does not edit the live Frontier storage policy and does
not submit jobs.  After a clean freeze, it binds that exact manifest,
executable and installed environment profile into four separate policy slices
and one fresh selected-case reviewed pre-submit config per invocation.
"""

from __future__ import annotations

import argparse
import copy
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import hashlib
import io
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import types
from typing import Any
import uuid


REPO_ROOT = Path(__file__).resolve().parents[2]
READINESS_ROOT = REPO_ROOT / "tst/publication/readiness"
MATERIALIZER_SOURCE = REPO_ROOT / "tst/publication/q011_section54_pressure_pilot_execution.py"
GENERATOR_SOURCE = REPO_ROOT / "src/pgen/tests/pic_parallel_shock.cpp"
PILOT_PREREGISTRATION = (
    READINESS_ROOT / "q011_section54_pressure_pilot_preregistration_2026-06-01.json"
)
EXECUTION_PREREGISTRATION = (
    READINESS_ROOT
    / "q011_section54_pressure_pilot_registered_execution_retry_successor_v2_2026-06-02.json"
)
JOB_SCRIPT = REPO_ROOT / "tst/publication/frontier_q011_section54_pressure_pilot_job.sh"
INPUT_DECK = (
    REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
ENVIRONMENT_PROFILE_SOURCE = (
    REPO_ROOT / "tst/publication/frontier_control_plane/frontier_pic_environment.sh"
)
OPERATOR_ATTESTATION_SOURCE = (
    REPO_ROOT / "tst/publication/frontier_control_plane/operator_attestation.py"
)
PUBLICATION_SCRIPTS = (
    REPO_ROOT / "tst/publication/analyze_q011_section54_pressure_pilot.py",
    REPO_ROOT / "tst/publication/publish_q011_section54_pressure_pilot_bundle.py",
)
ANALYSIS_SCRIPTS = (
    REPO_ROOT / "tst/publication/analyze_q011_section54_pressure_pilot_case.py",
    REPO_ROOT / "tst/publication/frontier_f1_structured_artifacts.py",
)
AUTHORIZED_PIC_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/PIC")
AUTHORIZED_CANONICAL_PROJECT_HOME_ROOT = Path(
    "/autofs/nccs-svm1_proj/ast207/proj-shared/PIC"
)
AUTHORIZED_PROJECT_HOME_ROOT = AUTHORIZED_CANONICAL_PROJECT_HOME_ROOT
AUTHORIZED_PROJECT_HOME_LEDGER_ROOT = Path("/ccs/proj/ast207/proj-shared/PIC")
EVIDENCE_CLASS = "engineering_calibration_only"
PHYSICAL_MODE = "paper_mhd_pic_vl2_tsc"
RUNTIME_PROFILE = "frontier_minimum_supported"
SELECTED_QOS = "debug"
WALLTIME_SECONDS = 900
SEED_ATHENA_WALLTIME_SECONDS = 600
SEED_TIMEOUT_MAX_VALIDITY_SECONDS = 4 * 60 * 60
SEED_TIMEOUT_FILENAME = "timeout_margin.json"
SEED_TIMEOUT_RATIONALE_FILENAME = "timeout_margin_seed_rationale.json"
QUEUE_SNAPSHOT_FORMAT = "%i|%P|%q|%T|%j|%k"
_SHA256 = re.compile(r"[0-9a-f]{64}")
_GIT_COMMIT = re.compile(r"[0-9a-f]{40}")
_TRUSTED_GIT = "/usr/bin/git"
_TRUSTED_GIT_OPTIONS = (
    "-c",
    "core.fsmonitor=false",
    "-c",
    "core.hooksPath=/dev/null",
)
_UUID = re.compile(
    r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}"
)
_STORAGE_PREFLIGHT_METHOD = "local_create_write_sync_remove_probe"
_REGISTERED_EXECUTION_CONTRACT_SHA256 = (
    "bbd58928a8c5310c9f109a77f53ffe0af245889d7dd8fe84686b82315faa1b5a"
)
_REPAIRED_GENERATOR_SOURCE_SHA256 = (
    "c972ea20eaf5d7e32dbc27879a261b01cffabe46c719dfb1e851b63c17e80f92"
)
_FAILED_V1_GIT_COMMIT = "4972f998589e6f16d1ccff4f6b9facbcb127909d"
_FAILED_V1_CLEAN_CANDIDATE_MANIFEST_SHA256 = (
    "e09b1ab8c7eb017a929a8c063fb53f176b43ca35050b13aa5a0db1f6432b7e97"
)
_FAILED_V1_EXECUTABLE_SHA256 = (
    "c7c3a986fdbd1bb68dd2e8a0aa95df13934acacb33714b38bdb77ef73a1579c5"
)
_CONSUMED_HISTORICAL_V2_SLICES_SHA256 = (
    "a839602e625ea97f301888a18467f6bf6f6ca5b18999fe726df38e68fd1b408a"
)
_HISTORICAL_V2_SCIENCE_FREEZE_SHA256 = (
    "045a80e0a646494c8644e875c7cc733c0f918c56daa981b28760a13401a72f62"
)
_CONSUMED_HISTORICAL_V2_CLEAN_CANDIDATE_MANIFEST_SHA256 = (
    "ef527ed467995bd60fda07b5a3b09b56ea871595ace12fd64a948e246720dbe3"
)
_TIMESTAMP = re.compile(
    r"[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"(?:\.[0-9]{1,6})?Z"
)
_TIMEOUT_MARGIN_KEYS = {
    "athena_walltime_seconds",
    "scheduler_walltime_seconds",
    "environment_profile_sha256",
    "measured_utc",
    "expires_utc",
}
FIXED_OVERRIDES = (
    "mesh/nx1=100",
    "mesh/x1max=1200",
    "mesh/nx2=20",
    "mesh/x2max=240",
    "mesh_refinement/refinement=none",
    "mesh_refinement/num_levels=1",
    "time/tlim=60",
    "time/nlim=4096",
    "time/ndiag=50",
    "problem/ps_enable_curvature_amr=false",
    "problem/ps_feedback_diag_dcycle=50",
    "output1/variable=mhd_w_bcc",
    "output1/id=mhd_w_bcc",
    "output1/dt=15",
    "output2/dt=15",
    "output3/dt=15",
    "output4/dt=15",
    "output5/dt=15",
    "output6/dt=15",
)


@dataclass(frozen=True)
class PressureCase:
    case_id: str
    problem_ps_p0: float
    argv_value: str

    @property
    def campaign(self) -> str:
        return f"q011_section54_pressure_{self.case_id}"

    @property
    def test_id(self) -> str:
        return f"pic_parallel_shock_section54_pressure_{self.case_id}"

    @property
    def authorization_id(self) -> str:
        return f"q011-section54-pressure-{self.case_id.replace('_', '-')}-v2"

    @property
    def launch_contract_path(self) -> Path:
        return READINESS_ROOT / (
            f"frontier_q011_section54_pressure_{self.case_id}_launch_contract.json"
        )


CASES = (
    PressureCase("ps_p0_1p00", 1.0, "1.0"),
    PressureCase("ps_p0_0p05", 0.05, "0.05"),
    PressureCase("ps_p0_0p10", 0.1, "0.10"),
    PressureCase("ps_p0_0p20", 0.2, "0.20"),
)


class ContractError(ValueError):
    """Raised when the registered-execution source tranche drifts."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ContractError(message)


def _expected_git_commit(value: object) -> str:
    _require(
        isinstance(value, str) and _GIT_COMMIT.fullmatch(value) is not None,
        "expected Git commit must be a full lowercase hexadecimal commit",
    )
    return value


def _trusted_git_environment() -> dict[str, str]:
    return {
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "HOME": "/",
        "LANG": "C",
        "LC_ALL": "C",
        "PATH": "/usr/bin:/bin",
    }


def _trusted_git_command(*arguments: str) -> list[str]:
    return [_TRUSTED_GIT, *_TRUSTED_GIT_OPTIONS, *arguments]


def _require_reviewed_git_commit(expected_git_commit: object) -> str:
    """Require this tracked source repository to be clean at one reviewed HEAD."""
    expected = _expected_git_commit(expected_git_commit)
    environment = _trusted_git_environment()
    try:
        repository = Path(
            subprocess.check_output(
                _trusted_git_command(
                    "-C", str(REPO_ROOT), "rev-parse", "--show-toplevel"
                ),
                text=True,
                env=environment,
            ).strip()
        ).resolve()
        _require(
            repository == REPO_ROOT,
            "materializer source repository differs from the expected repository",
        )
        source_path = str(MATERIALIZER_SOURCE.relative_to(repository))
        subprocess.run(
            _trusted_git_command(
                "-C",
                str(repository),
                "ls-files",
                "--error-unmatch",
                "--",
                source_path,
            ),
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            env=environment,
        )
        head_before = subprocess.check_output(
            _trusted_git_command("-C", str(repository), "rev-parse", "HEAD"),
            text=True,
            env=environment,
        ).strip()
        status = subprocess.check_output(
            _trusted_git_command(
                "-C",
                str(repository),
                "status",
                "--porcelain=v1",
                "--untracked-files=no",
                "--ignore-submodules=none",
            ),
            text=True,
            env=environment,
        )
        head_after = subprocess.check_output(
            _trusted_git_command("-C", str(repository), "rev-parse", "HEAD"),
            text=True,
            env=environment,
        ).strip()
    except ContractError:
        raise
    except (OSError, subprocess.CalledProcessError) as error:
        raise ContractError(
            "could not authenticate the reviewed materializer Git commit"
        ) from error
    _require(
        head_before == expected and head_after == expected,
        "materializer source HEAD differs from the expected Git commit",
    )
    _require(
        status == "",
        "materializer source repository must have a clean tracked HEAD",
    )
    return expected


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def _launch_contract_sha256(value: object) -> str:
    return _sha256_bytes(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    )


def _absolute(path: Path, *, label: str) -> Path:
    _require(path.is_absolute(), f"{label} must be absolute")
    return Path(os.path.abspath(path))


def _stable_regular_bytes(
    path: Path,
    *,
    label: str,
    require_read_only: bool = False,
    require_executable: bool = False,
) -> tuple[Path, bytes]:
    path = _absolute(path, label=label)
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError as error:
        raise ContractError(f"{label} is not an openable regular file: {path}") from error
    try:
        before = os.fstat(descriptor)
        _require(stat.S_ISREG(before.st_mode), f"{label} is not a regular file")
        if require_read_only:
            _require(not before.st_mode & 0o222, f"{label} must be read-only")
        if require_executable:
            _require(bool(before.st_mode & 0o111), f"{label} must be executable")
        with os.fdopen(os.dup(descriptor), "rb") as stream:
            payload = stream.read()
        after = os.fstat(descriptor)
        identity = lambda value: (
            value.st_dev,
            value.st_ino,
            value.st_mode,
            value.st_size,
            value.st_mtime_ns,
            value.st_ctime_ns,
        )
        _require(identity(before) == identity(after), f"{label} changed while reading")
        lexical = os.stat(path, follow_symlinks=False)
        _require(
            (lexical.st_dev, lexical.st_ino) == (after.st_dev, after.st_ino),
            f"{label} path changed while reading",
        )
        return path, payload
    finally:
        os.close(descriptor)


def _decode_json(payload: bytes, *, label: str) -> Any:
    def reject_constant(value: str) -> None:
        raise ContractError(f"{label} contains forbidden JSON constant {value}")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            _require(key not in result, f"{label} contains duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        return json.loads(
            payload.decode("utf-8"),
            parse_constant=reject_constant,
            object_pairs_hook=reject_duplicates,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ContractError(f"{label} is not valid UTF-8 JSON") from error


def _read_json(path: Path, *, label: str) -> tuple[Path, bytes, Any]:
    lexical, payload = _stable_regular_bytes(path, label=label)
    return lexical, payload, _decode_json(payload, label=label)


def _utc_datetime(value: object, *, label: str) -> datetime:
    _require(isinstance(value, str), f"{label} must be a canonical UTC timestamp")
    _require(_TIMESTAMP.fullmatch(value) is not None, f"{label} must be canonical UTC")
    try:
        parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    except ValueError as error:
        raise ContractError(f"{label} is invalid") from error
    _require(parsed.tzinfo == timezone.utc, f"{label} must use UTC")
    return parsed


def _timestamp(value: object, *, label: str) -> str:
    _utc_datetime(value, label=label)
    assert isinstance(value, str)
    return value


def _relative(path: Path) -> str:
    return str(path.relative_to(REPO_ROOT))


def expected_launch_contract(case: PressureCase) -> dict[str, object]:
    """Return one exact single-action pressure-pilot launch contract."""
    return {
        "schema_version": 1,
        "executor": "trusted_trampoline_athena_argv_v1",
        "pre_actions": [],
        "actions": [
            {
                "action_id": f"q011-pressure-{case.case_id.replace('_', '-')}",
                "kind": "athena",
                "resources": {
                    "nodes": 1,
                    "tasks": 1,
                    "cpus_per_task": 7,
                    "gpus_per_task": 1,
                    "gpu_bind": "closest",
                },
                "arguments": [
                    {"literal": "-i"},
                    {"snapshot_role": "input-deck"},
                    {"literal": "-d"},
                    {"artifact_directory": "output"},
                    {"literal": f"job/basename={case.case_id}"},
                    *({"literal": override} for override in FIXED_OVERRIDES),
                    {"literal": f"problem/ps_p0={case.argv_value}"},
                ],
                "stdout_artifact": "athena_stdout.txt",
                "stderr_artifact": "athena_stderr.txt",
            }
        ],
        "post_actions": [
            {
                "action_id": "require-pressure-stdout",
                "kind": "artifact_nonempty",
                "artifact": "athena_stdout.txt",
            },
            {
                "action_id": "sha-pressure-stdout",
                "kind": "artifact_sha256",
                "artifact": "athena_stdout.txt",
                "output_artifact": "athena_stdout.sha256",
            },
        ],
    }


def _load_launch_contract(case: PressureCase) -> tuple[bytes, dict[str, object]]:
    _, payload, contract = _read_json(
        case.launch_contract_path, label=f"{case.case_id} launch contract"
    )
    _require(
        contract == expected_launch_contract(case),
        f"{case.case_id} launch contract differs from its exact source contract",
    )
    return payload, contract


def _load_pilot_preregistration() -> tuple[bytes, dict[str, object]]:
    _, payload, preregistration = _read_json(
        PILOT_PREREGISTRATION, label="Q-011 pressure-pilot preregistration"
    )
    _require(isinstance(preregistration, dict), "pressure-pilot preregistration is malformed")
    active_deck = preregistration.get("active_deck_binding")
    pilot_contract = preregistration.get("pilot_contract")
    _require(isinstance(active_deck, dict), "pressure-pilot deck binding is malformed")
    _require(isinstance(pilot_contract, dict), "pressure-pilot contract is malformed")
    _require(
        active_deck
        == {
            "path": _relative(INPUT_DECK),
            "sha256": "0b1cbd62d54027ec81a5f4f5c88d5ee56b86b8cc0cb018c3fbebfb37a11be7b1",
        },
        "pressure-pilot canonical VL2/TSC deck binding drifted",
    )
    _require(
        pilot_contract.get("physical_mode") == PHYSICAL_MODE,
        "pressure-pilot physical mode drifted",
    )
    _require(
        pilot_contract.get("fixed_overrides") == list(FIXED_OVERRIDES),
        "pressure-pilot fixed overrides drifted",
    )
    _require(
        pilot_contract.get("cases")
        == [
            {
                "case_id": case.case_id,
                "problem_ps_p0": case.problem_ps_p0,
                "argv_value": case.argv_value,
            }
            for case in CASES
        ],
        "pressure-pilot case matrix drifted",
    )
    return payload, preregistration


def source_bindings() -> dict[str, object]:
    """Measure the exact local source tranche without binding post-freeze files."""
    pilot_payload, _ = _load_pilot_preregistration()
    _, generator_payload = _stable_regular_bytes(
        GENERATOR_SOURCE, label="repaired parallel-shock generator source"
    )
    _require(
        _sha256_bytes(generator_payload) == _REPAIRED_GENERATOR_SOURCE_SHA256,
        "repaired parallel-shock generator source drifted",
    )
    _, job_payload = _stable_regular_bytes(JOB_SCRIPT, label="pressure-pilot job template")
    _, deck_payload = _stable_regular_bytes(INPUT_DECK, label="canonical VL2/TSC input deck")
    _, environment_payload = _stable_regular_bytes(
        ENVIRONMENT_PROFILE_SOURCE, label="reviewed environment-profile source"
    )
    _, materializer_payload = _stable_regular_bytes(
        MATERIALIZER_SOURCE, label="registered-execution materializer source"
    )
    analyses = []
    for path in ANALYSIS_SCRIPTS:
        _, payload = _stable_regular_bytes(path, label=f"analysis source {_relative(path)}")
        analyses.append({"path": _relative(path), "sha256": _sha256_bytes(payload)})
    publications = []
    for path in PUBLICATION_SCRIPTS:
        _, payload = _stable_regular_bytes(path, label=f"publication source {_relative(path)}")
        publications.append({"path": _relative(path), "sha256": _sha256_bytes(payload)})
    contracts = []
    for case in CASES:
        payload, contract = _load_launch_contract(case)
        contracts.append(
            {
                "case_id": case.case_id,
                "path": _relative(case.launch_contract_path),
                "file_sha256": _sha256_bytes(payload),
                "launch_contract_sha256": _launch_contract_sha256(contract),
            }
        )
    return {
        "pressure_pilot_preregistration": {
            "path": _relative(PILOT_PREREGISTRATION),
            "sha256": _sha256_bytes(pilot_payload),
        },
        "generator_source": {
            "path": _relative(GENERATOR_SOURCE),
            "sha256": _sha256_bytes(generator_payload),
        },
        "job_script": {
            "path": _relative(JOB_SCRIPT),
            "sha256": _sha256_bytes(job_payload),
        },
        "input_deck": {
            "path": _relative(INPUT_DECK),
            "sha256": _sha256_bytes(deck_payload),
        },
        "environment_profile_source": {
            "path": _relative(ENVIRONMENT_PROFILE_SOURCE),
            "sha256": _sha256_bytes(environment_payload),
        },
        "materializer": {
            "path": _relative(MATERIALIZER_SOURCE),
            "sha256": _sha256_bytes(materializer_payload),
        },
        "analysis_scripts": analyses,
        "publication_scripts": publications,
        "launch_contracts": contracts,
    }


def _historical_v2_preregistration() -> dict[str, object]:
    """Load the immutable launch preregistration without reauthorizing it."""
    _, _, preregistration = _read_json(
        EXECUTION_PREREGISTRATION,
        label="Q-011 pressure-pilot registered-execution preregistration",
    )
    _require(isinstance(preregistration, dict), "execution preregistration is malformed")
    contract = dict(preregistration)
    contract.pop("source_bindings", None)
    _require(
        _sha256_bytes(_json_bytes(contract)) == _REGISTERED_EXECUTION_CONTRACT_SHA256,
        "registered-execution preregistration contract drifted",
    )
    return preregistration


def historical_v2_source_tranche_status() -> dict[str, object]:
    """Describe whether the consumed v2 launch tranche still matches checkout."""
    preregistration = _historical_v2_preregistration()
    matches = preregistration.get("source_bindings") == source_bindings()
    return {
        "record_type": "q011_section54_pressure_pilot_historical_v2_source_tranche_status",
        "schema_version": 1,
        "state": "historical_consumed_slice_non_authorizing",
        "source_bindings_match_current_checkout": matches,
        "launch_reauthorization_effect": "none",
        "consumed_slice_reauthorization_allowed": False,
        "required_postrun_boundary": (
            "reviewed_immutable_postrun_aggregate_source_authorization_successor"
        ),
    }


def validate_source_tranche() -> dict[str, object]:
    """Reject every fresh launch attempt through the consumed historical v2 tranche."""
    _historical_v2_preregistration()
    raise ContractError(
        "historical v2 registered-execution tranche is consumed and nonauthorizing; "
        "fresh launches require a separately reviewed successor"
    )


def _generator_bytes_from_source_archive(payload: bytes) -> bytes:
    """Read the exact repaired generator member from a frozen Git source archive."""
    path = _relative(GENERATOR_SOURCE)
    try:
        with tarfile.open(fileobj=io.BytesIO(payload), mode="r:") as archive:
            members = [member for member in archive.getmembers() if member.name == path]
            _require(
                len(members) == 1 and members[0].isfile(),
                "clean-candidate source archive lacks one regular generator source",
            )
            stream = archive.extractfile(members[0])
            _require(stream is not None, "clean-candidate generator source cannot be read")
            return stream.read()
    except tarfile.TarError as error:
        raise ContractError("clean-candidate source archive is not a readable tar file") from error


def _bound_candidate_artifacts(
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, str]:
    clean_candidate_manifest, manifest_payload = _stable_regular_bytes(
        clean_candidate_manifest,
        label="clean-candidate manifest",
        require_read_only=True,
    )
    executable, executable_payload = _stable_regular_bytes(
        executable,
        label="clean-candidate executable",
        require_read_only=True,
        require_executable=True,
    )
    environment_profile, environment_payload = _stable_regular_bytes(
        environment_profile,
        label="installed environment profile",
        require_read_only=True,
    )
    _require(
        clean_candidate_manifest.name == "clean_candidate_manifest.json",
        "clean-candidate manifest filename is not canonical",
    )
    _require(
        executable == clean_candidate_manifest.parent / "athena",
        "clean-candidate executable is not adjacent to the supplied manifest",
    )
    source_archive, source_archive_payload = _stable_regular_bytes(
        clean_candidate_manifest.parent / "source.tar",
        label="clean-candidate source archive",
        require_read_only=True,
    )
    manifest = _decode_json(manifest_payload, label="clean-candidate manifest")
    _require(isinstance(manifest, dict), "clean-candidate manifest is malformed")
    _require(
        type(manifest.get("schema_version")) is int
        and manifest["schema_version"] == 4,
        "clean-candidate manifest schema is unsupported",
    )
    source = manifest.get("source")
    build = manifest.get("build")
    _require(isinstance(source, dict), "clean-candidate source binding is malformed")
    _require(isinstance(build, dict), "clean-candidate build binding is malformed")
    clean_candidate_manifest_sha256 = _sha256_bytes(manifest_payload)
    _require(
        clean_candidate_manifest_sha256 != _FAILED_V1_CLEAN_CANDIDATE_MANIFEST_SHA256,
        "clean-candidate manifest is the failed v1 carrier",
    )
    git_commit = source.get("git_commit")
    _require(
        isinstance(git_commit, str) and _GIT_COMMIT.fullmatch(git_commit) is not None,
        "clean-candidate Git commit is malformed",
    )
    _require(git_commit != _FAILED_V1_GIT_COMMIT, "clean-candidate Git commit is the failed v1 carrier")
    source_archive_sha256 = _sha256_bytes(source_archive_payload)
    _require(
        source_archive == clean_candidate_manifest.parent / "source.tar"
        and source.get("archive_path") == str(source_archive)
        and source.get("archive_sha256") == source_archive_sha256
        and build.get("source_archive_sha256") == source_archive_sha256,
        "clean-candidate source archive differs from its manifest binding",
    )
    expected_generator_sha256 = source_bindings()["generator_source"]["sha256"]
    _require(
        _sha256_bytes(_generator_bytes_from_source_archive(source_archive_payload))
        == expected_generator_sha256,
        "clean-candidate source archive does not contain the repaired generator source",
    )
    executable_sha256 = _sha256_bytes(executable_payload)
    _require(
        executable_sha256 != _FAILED_V1_EXECUTABLE_SHA256,
        "clean-candidate executable is the failed v1 carrier",
    )
    _require(
        build.get("executable_path") == str(executable)
        and build.get("executable_sha256") == executable_sha256,
        "clean-candidate executable differs from its manifest binding",
    )
    expected_environment_sha256 = source_bindings()["environment_profile_source"]["sha256"]
    environment_sha256 = _sha256_bytes(environment_payload)
    _require(
        environment_sha256 == expected_environment_sha256,
        "installed environment profile differs from preregistered source bytes",
    )
    return {
        "clean_candidate_manifest_path": str(clean_candidate_manifest),
        "clean_candidate_manifest_sha256": clean_candidate_manifest_sha256,
        "source_archive_path": str(source_archive),
        "source_archive_sha256": source_archive_sha256,
        "generator_source_sha256": expected_generator_sha256,
        "executable_path": str(executable),
        "executable_sha256": executable_sha256,
        "environment_profile_path": str(environment_profile),
        "environment_profile_sha256": environment_sha256,
        "git_commit": git_commit,
    }


def _bound_final_artifacts(
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, str]:
    """Bind candidate artifacts only through the historical v2 launch guard."""
    validate_source_tranche()
    return _bound_candidate_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )


def materialize_registered_science_slices(
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> list[dict[str, object]]:
    """Bind four separate one-attempt policy slices after the clean freeze."""
    final = _bound_final_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    sources = source_bindings()
    contract_digests = {
        record["case_id"]: record["launch_contract_sha256"]
        for record in sources["launch_contracts"]
    }
    return [
        {
            "authorization_id": case.authorization_id,
            "status": "authorized",
            "campaign": case.campaign,
            "test_id": case.test_id,
            "evidence_class": EVIDENCE_CLASS,
            "physical_mode": PHYSICAL_MODE,
            "runtime_profile": RUNTIME_PROFILE,
            "selected_qos": SELECTED_QOS,
            "registered_short_nonproduction": True,
            "maximum_nodes": 1,
            "maximum_walltime_seconds": WALLTIME_SECONDS,
            "maximum_attempts": 1,
            "job_script_sha256": sources["job_script"]["sha256"],
            "input_deck_sha256": sources["input_deck"]["sha256"],
            "environment_profile_sha256": final["environment_profile_sha256"],
            "analysis_script_sha256": [
                record["sha256"] for record in sources["analysis_scripts"]
            ],
            "executable_sha256": final["executable_sha256"],
            "launch_contract_sha256": contract_digests[case.case_id],
            "clean_candidate_manifest_sha256": final[
                "clean_candidate_manifest_sha256"
            ],
        }
        for case in CASES
    ]


def materialize_policy_fragment(
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, object]:
    """Return a reviewable additive fragment; never edit storage_policy.json."""
    validate_source_tranche()
    _, preregistration_payload = _stable_regular_bytes(
        EXECUTION_PREREGISTRATION,
        label="Q-011 pressure-pilot registered-execution preregistration",
    )
    final = _bound_final_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    return {
        "record_type": "q011_section54_pressure_pilot_post_freeze_policy_fragment",
        "schema_version": 1,
        "source_preregistration": {
            "path": _relative(EXECUTION_PREREGISTRATION),
            "sha256": _sha256_bytes(preregistration_payload),
        },
        "final_clean_binding": final,
        "registered_science_slices": materialize_registered_science_slices(
            clean_candidate_manifest=clean_candidate_manifest,
            executable=executable,
            environment_profile=environment_profile,
        ),
        "integration_policy": (
            "review_then_add_registered_science_slices_via_separate_storage_policy_"
            "successor_and_installed_control_plane_promotion"
        ),
    }


def _control_plane_version(value: object) -> str:
    _require(
        isinstance(value, str) and _SHA256.fullmatch(value) is not None,
        "control-plane version must be one lowercase SHA-256 digest",
    )
    return value


def _project_home_ledger_root() -> Path:
    """Retain the production ledger chain's historical lexical destination."""
    if AUTHORIZED_PROJECT_HOME_ROOT == AUTHORIZED_CANONICAL_PROJECT_HOME_ROOT:
        return AUTHORIZED_PROJECT_HOME_LEDGER_ROOT
    return AUTHORIZED_PROJECT_HOME_ROOT


def _storage_preflight_binding(path: Path) -> dict[str, object]:
    """Read one immutable canonical storage-probe policy fragment."""
    _, payload = _stable_regular_bytes(
        path,
        label="storage-preflight binding",
        require_read_only=True,
    )
    value = _decode_json(payload, label="storage-preflight binding")
    _require(isinstance(value, dict), "storage-preflight binding must be one object")
    _require(
        set(value)
        == {
            "last_preflight_utc",
            "orion_simulation_root_preflight",
            "project_home_preflight",
            "storage_preflight_evidence",
        },
        "storage-preflight binding schema is malformed",
    )
    completed_utc = _timestamp(
        value["last_preflight_utc"],
        label="storage-preflight binding completion time",
    )
    expected_records = {
        "orion_simulation_root_preflight": {
            "method": _STORAGE_PREFLIGHT_METHOD,
            "path": str(AUTHORIZED_PIC_ROOT),
            "status": "passed",
        },
        "project_home_preflight": {
            "method": _STORAGE_PREFLIGHT_METHOD,
            "path": str(AUTHORIZED_PROJECT_HOME_ROOT),
            "status": "passed",
        },
    }
    for key, expected in expected_records.items():
        _require(
            value.get(key) == expected,
            f"storage-preflight binding {key} is malformed",
        )
    evidence = value.get("storage_preflight_evidence")
    _require(
        isinstance(evidence, dict)
        and set(evidence) == {"orion_path", "probe_id", "project_home_path", "sha256"},
        "storage-preflight evidence binding is malformed",
    )
    probe_id = evidence.get("probe_id")
    _require(
        isinstance(probe_id, str) and _UUID.fullmatch(probe_id) is not None,
        "storage-preflight evidence probe ID is malformed",
    )
    _require(
        evidence
        == {
            "orion_path": str(
                AUTHORIZED_PIC_ROOT
                / "policy"
                / "storage_preflight_evidence"
                / f"{probe_id}.json"
            ),
            "probe_id": probe_id,
            "project_home_path": str(
                AUTHORIZED_PROJECT_HOME_ROOT
                / "policy"
                / "storage_preflight_evidence"
                / f"{probe_id}.json"
            ),
            "sha256": evidence.get("sha256"),
        }
        and isinstance(evidence.get("sha256"), str)
        and _SHA256.fullmatch(evidence["sha256"]) is not None,
        "storage-preflight evidence binding is malformed",
    )
    return {
        "last_preflight_utc": completed_utc,
        **expected_records,
        "storage_preflight_evidence": dict(evidence),
    }


def _canonicalize_project_home_successor_storage(storage: dict[str, object]) -> None:
    """Move inherited policy mirrors to the strict physical Project Home root."""
    mirror_root = storage.get("project_home_mirror_root")
    if mirror_root is not None:
        _require(
            mirror_root
            in {
                str(AUTHORIZED_PROJECT_HOME_ROOT),
                str(AUTHORIZED_CANONICAL_PROJECT_HOME_ROOT),
                str(AUTHORIZED_PROJECT_HOME_LEDGER_ROOT),
            },
            "baseline Project Home policy root is not an expected reviewed spelling",
        )
    storage["project_home_mirror_root"] = str(AUTHORIZED_PROJECT_HOME_ROOT)
    authorizations = storage.get("manual_accounting_authorizations")
    if authorizations is None:
        return
    _require(
        isinstance(authorizations, list),
        "baseline manual-accounting authorizations are malformed",
    )
    for authorization in authorizations:
        _require(
            isinstance(authorization, dict)
            and isinstance(authorization.get("project_home_path"), str),
            "baseline manual-accounting authorization mirror path is malformed",
        )
        lexical = Path(os.path.abspath(authorization["project_home_path"]))
        relative: Path | None = None
        for root in [
            AUTHORIZED_PROJECT_HOME_LEDGER_ROOT,
            AUTHORIZED_CANONICAL_PROJECT_HOME_ROOT,
            AUTHORIZED_PROJECT_HOME_ROOT,
        ]:
            try:
                relative = lexical.relative_to(root)
                break
            except ValueError:
                continue
        _require(
            relative is not None
            and len(relative.parts) == 3
            and relative.parts[:2]
            == ("policy", "manual_accounting_authorizations")
            and relative.name.endswith(".json"),
            "baseline manual-accounting authorization mirror path is not reviewed",
        )
        authorization["project_home_path"] = str(
            AUTHORIZED_PROJECT_HOME_ROOT / relative
        )


def _advance_control_plane_fields(
    successor: dict[str, object],
    *,
    control_plane_version: str,
    storage_preflight_binding: Path,
    require_fresh_preflight: bool = False,
    require_new_control_plane: bool = False,
    require_same_control_plane: bool = False,
    allow_equal_preflight: bool = False,
) -> dict[str, object]:
    storage = successor.get("olcf_side_storage")
    _require(isinstance(storage, dict), "baseline OLCF-side storage policy is malformed")
    version = _control_plane_version(control_plane_version)
    _require(
        not (require_new_control_plane and require_same_control_plane),
        "control-plane successor requirement is contradictory",
    )
    if require_new_control_plane or require_same_control_plane:
        installed = storage.get("installed_control_plane_version")
        staged = storage.get("staged_control_plane_candidate_version")
        _require(
            installed == staged
            and isinstance(installed, str)
            and _SHA256.fullmatch(installed) is not None,
            "baseline installed/staged control-plane versions are malformed",
        )
    if require_new_control_plane:
        _require(version != installed, "successor control-plane version must be new")
    if require_same_control_plane:
        _require(
            version == installed,
            "candidate-only control-plane version must match baseline "
            "installed/staged version",
        )
    binding = _storage_preflight_binding(storage_preflight_binding)
    _canonicalize_project_home_successor_storage(storage)
    checked = str(binding["last_preflight_utc"])
    if require_fresh_preflight:
        previous = storage.get("last_preflight_utc")
        checked_time = _utc_datetime(
            checked, label="OLCF-side storage preflight time"
        )
        previous_time = _utc_datetime(
            previous, label="baseline OLCF-side storage preflight time"
        )
        _require(
            checked_time >= previous_time if allow_equal_preflight
            else checked_time > previous_time,
            "successor OLCF-side storage preflight time must be equal to or newer "
            "than baseline"
            if allow_equal_preflight
            else "successor OLCF-side storage preflight time must be newer than baseline",
        )
    storage["installed_control_plane_version"] = version
    storage["staged_control_plane_candidate_version"] = version
    storage.update(binding)
    return successor


def materialize_baseline_policy_successor(
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Bind a full launch-prohibited policy copy to one installed successor."""
    _, _, baseline = _read_json(baseline_policy, label="baseline storage policy")
    _require(isinstance(baseline, dict), "baseline storage policy is malformed")
    slices = baseline.get("registered_science_slices")
    _require(
        slices == [],
        "baseline storage policy must have an empty registered-science allowlist",
    )
    successor = copy.deepcopy(baseline)
    return _advance_control_plane_fields(
        successor,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )


def materialize_exact_reviewed_preflight_predecessor_policy_successor(
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Advance only the controller and preflight of an empty-allowlist policy."""
    _, _, baseline = _read_json(baseline_policy, label="baseline storage policy")
    _require(isinstance(baseline, dict), "baseline storage policy is malformed")
    _require(
        baseline.get("registered_science_slices") == [],
        "baseline storage policy must have an empty registered-science allowlist",
    )
    baseline_storage = baseline.get("olcf_side_storage")
    _require(
        isinstance(baseline_storage, dict),
        "baseline OLCF-side storage policy is malformed",
    )
    baseline_binding = baseline_storage.get("storage_preflight_evidence")
    _require(
        isinstance(baseline_binding, dict),
        "baseline storage-preflight evidence binding is malformed",
    )
    successor = _advance_control_plane_fields(
        copy.deepcopy(baseline),
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
        require_fresh_preflight=True,
        require_new_control_plane=True,
    )
    successor_storage = successor["olcf_side_storage"]
    _require(
        isinstance(successor_storage, dict),
        "successor OLCF-side storage policy is malformed",
    )
    successor_binding = successor_storage.get("storage_preflight_evidence")
    _require(
        isinstance(successor_binding, dict)
        and successor_binding != baseline_binding
        and successor_binding.get("probe_id") != baseline_binding.get("probe_id")
        and successor_binding.get("sha256") != baseline_binding.get("sha256"),
        "successor storage-preflight evidence binding must be different",
    )
    mutable_storage_fields = {
        "installed_control_plane_version",
        "staged_control_plane_candidate_version",
        "last_preflight_utc",
        "storage_preflight_evidence",
    }
    baseline_comparable = copy.deepcopy(baseline)
    successor_comparable = copy.deepcopy(successor)
    for value in [baseline_comparable, successor_comparable]:
        storage = value["olcf_side_storage"]
        for field in mutable_storage_fields:
            storage[field] = None
    _require(
        _json_bytes(successor_comparable) == _json_bytes(baseline_comparable),
        "exact reviewed predecessor successor changed unrelated policy fields",
    )
    return successor


def materialize_retire_consumed_slices_baseline_policy_successor(
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Retire only the exact consumed historical Q011 pressure-v2 allowlist."""
    _, _, baseline = _read_json(baseline_policy, label="historical storage policy")
    _require(isinstance(baseline, dict), "historical storage policy is malformed")
    slices = baseline.get("registered_science_slices")
    _require(
        isinstance(slices, list)
        and _sha256_bytes(_json_bytes(slices))
        == _CONSUMED_HISTORICAL_V2_SLICES_SHA256,
        "historical storage policy must contain exactly the consumed Q011 "
        "pressure-v2 slices",
    )
    science_freeze = baseline.get("science_submission_freeze")
    _require(
        isinstance(science_freeze, dict)
        and _sha256_bytes(_json_bytes(science_freeze))
        == _HISTORICAL_V2_SCIENCE_FREEZE_SHA256,
        "historical storage policy must retain the exact consumed Q011 frozen candidate",
    )
    successor = copy.deepcopy(baseline)
    successor["registered_science_slices"] = []
    successor["science_submission_freeze"] = {
        "status": "pending_clean_candidate_freeze"
    }
    return _advance_control_plane_fields(
        successor,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
        require_fresh_preflight=True,
        require_new_control_plane=True,
    )


def materialize_candidate_only_policy_successor(
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    expected_git_commit: str,
) -> dict[str, object]:
    """Authorize one exact clean candidate while keeping all launches prohibited."""
    expected_git_commit = _expected_git_commit(expected_git_commit)
    _, _, baseline = _read_json(baseline_policy, label="baseline storage policy")
    _require(isinstance(baseline, dict), "baseline storage policy is malformed")
    _require(
        baseline.get("registered_science_slices") == [],
        "baseline storage policy must have an empty registered-science allowlist",
    )
    baseline_freeze = baseline.get("science_submission_freeze")
    pending_freeze = baseline_freeze == {"status": "pending_clean_candidate_freeze"}
    authorized_freeze = (
        isinstance(baseline_freeze, dict)
        and set(baseline_freeze)
        == {
            "status",
            "manifest_path",
            "manifest_sha256",
            "build_profile_control_plane_version",
        }
        and baseline_freeze.get("status") == "authorized"
        and isinstance(baseline_freeze.get("manifest_path"), str)
        and _SHA256.fullmatch(str(baseline_freeze.get("manifest_sha256"))) is not None
        and _SHA256.fullmatch(
            str(baseline_freeze.get("build_profile_control_plane_version"))
        )
        is not None
    )
    _require(
        pending_freeze or authorized_freeze,
        "baseline storage policy must have a pending or exact authorized "
        "clean-candidate freeze",
    )
    successor = _advance_control_plane_fields(
        copy.deepcopy(baseline),
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
        require_fresh_preflight=True,
        require_same_control_plane=True,
        allow_equal_preflight=True,
    )
    final = _bound_candidate_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    _require(
        final["clean_candidate_manifest_sha256"]
        != _CONSUMED_HISTORICAL_V2_CLEAN_CANDIDATE_MANIFEST_SHA256,
        "candidate-only clean-candidate manifest must be fresh",
    )
    _require(
        final["git_commit"] == expected_git_commit,
        "candidate-only clean-candidate Git commit differs from the expected Git commit",
    )
    if authorized_freeze:
        assert isinstance(baseline_freeze, dict)
        _require(
            final["clean_candidate_manifest_sha256"]
            != baseline_freeze["manifest_sha256"],
            "candidate-only clean-candidate manifest must differ from the currently "
            "authorized freeze",
        )
    successor["science_submission_freeze"] = {
        "status": "authorized",
        "manifest_path": final["clean_candidate_manifest_path"],
        "manifest_sha256": final["clean_candidate_manifest_sha256"],
        "build_profile_control_plane_version": _control_plane_version(
            control_plane_version
        ),
    }
    return successor


def materialize_pilot_policy_successor(
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, object]:
    """Bind a complete reviewed policy successor to the fresh freeze and four pilots."""
    successor = materialize_baseline_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )
    final = _bound_final_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    successor["science_submission_freeze"] = {
        "status": "authorized",
        "manifest_path": final["clean_candidate_manifest_path"],
        "manifest_sha256": final["clean_candidate_manifest_sha256"],
        "build_profile_control_plane_version": _control_plane_version(
            control_plane_version
        ),
    }
    successor["registered_science_slices"] = materialize_registered_science_slices(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    return successor


def _validate_timeout_margin(
    path: Path,
    *,
    environment_profile_sha256: str,
    now: datetime | None = None,
) -> str:
    path, _, timeout = _read_json(path, label="timeout-margin artifact")
    _require(
        isinstance(timeout, dict) and set(timeout) == _TIMEOUT_MARGIN_KEYS,
        "timeout-margin artifact schema drifted",
    )
    _require(
        type(timeout["scheduler_walltime_seconds"]) is int
        and timeout["scheduler_walltime_seconds"] == WALLTIME_SECONDS,
        "timeout-margin scheduler walltime must be exactly 900 seconds",
    )
    _require(
        type(timeout["athena_walltime_seconds"]) is int
        and 0 < timeout["athena_walltime_seconds"] < WALLTIME_SECONDS,
        "timeout-margin Athena walltime must be a bounded exact integer",
    )
    _require(
        timeout["environment_profile_sha256"] == environment_profile_sha256,
        "timeout-margin artifact belongs to another environment profile",
    )
    measured = _utc_datetime(timeout["measured_utc"], label="timeout-margin measured_utc")
    expires = _utc_datetime(timeout["expires_utc"], label="timeout-margin expires_utc")
    _require(measured < expires, "timeout-margin validity interval is empty")
    checked = datetime.now(timezone.utc) if now is None else now
    _require(
        measured <= checked < expires,
        "timeout-margin artifact is stale or not yet valid",
    )
    return str(path)


def _selected_case(case_id: object) -> PressureCase:
    _require(
        isinstance(case_id, str) and bool(case_id),
        "one pressure-pilot case ID is required",
    )
    matches = [case for case in CASES if case.case_id == case_id]
    _require(len(matches) == 1, f"unknown pressure-pilot case ID: {case_id}")
    return matches[0]


def _submission_id(value: object, *, case_id: str) -> str:
    _require(isinstance(value, str), f"{case_id} submission ID is malformed")
    try:
        parsed = uuid.UUID(value)
    except ValueError as error:
        raise ContractError(f"{case_id} submission ID is malformed") from error
    _require(str(parsed) == value, f"{case_id} submission ID is not canonical")
    return value


def _validate_queue_snapshot(path: Path) -> tuple[Path, dict[str, str]]:
    path, payload = _stable_regular_bytes(
        path, label="six-field queue snapshot", require_read_only=True
    )
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise ContractError("six-field queue snapshot is not UTF-8") from error
    _require(not payload or text.endswith("\n"), "six-field queue snapshot lacks final LF")
    for index, line in enumerate(text.splitlines()):
        _require(
            len(line.split("|")) == 6,
            f"six-field queue snapshot line {index + 1} is malformed",
        )
    return path, {
        "path": str(path),
        "sha256": _sha256_bytes(payload),
        "format": QUEUE_SNAPSHOT_FORMAT,
    }


def _validate_pre_manifest_attestation(
    path: Path, *, case: PressureCase, control_plane_version: str
) -> tuple[Path, dict[str, str]]:
    specification = importlib.util.spec_from_file_location(
        "_q011_operator_attestation", OPERATOR_ATTESTATION_SOURCE
    )
    _require(
        specification is not None and specification.loader is not None,
        "operator-attestation validator cannot be loaded",
    )
    module = importlib.util.module_from_spec(specification)
    sys.modules[specification.name] = module
    try:
        specification.loader.exec_module(module)
        binding = module.validate_sealed_operator_attestation(
            path,
            authorization_id=case.authorization_id,
            phase="pre_manifest",
            control_plane_version=_control_plane_version(control_plane_version),
            authorized_pic_root=AUTHORIZED_PIC_ROOT,
            authorized_project_home_root=_project_home_ledger_root(),
        )
    except ValueError as error:
        raise ContractError(f"pre-manifest attestation is invalid: {error}") from error
    finally:
        sys.modules.pop(specification.name, None)
    return Path(binding["path"]), binding


def _mirrored_ledger_records() -> list[dict[str, object]]:
    local = AUTHORIZED_PIC_ROOT / "ledger/node_hours.jsonl"
    mirror = _project_home_ledger_root() / "ledger/node_hours.jsonl"
    receipts = AUTHORIZED_PIC_ROOT / "ledger/mirror_receipts.jsonl"
    control_plane = REPO_ROOT / "tst/publication/frontier_control_plane"
    specification = importlib.util.spec_from_file_location(
        "_q011_frontier_ledger", control_plane / "ledger.py"
    )
    _require(
        specification is not None and specification.loader is not None,
        "Frontier ledger validator cannot be loaded",
    )
    module = importlib.util.module_from_spec(specification)
    sys.path.insert(0, str(control_plane))
    try:
        specification.loader.exec_module(module)
        records = module.validate_mirrored_state(local, receipts, mirror)
        module.require_explicit_genesis(records)
        return records
    except ValueError as error:
        raise ContractError(f"mirrored node-hours ledger is invalid: {error}") from error
    finally:
        sys.path.pop(0)


def _load_raw_case_analyzer() -> object:
    _require(
        os.environ.get("PIC_F1_ANALYSIS_HELPER_FD") is None,
        "raw-case analyzer rejects inherited helper overrides",
    )
    path = ANALYSIS_SCRIPTS[0]
    path, payload = _stable_regular_bytes(path, label="raw-case analyzer source")
    module = types.ModuleType("_q011_raw_case_analyzer")
    module.__file__ = str(path)
    code = compile(payload, str(path), "exec", dont_inherit=True)
    exec(code, module.__dict__)
    return module


def _verified_prior_case_closure(
    case: PressureCase,
    submission_id: str,
    descriptor_sha256: str,
    final: dict[str, str],
) -> dict[str, str]:
    identifier = _submission_id(submission_id, case_id=case.case_id)
    _require(
        _SHA256.fullmatch(descriptor_sha256) is not None,
        f"{case.case_id} descriptor SHA-256 is malformed",
    )
    artifact_dir = AUTHORIZED_PIC_ROOT / "runs" / case.campaign / identifier
    matches = [
        record
        for record in _mirrored_ledger_records()
        if record.get("event_type") == "reconciliation"
        and record.get("submission_scope") == "registered_science"
        and record.get("submission_id") == identifier
        and record.get("campaign") == case.campaign
        and record.get("registered_science_authorization_id") == case.authorization_id
        and record.get("artifact_dir") == str(artifact_dir)
        and record.get("clean_candidate_manifest_sha256")
        == final["clean_candidate_manifest_sha256"]
        and record.get("git_commit") == final["git_commit"]
        and record.get("executable_sha256") == final["executable_sha256"]
        and record.get("reconciled") is True
        and record.get("state") == "COMPLETED"
    ]
    _require(
        len(matches) == 1,
        f"{case.case_id} lacks one completed registered-science reconciliation",
    )
    module = _load_raw_case_analyzer()
    with module.StructuredArtifactTree(artifact_dir) as tree:
        module.verify_published_case_descriptor(tree, case.case_id, descriptor_sha256)
    return {
        "case_id": case.case_id,
        "submission_id": identifier,
        "artifact_dir": str(artifact_dir),
        "descriptor_path": str(artifact_dir / "analysis/analysis.json"),
        "descriptor_sha256": descriptor_sha256,
        "reconciliation_event_sha256": str(matches[0]["event_sha256"]),
    }


def _validate_prior_case_closures(
    case: PressureCase, values: tuple[str, ...], final: dict[str, str]
) -> list[dict[str, str]]:
    selected_index = CASES.index(case)
    expected = CASES[:selected_index]
    _require(
        len(values) == len(expected),
        f"{case.case_id} requires exactly {len(expected)} prior-case closure bindings",
    )
    closures = []
    for prior, value in zip(expected, values):
        fields = value.split("=")
        _require(
            len(fields) == 3 and fields[0] == prior.case_id,
            f"{case.case_id} prior-case closures are not in exact preregistered order",
        )
        closures.append(_verified_prior_case_closure(prior, fields[1], fields[2], final))
    return closures


def materialize_reviewed_pre_submit_config(
    *,
    case_id: str,
    submission_id: str,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    control_plane_version: str,
    pre_manifest_attestation: Path,
    timeout_margin_artifact: Path,
    queue_snapshot: Path,
    site_policy_checked_utc: str,
    prior_case_closures: tuple[str, ...] = (),
) -> dict[str, object]:
    """Build one exact config after a fresh selected-case pre-manifest capture."""
    case = _selected_case(case_id)
    identifier = _submission_id(submission_id, case_id=case.case_id)
    final = _bound_final_artifacts(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    _validate_pre_manifest_attestation(
        pre_manifest_attestation,
        case=case,
        control_plane_version=control_plane_version,
    )
    timeout_margin_artifact = Path(
        _validate_timeout_margin(
            timeout_margin_artifact,
            environment_profile_sha256=final["environment_profile_sha256"],
        )
    )
    queue_snapshot, _ = _validate_queue_snapshot(queue_snapshot)
    closures = _validate_prior_case_closures(case, prior_case_closures, final)
    checked = _timestamp(site_policy_checked_utc, label="site-policy checked time")
    return {
        "pic_root": str(AUTHORIZED_PIC_ROOT),
        "campaign": case.campaign,
        "test_id": case.test_id,
        "submission_scope": "registered_science",
        "registered_science_authorization_id": case.authorization_id,
        "pre_manifest_attestation": str(pre_manifest_attestation),
        "submission_id": identifier,
        "git_commit": final["git_commit"],
        "evidence_class": EVIDENCE_CLASS,
        "physical_mode": PHYSICAL_MODE,
        "selected_qos": SELECTED_QOS,
        "qos_selection_reason": "debug_available",
        "site_policy_checked_utc": checked,
        "registered_short_nonproduction": True,
        "artifact_dir": str(AUTHORIZED_PIC_ROOT / "runs" / case.campaign / identifier),
        "job_script_executable_env": "PIC_EXECUTABLE",
        "job_script": str(JOB_SCRIPT),
        "executable": final["executable_path"],
        "input_deck": str(INPUT_DECK),
        "environment_profile": final["environment_profile_path"],
        "timeout_margin_artifact": str(timeout_margin_artifact),
        "analysis_scripts": [str(path) for path in ANALYSIS_SCRIPTS],
        "queue_snapshot": str(queue_snapshot),
        "prior_case_closures": closures,
        "clean_candidate_manifest": final["clean_candidate_manifest_path"],
        "launch_contract": _load_launch_contract(case)[1],
    }


def _write_new_file(path: Path, payload: bytes) -> None:
    path = _absolute(path, label="output file")
    _require(path.parent.is_dir(), f"output parent is not a directory: {path.parent}")
    directory_flags = (
        os.O_RDONLY | os.O_DIRECTORY | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        directory = os.open(path.parent, directory_flags)
    except OSError as error:
        raise ContractError(f"output parent is not a stable directory: {path.parent}") from error

    def require_same_parent() -> None:
        lexical = os.open(path.parent, directory_flags)
        try:
            expected = os.fstat(directory)
            actual = os.fstat(lexical)
            _require(
                (actual.st_dev, actual.st_ino) == (expected.st_dev, expected.st_ino),
                f"output parent changed while writing: {path.parent}",
            )
        finally:
            os.close(lexical)

    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(path.name, flags, 0o444, dir_fd=directory)
        except OSError as error:
            raise ContractError(f"refusing to overwrite output file: {path}") from error
        try:
            offset = 0
            while offset < len(payload):
                offset += os.write(descriptor, payload[offset:])
            os.fsync(descriptor)
            os.fchmod(descriptor, 0o444)
            os.fsync(descriptor)
        finally:
            os.close(descriptor)
        require_same_parent()
        os.fsync(directory)
        require_same_parent()
    finally:
        os.close(directory)


def write_seed_timeout_margin(
    output_root: Path,
    *,
    case_id: str,
    environment_profile: Path,
    materialized_utc: str,
    expires_utc: str,
) -> dict[str, object]:
    """Write an explicit non-measurement timeout seed for one engineering pilot."""
    output_root = _absolute(output_root, label="seed-timeout output root")
    _require(output_root.parent.is_dir(), "seed-timeout output parent does not exist")
    case = _selected_case(case_id)
    environment_profile, environment_payload = _stable_regular_bytes(
        environment_profile,
        label="installed environment profile",
        require_read_only=True,
    )
    materialized = _utc_datetime(materialized_utc, label="seed-timeout materialized_utc")
    expires = _utc_datetime(expires_utc, label="seed-timeout expires_utc")
    _require(materialized < expires, "seed-timeout validity interval is empty")
    _require(
        expires - materialized
        <= timedelta(seconds=SEED_TIMEOUT_MAX_VALIDITY_SECONDS),
        "seed-timeout validity interval exceeds four hours",
    )
    environment_sha256 = _sha256_bytes(environment_payload)
    timeout = {
        "athena_walltime_seconds": SEED_ATHENA_WALLTIME_SECONDS,
        "scheduler_walltime_seconds": WALLTIME_SECONDS,
        "environment_profile_sha256": environment_sha256,
        "measured_utc": materialized_utc,
        "expires_utc": expires_utc,
    }
    timeout_payload = _json_bytes(timeout)
    rationale = {
        "record_type": "q011_section54_pressure_pilot_seed_timeout_margin_rationale",
        "schema_version": 1,
        "classification": "engineering_seed_margin_bootstrap_only_not_measurement",
        "qualification_effect": (
            "none_not_empirical_runtime_evidence_not_production_sizing_not_"
            "sun_bai_qualification"
        ),
        "case_id": case.case_id,
        "authorization_id": case.authorization_id,
        "environment_profile": {
            "path": str(environment_profile),
            "sha256": environment_sha256,
        },
        "timeout_margin_artifact": {
            "path": SEED_TIMEOUT_FILENAME,
            "sha256": _sha256_bytes(timeout_payload),
        },
        "controller_field_semantics": {
            "measured_utc": (
                "seed artifact materialization timestamp only; this bootstrap "
                "record is explicitly not an empirical runtime measurement"
            )
        },
        "selected_margin": {
            "scheduler_walltime_seconds": WALLTIME_SECONDS,
            "athena_walltime_seconds": SEED_ATHENA_WALLTIME_SECONDS,
            "shutdown_margin_seconds": (
                WALLTIME_SECONDS - SEED_ATHENA_WALLTIME_SECONDS
            ),
        },
        "invalidated_by": [
            "code_or_toolchain_change",
            "environment_profile_change",
            "mesh_or_ppc_change",
            "rank_layout_change",
            "restart_or_checkpoint_cadence_change",
            "output_mode_change",
            "case_change",
            "expiry",
        ],
        "required_followup": [
            "retain_scheduler_elapsed_seconds_for_this_case",
            "retain_observed_cycle_and_checkpoint_timing_if_emitted",
            "replace_seed_rationale_with_case_specific_empirical_margin_before_reuse",
            "do_not_use_seed_rationale_for_production_sizing_or_qualification",
        ],
    }
    try:
        output_root.mkdir(mode=0o700)
    except OSError as error:
        raise ContractError("seed-timeout output root already exists") from error
    try:
        _write_new_file(output_root / SEED_TIMEOUT_FILENAME, timeout_payload)
        _write_new_file(
            output_root / SEED_TIMEOUT_RATIONALE_FILENAME, _json_bytes(rationale)
        )
        os.chmod(output_root, 0o555)
        return rationale
    except BaseException:
        shutil.rmtree(output_root)
        raise


def write_policy_fragment(
    output: Path,
    *,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, object]:
    """Write one new read-only additive policy fragment."""
    fragment = materialize_policy_fragment(
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    _write_new_file(output, _json_bytes(fragment))
    return fragment


def write_baseline_policy_successor(
    output: Path,
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Write one complete launch-prohibited reviewed policy successor."""
    successor = materialize_baseline_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )
    _write_new_file(output, _json_bytes(successor))
    return successor


def write_retire_consumed_slices_baseline_policy_successor(
    output: Path,
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
) -> dict[str, object]:
    """Write one complete successor retiring the exact consumed Q011 slices."""
    successor = materialize_retire_consumed_slices_baseline_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )
    _write_new_file(output, _json_bytes(successor))
    return successor


def write_exact_reviewed_preflight_predecessor_policy_successor(
    output: Path,
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
    expected_git_commit: str,
) -> dict[str, object]:
    """Write one exact empty-allowlist predecessor migration successor."""
    expected_git_commit = _require_reviewed_git_commit(expected_git_commit)
    successor = materialize_exact_reviewed_preflight_predecessor_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
    )
    payload = _json_bytes(successor)
    _require_reviewed_git_commit(expected_git_commit)
    _write_new_file(output, payload)
    return successor


def write_candidate_only_policy_successor(
    output: Path,
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    expected_git_commit: str,
) -> dict[str, object]:
    """Write one exact-candidate successor without adding launch slices."""
    expected_git_commit = _require_reviewed_git_commit(expected_git_commit)
    successor = materialize_candidate_only_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
        expected_git_commit=expected_git_commit,
    )
    payload = _json_bytes(successor)
    _require_reviewed_git_commit(expected_git_commit)
    _write_new_file(output, payload)
    return successor


def write_pilot_policy_successor(
    output: Path,
    *,
    baseline_policy: Path,
    control_plane_version: str,
    storage_preflight_binding: Path,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
) -> dict[str, object]:
    """Write one complete reviewed post-freeze four-slice policy successor."""
    successor = materialize_pilot_policy_successor(
        baseline_policy=baseline_policy,
        control_plane_version=control_plane_version,
        storage_preflight_binding=storage_preflight_binding,
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
    )
    _write_new_file(output, _json_bytes(successor))
    return successor


def write_reviewed_pre_submit_config(
    output_root: Path,
    *,
    case_id: str,
    submission_id: str,
    clean_candidate_manifest: Path,
    executable: Path,
    environment_profile: Path,
    control_plane_version: str,
    pre_manifest_attestation: Path,
    timeout_margin_artifact: Path,
    queue_snapshot: Path,
    site_policy_checked_utc: str,
    prior_case_closures: tuple[str, ...] = (),
) -> dict[str, object]:
    """Create one new read-only selected-case config handoff directory."""
    output_root = _absolute(output_root, label="reviewed-config output root")
    _require(output_root.parent.is_dir(), "reviewed-config output parent does not exist")
    case = _selected_case(case_id)
    config = materialize_reviewed_pre_submit_config(
        case_id=case.case_id,
        submission_id=submission_id,
        clean_candidate_manifest=clean_candidate_manifest,
        executable=executable,
        environment_profile=environment_profile,
        control_plane_version=control_plane_version,
        pre_manifest_attestation=pre_manifest_attestation,
        timeout_margin_artifact=timeout_margin_artifact,
        queue_snapshot=queue_snapshot,
        site_policy_checked_utc=site_policy_checked_utc,
        prior_case_closures=prior_case_closures,
    )
    _, pre_manifest_binding = _validate_pre_manifest_attestation(
        pre_manifest_attestation,
        case=case,
        control_plane_version=control_plane_version,
    )
    _, queue_binding = _validate_queue_snapshot(queue_snapshot)
    timeout_path, timeout_payload = _stable_regular_bytes(
        timeout_margin_artifact, label="timeout-margin artifact"
    )
    try:
        output_root.mkdir(mode=0o700)
    except OSError as error:
        raise ContractError("reviewed-config output root already exists") from error
    try:
        filename = "pre_submit_config.json"
        payload = _json_bytes(config)
        _write_new_file(output_root / filename, payload)
        manifest = {
            "record_type": "q011_section54_pressure_pilot_reviewed_pre_submit_config",
            "schema_version": 1,
            "source_preregistration": _relative(EXECUTION_PREREGISTRATION),
            "case_id": case.case_id,
            "authorization_id": case.authorization_id,
            "submission_id": submission_id,
            "captured_inputs": {
                "pre_manifest_attestation": pre_manifest_binding,
                "six_field_queue_snapshot": queue_binding,
                "timeout_margin_artifact": {
                    "path": str(timeout_path),
                    "sha256": _sha256_bytes(timeout_payload),
                },
                "verified_prior_case_closures": config["prior_case_closures"],
            },
            "config": {
                "path": filename,
                "sha256": _sha256_bytes(payload),
            },
            "required_next_steps": [
                "create_pre_submit_manifest_for_this_selected_case_only",
                "capture_fresh_pre_submit_wrapper_attestation_for_this_selected_case",
                "invoke_registered_submission_wrapper_for_this_selected_case",
                "wait_for_terminal_scheduler_state_for_this_selected_case",
                "reconcile_terminal_job_through_installed_control_plane",
                "retain_elapsed_seconds_and_checkpoint_timing_evidence",
                "verify_absent_pending_marker_zero_active_reservations_and_empty_same_account_pic_queue",
                "run_snapshotted_raw_case_analyzer_and_retain_descriptor_sha256",
                "materialize_no_other_case_until_this_terminal_boundary_is_complete",
            ],
        }
        _write_new_file(output_root / "materialization_manifest.json", _json_bytes(manifest))
        os.chmod(output_root, 0o555)
        return manifest
    except BaseException:
        shutil.rmtree(output_root)
        raise


def _add_final_binding_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--clean-candidate-manifest", required=True, type=Path)
    parser.add_argument("--executable", required=True, type=Path)
    parser.add_argument("--environment-profile", required=True, type=Path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    policy = subparsers.add_parser("policy-slices")
    _add_final_binding_arguments(policy)
    policy.add_argument("--output", required=True, type=Path)
    baseline_policy = subparsers.add_parser("baseline-policy-successor")
    baseline_policy.add_argument("--baseline-policy", required=True, type=Path)
    baseline_policy.add_argument("--control-plane-version", required=True)
    baseline_policy.add_argument("--storage-preflight-binding", required=True, type=Path)
    baseline_policy.add_argument("--output", required=True, type=Path)
    exact_predecessor_policy = subparsers.add_parser(
        "exact-reviewed-preflight-predecessor-policy-successor"
    )
    exact_predecessor_policy.add_argument("--baseline-policy", required=True, type=Path)
    exact_predecessor_policy.add_argument("--control-plane-version", required=True)
    exact_predecessor_policy.add_argument(
        "--storage-preflight-binding", required=True, type=Path
    )
    exact_predecessor_policy.add_argument("--expected-git-commit", required=True)
    exact_predecessor_policy.add_argument("--output", required=True, type=Path)
    retirement_policy = subparsers.add_parser(
        "retire-consumed-slices-baseline-policy-successor"
    )
    retirement_policy.add_argument("--baseline-policy", required=True, type=Path)
    retirement_policy.add_argument("--control-plane-version", required=True)
    retirement_policy.add_argument(
        "--storage-preflight-binding", required=True, type=Path
    )
    retirement_policy.add_argument("--output", required=True, type=Path)
    candidate_policy = subparsers.add_parser("candidate-only-policy-successor")
    candidate_policy.add_argument("--baseline-policy", required=True, type=Path)
    candidate_policy.add_argument("--control-plane-version", required=True)
    candidate_policy.add_argument("--storage-preflight-binding", required=True, type=Path)
    _add_final_binding_arguments(candidate_policy)
    candidate_policy.add_argument("--expected-git-commit", required=True)
    candidate_policy.add_argument("--output", required=True, type=Path)
    pilot_policy = subparsers.add_parser("pilot-policy-successor")
    pilot_policy.add_argument("--baseline-policy", required=True, type=Path)
    pilot_policy.add_argument("--control-plane-version", required=True)
    pilot_policy.add_argument("--storage-preflight-binding", required=True, type=Path)
    _add_final_binding_arguments(pilot_policy)
    pilot_policy.add_argument("--output", required=True, type=Path)
    seed = subparsers.add_parser("seed-timeout-margin")
    seed.add_argument("--case-id", required=True, choices=[case.case_id for case in CASES])
    seed.add_argument("--environment-profile", required=True, type=Path)
    seed.add_argument("--materialized-utc", required=True)
    seed.add_argument("--expires-utc", required=True)
    seed.add_argument("--output-root", required=True, type=Path)
    configs = subparsers.add_parser("pre-submit-config")
    _add_final_binding_arguments(configs)
    configs.add_argument("--case-id", required=True, choices=[case.case_id for case in CASES])
    configs.add_argument("--submission-id", required=True)
    configs.add_argument("--control-plane-version", required=True)
    configs.add_argument("--output-root", required=True, type=Path)
    configs.add_argument("--pre-manifest-attestation", required=True, type=Path)
    configs.add_argument("--timeout-margin-artifact", required=True, type=Path)
    configs.add_argument("--queue-snapshot", required=True, type=Path)
    configs.add_argument("--site-policy-checked-utc", required=True)
    configs.add_argument("--prior-case-closure", action="append", default=[])
    return parser


def main() -> None:
    arguments = build_parser().parse_args()
    if arguments.command == "seed-timeout-margin":
        rationale = write_seed_timeout_margin(
            arguments.output_root,
            case_id=arguments.case_id,
            environment_profile=arguments.environment_profile,
            materialized_utc=arguments.materialized_utc,
            expires_utc=arguments.expires_utc,
        )
        print(json.dumps(rationale, sort_keys=True, separators=(",", ":")))
        return
    if arguments.command == "baseline-policy-successor":
        successor = write_baseline_policy_successor(
            arguments.output,
            baseline_policy=arguments.baseline_policy,
            control_plane_version=arguments.control_plane_version,
            storage_preflight_binding=arguments.storage_preflight_binding,
        )
        print(json.dumps(successor, sort_keys=True, separators=(",", ":")))
        return
    if arguments.command == "retire-consumed-slices-baseline-policy-successor":
        successor = write_retire_consumed_slices_baseline_policy_successor(
            arguments.output,
            baseline_policy=arguments.baseline_policy,
            control_plane_version=arguments.control_plane_version,
            storage_preflight_binding=arguments.storage_preflight_binding,
        )
        print(json.dumps(successor, sort_keys=True, separators=(",", ":")))
        return
    if arguments.command == "exact-reviewed-preflight-predecessor-policy-successor":
        successor = write_exact_reviewed_preflight_predecessor_policy_successor(
            arguments.output,
            baseline_policy=arguments.baseline_policy,
            control_plane_version=arguments.control_plane_version,
            storage_preflight_binding=arguments.storage_preflight_binding,
            expected_git_commit=arguments.expected_git_commit,
        )
        print(json.dumps(successor, sort_keys=True, separators=(",", ":")))
        return
    common = {
        "clean_candidate_manifest": arguments.clean_candidate_manifest,
        "executable": arguments.executable,
        "environment_profile": arguments.environment_profile,
    }
    if arguments.command == "policy-slices":
        fragment = write_policy_fragment(arguments.output, **common)
        print(json.dumps(fragment, sort_keys=True, separators=(",", ":")))
        return
    if arguments.command == "candidate-only-policy-successor":
        successor = write_candidate_only_policy_successor(
            arguments.output,
            baseline_policy=arguments.baseline_policy,
            control_plane_version=arguments.control_plane_version,
            storage_preflight_binding=arguments.storage_preflight_binding,
            expected_git_commit=arguments.expected_git_commit,
            **common,
        )
        print(json.dumps(successor, sort_keys=True, separators=(",", ":")))
        return
    if arguments.command == "pilot-policy-successor":
        successor = write_pilot_policy_successor(
            arguments.output,
            baseline_policy=arguments.baseline_policy,
            control_plane_version=arguments.control_plane_version,
            storage_preflight_binding=arguments.storage_preflight_binding,
            **common,
        )
        print(json.dumps(successor, sort_keys=True, separators=(",", ":")))
        return
    manifest = write_reviewed_pre_submit_config(
        arguments.output_root,
        case_id=arguments.case_id,
        submission_id=arguments.submission_id,
        control_plane_version=arguments.control_plane_version,
        pre_manifest_attestation=arguments.pre_manifest_attestation,
        timeout_margin_artifact=arguments.timeout_margin_artifact,
        queue_snapshot=arguments.queue_snapshot,
        site_policy_checked_utc=arguments.site_policy_checked_utc,
        prior_case_closures=tuple(arguments.prior_case_closure),
        **common,
    )
    print(json.dumps(manifest, sort_keys=True, separators=(",", ":")))


if __name__ == "__main__":
    main()
