#!/usr/bin/env python3
"""Replay R14/R15 after the known zero-forcing normalization guard.

The corrected bounded launcher intentionally leaves failed attempts in place.
This adapter recognizes only the exact turbulence-driver zero-field failure,
then starts a fresh attempt from the nearest earlier segment whose analysis
authenticated a complete restart.  Failed attempts are retained and bound in
the recovery manifest; simulation inputs, executable, ranks, and physics
overrides are unchanged.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve()
BOUNDED_LAUNCHER = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_fast_corrected_bounded.py"
)
BOUNDED_LAUNCHER_SHA256 = (
    "90104d65f3b207634d300d8d277d61991c8e72d85b4568e68a37179541ac03ea"
)
RECOVERY_POLICY = (
    "replay_last_authenticated_checkpoint_after_zero_forcing_guard_v1"
)
FAILURE_REASON = "turbulence_driver_zero_forcing_field_guard"
FAILURE_SIGNATURE = "cannot inject non-zero dedt with a zero forcing field"
FAILED_STATES = frozenset(
    {
        "BOOT_FAIL",
        "CANCELLED",
        "DEADLINE",
        "FAILED",
        "NODE_FAIL",
        "OUT_OF_MEMORY",
        "PREEMPTED",
        "REVOKED",
        "TIMEOUT",
    }
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_bounded_launcher():
    if sha256_file(BOUNDED_LAUNCHER) != BOUNDED_LAUNCHER_SHA256:
        raise RuntimeError(f"bounded launcher SHA-256 differs: {BOUNDED_LAUNCHER}")
    name = "cgl_lf_stage_i_fast_corrected_bounded_recovery_base"
    spec = importlib.util.spec_from_file_location(name, BOUNDED_LAUNCHER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import bounded launcher: {BOUNDED_LAUNCHER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


bounded = _load_bounded_launcher()


class BoundedRecoveryError(RuntimeError):
    """A bounded zero-field recovery contract failed."""


def normalized_state(value: str) -> str:
    stripped = value.strip()
    return stripped.split("+", 1)[0].split()[0].upper() if stripped else "UNKNOWN"


def require_inactive_case(case_id: str) -> None:
    completed = subprocess.run(
        ["/usr/bin/squeue", "-h", "-u", subprocess.check_output(
            ["/usr/bin/id", "-un"], text=True
        ).strip(), "-o", "%j"],
        check=True,
        text=True,
        capture_output=True,
    )
    prefix = f"cglc_{case_id}_"
    active = [line.strip() for line in completed.stdout.splitlines()
              if line.strip().startswith(prefix)]
    if active:
        raise BoundedRecoveryError(
            f"{case_id} already has an active bounded job: {', '.join(active)}"
        )


def slurm_log_path(manifest: dict[str, object]) -> Path:
    pattern = str(manifest.get("slurm_log") or "")
    job_id = str(manifest.get("job_id") or "")
    case_id = str(manifest.get("case_id") or "")
    sequence = manifest.get("sequence")
    if (
        not pattern
        or not job_id
        or not isinstance(sequence, int)
        or isinstance(sequence, bool)
    ):
        raise BoundedRecoveryError("failed segment lacks Slurm log metadata")
    job_name = f"cglc_{case_id}_s{sequence:03d}"
    return Path(pattern.replace("%x", job_name).replace("%j", job_id))


def failed_attempt_record(segment: Path) -> dict[str, object]:
    manifest = bounded.validate_bounded_segment(segment)
    job_id = str(manifest.get("job_id") or "")
    if not re.fullmatch(r"[1-9][0-9]*", job_id):
        raise BoundedRecoveryError(f"failed segment has no valid job ID: {segment}")
    state = normalized_state(bounded.corrected.fast.job_state(job_id))
    if state not in FAILED_STATES:
        raise BoundedRecoveryError(
            f"segment is not in a terminal failed state ({state}): {segment}"
        )
    exit_path = segment / "manifest/run_exit_code"
    try:
        exit_code = int(exit_path.read_text(encoding="utf-8").strip())
    except (OSError, ValueError) as error:
        raise BoundedRecoveryError(
            f"failed segment lacks an integer exit artifact: {segment}"
        ) from error
    if exit_code == 0:
        raise BoundedRecoveryError(f"failed segment records exit zero: {segment}")
    analysis_path = segment / "manifest/fast_analysis.json"
    if analysis_path.is_file():
        analysis = bounded.corrected.fast.load_json(analysis_path)
        if analysis.get("passed") is True:
            raise BoundedRecoveryError(
                f"failed segment already has passing analysis: {segment}"
            )
    log = slurm_log_path(manifest)
    try:
        log_text = log.read_text(encoding="utf-8", errors="replace")
    except OSError as error:
        raise BoundedRecoveryError(f"failed segment log is unavailable: {log}") from error
    if FAILURE_SIGNATURE not in log_text:
        raise BoundedRecoveryError(
            f"failed segment lacks the zero-field signature: {segment}"
        )
    manifest_path = bounded.corrected.fast.segment_manifest(segment)
    return {
        "segment": str(segment.resolve()),
        "sequence": int(manifest["sequence"]),
        "job_id": job_id,
        "job_state": state,
        "exit_code": exit_code,
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": sha256_file(manifest_path),
        "slurm_log": str(log.resolve()),
        "slurm_log_sha256": sha256_file(log),
        "failure_signature": FAILURE_SIGNATURE,
    }


def authenticated_predecessor(
    segments: list[Path],
) -> tuple[Path, dict[str, object], dict[str, object]]:
    for segment in reversed(segments):
        analysis_path = segment / "manifest/fast_analysis.json"
        if not analysis_path.is_file():
            continue
        manifest = bounded.validate_bounded_segment(segment)
        analysis = bounded.corrected.fast.load_json(analysis_path)
        restart = analysis.get("terminal_restart")
        if (
            analysis.get("passed") is True
            and analysis.get("complete") is not True
            and isinstance(restart, dict)
            and restart.get("rank_zero")
            and restart.get("sha256")
            and float(analysis.get("final_time", -1.0))
            == float(restart.get("physical_time", -2.0))
        ):
            return segment, manifest, analysis
    raise BoundedRecoveryError("no earlier authenticated incomplete checkpoint exists")


def validate_prepared_recovery(segment: Path) -> dict[str, object]:
    manifest = bounded.validate_bounded_segment(segment)
    if (
        manifest.get("recovery_policy") != RECOVERY_POLICY
        or manifest.get("recovery_reason") != FAILURE_REASON
        or manifest.get("recovery_changes_physics") is not False
    ):
        raise BoundedRecoveryError(
            f"prepared segment is not a zero-field recovery: {segment}"
        )
    failed = manifest.get("supersedes_failed_segments")
    if not isinstance(failed, list) or not failed:
        raise BoundedRecoveryError(
            f"prepared recovery lacks failed-attempt evidence: {segment}"
        )
    return manifest


def recover_case(root: Path, case_id: str, *, submit: bool) -> Path | None:
    bounded.require_bounded_case(case_id)
    root = bounded.configure(root)
    bounded.corrected.validate_provenance(root, [case_id])
    require_inactive_case(case_id)
    segments = bounded.corrected.fast.segment_directories(root, case_id)
    if not segments:
        raise BoundedRecoveryError(f"{case_id} has no bounded lineage")

    latest = segments[-1]
    latest_manifest = bounded.validate_bounded_segment(latest)
    if latest_manifest.get("job_id") is None:
        validate_prepared_recovery(latest)
        if submit:
            bounded.corrected.fast.submit_segment(latest)
        else:
            print(latest)
        return latest

    state = normalized_state(
        bounded.corrected.fast.job_state(str(latest_manifest["job_id"]))
    )
    if state not in FAILED_STATES:
        print(f"no recovery required for {case_id}: latest state={state}")
        return None

    predecessor, predecessor_manifest, analysis = authenticated_predecessor(
        segments[:-1]
    )
    predecessor_sequence = int(predecessor_manifest["sequence"])
    failed_segments = [
        segment
        for segment in segments
        if int(bounded.corrected.fast.load_json(
            bounded.corrected.fast.segment_manifest(segment)
        )["sequence"]) > predecessor_sequence
    ]
    failed_records = [failed_attempt_record(segment) for segment in failed_segments]
    sequences = [int(record["sequence"]) for record in failed_records]
    expected = list(range(predecessor_sequence + 1, max(sequences) + 1))
    if sequences != expected:
        raise BoundedRecoveryError(
            f"{case_id} failed-attempt sequence is not contiguous: {sequences}"
        )

    restart = analysis["terminal_restart"]
    assert isinstance(restart, dict)
    next_sequence = max(sequences) + 1
    new_segment = bounded.prepare_segment(
        root,
        bounded.corrected.fast.campaign_cases()[case_id],
        next_sequence,
        float(analysis["final_time"]),
        Path(str(restart["rank_zero"])),
        str(restart["sha256"]),
    )
    manifest_path = bounded.corrected.fast.segment_manifest(new_segment)
    manifest = bounded.corrected.fast.load_json(manifest_path)
    analysis_path = predecessor / "manifest/fast_analysis.json"
    manifest.update(
        {
            "recovery_policy": RECOVERY_POLICY,
            "recovery_reason": FAILURE_REASON,
            "recovery_changes_physics": False,
            "restarts_from_authenticated_segment": str(predecessor.resolve()),
            "restarts_from_authenticated_analysis": str(analysis_path.resolve()),
            "restarts_from_authenticated_analysis_sha256": sha256_file(analysis_path),
            "supersedes_failed_segments": failed_records,
        }
    )
    bounded.corrected.fast.write_json(manifest_path, manifest)
    validate_prepared_recovery(new_segment)
    if submit:
        bounded.corrected.fast.submit_segment(new_segment)
    else:
        print(new_segment)
    return new_segment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=bounded.corrected.CAMPAIGN_ROOT,
        help="fixed corrected-production campaign root",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    recover = commands.add_parser(
        "recover", help="replay a failed R14/R15 attempt if and only if eligible"
    )
    recover.add_argument("case", choices=bounded.BOUNDED_CASES)
    recover.add_argument("--submit", action="store_true")

    args = parser.parse_args()
    try:
        recover_case(args.root, args.case, submit=args.submit)
        return 0
    except (
        BoundedRecoveryError,
        bounded.BoundedContinuationError,
        bounded.corrected.CorrectedFastError,
        bounded.corrected.fast.FastRunError,
        OSError,
        subprocess.SubprocessError,
        ValueError,
        KeyError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
