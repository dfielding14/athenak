#!/usr/bin/env python3
"""Launch independent final-retained-state CGL hyperbolicity audits.

This is a lean Slurm launcher, not a campaign controller.  It authenticates an
externally SHA-bound direct-fast assembled inventory, selects the latest
complete retained snapshot for active CGL cases, and prepares or submits one
one-node audit job per case.  R10 is admitted only when it is both explicitly
selected and marked exploratory.  Target-complete cases are the default;
``--include-partial`` explicitly admits current retained states.

Each attempt and its result live outside simulation roots.  Explicit ``retry``
creates a new attempt for failed audits or when a partial campaign has advanced
to a different retained state.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import stat
import subprocess
import sys
import tempfile
from typing import Iterable


DEFAULT_ACCOUNT = "ast207"
DEFAULT_PARTITION = "batch"
DEFAULT_WALLTIME = "02:00:00"
DEFAULT_CPUS_PER_TASK = 56
DEFAULT_AUDIT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/"
    "analysis/audit_cgl_hyperbolicity.py"
)
CASE_ID = re.compile(r"R(?:0[2-9]|1[0-7])")
SHA256 = re.compile(r"[0-9a-f]{64}")
RANK_DIRECTORY = re.compile(r"rank_(\d{8})")
JOB_ID = re.compile(r"[1-9]\d*(?:;[A-Za-z0-9_.-]+)?")
FAILED_STATES = {
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


class HyperbolicityLaunchError(RuntimeError):
    """A retained-state audit launch precondition failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def canonical_json(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return hashlib.sha256(canonical_json(value)).hexdigest()


def unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise HyperbolicityLaunchError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def reject_constant(value: str) -> object:
    raise HyperbolicityLaunchError(f"invalid JSON numeric constant: {value}")


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise HyperbolicityLaunchError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise HyperbolicityLaunchError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise HyperbolicityLaunchError(f"{label} must be a nonempty string")
    return value


def require_int(value: object, label: str, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise HyperbolicityLaunchError(f"{label} must be an integer >= {minimum}")
    return value


def require_finite(value: object, label: str) -> float:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(float(value))
    ):
        raise HyperbolicityLaunchError(f"{label} must be finite")
    return float(value)


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise HyperbolicityLaunchError(f"{label} must be a lowercase SHA-256")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_binding(path: Path) -> dict[str, object]:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise HyperbolicityLaunchError(f"artifact does not resolve: {path}") from error
    before = resolved.stat()
    if not stat.S_ISREG(before.st_mode):
        raise HyperbolicityLaunchError(f"artifact is not a regular file: {resolved}")
    digest = sha256_file(resolved)
    after = resolved.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise HyperbolicityLaunchError(f"artifact changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
        "sha256": digest,
    }


def load_bound_json(
    path: Path, label: str, expected_sha256: str | None = None
) -> tuple[dict[str, object], dict[str, object]]:
    binding = artifact_binding(path)
    if expected_sha256 is not None and binding["sha256"] != expected_sha256:
        raise HyperbolicityLaunchError(f"{label} SHA-256 differs from authority")
    payload = Path(str(binding["path"])).read_bytes()
    if (
        len(payload) != binding["size_bytes"]
        or hashlib.sha256(payload).hexdigest() != binding["sha256"]
    ):
        raise HyperbolicityLaunchError(f"{label} changed between binding and read")
    try:
        value = json.loads(
            payload.decode("utf-8"),
            object_pairs_hook=unique_object,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise HyperbolicityLaunchError(f"{label} is invalid JSON") from error
    return require_dict(value, label), binding


def verify_declared_binding(value: object, label: str) -> dict[str, object]:
    record = require_dict(value, label)
    expected_path = Path(require_text(record.get("path"), f"{label} path")).resolve()
    expected_size = require_int(record.get("size_bytes"), f"{label} size", 1)
    expected_sha = require_sha256(record.get("sha256"), f"{label} SHA-256")
    observed = artifact_binding(expected_path)
    if (
        observed["path"] != str(expected_path)
        or observed["size_bytes"] != expected_size
        or observed["sha256"] != expected_sha
    ):
        raise HyperbolicityLaunchError(f"{label} differs from its declared binding")
    return observed


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as stream:
        stream.write(payload)
        temporary = Path(stream.name)
    temporary.replace(path)


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except ValueError:
        return False


def is_lexically_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.expanduser().absolute().relative_to(parent.expanduser().absolute())
        return True
    except ValueError:
        return False


def inventory_context(
    analysis_argument: Path,
    expected_inventory_sha256: str,
    jobs_argument: Path | None,
    audit_script: Path,
) -> dict[str, object]:
    expected_sha = require_sha256(
        expected_inventory_sha256, "expected inventory SHA-256"
    )
    supplied = analysis_argument.expanduser().resolve()
    inventory_path = supplied if supplied.is_file() else supplied / "inventory.json"
    inventory, inventory_binding = load_bound_json(
        inventory_path, "direct-fast assembled inventory", expected_sha
    )
    if inventory.get("schema_version") != 1:
        raise HyperbolicityLaunchError("direct-fast assembled inventory schema differs")
    output = Path(
        require_text(inventory.get("output"), "assembled output")
    ).expanduser().resolve()
    if output != inventory_path.parent:
        raise HyperbolicityLaunchError(
            "assembled inventory path differs from declared output"
        )
    root = Path(
        require_text(inventory.get("root"), "campaign root")
    ).expanduser().resolve()
    simulation_root = root / "runs"
    jobs = (
        jobs_argument.expanduser().resolve()
        if jobs_argument is not None
        else output / "slurm-hyperbolicity"
    )
    if is_relative_to(output, simulation_root) or is_relative_to(jobs, simulation_root):
        raise HyperbolicityLaunchError(
            "assembled inventory, job manifests, and results must be outside runs/"
        )
    cases = require_dict(inventory.get("cases"), "assembled cases")
    matrix = verify_declared_binding(inventory.get("matrix"), "assembled matrix")
    adapter = verify_declared_binding(
        inventory.get("adapter"), "direct-fast reporter adapter"
    )
    audit = artifact_binding(audit_script)
    return {
        "root": root,
        "simulation_root": simulation_root,
        "analysis": output,
        "inventory_path": inventory_path,
        "inventory": inventory,
        "inventory_binding": inventory_binding,
        "cases": cases,
        "matrix": matrix,
        "adapter": adapter,
        "audit": audit,
        "jobs": jobs,
    }


def expand_case_tokens(tokens: Iterable[str], available: Iterable[str]) -> list[str]:
    requested: list[str] = []
    for raw_token in tokens:
        for token in raw_token.split(","):
            token = token.strip().upper()
            if not token:
                continue
            match = re.fullmatch(r"R(\d{2})-R(\d{2})", token)
            if match:
                start, end = (int(item) for item in match.groups())
                if start > end:
                    raise HyperbolicityLaunchError(
                        f"descending case range is invalid: {token}"
                    )
                requested.extend(f"R{number:02d}" for number in range(start, end + 1))
            else:
                requested.append(token)
    if not requested:
        requested = sorted(available)
    invalid = [case for case in requested if CASE_ID.fullmatch(case) is None]
    if invalid:
        raise HyperbolicityLaunchError(
            f"unsupported cases: {', '.join(sorted(set(invalid)))}"
        )
    return list(dict.fromkeys(requested))


def rank_file_record(path: Path, expected_size: int, label: str) -> dict[str, object]:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise HyperbolicityLaunchError(f"{label} does not resolve: {path}") from error
    before = resolved.stat()
    if not stat.S_ISREG(before.st_mode) or before.st_size != expected_size:
        raise HyperbolicityLaunchError(f"{label} size or file type differs")
    after = resolved.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise HyperbolicityLaunchError(f"{label} changed while it was inspected")
    return {
        "path": str(resolved),
        "size_bytes": after.st_size,
        "mtime_ns": after.st_mtime_ns,
    }


def authenticate_latest_snapshot(
    context: dict[str, object], case_id: str, case: dict[str, object]
) -> dict[str, object]:
    analysis = Path(str(context["analysis"]))
    case_dir = analysis / "cases" / case_id
    lineage_path = case_dir / "lineage.json"
    lineage, lineage_binding = load_bound_json(lineage_path, f"{case_id} lineage")
    if lineage != case:
        raise HyperbolicityLaunchError(
            f"{case_id} inventory record differs from bound lineage.json"
        )
    if lineage.get("case_id") != case_id:
        raise HyperbolicityLaunchError(f"{case_id} lineage identity differs")
    retained_errors = require_list(lineage.get("errors"), f"{case_id} lineage errors")
    if retained_errors:
        raise HyperbolicityLaunchError(f"{case_id} lineage retains assembly errors")
    input_binding = verify_declared_binding(lineage.get("input"), f"{case_id} input")
    identities = lineage.get("lineage_identities")
    if identities is not None:
        identity_record = require_dict(identities, f"{case_id} lineage identities")
        if identity_record.get("matrix_sha256") != [context["matrix"]["sha256"]]:
            raise HyperbolicityLaunchError(f"{case_id} matrix lineage identity differs")
        if identity_record.get("input_sha256") != [input_binding["sha256"]]:
            raise HyperbolicityLaunchError(f"{case_id} input lineage identity differs")

    output_roots: list[Path] = []
    segments = require_list(lineage.get("lineage"), f"{case_id} selected lineage")
    if not segments:
        raise HyperbolicityLaunchError(f"{case_id} selected lineage is empty")
    manifests: list[dict[str, object]] = []
    for order, value in enumerate(segments):
        segment = require_dict(value, f"{case_id} lineage segment {order}")
        if segment.get("order") != order:
            raise HyperbolicityLaunchError(f"{case_id} lineage segment order differs")
        output_roots.append(
            Path(
                require_text(segment.get("output"), f"{case_id} segment output")
            ).absolute()
        )
        manifests.append(
            verify_declared_binding(
                segment.get("manifest"), f"{case_id} lineage segment {order} manifest"
            )
        )

    snapshot_record = require_dict(lineage.get("snapshots"), f"{case_id} snapshots")
    snapshot_path = Path(
        require_text(snapshot_record.get("path"), f"{case_id} snapshot index path")
    ).resolve()
    expected_snapshot_path = case_dir / "snapshots.json"
    if snapshot_path != expected_snapshot_path:
        raise HyperbolicityLaunchError(f"{case_id} snapshot index path differs")
    index, snapshot_binding = load_bound_json(
        snapshot_path, f"{case_id} snapshot index"
    )
    if index.get("schema_version") != 1:
        raise HyperbolicityLaunchError(f"{case_id} snapshot index schema differs")
    snapshots = require_list(index.get("snapshots"), f"{case_id} snapshots")
    if require_int(index.get("snapshot_count"), f"{case_id} snapshot count") != len(
        snapshots
    ):
        raise HyperbolicityLaunchError(f"{case_id} snapshot count differs")
    if snapshot_record.get("snapshot_count") != len(snapshots):
        raise HyperbolicityLaunchError(f"{case_id} assembled snapshot count differs")

    complete: list[dict[str, object]] = []
    times: list[float] = []
    for position, value in enumerate(snapshots):
        group = require_dict(value, f"{case_id} snapshot {position}")
        times.append(require_finite(group.get("time"), f"{case_id} snapshot time"))
        flag = group.get("complete")
        if not isinstance(flag, bool):
            raise HyperbolicityLaunchError(
                f"{case_id} snapshot {position} complete flag is malformed"
            )
        if flag:
            complete.append(group)
    if any(right <= left for left, right in zip(times, times[1:])):
        raise HyperbolicityLaunchError(
            f"{case_id} snapshot times are not strictly increasing"
        )
    declared_complete = require_int(
        index.get("complete_snapshot_count"), f"{case_id} complete snapshot count"
    )
    if declared_complete != len(complete) or snapshot_record.get(
        "complete_snapshot_count"
    ) != len(complete):
        raise HyperbolicityLaunchError(f"{case_id} complete snapshot count differs")
    if not complete:
        raise HyperbolicityLaunchError(f"{case_id} has no complete retained snapshot")

    selected = complete[-1]
    expected_ranks = require_int(
        selected.get("expected_ranks"), f"{case_id} selected snapshot ranks", 1
    )
    members = require_list(
        selected.get("rank_files"), f"{case_id} selected snapshot rank files"
    )
    if len(members) != expected_ranks:
        raise HyperbolicityLaunchError(
            f"{case_id} selected snapshot rank count differs"
        )
    rank_inventory: list[dict[str, object]] = []
    rank_ids: list[int] = []
    for position, value in enumerate(members):
        member = require_dict(value, f"{case_id} selected snapshot rank {position}")
        declared_path = Path(
            require_text(member.get("path"), f"{case_id} snapshot rank path")
        )
        record = rank_file_record(
            declared_path,
            require_int(member.get("size_bytes"), f"{case_id} snapshot rank size", 1),
            f"{case_id} selected snapshot rank {position}",
        )
        resolved = Path(str(record["path"]))
        match = RANK_DIRECTORY.fullmatch(resolved.parent.name)
        if match is None:
            raise HyperbolicityLaunchError(
                f"{case_id} selected snapshot rank path is noncanonical"
            )
        rank_ids.append(int(match.group(1)))
        if not any(
            is_lexically_relative_to(declared_path, output) for output in output_roots
        ) or not is_relative_to(resolved, Path(str(context["simulation_root"]))):
            raise HyperbolicityLaunchError(
                f"{case_id} selected snapshot escapes selected lineage outputs"
            )
        rank_inventory.append(record)
    if rank_ids != list(range(expected_ranks)):
        raise HyperbolicityLaunchError(
            f"{case_id} selected snapshot rank inventory is not contiguous"
        )
    if len({record["path"] for record in rank_inventory}) != expected_ranks:
        raise HyperbolicityLaunchError(
            f"{case_id} selected snapshot rank inventory contains duplicates"
        )
    representative = Path(
        require_text(
            selected.get("representative"), f"{case_id} snapshot representative"
        )
    ).resolve()
    if representative != Path(str(rank_inventory[0]["path"])):
        raise HyperbolicityLaunchError(
            f"{case_id} selected snapshot representative is not rank zero"
        )

    model = require_dict(lineage.get("model_choices"), f"{case_id} model choices")
    passive = require_text(
        model.get("passive_delta"), f"{case_id} passive_delta"
    ).lower()
    if passive not in {"true", "false"}:
        raise HyperbolicityLaunchError(f"{case_id} passive_delta is not boolean")
    snapshot = {
        "time": require_finite(
            selected.get("time"), f"{case_id} selected snapshot time"
        ),
        "lineage_order": selected.get("lineage_order"),
        "representative": str(representative),
        "expected_ranks": expected_ranks,
        "rank_files": rank_inventory,
        "rank_inventory_sha256": canonical_sha256(rank_inventory),
    }
    selection_identity = {
        "case_id": case_id,
        "case_name": lineage.get("case_name"),
        "matrix_sha256": context["matrix"]["sha256"],
        "input_sha256": input_binding["sha256"],
        "passive_delta": passive,
        "snapshot": snapshot,
        "audit_script_sha256": context["audit"]["sha256"],
    }
    return {
        "case_id": case_id,
        "case_name": lineage.get("case_name"),
        "assembled_status": lineage.get("status"),
        "active": passive == "false",
        "input": input_binding,
        "lineage": lineage_binding,
        "segment_manifests": manifests,
        "snapshot_index": snapshot_binding,
        "snapshot": snapshot,
        "selection_sha256": canonical_sha256(selection_identity),
    }


def eligible_selections(
    context: dict[str, object],
    selectors: Iterable[str],
    exploratory: bool,
    include_partial: bool,
) -> list[dict[str, object]]:
    cases = require_dict(context["cases"], "assembled cases")
    selector_list = list(selectors)
    selected = expand_case_tokens(selector_list, cases)
    explicit = bool(selector_list)
    eligible: list[dict[str, object]] = []
    for case_id in selected:
        value = cases.get(case_id)
        if value is None:
            print(f"skip {case_id}: absent from partial assembled inventory")
            continue
        case = require_dict(value, f"{case_id} assembled case")
        if case.get("status") in {"not_started", "assembly_error"}:
            print(f"skip {case_id}: assembled status is {case.get('status')}")
            continue
        try:
            selection = authenticate_latest_snapshot(context, case_id, case)
        except HyperbolicityLaunchError as error:
            if "has no complete retained snapshot" in str(error):
                print(f"skip {case_id}: {error}")
                continue
            raise
        if not selection["active"]:
            print(f"skip {case_id}: passive_delta=true")
            continue
        if case_id == "R10" and not (explicit and exploratory):
            print(
                "skip R10: requires an explicit case selector and --exploratory"
            )
            continue
        if not include_partial and selection["assembled_status"] != "complete":
            print(
                f"skip {case_id}: assembled lineage is "
                f"{selection['assembled_status']}, not complete"
            )
            continue
        eligible.append(selection)
    return eligible


def attempt_directories(jobs: Path, case_id: str) -> list[Path]:
    parent = jobs / case_id
    if not parent.is_dir():
        return []
    return sorted(
        path
        for path in parent.glob("attempt-[0-9][0-9][0-9]")
        if path.is_dir() and (path / "manifest.json").is_file()
    )


def batch_script(manifest: dict[str, object]) -> str:
    command = manifest.get("command")
    if not isinstance(command, list) or not all(
        isinstance(item, str) for item in command
    ):
        raise HyperbolicityLaunchError("internal error: malformed audit command")
    attempt_dir = Path(str(manifest["attempt_dir"]))
    result = require_dict(manifest.get("result"), "result paths")
    result_path = Path(str(result["path"]))
    result_sha = Path(str(result["sha256_path"]))
    temporary = attempt_dir / ".result.json.tmp"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --account={manifest['account']}",
        f"#SBATCH --partition={manifest['partition']}",
        f"#SBATCH --time={manifest['walltime']}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={manifest['cpus_per_task']}",
        f"#SBATCH --job-name={manifest['job_name']}",
        f"#SBATCH --output={manifest['slurm_log']}",
        "",
        "set -uo pipefail",
        "umask 027",
        f"cd {shlex.quote(str(attempt_dir))}",
        f"export OMP_NUM_THREADS={int(manifest['cpus_per_task'])}",
        f"rm -f {shlex.quote(str(temporary))}",
        "set +e",
        (
            "srun --nodes=1 --ntasks=1 "
            f"--cpus-per-task={int(manifest['cpus_per_task'])} "
            + shlex.join(command)
            + f" > {shlex.quote(str(temporary))}"
        ),
        "status=$?",
        "set -e",
        f"if [[ -s {shlex.quote(str(temporary))} ]]; then",
        f"  mv {shlex.quote(str(temporary))} {shlex.quote(str(result_path))}",
        (
            f"  sha256sum {shlex.quote(str(result_path))} "
            f"> {shlex.quote(str(result_sha))}"
        ),
        "else",
        f"  rm -f {shlex.quote(str(temporary))}",
        "fi",
        (
            f"printf '%s\\n' \"$status\" > "
            f"{shlex.quote(str(attempt_dir / 'exit_code.txt'))}"
        ),
        "exit \"$status\"",
        "",
    ]
    return "\n".join(lines)


def prepare_attempt(
    *,
    context: dict[str, object],
    selection: dict[str, object],
    account: str,
    partition: str,
    walltime: str,
    cpus_per_task: int,
    python: Path,
    retry_of: str | None = None,
) -> Path:
    jobs = Path(str(context["jobs"]))
    case_id = str(selection["case_id"])
    existing = attempt_directories(jobs, case_id)
    number = len(existing)
    attempt_dir = jobs / case_id / f"attempt-{number:03d}"
    if attempt_dir.exists():
        raise HyperbolicityLaunchError(f"attempt path already exists: {attempt_dir}")
    attempt_dir.mkdir(parents=True)
    log_dir = jobs / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    result_path = attempt_dir / "result.json"
    result_sha = attempt_dir / "result.sha256"
    snapshot = require_dict(selection["snapshot"], "selected snapshot")
    rank_files = require_list(snapshot["rank_files"], "selected rank inventory")
    rank_paths = [
        require_text(require_dict(value, "selected rank file").get("path"), "rank path")
        for value in rank_files
    ]
    audit_path = require_text(
        require_dict(context["audit"], "audit binding").get("path"), "audit path"
    )
    command = [
        str(python.expanduser().absolute()),
        audit_path,
        *rank_paths,
        "--format",
        "json",
        "--hash-inputs",
    ]
    manifest: dict[str, object] = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_direct_fast_hyperbolicity_job",
        "prepared_utc": utc_now(),
        "case_id": case_id,
        "case_name": selection.get("case_name"),
        "attempt": number,
        "attempt_dir": str(attempt_dir),
        "assembled_status": selection.get("assembled_status"),
        "inventory": context["inventory_binding"],
        "matrix": context["matrix"],
        "reporter_adapter": context["adapter"],
        "case_input": selection["input"],
        "case_lineage": selection["lineage"],
        "segment_manifests": selection["segment_manifests"],
        "snapshot_index": selection["snapshot_index"],
        "selected_snapshot": snapshot,
        "selection_sha256": selection["selection_sha256"],
        "audit_script": context["audit"],
        "launcher": artifact_binding(Path(__file__)),
        "python": str(python.expanduser().absolute()),
        "command": command,
        "result": {
            "path": str(result_path),
            "sha256_path": str(result_sha),
        },
        "account": account,
        "partition": partition,
        "walltime": walltime,
        "nodes": 1,
        "ntasks": 1,
        "cpus_per_task": cpus_per_task,
        "job_name": f"cglfh_{case_id}_a{number:03d}",
        "slurm_log": str(log_dir / f"{case_id}.attempt-{number:03d}.%j.log"),
        "retry_of": retry_of,
        "job_id": None,
    }
    write_json(attempt_dir / "manifest.json", manifest)
    script = attempt_dir / "run.sbatch"
    script.write_text(batch_script(manifest), encoding="utf-8")
    script.chmod(0o755)
    return attempt_dir


def submit_attempt(attempt_dir: Path) -> str:
    manifest_path = attempt_dir / "manifest.json"
    manifest, _ = load_bound_json(manifest_path, "attempt manifest")
    if manifest.get("job_id") is not None:
        raise HyperbolicityLaunchError(f"attempt already submitted: {attempt_dir}")
    completed = subprocess.run(
        ["/usr/bin/sbatch", "--parsable", str(attempt_dir / "run.sbatch")],
        check=True,
        text=True,
        capture_output=True,
    )
    response = completed.stdout.strip()
    if JOB_ID.fullmatch(response) is None:
        raise HyperbolicityLaunchError(f"unexpected sbatch response: {response!r}")
    job_id = response.split(";", 1)[0]
    manifest["job_id"] = job_id
    manifest["submitted_utc"] = utc_now()
    write_json(manifest_path, manifest)
    return job_id


def normalized_state(value: str) -> str:
    return value.strip().split("+", 1)[0].split()[0].upper() if value.strip() else ""


def validate_result(attempt_dir: Path, manifest: dict[str, object]) -> str | None:
    result_record = require_dict(manifest.get("result"), "attempt result")
    result_path = Path(require_text(result_record.get("path"), "result path"))
    result_sha_path = Path(
        require_text(result_record.get("sha256_path"), "result SHA-256 path")
    )
    if (
        result_path != attempt_dir / "result.json"
        or result_sha_path != attempt_dir / "result.sha256"
    ):
        return "result paths differ from attempt directory"
    if not result_path.is_file() or not result_sha_path.is_file():
        return "result or result SHA-256 is missing"
    try:
        expected_sha = result_sha_path.read_text(encoding="utf-8").split()[0]
    except (OSError, IndexError):
        return "result SHA-256 record is unreadable"
    if (
        SHA256.fullmatch(expected_sha) is None
        or sha256_file(result_path) != expected_sha
    ):
        return "result SHA-256 differs"
    try:
        result, _ = load_bound_json(result_path, "hyperbolicity result", expected_sha)
        provenance = require_dict(result.get("provenance"), "result provenance")
        audit = require_dict(manifest.get("audit_script"), "audit script binding")
        snapshot = require_dict(manifest.get("selected_snapshot"), "selected snapshot")
        rank_files = [
            require_dict(value, "selected rank file")
            for value in require_list(snapshot.get("rank_files"), "selected rank files")
        ]
        expected_paths = [str(value["path"]) for value in rank_files]
        if (
            provenance.get("script_path") != audit.get("path")
            or provenance.get("script_sha256") != audit.get("sha256")
            or provenance.get("input_patterns") != expected_paths
            or provenance.get("hash_inputs") is not True
        ):
            return "result provenance differs from attempt manifest"
        snapshots = require_list(result.get("snapshots"), "result snapshots")
        if len(snapshots) != 1:
            return "result does not contain exactly one retained snapshot"
        observed = require_dict(snapshots[0], "result retained snapshot")
        if observed.get("active_cgl_signal_speed") is not True:
            return "audited snapshot is not active CGL signal-speed state"
        if not math.isclose(
            require_finite(observed.get("time"), "result snapshot time"),
            require_finite(snapshot.get("time"), "selected snapshot time"),
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ):
            return "result snapshot time differs"
        observed_files = [
            require_dict(value, "result rank file")
            for value in require_list(observed.get("rank_files"), "result rank files")
        ]
        expected_profiles = [
            {
                "path": value["path"],
                "size_bytes": value["size_bytes"],
                "mtime_ns": value["mtime_ns"],
            }
            for value in rank_files
        ]
        observed_profiles = [
            {
                "path": value.get("path"),
                "size_bytes": value.get("size_bytes"),
                "mtime_ns": value.get("mtime_ns"),
            }
            for value in observed_files
        ]
        if observed_profiles != expected_profiles:
            return "result rank inventory differs from selected retained state"
        if any(
            SHA256.fullmatch(str(value.get("sha256"))) is None
            for value in observed_files
        ):
            return "result rank inventory lacks content SHA-256"
        if observed.get("input_inventory_sha256") != canonical_sha256(observed_files):
            return "result rank inventory digest differs"
        if observed.get("ranks_contiguous_from_zero") is not True:
            return "result rank inventory is not contiguous"
    except (HyperbolicityLaunchError, OSError, ValueError) as error:
        return str(error)
    return None


def audit_disposition(attempt_dir: Path, manifest: dict[str, object]) -> str:
    if validate_result(attempt_dir, manifest) is not None:
        return "-"
    result_record = require_dict(manifest["result"], "result")
    result, _ = load_bound_json(Path(str(result_record["path"])), "result")
    snapshot = require_dict(
        require_list(result["snapshots"], "snapshots")[0], "snapshot"
    )
    aggregate = require_dict(snapshot.get("aggregate"), "aggregate")
    if require_int(
        aggregate.get("nonfinite_discriminant"), "nonfinite discriminants"
    ) > 0:
        return "NONFINITE"
    if require_int(aggregate.get("negative"), "negative discriminants") > 0:
        return "NEGATIVE"
    return "HYPERBOLIC"


def scheduler_state(manifest: dict[str, object]) -> str:
    job_id = manifest.get("job_id")
    if not isinstance(job_id, str):
        return "PREPARED"
    queued = subprocess.run(
        ["/usr/bin/squeue", "-h", "-j", job_id, "-o", "%T"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if queued:
        return normalized_state(queued.splitlines()[0])
    accounted = subprocess.run(
        ["/usr/bin/sacct", "-X", "-n", "-P", "-j", job_id, "-o", "State"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if accounted:
        return normalized_state(accounted.splitlines()[0].split("|", 1)[0])
    return "UNKNOWN"


def attempt_state(attempt_dir: Path, manifest: dict[str, object]) -> str:
    exit_path = attempt_dir / "exit_code.txt"
    if exit_path.is_file():
        try:
            exit_code = int(exit_path.read_text(encoding="utf-8").strip())
        except (OSError, ValueError):
            return "INVALID_EXIT_RECORD"
        if exit_code != 0:
            return f"FAILED_EXIT_{exit_code}"
        invalid = validate_result(attempt_dir, manifest)
        return "COMPLETED" if invalid is None else "INVALID_RESULT"
    state = scheduler_state(manifest)
    return "INVALID_EXIT_RECORD" if state == "COMPLETED" else state


def current_selections(context: dict[str, object]) -> dict[str, dict[str, object]]:
    cases = require_dict(context["cases"], "assembled cases")
    result: dict[str, dict[str, object]] = {}
    for case_id, value in cases.items():
        if CASE_ID.fullmatch(case_id) is None or not isinstance(value, dict):
            continue
        if value.get("status") in {"not_started", "assembly_error"}:
            continue
        try:
            result[case_id] = authenticate_latest_snapshot(context, case_id, value)
        except HyperbolicityLaunchError as error:
            if "has no complete retained snapshot" not in str(error):
                raise
    return result


def state_with_staleness(
    attempt_dir: Path,
    manifest: dict[str, object],
    current: dict[str, dict[str, object]],
) -> str:
    state = attempt_state(attempt_dir, manifest)
    case_id = str(manifest.get("case_id"))
    selection = current.get(case_id)
    if selection is None:
        return f"STALE_{state}"
    if manifest.get("selection_sha256") != selection.get("selection_sha256"):
        return f"STALE_{state}"
    return state


def print_attempt(
    case_id: str,
    attempt_dir: Path,
    current: dict[str, dict[str, object]],
) -> str:
    manifest, _ = load_bound_json(attempt_dir / "manifest.json", "attempt manifest")
    state = state_with_staleness(attempt_dir, manifest, current)
    disposition = (
        audit_disposition(attempt_dir, manifest)
        if state == "COMPLETED"
        else "-"
    )
    selected = require_dict(manifest.get("selected_snapshot"), "selected snapshot")
    print(
        f"{case_id}\t{attempt_dir.name}\t{manifest.get('job_id') or '-'}\t"
        f"{state}\t{selected.get('time')}\t{disposition}\t"
        f"{manifest.get('slurm_log')}"
    )
    return state


def launch(args: argparse.Namespace) -> int:
    context = inventory_context(
        args.analysis, args.inventory_sha256, args.jobs_dir, args.audit_script
    )
    selections = eligible_selections(
        context, args.cases, args.exploratory, args.include_partial
    )
    jobs = Path(str(context["jobs"]))
    for selection in selections:
        case_id = str(selection["case_id"])
        if attempt_directories(jobs, case_id):
            print(f"skip {case_id}: attempt exists; use retry")
            continue
        attempt = prepare_attempt(
            context=context,
            selection=selection,
            account=args.account,
            partition=args.partition,
            walltime=args.walltime,
            cpus_per_task=args.cpus_per_task,
            python=args.python,
        )
        if args.submit:
            print(f"submitted {case_id}: {submit_attempt(attempt)}")
        else:
            print(f"prepared {case_id}: {attempt}")
    return 0


def status(args: argparse.Namespace) -> int:
    context = inventory_context(
        args.analysis, args.inventory_sha256, args.jobs_dir, args.audit_script
    )
    jobs = Path(str(context["jobs"]))
    cases = require_dict(context["cases"], "assembled cases")
    selected = expand_case_tokens(args.cases, cases)
    current = current_selections(context)
    print("CASE\tATTEMPT\tJOB_ID\tSTATE\tRETAINED_T\tAUDIT\tLOG")
    for case_id in selected:
        attempts = attempt_directories(jobs, case_id)
        attempts = attempts if args.all_attempts else attempts[-1:]
        if not attempts:
            print(f"{case_id}\t-\t-\tNOT_PREPARED\t-\t-\t-")
        for attempt in attempts:
            print_attempt(case_id, attempt, current)
    return 0


def retryable_state(state: str) -> bool:
    base = state.removeprefix("STALE_")
    return (
        state.startswith("STALE_")
        or base.startswith("FAILED_EXIT_")
        or base.startswith("INVALID_")
        or base in FAILED_STATES
    )


def retry(args: argparse.Namespace) -> int:
    context = inventory_context(
        args.analysis, args.inventory_sha256, args.jobs_dir, args.audit_script
    )
    selections = eligible_selections(
        context, args.cases, args.exploratory, args.include_partial
    )
    current = {str(value["case_id"]): value for value in selections}
    jobs = Path(str(context["jobs"]))
    for selection in selections:
        case_id = str(selection["case_id"])
        attempts = attempt_directories(jobs, case_id)
        retry_of: str | None = None
        settings = {
            "account": args.account,
            "partition": args.partition,
            "walltime": args.walltime,
            "cpus_per_task": args.cpus_per_task,
            "python": args.python,
        }
        if attempts:
            latest = attempts[-1]
            manifest, _ = load_bound_json(latest / "manifest.json", "attempt manifest")
            state = state_with_staleness(latest, manifest, current)
            if state == "PREPARED":
                attempt = latest
            elif retryable_state(state):
                retry_of = str(latest)
                settings = {
                    "account": str(manifest.get("account", args.account)),
                    "partition": str(manifest.get("partition", args.partition)),
                    "walltime": str(manifest.get("walltime", args.walltime)),
                    "cpus_per_task": int(
                        manifest.get("cpus_per_task", args.cpus_per_task)
                    ),
                    "python": Path(str(manifest.get("python", args.python))),
                }
                attempt = prepare_attempt(
                    context=context,
                    selection=selection,
                    retry_of=retry_of,
                    **settings,
                )
            else:
                print(f"skip {case_id}: latest attempt is {state}")
                continue
        else:
            attempt = prepare_attempt(
                context=context,
                selection=selection,
                **settings,
            )
        if args.submit:
            print(f"submitted {case_id}: {submit_attempt(attempt)}")
        else:
            print(f"retry-ready {case_id}: {attempt}")
    return 0


def add_selection_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("cases", nargs="*", help="case IDs or ranges; default: all")
    parser.add_argument(
        "--exploratory",
        action="store_true",
        help="admit explicitly selected R10 as exploratory-only",
    )
    parser.add_argument(
        "--include-partial",
        action="store_true",
        help="admit latest complete retained states from in-progress/partial cases",
    )


def add_slurm_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--account", default=DEFAULT_ACCOUNT)
    parser.add_argument("--partition", default=DEFAULT_PARTITION)
    parser.add_argument("--walltime", default=DEFAULT_WALLTIME)
    parser.add_argument("--cpus-per-task", type=int, default=DEFAULT_CPUS_PER_TASK)
    parser.add_argument(
        "--python",
        type=Path,
        default=Path(sys.executable),
        help="Python used inside audit jobs; default: current Python",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="submit with sbatch; otherwise only prepare provenance records",
    )


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "analysis",
        type=Path,
        help="direct-fast assembled output directory or its inventory.json",
    )
    command.add_argument(
        "--inventory-sha256",
        required=True,
        help="externally recorded SHA-256 of the assembled inventory.json",
    )
    command.add_argument(
        "--jobs-dir",
        type=Path,
        help="external job/result directory; default: ANALYSIS/slurm-hyperbolicity",
    )
    command.set_defaults(audit_script=DEFAULT_AUDIT)
    subcommands = command.add_subparsers(dest="command", required=True)

    launch_command = subcommands.add_parser(
        "launch", help="prepare or submit one latest-state audit per eligible case"
    )
    add_selection_options(launch_command)
    add_slurm_options(launch_command)

    status_command = subcommands.add_parser(
        "status", help="show current state of prepared/submitted audits"
    )
    status_command.add_argument("cases", nargs="*", help="case IDs or ranges")
    status_command.add_argument(
        "--all-attempts", action="store_true", help="show every attempt"
    )

    retry_command = subcommands.add_parser(
        "retry", help="prepare or submit missing, failed, invalid, or stale audits"
    )
    add_selection_options(retry_command)
    add_slurm_options(retry_command)
    return command


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "launch":
            return launch(args)
        if args.command == "status":
            return status(args)
        return retry(args)
    except (
        HyperbolicityLaunchError,
        OSError,
        subprocess.CalledProcessError,
        ValueError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
