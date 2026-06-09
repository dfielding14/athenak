#!/usr/bin/env python3
"""Prepare and authenticate exact-t=9 CT-only replays for Stage I fast lineages.

The replay lane is deliberately separate from the accepted science histories.
Each replay starts from the latest complete native restart group before t=9
within the selected assembled lineage, changes only ``time/tlim`` to 9, and
retains a complete exact-t=9 native restart group for the reviewed CT audit.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
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


SCRIPT_PATH = Path(__file__).resolve()
FAST_PATH = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast.py")
FAST_SPEC = importlib.util.spec_from_file_location("_cgl_ct_replay_fast", FAST_PATH)
if FAST_SPEC is None or FAST_SPEC.loader is None:
    raise RuntimeError(f"cannot load direct-fast utility: {FAST_PATH}")
fast = importlib.util.module_from_spec(FAST_SPEC)
sys.modules[FAST_SPEC.name] = fast
FAST_SPEC.loader.exec_module(fast)

SCHEMA_VERSION = 1
RUN_RECORD_TYPE = "stage-i-direct-fast-ct-exact-state-replay-run"
COMPLETION_RECORD_TYPE = "stage-i-direct-fast-ct-exact-state-replay-completion"
INVENTORY_RECORD_TYPE = "stage-i-direct-fast-ct-exact-state-replay-inventory"
TARGET_TIME = 9.0
CASE_ID = re.compile(r"R(?:0[2-9]|1[0-7])")
RANK_DIRECTORY = re.compile(r"rank_(\d{8})")
JOB_ID = re.compile(r"[1-9][0-9]*")
SHA256 = re.compile(r"[0-9a-f]{64}")
TERMINAL_FAILURE_CASES = {"R14", "R15"}
ACTIVE_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING"}


class CtReplayError(RuntimeError):
    """Raised when replay preparation or authentication fails."""


def utc_now() -> str:
    return fast.utc_now()


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise CtReplayError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise CtReplayError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise CtReplayError(f"{label} must be nonempty text")
    return value


def require_int(value: object, label: str, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise CtReplayError(f"{label} must be an integer >= {minimum}")
    return value


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256.fullmatch(value) is None:
        raise CtReplayError(f"{label} must be a lowercase SHA-256")
    return value


def require_finite(value: object, label: str) -> float:
    if (
        not isinstance(value, (int, float))
        or isinstance(value, bool)
        or not math.isfinite(float(value))
    ):
        raise CtReplayError(f"{label} must be finite")
    return float(value)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while block := stream.read(1024 * 1024):
            digest.update(block)
    return digest.hexdigest()


def stable_profile(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def binding(path: Path, label: str, expected_sha256: str | None = None) -> dict[str, object]:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise CtReplayError(f"{label} does not resolve: {path}") from error
    before = resolved.stat()
    if not stat.S_ISREG(before.st_mode):
        raise CtReplayError(f"{label} is not a regular file: {resolved}")
    digest = sha256_file(resolved)
    after = resolved.stat()
    if stable_profile(before) != stable_profile(after):
        raise CtReplayError(f"{label} changed while hashing: {resolved}")
    if expected_sha256 is not None and digest != require_sha256(
        expected_sha256, f"{label} expected SHA-256"
    ):
        raise CtReplayError(f"{label} SHA-256 differs")
    return {
        "path": str(resolved),
        "size_bytes": after.st_size,
        "sha256": digest,
    }


def verify_binding(value: object, label: str) -> dict[str, object]:
    declared = require_dict(value, f"{label} binding")
    current = binding(
        Path(require_text(declared.get("path"), f"{label} path")),
        label,
        require_sha256(declared.get("sha256"), f"{label} SHA-256"),
    )
    if current["size_bytes"] != require_int(
        declared.get("size_bytes"), f"{label} size"
    ):
        raise CtReplayError(f"{label} size differs")
    return current


def load_json(
    path: Path, label: str, expected_sha256: str | None = None
) -> tuple[dict[str, object], dict[str, object]]:
    observed = binding(path, label, expected_sha256)
    try:
        value = json.loads(Path(str(observed["path"])).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CtReplayError(f"{label} is not valid UTF-8 JSON") from error
    return require_dict(value, label), observed


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def write_json(path: Path, value: object, *, replace: bool = False) -> None:
    payload = (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False).encode("utf-8")
        + b"\n"
    )
    if path.exists() and not replace:
        if path.read_bytes() != payload:
            raise CtReplayError(f"refusing to replace differing artifact: {path}")
        return
    atomic_write(path, payload)


def write_timestamped_json(
    path: Path,
    value: dict[str, object],
    timestamp_key: str,
) -> dict[str, object]:
    """Write once while allowing exact idempotent finalization retries."""

    if path.exists():
        retained, _ = load_json(path, f"retained {path.name}")
        comparison = dict(retained)
        timestamp = comparison.pop(timestamp_key, None)
        if not isinstance(timestamp, str) or comparison != value:
            raise CtReplayError(f"retained artifact differs: {path}")
        return retained
    result = {timestamp_key: utc_now(), **value}
    write_json(path, result)
    return result


def expand_cases(values: Iterable[str], present: dict[str, object]) -> list[str]:
    requested: list[str] = []
    for value in values:
        for token in value.split(","):
            token = token.strip()
            if not token:
                continue
            if "-" in token:
                start, end = token.split("-", 1)
                if CASE_ID.fullmatch(start) is None or CASE_ID.fullmatch(end) is None:
                    raise CtReplayError(f"invalid case range: {token}")
                first = int(start[1:])
                last = int(end[1:])
                if first > last:
                    raise CtReplayError(f"descending case range: {token}")
                requested.extend(f"R{number:02d}" for number in range(first, last + 1))
            else:
                if CASE_ID.fullmatch(token) is None:
                    raise CtReplayError(f"invalid case ID: {token}")
                requested.append(token)
    selected = requested or sorted(present)
    result: list[str] = []
    for case_id in selected:
        if case_id in TERMINAL_FAILURE_CASES:
            continue
        if case_id not in present:
            raise CtReplayError(f"case is absent from primary inventory: {case_id}")
        if case_id not in result:
            result.append(case_id)
    return result


def rank_group(root: Path, name: str, rank_count: int, label: str) -> list[dict[str, object]]:
    paths = [root / f"rank_{rank:08d}" / name for rank in range(rank_count)]
    if any(not path.is_file() or path.stat().st_size == 0 for path in paths):
        raise CtReplayError(f"{label} rank group is incomplete")
    return [
        {**binding(path, f"{label} rank {rank}"), "rank": rank}
        for rank, path in enumerate(paths)
    ]


def candidate_parent_groups(
    case_id: str, case: dict[str, object], matrix_sha256: str
) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    for order, value in enumerate(require_list(case.get("lineage"), f"{case_id} lineage")):
        segment = require_dict(value, f"{case_id} lineage segment {order}")
        if segment.get("kind") != "fast" or segment.get("order") != order:
            raise CtReplayError(f"{case_id} selected lineage segment {order} is malformed")
        if segment.get("case_id") != case_id or segment.get("case_name") != case.get(
            "case_name"
        ):
            raise CtReplayError(f"{case_id} selected lineage identity differs")
        if segment.get("matrix_sha256") != matrix_sha256:
            raise CtReplayError(f"{case_id} selected lineage matrix differs")
        segment_binding = verify_binding(
            segment.get("manifest"), f"{case_id} segment {order} manifest"
        )
        manifest, _ = load_json(
            Path(str(segment_binding["path"])),
            f"{case_id} segment {order} manifest",
            str(segment_binding["sha256"]),
        )
        output = Path(require_text(segment.get("output"), f"{case_id} segment output"))
        if Path(require_text(manifest.get("output_dir"), "segment manifest output")) != output:
            raise CtReplayError(f"{case_id} segment {order} output differs")
        rank_count = require_int(segment.get("ranks"), f"{case_id} segment ranks", 1)
        if require_int(manifest.get("ranks"), "segment manifest ranks", 1) != rank_count:
            raise CtReplayError(f"{case_id} segment {order} rank count differs")
        rank_zero = output / "rst/rank_00000000"
        if not rank_zero.is_dir():
            continue
        for path in sorted(rank_zero.glob("*.rst")):
            time_value = fast.restart_time(path)
            if time_value >= TARGET_TIME:
                continue
            siblings = [
                output / "rst" / f"rank_{rank:08d}" / path.name
                for rank in range(rank_count)
            ]
            if any(not sibling.is_file() or sibling.stat().st_size == 0 for sibling in siblings):
                continue
            candidates.append({
                "time": time_value,
                "name": path.name,
                "rank_count": rank_count,
                "root": str((output / "rst").resolve()),
                "source_segment_manifest": segment_binding,
                "source_lineage_order": order,
                "source_segment": segment.get("segment"),
            })
    return candidates


def select_parent_group(
    case_id: str, case: dict[str, object], matrix_sha256: str
) -> dict[str, object]:
    candidates = candidate_parent_groups(case_id, case, matrix_sha256)
    if not candidates:
        raise CtReplayError(f"{case_id} has no complete selected-lineage restart before t=9")
    selected = max(
        candidates,
        key=lambda value: (
            float(value["time"]),
            int(value["source_lineage_order"]),
            str(value["name"]),
        ),
    )
    selected["rank_files"] = rank_group(
        Path(str(selected["root"])),
        str(selected["name"]),
        int(selected["rank_count"]),
        f"{case_id} parent restart",
    )
    return selected


def case_identity(
    inventory: dict[str, object], case_id: str, case: dict[str, object]
) -> dict[str, object]:
    matrix = verify_binding(inventory.get("matrix"), "primary matrix")
    executable = verify_binding(case.get("execution_executable"), f"{case_id} executable")
    authority = require_dict(case.get("execution_authority"), f"{case_id} authority")
    input_binding = verify_binding(authority.get("input"), f"{case_id} input")
    identities = require_dict(case.get("lineage_identities"), f"{case_id} identities")
    if identities.get("matrix_sha256") != [matrix["sha256"]]:
        raise CtReplayError(f"{case_id} matrix identity differs")
    if identities.get("input_sha256") != [input_binding["sha256"]]:
        raise CtReplayError(f"{case_id} input identity differs")
    if identities.get("executable_sha256") != [executable["sha256"]]:
        raise CtReplayError(f"{case_id} executable identity differs")
    return {
        "matrix": matrix,
        "executable": executable,
        "input": input_binding,
    }


def batch_script_text(manifest: dict[str, object]) -> str:
    allocation = require_dict(manifest["allocation"], "allocation")
    identity = require_dict(manifest["execution_identity"], "execution identity")
    parent = require_dict(manifest["parent_restart"], "parent restart")
    paths = require_dict(manifest["paths"], "paths")
    tools = require_dict(manifest["tools"], "tools")
    rank_files = require_list(parent["rank_files"], "parent rank files")
    representative = require_dict(rank_files[0], "parent rank zero")
    job = require_dict(manifest["job"], "job")
    return f"""#!/bin/bash
#SBATCH -J cglct9_{manifest['case_id']}
#SBATCH -A {job['account']}
#SBATCH -p {job['partition']}
#SBATCH -t {job['walltime']}
#SBATCH -N {allocation['nodes']}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1
#SBATCH -o {shlex.quote(str(paths['slurm_log']))}

set -euo pipefail
ATHENA={shlex.quote(str(require_dict(identity['executable'], 'executable')['path']))}
INPUT={shlex.quote(str(require_dict(identity['input'], 'input')['path']))}
PRIMARY={shlex.quote(str(require_dict(manifest['primary_inventory'], 'primary inventory')['path']))}
REPLAY_TOOL={shlex.quote(str(require_dict(tools['replay_tool'], 'replay tool')['path']))}
RESTART={shlex.quote(str(representative['path']))}
RUN_DIR={shlex.quote(str(paths['run_dir']))}
OUT_DIR={shlex.quote(str(paths['output_dir']))}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * 8))"

require_sha() {{
  local expected="$1" path="$2" label="$3" actual
  test -f "$path" || {{ echo "missing $label: $path" >&2; exit 1; }}
  actual="$(sha256sum "$path" | awk '{{print $1}}')"
  test "$actual" = "$expected" || {{ echo "$label checksum mismatch: $path" >&2; exit 1; }}
}}

require_sha {require_dict(identity['executable'], 'executable')['sha256']} "$ATHENA" executable
require_sha {require_dict(identity['input'], 'input')['sha256']} "$INPUT" input
require_sha {require_dict(manifest['primary_inventory'], 'primary inventory')['sha256']} "$PRIMARY" primary_inventory
require_sha {require_dict(tools['replay_tool'], 'replay tool')['sha256']} "$REPLAY_TOOL" replay_tool
require_sha {representative['sha256']} "$RESTART" parent_restart_rank_zero
test "$NNODES" -eq {allocation['nodes']}
test "$NRANKS" -eq {allocation['ranks']}
mkdir "$OUT_DIR"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1

date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ" > "$RUN_DIR/run_environment.txt"
module -t list 2>&1 >> "$RUN_DIR/run_environment.txt"
set +e
srun -N "$NNODES" -n "$NRANKS" --ntasks-per-node=8 \\
  -c 7 --threads-per-core=1 --cpu-bind=threads \\
  --gpus-per-task=1 --gpu-bind=closest \\
  "$ATHENA" -r "$RESTART" -d "$OUT_DIR" -t 03:50:00 \\
  job/basename={shlex.quote(str(manifest['run_basename']))} time/tlim=9.0
run_rc=$?
set -e
echo "$run_rc" > "$RUN_DIR/run_exit_code"
date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "$RUN_DIR/run_environment.txt"
exit "$run_rc"
"""


def prepare_case(
    output_root: Path,
    primary_binding: dict[str, object],
    inventory: dict[str, object],
    case_id: str,
    case: dict[str, object],
    account: str,
    partition: str,
    walltime: str,
) -> Path:
    if case.get("status") != "complete":
        raise CtReplayError(f"{case_id} is not complete")
    identity = case_identity(inventory, case_id, case)
    parent = select_parent_group(case_id, case, str(identity["matrix"]["sha256"]))
    rank_count = int(parent["rank_count"])
    if rank_count % 8 != 0:
        raise CtReplayError(f"{case_id} replay rank count is not divisible by 8")
    run_dir = output_root / "cases" / case_id / "t9"
    run_dir.mkdir(parents=True, exist_ok=False)
    output_dir = run_dir / "output"
    lineage_path = Path(str(inventory["output"])) / "cases" / case_id / "lineage.json"
    lineage_binding = binding(lineage_path, f"{case_id} lineage")
    if json.loads(lineage_path.read_text(encoding="utf-8")) != case:
        raise CtReplayError(f"{case_id} lineage differs from primary inventory")
    run_basename = f"E03_ct_exact_t9_{case['case_name']}"
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "record_type": RUN_RECORD_TYPE,
        "prepared_utc": utc_now(),
        "case_id": case_id,
        "case_name": case.get("case_name"),
        "purpose": "CT-only exact-state replay; excluded from accepted science histories",
        "primary_inventory": primary_binding,
        "primary_case_lineage": lineage_binding,
        "execution_identity": identity,
        "allocation": {
            "nodes": rank_count // 8,
            "ranks_per_node": 8,
            "ranks": rank_count,
        },
        "parent_restart": parent,
        "target_time": TARGET_TIME,
        "command_line_overrides": ["time/tlim=9.0"],
        "run_basename": run_basename,
        "paths": {
            "run_dir": str(run_dir.resolve()),
            "output_dir": str(output_dir.resolve()),
            "slurm_log": str((run_dir / "cglct9.%j.log").resolve()),
        },
        "tools": {
            "replay_tool": binding(SCRIPT_PATH, "CT replay tool"),
            "restart_parser": binding(FAST_PATH, "direct-fast restart parser"),
        },
        "job": {
            "account": account,
            "partition": partition,
            "walltime": walltime,
            "job_id": None,
        },
    }
    manifest_path = run_dir / "replay_run.json"
    write_json(manifest_path, manifest)
    script_path = run_dir / "run.sbatch"
    atomic_write(script_path, batch_script_text(manifest).encode("utf-8"))
    script_path.chmod(0o755)
    return run_dir


def prepare(args: argparse.Namespace) -> Path:
    primary, primary_binding = load_json(
        args.inventory, "primary assembled inventory", args.inventory_sha256
    )
    if primary.get("schema_version") != 1:
        raise CtReplayError("primary assembled inventory schema differs")
    output = Path(require_text(primary.get("output"), "primary output")).resolve()
    if Path(str(primary_binding["path"])) != output / "inventory.json":
        raise CtReplayError("primary assembled inventory path differs")
    cases = require_dict(primary.get("cases"), "primary cases")
    selected = expand_cases(args.cases, cases)
    root = args.output.expanduser().resolve(strict=False)
    if root.exists():
        raise CtReplayError(f"replay output already exists: {root}")
    root.mkdir(parents=True)
    for case_id in selected:
        prepare_case(
            root,
            primary_binding,
            primary,
            case_id,
            require_dict(cases[case_id], f"{case_id} case"),
            args.account,
            args.partition,
            args.walltime,
        )
        print(f"prepared {case_id}")
    plan = {
        "schema_version": SCHEMA_VERSION,
        "record_type": "stage-i-direct-fast-ct-exact-state-replay-plan",
        "prepared_utc": utc_now(),
        "primary_inventory": primary_binding,
        "cases": selected,
        "target_time": TARGET_TIME,
        "replay_tool": binding(SCRIPT_PATH, "CT replay tool"),
    }
    write_json(root / "plan.json", plan)
    print(root / "plan.json")
    return root


def manifest_paths(root: Path) -> list[Path]:
    plan, _ = load_json(root / "plan.json", "CT replay plan")
    return [
        root / "cases" / require_text(case_id, "replay case") / "t9/replay_run.json"
        for case_id in require_list(plan.get("cases"), "replay cases")
    ]


def submit(args: argparse.Namespace) -> None:
    root = args.output.expanduser().resolve(strict=True)
    for path in manifest_paths(root):
        manifest, _ = load_json(path, "replay run manifest")
        job = require_dict(manifest.get("job"), "replay job")
        if job.get("job_id") is not None:
            print(f"already submitted {manifest['case_id']}: {job['job_id']}")
            continue
        completed = subprocess.run(
            ["/usr/bin/sbatch", "--parsable", str(path.parent / "run.sbatch")],
            check=True,
            text=True,
            capture_output=True,
        )
        response = completed.stdout.strip()
        job_id = response.split(";", 1)[0]
        if JOB_ID.fullmatch(job_id) is None:
            raise CtReplayError(f"unexpected sbatch response: {response!r}")
        job["job_id"] = job_id
        job["submitted_utc"] = utc_now()
        manifest["job"] = job
        write_json(path, manifest, replace=True)
        print(f"submitted {manifest['case_id']}: {job_id}")


def job_observation(job_id: str) -> dict[str, object]:
    queued = subprocess.run(
        ["/usr/bin/squeue", "-h", "-j", job_id, "-o", "%T|%S"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if queued:
        state, start = (queued.split("|", 1) + [""])[:2]
        return {"state": state, "exit_code": None, "start": start, "end": None}
    accounted = subprocess.run(
        [
            "/usr/bin/sacct", "-X", "-n", "-P", "-j", job_id,
            "-o", "State,ExitCode,Start,End",
        ],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if not accounted:
        return {"state": "UNKNOWN", "exit_code": None, "start": None, "end": None}
    state, exit_code, start, end = (accounted.splitlines()[0].split("|") + [""] * 4)[:4]
    return {
        "state": state.split("+", 1)[0],
        "exit_code": exit_code,
        "start": start or None,
        "end": end or None,
    }


def terminal_restart_group(
    output: Path, rank_count: int, required_time: float
) -> dict[str, object]:
    rank_zero = output / "rst/rank_00000000"
    matches = [
        path for path in sorted(rank_zero.glob("*.rst"))
        if fast.restart_time(path) == required_time
    ]
    if len(matches) != 1:
        raise CtReplayError(
            f"replay must retain exactly one t={required_time:g} rank-zero restart"
        )
    path = matches[0]
    return {
        "time": required_time,
        "name": path.name,
        "rank_count": rank_count,
        "rank_files": rank_group(
            output / "rst", path.name, rank_count, "exact replay terminal restart"
        ),
    }


def replay_sanity(manifest: dict[str, object]) -> dict[str, object]:
    paths = require_dict(manifest.get("paths"), "replay paths")
    output = Path(require_text(paths.get("output_dir"), "replay output"))
    mhd_paths = sorted(output.glob("*.mhd.hst"))
    user_paths = sorted(output.glob("*.user.hst"))
    if len(mhd_paths) != 1 or len(user_paths) != 1:
        raise CtReplayError("replay must retain one MHD and one user history")
    mhd = fast.parse_history(mhd_paths[0])
    user = fast.parse_history(user_paths[0])
    if mhd.get("time") != user.get("time"):
        raise CtReplayError("replay histories are not synchronized")
    final_time = float(mhd["time"][-1])
    if final_time != TARGET_TIME:
        raise CtReplayError(f"replay final time is not exact t=9: {final_time:.17g}")
    strict = {
        name: max(abs(value) for value in mhd.get(name, [math.inf]))
        for name in fast.STRICT_FAILURE_COLUMNS
    }
    if any(value != 0.0 for value in strict.values()):
        raise CtReplayError("replay strict failure counters are nonzero")
    initial_mass = float(mhd["mass"][0])
    mass_drift = max(abs(float(value) - initial_mass) for value in mhd["mass"]) / max(
        abs(initial_mass), 1.0
    )
    mismatch = max(
        abs(float(left) - float(right))
        for left, right in zip(mhd["mass"], user["mass"])
    ) / max(abs(initial_mass), 1.0)
    if mass_drift > 1.0e-8 or mismatch > 1.0e-8:
        raise CtReplayError("replay mass sanity check failed")
    allocation = require_dict(manifest.get("allocation"), "replay allocation")
    terminal = terminal_restart_group(
        output, require_int(allocation.get("ranks"), "replay ranks", 1), TARGET_TIME
    )
    return {
        "final_time": final_time,
        "strict_lf_failure_maxima": strict,
        "mass_relative_drift": mass_drift,
        "mhd_user_mass_relative_mismatch": mismatch,
        "mhd_history": binding(mhd_paths[0], "replay MHD history"),
        "user_history": binding(user_paths[0], "replay user history"),
        "terminal_restart": terminal,
    }


def finalize(args: argparse.Namespace) -> Path:
    root = args.output.expanduser().resolve(strict=True)
    plan, plan_binding = load_json(root / "plan.json", "CT replay plan")
    primary = verify_binding(plan.get("primary_inventory"), "primary inventory")
    cases: dict[str, object] = {}
    for manifest_path in manifest_paths(root):
        manifest, manifest_binding = load_json(manifest_path, "replay run manifest")
        verify_binding(manifest.get("primary_inventory"), "run primary inventory")
        verify_binding(require_dict(manifest.get("tools"), "run tools").get("replay_tool"), "run replay tool")
        job = require_dict(manifest.get("job"), "replay job")
        job_id = require_text(job.get("job_id"), "replay job ID")
        if JOB_ID.fullmatch(job_id) is None:
            raise CtReplayError(f"invalid replay job ID: {job_id}")
        observation = job_observation(job_id)
        if observation["state"] in ACTIVE_STATES:
            raise CtReplayError(f"{manifest['case_id']} replay job is still active")
        if observation["state"] != "COMPLETED" or observation["exit_code"] != "0:0":
            raise CtReplayError(
                f"{manifest['case_id']} replay job did not complete cleanly: {observation}"
            )
        run_dir = manifest_path.parent
        exit_path = run_dir / "run_exit_code"
        if exit_path.read_text(encoding="utf-8").strip() != "0":
            raise CtReplayError(f"{manifest['case_id']} replay exit artifact is nonzero")
        sanity = replay_sanity(manifest)
        completion_body = {
            "schema_version": SCHEMA_VERSION,
            "record_type": COMPLETION_RECORD_TYPE,
            "case_id": manifest.get("case_id"),
            "case_name": manifest.get("case_name"),
            "purpose": manifest.get("purpose"),
            "target_time": TARGET_TIME,
            "command_line_overrides": ["time/tlim=9.0"],
            "primary_inventory": primary,
            "primary_case_lineage": verify_binding(
                manifest.get("primary_case_lineage"), "primary case lineage"
            ),
            "run_manifest": manifest_binding,
            "source_segment_manifest": verify_binding(
                require_dict(manifest.get("parent_restart"), "parent restart").get(
                    "source_segment_manifest"
                ),
                "source segment manifest",
            ),
            "execution_identity": {
                key: verify_binding(value, f"execution {key}")
                for key, value in require_dict(
                    manifest.get("execution_identity"), "execution identity"
                ).items()
            },
            "allocation": manifest.get("allocation"),
            "parent_restart": manifest.get("parent_restart"),
            "scheduler": {"job_id": job_id, **observation},
            "run_exit_code": binding(exit_path, "replay exit code"),
            "run_environment": binding(
                run_dir / "run_environment.txt", "replay environment"
            ),
            "scientific_sanity": sanity,
        }
        completion_path = run_dir / "completion.json"
        write_timestamped_json(completion_path, completion_body, "completed_utc")
        cases[str(manifest["case_id"])] = binding(
            completion_path, f"{manifest['case_id']} replay completion"
        )
        print(f"authenticated {manifest['case_id']}")
    inventory_body = {
        "schema_version": SCHEMA_VERSION,
        "record_type": INVENTORY_RECORD_TYPE,
        "purpose": "supplemental CT-only exact-t=9 evidence",
        "target_time": TARGET_TIME,
        "primary_inventory": primary,
        "plan": plan_binding,
        "replay_tool": binding(SCRIPT_PATH, "CT replay tool"),
        "cases": cases,
    }
    inventory_path = root / "inventory.json"
    write_timestamped_json(inventory_path, inventory_body, "completed_utc")
    print(inventory_path)
    return inventory_path


def status(args: argparse.Namespace) -> None:
    root = args.output.expanduser().resolve(strict=True)
    for path in manifest_paths(root):
        manifest, _ = load_json(path, "replay run manifest")
        job = require_dict(manifest.get("job"), "replay job")
        job_id = job.get("job_id")
        observed = job_observation(str(job_id)) if job_id is not None else {"state": "PREPARED"}
        print(f"{manifest['case_id']}\t{job_id}\t{observed['state']}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)
    prepare_parser = subcommands.add_parser("prepare")
    prepare_parser.add_argument("--inventory", type=Path, required=True)
    prepare_parser.add_argument("--inventory-sha256", required=True)
    prepare_parser.add_argument("--output", type=Path, required=True)
    prepare_parser.add_argument("--cases", action="append", default=[])
    prepare_parser.add_argument("--account", default="AST207")
    prepare_parser.add_argument("--partition", default="extended")
    prepare_parser.add_argument("--walltime", default="04:00:00")
    for name in ("submit", "status", "finalize"):
        command = subcommands.add_parser(name)
        command.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "prepare":
            prepare(args)
        elif args.command == "submit":
            submit(args)
        elif args.command == "status":
            status(args)
        elif args.command == "finalize":
            finalize(args)
        return 0
    except (CtReplayError, OSError, subprocess.SubprocessError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
