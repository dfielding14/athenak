#!/usr/bin/env python3
"""Prepare and account Frontier MKS24 Stage I production segments.

This utility is separate from ``cgl_lf_frontier.py`` because that launcher is
restricted to non-production ``debug`` QOS work.  Stage I jobs use Frontier's
``batch`` partition with its default production QOS, are prepared and
submitted one segment at a time, and must belong to the tracked mapped-case
manifest.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
from datetime import datetime, timedelta, timezone
import errno
import fcntl
from functools import wraps
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import pwd
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import uuid


ROOT_DIR = Path(__file__).resolve().parents[2]
PRODUCTION_UTILITY_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i.py")
DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
SACCT = Path("/usr/bin/sacct")
SBATCH = Path("/usr/bin/sbatch")
SCONTROL = Path("/usr/bin/scontrol")
SQUEUE = Path("/usr/bin/squeue")
DEFAULT_MATRIX = ROOT_DIR / "inputs/cgl_lf_paper/mks24_stage_i_manifest.json"
ACCOUNT = "AST207"
PARTITION = "batch"
PRODUCTION_QOS = "normal (Frontier default; no -q directive)"
PROJECT_BUDGET_NODE_HOURS = 4000.0
HISTORICAL_DEBUG_NODE_HOURS = 0.851670
HISTORICAL_E01_STAGE_I_NODE_HOURS = 9.962778
HISTORICAL_E02_PIPELINE_NODE_HOURS = 15.628610
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
AUTHORIZED_CASE_IDS = frozenset(f"R{number:02d}" for number in range(2, 18))
R17_CASE_ID = "R17"
R17_PREDECESSOR_CASE_IDS = tuple(
    f"R{number:02d}" for number in range(2, 17)
)
REQUIRED_CASE_FINAL_TIME = 10.0
COMPLETED_R16_NODE_HOURS = 6.145556
COMPLETED_R02_STANDARD_LAYOUT_PILOT_NODE_HOURS = 0.473333
COMPLETED_R17_HIGH_RESOLUTION_PILOT_NODE_HOURS = 4.235556
MEASURED_STAGE_I_RESERVED_NODE_HOURS = 900.0
CURRENT_STAGE_I_RESERVED_NODE_HOURS = MEASURED_STAGE_I_RESERVED_NODE_HOURS
MAX_SEGMENT_SECONDS = 2 * 60 * 60
MAX_RESTART_PARAMETER_DUMP_BYTES = 16 * 1024 * 1024
# Allow scheduler timestamp formatting and host-clock skew around the persisted
# pre-sbatch ambiguity barrier, but never an unrelated later submission.
SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS = 5 * 60
EXPECTED_RANKS_PER_NODE = 8
EXPECTED_CPUS_PER_TASK = 7
SEGMENT_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,28}")
JOB_ID_PATTERN = re.compile(r"[1-9][0-9]*")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
GIT_REVISION_PATTERN = re.compile(r"[0-9a-f]{40}")
BATCH_SCRIPT_DIGEST_PLACEHOLDER = "0" * 64
BATCH_SCRIPT_DIGEST_PATTERN = re.compile(
    r"(?m)^BATCH_SCRIPT_SHA256=([0-9a-f]{64})$"
)
LEDGER_NODE_HOUR_TOLERANCE = 5.0e-7 + 1.0e-12
LEDGER_CUMULATIVE_NODE_HOUR_TOLERANCE = 5.0e-7 + 1.0e-12
LEDGER_COLUMNS = (
    "execution_epoch",
    "job_id",
    "submitted_utc",
    "completed_utc",
    "case_id",
    "case_name",
    "segment",
    "state",
    "exit_code",
    "nodes",
    "requested_walltime",
    "elapsed_seconds",
    "reserved_node_hours",
    "actual_node_hours",
    "cumulative_stage_i_node_hours",
    "executable_revision",
    "executable_sha256",
    "input_revision",
    "input_file",
    "output_dir",
    "result",
    "notes",
)
RESERVATION_REQUIRED_COLUMNS = frozenset({
    "execution_epoch",
    "manifest",
    "case_id",
    "case_name",
    "segment",
    "nodes",
    "requested_walltime",
    "reserved_node_hours",
    "state",
    "prepared_utc",
})
RESERVATION_OPTIONAL_COLUMNS = frozenset({
    "execution_intent_sha256",
    "job_id",
    "actual_node_hours",
    "result",
    "notes",
})
TRANSACTION_KINDS = frozenset({
    "prepared", "submit_pending", "submitted", "submit_cleared",
    "recorded", "cancelled",
})
TRANSACTION_COMMON_COLUMNS = frozenset({
    "schema_version",
    "execution_epoch",
    "transaction_id",
    "kind",
    "created_utc",
    "manifest_path",
    "prior_reservations",
    "prior_reservations_sha256",
})
TRANSACTION_PAYLOAD_COLUMNS = frozenset({
    "manifest",
    "reservations",
    "ledger_row",
})
NONTERMINAL_STATES = {
    "PENDING",
    "RUNNING",
    "CONFIGURING",
    "COMPLETING",
    "SUSPENDED",
}
SHARED_ROOT_ACTIVE_STATES = {
    "prepared",
    "submitted",
    "pending",
    "running",
}
STRICT_LF_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
_ACTIVE_ROOT_LOCKS: dict[Path, tuple[object, int]] = {}


def utc_now() -> str:
    """Return a deterministic UTC timestamp."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def fsync_directory(path: Path) -> None:
    """Persist directory-entry updates after atomic replacement or removal."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def fsync_file(path: Path) -> None:
    """Persist one retained file after a copied or appended payload."""

    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def copy_file(source: Path, destination: Path) -> None:
    """Copy one retained artifact and persist its file and directory entries."""

    shutil.copy2(source, destination)
    fsync_file(destination)
    fsync_directory(destination.parent)


def mkdir_durable(path: Path) -> None:
    """Create a directory tree and persist each new directory entry."""

    missing = []
    current = path
    while not current.exists():
        missing.append(current)
        current = current.parent
    path.mkdir(parents=True, exist_ok=True)
    for directory in reversed(missing):
        fsync_directory(directory)
        fsync_directory(directory.parent)


def write_text(path: Path, value: str, mode: int | None = None) -> None:
    """Atomically and durably write retained text metadata."""

    temporary = path.with_name(
        f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp"
    )
    descriptor = os.open(
        temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o666
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            if mode is not None:
                os.fchmod(stream.fileno(), mode)
            elif path.exists():
                os.fchmod(stream.fileno(), stat.S_IMODE(path.stat().st_mode))
            stream.write(value)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        fsync_directory(path.parent)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def write_json(path: Path, value: object) -> None:
    """Atomically and durably write stable JSON metadata."""

    write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def stable_json_sha256(value: object) -> str:
    """Return the checksum produced by ``write_json`` for one value."""

    return hashlib.sha256(
        (json.dumps(value, indent=2, sort_keys=True) + "\n").encode("utf-8")
    ).hexdigest()


def unlink_durable(path: Path) -> None:
    """Remove one retained entry and persist the directory update."""

    path.unlink()
    fsync_directory(path.parent)


def append_ledger_row(path: Path, row: dict[str, object],
                      prior_rows: list[dict[str, str]]) -> None:
    """Durably append one allocation row."""

    if frozenset(row) != frozenset(LEDGER_COLUMNS):
        raise ValueError("ledger append row has invalid columns")
    prior_actual = sum(float(item["actual_node_hours"]) for item in prior_rows)
    validate_ledger_cumulative_fields(row, prior_actual, "ledger append row")
    with path.open("a", newline="", encoding="utf-8") as stream:
        csv.DictWriter(stream, fieldnames=LEDGER_COLUMNS).writerow(row)
        stream.flush()
        os.fsync(stream.fileno())


def canonical_root_lock_path(root: Path) -> Path:
    """Return the cooperative Stage I mutation lock location."""

    return root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock"


def require_canonical_root_lock_profile(profile: os.stat_result,
                                        lock_path: Path) -> None:
    """Require one cooperative canonical-root lock profile."""

    mode = stat.S_IMODE(profile.st_mode)
    if not stat.S_ISREG(profile.st_mode):
        raise ValueError(f"Stage I lock is not a regular file: {lock_path}")
    if profile.st_uid != os.geteuid():
        raise ValueError(f"Stage I lock owner differs: {lock_path}")
    if profile.st_nlink != 1:
        raise ValueError(f"Stage I lock link count differs: {lock_path}")
    if mode & ~0o644:
        raise ValueError(f"Stage I lock mode is too permissive: {lock_path}")


@contextmanager
def canonical_root_lock(root: Path):
    """Take a nonblocking reentrant lock for canonical-root mutations."""

    resolved = root.expanduser().resolve()
    if resolved != DEFAULT_ROOT.expanduser().resolve():
        yield
        return
    active = _ACTIVE_ROOT_LOCKS.get(resolved)
    if active is not None:
        stream, depth = active
        _ACTIVE_ROOT_LOCKS[resolved] = (stream, depth + 1)
        try:
            yield
        finally:
            _ACTIVE_ROOT_LOCKS[resolved] = (stream, depth)
        return
    if not resolved.is_dir():
        raise ValueError(f"canonical Stage I root is unavailable: {resolved}")
    lock_path = canonical_root_lock_path(resolved)
    created = False
    previous = os.umask(0)
    try:
        try:
            descriptor = os.open(
                lock_path,
                os.O_RDWR | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW,
                0o644,
            )
            created = True
        except FileExistsError:
            try:
                descriptor = os.open(lock_path, os.O_RDWR | os.O_NOFOLLOW)
            except OSError as error:
                if error.errno == errno.ELOOP:
                    raise ValueError(
                        f"Stage I lock must not be a symlink: {lock_path}"
                    ) from error
                raise
    finally:
        os.umask(previous)
    try:
        profile = os.fstat(descriptor)
        require_canonical_root_lock_profile(profile, lock_path)
    except BaseException:
        os.close(descriptor)
        raise
    if created:
        fsync_directory(resolved)
    stream = os.fdopen(descriptor, "a+", encoding="utf-8")
    try:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            raise ValueError(
                f"another Stage I mutation holds {lock_path}"
            ) from error
        profile = os.fstat(stream.fileno())
        try:
            named_profile = os.stat(lock_path, follow_symlinks=False)
        except FileNotFoundError as error:
            raise ValueError(
                f"Stage I lock path changed while locking: {lock_path}"
            ) from error
        if (
            not stat.S_ISREG(named_profile.st_mode)
            or (named_profile.st_dev, named_profile.st_ino)
            != (profile.st_dev, profile.st_ino)
        ):
            raise ValueError(
                f"Stage I lock path changed while locking: {lock_path}"
            )
        require_canonical_root_lock_profile(profile, lock_path)
        if stat.S_IMODE(os.fstat(stream.fileno()).st_mode) != 0o644:
            os.fchmod(stream.fileno(), 0o644)
            os.fsync(stream.fileno())
            fsync_directory(resolved)
        _ACTIVE_ROOT_LOCKS[resolved] = (stream, 1)
        try:
            yield
        finally:
            del _ACTIVE_ROOT_LOCKS[resolved]
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    finally:
        stream.close()


def locked_root_action(function):
    """Lock a mutating action whose root is supplied directly."""

    @wraps(function)
    def wrapped(args: argparse.Namespace):
        root = require_root(Path(args.root), args.allow_local_root)
        with canonical_root_lock(root):
            return function(args)

    return wrapped


def locked_manifest_action(function):
    """Lock a mutating action whose root is retained in its manifest."""

    @wraps(function)
    def wrapped(args: argparse.Namespace):
        manifest = read_manifest(Path(args.manifest).expanduser().resolve())
        root = require_root(
            Path(str(manifest["project_root"])),
            getattr(args, "allow_local_root", False),
        )
        with canonical_root_lock(root):
            retained = read_manifest(Path(args.manifest).expanduser().resolve())
            if Path(str(retained.get("project_root", ""))).resolve() != root:
                raise ValueError("manifest project root changed while acquiring its lock")
            return function(args)

    return wrapped


def sha256(path: Path) -> str:
    """Return a file SHA-256 digest."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_walltime(value: str) -> int:
    """Parse an HH:MM:SS allocation walltime."""

    match = re.fullmatch(r"(\d{2}):([0-5]\d):([0-5]\d)", value)
    if match is None:
        raise ValueError(f"walltime must use HH:MM:SS: {value}")
    hours, minutes, seconds = (int(item) for item in match.groups())
    return hours * 3600 + minutes * 60 + seconds


def node_hours(nodes: int, seconds: int) -> float:
    """Compute allocated node hours."""

    return nodes * seconds / 3600.0


def validate_ledger_numeric_fields(row: dict[str, object],
                                   label: str) -> tuple[float, float, float]:
    """Require finite, non-negative retained Stage I accounting values."""

    try:
        nodes_text = str(row["nodes"])
        elapsed_text = str(row["elapsed_seconds"])
        if (
            re.fullmatch(r"[0-9]+", nodes_text) is None
            or re.fullmatch(r"[0-9]+", elapsed_text) is None
        ):
            raise ValueError
        nodes = int(nodes_text)
        elapsed_seconds = int(elapsed_text)
        requested_seconds = parse_walltime(str(row["requested_walltime"]))
        reserved_node_hours = float(row["reserved_node_hours"])
        actual_node_hours = float(row["actual_node_hours"])
        cumulative_node_hours = float(row["cumulative_stage_i_node_hours"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{label} has invalid numeric fields") from error
    if (
        nodes < 1
        or elapsed_seconds < 0
        or not math.isfinite(reserved_node_hours)
        or reserved_node_hours <= 0.0
        or not math.isfinite(actual_node_hours)
        or actual_node_hours < 0.0
        or not math.isfinite(cumulative_node_hours)
        or cumulative_node_hours < 0.0
        or abs(
            reserved_node_hours - node_hours(nodes, requested_seconds)
        ) > LEDGER_NODE_HOUR_TOLERANCE
        or abs(actual_node_hours - node_hours(nodes, elapsed_seconds))
        > LEDGER_NODE_HOUR_TOLERANCE
    ):
        raise ValueError(f"{label} has invalid numeric fields")
    return (
        actual_node_hours,
        cumulative_node_hours,
        node_hours(nodes, elapsed_seconds),
    )


def validate_ledger_cumulative_fields(row: dict[str, object],
                                      prior_actual_node_hours: float,
                                      label: str) -> float:
    """Require one retained cumulative total to follow its preceding rows."""

    actual, cumulative, derived_actual = validate_ledger_numeric_fields(row, label)
    expected = prior_actual_node_hours + derived_actual
    if (
        not math.isfinite(expected)
        or abs(cumulative - expected) > LEDGER_CUMULATIVE_NODE_HOUR_TOLERANCE
    ):
        raise ValueError(f"{label} has inconsistent cumulative node-hours")
    if (
        expected > CURRENT_STAGE_I_RESERVED_NODE_HOURS
        or expected > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError(f"{label} exceeds an accounting ceiling")
    return prior_actual_node_hours + actual


def require_safe_segment(value: str) -> str:
    """Require a path- and Slurm-safe retained segment identifier."""

    if value == "analysis" or SEGMENT_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "--segment must contain 1-29 ASCII letters, digits, underscores, "
            "or hyphens, must begin with a letter or digit, and must not use "
            "the reserved analysis namespace"
        )
    return value


def require_numeric_job_id(value: str) -> str:
    """Require a top-level numeric Slurm allocation identifier."""

    if JOB_ID_PATTERN.fullmatch(value) is None:
        raise ValueError("--job-id must be a positive numeric Slurm job ID")
    return value


def require_root(root: Path, allow_local_root: bool) -> Path:
    """Require the declared project run root outside offline validation."""

    resolved = root.expanduser().resolve()
    if resolved != DEFAULT_ROOT.expanduser().resolve() and not allow_local_root:
        raise ValueError(
            f"Stage I root must be {DEFAULT_ROOT}; "
            "use --allow-local-root only for offline validation"
        )
    return resolved


def is_offline_local_root(root: Path, allow_local_root: bool) -> bool:
    """Return whether relaxed validation is permitted for a nonproject fixture."""

    return (
        allow_local_root
        and root.expanduser().resolve() != DEFAULT_ROOT.expanduser().resolve()
    )


def require_beneath_root(path: Path, root: Path, label: str,
                         allow_local_root: bool) -> None:
    """Keep submitted products in the project filesystem."""

    if is_offline_local_root(root, allow_local_root):
        return
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} must be beneath {root}: {path}") from error


def require_current_epoch(manifest: dict[str, object], label: str) -> None:
    """Reject archival or untagged manifests in the current production path."""

    if manifest.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(
            f"{label} execution epoch is not {EXECUTION_EPOCH}: "
            f"{manifest.get('execution_epoch')!r}"
        )


def require_authorized_case(case_id: str) -> None:
    """Limit execution to the frozen mapped Stage I matrix."""

    if case_id not in AUTHORIZED_CASE_IDS:
        raise ValueError(
            f"{EXECUTION_EPOCH} Stage I is authorized only for mapped matrix cases "
            "R02-R17 under sequential inspection"
        )


def require_case_node_count(case_id: str, nodes: int) -> None:
    """Bind each mapped case to its reviewed Frontier allocation shape."""

    require_authorized_case(case_id)
    expected = 8 if case_id == R17_CASE_ID else 1
    if nodes != expected:
        raise ValueError(
            f"{case_id} canonical Stage I preparation requires --nodes={expected}"
        )


def retained_case_has_started(paths: dict[str, Path], case_id: str) -> bool:
    """Return whether one mapped case has retained any segment manifest."""

    require_authorized_case(case_id)
    return any(
        (paths["runs"] / case_id).glob("*/manifest/prepared_run.json")
    )


def require_r17_last(paths: dict[str, Path], case_id: str) -> None:
    """Require every lower-cost mapped case to complete before R17 starts."""

    if case_id != R17_CASE_ID:
        if retained_case_has_started(paths, R17_CASE_ID):
            raise ValueError(
                f"{R17_CASE_ID} has started; later {case_id} preparation is forbidden"
            )
        return
    incomplete = []
    for predecessor in R17_PREDECESSOR_CASE_IDS:
        lineage = accepted_case_lineage(paths, predecessor)
        if not lineage:
            incomplete.append(predecessor)
            continue
        final_time = float(lineage[-1]["scientific_inspection"]["final_time"])
        if (
            not math.isfinite(final_time)
            or final_time < REQUIRED_CASE_FINAL_TIME - 1.0e-10
        ):
            incomplete.append(predecessor)
    if incomplete:
        raise ValueError(
            "R17 must remain last; accepted t=10 lineages are incomplete for: "
            + ", ".join(incomplete)
        )


def require_prepare_case_policy(paths: dict[str, Path], case_id: str,
                                nodes: int, offline_local_root: bool) -> None:
    """Apply mapped allocation and ordering gates only to canonical production."""

    if offline_local_root:
        return
    require_case_node_count(case_id, nodes)
    require_r17_last(paths, case_id)


def parse_utc_timestamp(value: object, label: str) -> datetime:
    """Parse one retained timestamp and normalize it to UTC."""

    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} is not an ISO-8601 timestamp: {value!r}") from error
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=datetime.now().astimezone().tzinfo)
    return parsed.astimezone(timezone.utc)


def require_positive_finite_float(value: object, label: str) -> float:
    """Require one positive finite command or API threshold."""

    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError(f"{label} must be positive and finite")
    return result


def positive_finite_float_arg(value: str) -> float:
    """Parse a positive finite command-line floating-point value."""

    try:
        return require_positive_finite_float(value, "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def layout(root: Path) -> dict[str, Path]:
    """Return retained current-epoch Stage I production locations."""

    accounting = root / "accounting"
    return {
        "root": root,
        "accounting": accounting,
        "ledger": accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_node_hours.csv",
        "reservations": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_reservations.json"
        ),
        "summary": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_budget_summary.md"
        ),
        "transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions"
        ),
        "qualification": (
            accounting
            / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_qualification_approval.json"
        ),
        "runs": root / "runs" / "mks24-stage-i" / EXECUTION_EPOCH,
        "logs_slurm": root / "logs" / "slurm",
    }


def require_safe_manifest_path(paths: dict[str, Path],
                               value: object) -> Path:
    """Require the exact E03 manifest location for one authorized segment."""

    declared = Path(str(value)).expanduser()
    if not declared.is_absolute():
        raise ValueError(f"manifest target must be absolute: {declared}")
    manifest_path = declared.resolve()
    try:
        relative = manifest_path.relative_to(paths["runs"].resolve())
    except ValueError as error:
        raise ValueError(
            f"manifest target is outside the E03 run store: {manifest_path}"
        ) from error
    if len(relative.parts) != 4 or relative.parts[2:] != (
        "manifest", "prepared_run.json"
    ):
        raise ValueError(f"manifest target is not an exact E03 segment path: {manifest_path}")
    case_id, segment = relative.parts[:2]
    require_authorized_case(case_id)
    require_safe_segment(segment)
    expected = paths["runs"] / case_id / segment / "manifest" / "prepared_run.json"
    if declared != expected or manifest_path != expected.resolve():
        raise ValueError(f"manifest target does not resolve exactly beneath E03: {manifest_path}")
    return manifest_path


def orphaned_segment_run_directories(paths: dict[str, Path]) -> list[Path]:
    """Return segment directories left behind before a manifest was retained."""

    orphans = []
    if not paths["runs"].is_dir():
        return orphans
    for case_dir in sorted(paths["runs"].iterdir()):
        if not case_dir.is_dir() or case_dir.name not in AUTHORIZED_CASE_IDS:
            continue
        for run_dir in sorted(case_dir.iterdir()):
            if run_dir.name == "analysis":
                continue
            if (
                run_dir.is_dir()
                and not (run_dir / "manifest" / "prepared_run.json").is_file()
            ):
                orphans.append(run_dir)
    return orphans


def require_no_orphaned_segment_runs(paths: dict[str, Path]) -> None:
    """Fail closed after an interrupted prepare leaves unaudited run content."""

    orphans = orphaned_segment_run_directories(paths)
    if orphans:
        raise ValueError(
            "Stage I interrupted-prepare cleanup is required before mutation: "
            + ", ".join(str(path) for path in orphans)
        )


def read_qualification_approval(path: Path) -> dict[str, object]:
    """Read and validate one corrected-build Frontier qualification token."""

    if not path.is_file():
        raise ValueError(f"E03 qualification approval token is absent: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"E03 qualification approval token is invalid: {path}")
    if value.get("schema_version") != 1:
        raise ValueError(f"E03 qualification approval token has wrong schema: {path}")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(f"E03 qualification approval token has wrong epoch: {path}")
    executable_sha256 = value.get("approved_executable_sha256")
    revision = value.get("approved_executable_revision")
    if (
        not isinstance(executable_sha256, str)
        or SHA256_PATTERN.fullmatch(executable_sha256) is None
    ):
        raise ValueError(
            f"E03 qualification approval token has invalid executable digest: {path}"
        )
    if (
        not isinstance(revision, str)
        or GIT_REVISION_PATTERN.fullmatch(revision) is None
    ):
        raise ValueError(
            f"E03 qualification approval token has invalid git revision: {path}"
        )
    for key in ("approved_utc", "approved_by", "review_notes"):
        if not isinstance(value.get(key), str) or not str(value[key]).strip():
            raise ValueError(
                f"E03 qualification approval token lacks nonempty {key}: {path}"
            )
    return value


def qualification_approval_status(paths: dict[str, Path]) -> dict[str, object]:
    """Return a summary-safe view of the current E03 qualification token."""

    path = paths["qualification"]
    if not path.is_file():
        return {
            "state": "pending",
            "path": str(path),
            "reason": "approval token is absent",
        }
    try:
        approval = read_qualification_approval(path)
    except (OSError, ValueError, json.JSONDecodeError) as error:
        return {
            "state": "invalid",
            "path": str(path),
            "reason": str(error),
        }
    return {
        "state": "approved",
        "path": str(path),
        "sha256": sha256(path),
        "approved_executable_sha256": approval["approved_executable_sha256"],
        "approved_executable_revision": approval["approved_executable_revision"],
        "approved_utc": approval["approved_utc"],
        "approved_by": approval["approved_by"],
    }


def require_qualification_approval(
    paths: dict[str, Path],
    executable_sha256: str,
    executable_revision: str,
    offline_local_root: bool,
) -> dict[str, object] | None:
    """Require the fixed-location token to approve this exact corrected build."""

    path = paths["qualification"]
    if offline_local_root and not path.is_file():
        return None
    approval = read_qualification_approval(path)
    if approval["approved_executable_sha256"] != executable_sha256:
        raise ValueError("E03 qualification token does not approve this executable")
    if approval["approved_executable_revision"] != executable_revision:
        raise ValueError("E03 qualification token does not approve this git revision")
    return {
        "path": str(path),
        "sha256": sha256(path),
        "execution_epoch": EXECUTION_EPOCH,
        "approved_executable_sha256": approval["approved_executable_sha256"],
        "approved_executable_revision": approval["approved_executable_revision"],
        "token": approval,
    }


def initialize(root: Path) -> dict[str, Path]:
    """Create the production layout and accounting stores."""

    with canonical_root_lock(root):
        paths = layout(root)
        for key in ("accounting", "runs", "logs_slurm", "transactions"):
            mkdir_durable(paths[key])
        if not paths["ledger"].exists():
            with paths["ledger"].open("w", newline="", encoding="utf-8") as stream:
                csv.writer(stream).writerow(LEDGER_COLUMNS)
                stream.flush()
                os.fsync(stream.fileno())
            fsync_directory(paths["ledger"].parent)
        if not paths["reservations"].exists():
            write_json(paths["reservations"], [])
        read_ledger(paths)
        read_reservations(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        refresh_summary(paths)
        return paths


def read_ledger(paths: dict[str, Path]) -> list[dict[str, str]]:
    """Read retained Stage I allocation records."""

    with paths["ledger"].open(newline="", encoding="utf-8") as stream:
        reader = csv.reader(stream)
        try:
            header = next(reader)
        except StopIteration as error:
            raise ValueError(f"Stage I ledger is empty: {paths['ledger']}") from error
        if tuple(header) != LEDGER_COLUMNS:
            raise ValueError(
                f"Stage I ledger header is invalid: {paths['ledger']}: {header!r}"
            )
        rows = []
        cumulative = 0.0
        for index, row in enumerate(reader, start=2):
            if len(row) != len(LEDGER_COLUMNS):
                raise ValueError(
                    f"Stage I ledger row {index} has {len(row)} columns; "
                    f"expected {len(LEDGER_COLUMNS)}"
                )
            retained = dict(zip(LEDGER_COLUMNS, row))
            cumulative = validate_ledger_cumulative_fields(
                retained, cumulative, f"Stage I ledger row {index}"
            )
            rows.append(retained)
        return rows


def read_reservations(paths: dict[str, Path]) -> list[dict[str, object]]:
    """Read all Stage I preparation records."""

    value = json.loads(paths["reservations"].read_text(encoding="utf-8"))
    if not isinstance(value, list):
        raise ValueError("Stage I reservation store must contain a list")
    return value


def pending_transaction_paths(paths: dict[str, Path]) -> list[Path]:
    """Return durable metadata transitions awaiting completion."""

    directory = paths["transactions"]
    if not directory.is_dir():
        return []
    return sorted(directory.glob("*.json"))


def require_no_pending_transactions(paths: dict[str, Path]) -> None:
    """Fail closed while a prior cross-file transition awaits recovery."""

    pending = pending_transaction_paths(paths)
    if pending:
        raise ValueError(
            "Stage I metadata recovery is required before mutation: "
            + ", ".join(str(path) for path in pending)
        )


def validate_reservation_record(paths: dict[str, Path],
                                reservation: object) -> dict[str, object]:
    """Require one journal-controlled reservation to use the retained schema."""

    if not isinstance(reservation, dict):
        raise ValueError("transaction reservation record is not an object")
    columns = frozenset(reservation)
    if not RESERVATION_REQUIRED_COLUMNS.issubset(columns):
        missing = sorted(RESERVATION_REQUIRED_COLUMNS - columns)
        raise ValueError(f"transaction reservation lacks columns: {missing}")
    supported = RESERVATION_REQUIRED_COLUMNS | RESERVATION_OPTIONAL_COLUMNS
    if not columns.issubset(supported):
        raise ValueError(
            "transaction reservation has unsupported columns: "
            f"{sorted(columns - supported)}"
        )
    if reservation.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("transaction reservation has wrong execution epoch")
    manifest_path = require_safe_manifest_path(paths, reservation["manifest"])
    case_id = str(reservation["case_id"])
    segment = require_safe_segment(str(reservation["segment"]))
    if (
        case_id != manifest_path.parents[2].name
        or segment != manifest_path.parents[1].name
    ):
        raise ValueError("transaction reservation identity differs from manifest path")
    require_authorized_case(case_id)
    try:
        nodes = int(reservation["nodes"])
        requested_seconds = parse_walltime(str(reservation["requested_walltime"]))
        reserved = float(reservation["reserved_node_hours"])
    except (TypeError, ValueError) as error:
        raise ValueError("transaction reservation allocation is invalid") from error
    if (
        nodes < 1
        or requested_seconds <= 0
        or not math.isfinite(reserved)
        or reserved <= 0.0
    ):
        raise ValueError("transaction reservation allocation must be positive")
    if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
        require_case_node_count(case_id, nodes)
    if abs(reserved - node_hours(nodes, requested_seconds)) > 5.0e-12:
        raise ValueError("transaction reservation node-hours differ from allocation")
    parse_utc_timestamp(reservation["prepared_utc"], "reservation prepared_utc")
    state = reservation.get("state")
    if state not in {"prepared", "submitted", "recorded", "cancelled"}:
        raise ValueError(f"transaction reservation has invalid state: {state!r}")
    if state in {"submitted", "recorded"}:
        require_numeric_job_id(str(reservation.get("job_id", "")))
    if state == "recorded":
        try:
            actual = float(reservation["actual_node_hours"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("recorded transaction reservation lacks actual use") from error
        if not math.isfinite(actual) or actual < 0.0:
            raise ValueError("recorded transaction reservation has invalid actual use")
        if reservation.get("result") not in {
            "accepted", "clean_partial", "rejected", "failed", "aborted"
        }:
            raise ValueError("recorded transaction reservation has invalid result")
    if state == "cancelled" and not isinstance(reservation.get("notes"), str):
        raise ValueError("cancelled transaction reservation lacks notes")
    intent = reservation.get("execution_intent_sha256")
    if (
        paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve()
        and intent is None
    ):
        raise ValueError(
            "canonical transaction reservation lacks execution intent digest"
        )
    if intent is not None and (
        not isinstance(intent, str) or SHA256_PATTERN.fullmatch(intent) is None
    ):
        raise ValueError("transaction reservation execution intent is invalid")
    return reservation


def validate_transaction_ledger_row(row: object,
                                    manifest: dict[str, object]) -> None:
    """Require one journal-controlled ledger row to match its target manifest."""

    if not isinstance(row, dict) or frozenset(row) != frozenset(LEDGER_COLUMNS):
        raise ValueError("recorded transaction ledger row has invalid columns")
    validate_ledger_numeric_fields(row, "recorded transaction ledger row")
    require_numeric_job_id(str(row["job_id"]))
    if row["execution_epoch"] != EXECUTION_EPOCH:
        raise ValueError("recorded transaction ledger row has wrong execution epoch")
    run = manifest.get("run")
    command = manifest.get("command")
    manifest_paths = manifest.get("paths")
    allocation = manifest.get("allocation")
    if (
        not isinstance(run, dict)
        or not isinstance(command, dict)
        or not isinstance(manifest_paths, dict)
        or not isinstance(allocation, dict)
    ):
        raise ValueError("recorded transaction manifest lacks ledger provenance")
    expected = {
        "job_id": manifest.get("job_id"),
        "case_id": run.get("case_id"),
        "case_name": run.get("case_name"),
        "segment": run.get("segment"),
        "nodes": str(allocation.get("nodes")),
        "requested_walltime": allocation.get("requested_walltime"),
        "executable_revision": command.get("executable_revision"),
        "executable_sha256": command.get("executable_sha256"),
        "input_revision": command.get("input_revision"),
        "input_file": command.get("input_file"),
        "output_dir": manifest_paths.get("output_dir"),
    }
    for key, value in expected.items():
        if row.get(key) != value:
            raise ValueError(f"recorded transaction ledger {key} differs from manifest")
    if manifest.get("accounting") != row:
        raise ValueError("recorded transaction manifest accounting differs from ledger")


def validate_submission_audit(paths: dict[str, Path], value: object) -> None:
    """Require fixed machine-readable submit-policy evidence."""

    if not isinstance(value, dict):
        raise ValueError("submission journal lacks submission audit")
    required = {
        "created_utc",
        "offline_local_root",
        "skip_slurm_test",
        "slurm_test_only",
        "acknowledged_shared_root_campaigns",
    }
    optional = {"legacy_mark_submitted"}
    if not required.issubset(value) or not frozenset(value).issubset(required | optional):
        raise ValueError("submission journal has invalid submission audit columns")
    parse_utc_timestamp(value["created_utc"], "submission audit created_utc")
    if (
        not isinstance(value["offline_local_root"], bool)
        or not isinstance(value["skip_slurm_test"], bool)
        or not isinstance(value["slurm_test_only"], str)
        or not isinstance(value["acknowledged_shared_root_campaigns"], list)
        or not all(
            isinstance(item, str)
            for item in value["acknowledged_shared_root_campaigns"]
        )
        or (
            "legacy_mark_submitted" in value
            and value["legacy_mark_submitted"] is not True
        )
    ):
        raise ValueError("submission journal has invalid submission audit values")
    offline_local_root = (
        paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve()
    )
    if value["offline_local_root"] is not offline_local_root:
        raise ValueError("submission journal root mode differs from submission audit")


def scheduler_output_contains_job(output: str, expected_job_name: str) -> bool:
    """Return whether retained scheduler rows name the prepared allocation."""

    for row in csv.reader(output.splitlines(), delimiter="|"):
        if not row:
            continue
        if len(row) < 2:
            raise ValueError("retained scheduler evidence contains a malformed row")
        if row[1] == expected_job_name:
            return True
    return False


def validate_scheduler_absence_evidence(paths: dict[str, Path],
                                        value: object,
                                        transaction: dict[str, object]) -> None:
    """Require retained scheduler-clear evidence to use one fixed schema."""

    if not isinstance(value, dict):
        raise ValueError("cleared submission transaction lacks scheduler evidence")
    common = {
        "mode", "checked_utc", "expected_job_name", "ambiguity_created_utc",
    }
    mode = value.get("mode")
    optional = {
        "live scheduler absence query": {
            "squeue_command", "squeue_output", "sacct_command", "sacct_output",
        },
        "offline-local fixture": {"fixture"},
        "offline-local operator confirmation": set(),
        "break-glass after scheduler query failure": {
            "operator_evidence", "query_error",
        },
    }
    if mode not in optional or frozenset(value) != frozenset(common | optional[mode]):
        raise ValueError("cleared submission scheduler evidence has invalid schema")
    parse_utc_timestamp(value["checked_utc"], "scheduler absence checked_utc")
    manifest = transaction.get("manifest")
    if not isinstance(manifest, dict):
        raise ValueError("cleared submission transaction lacks manifest evidence")
    expected_name = expected_job_name(manifest)
    if (
        value["ambiguity_created_utc"] != transaction.get("created_utc")
        or value["expected_job_name"] != expected_name
    ):
        raise ValueError("cleared submission scheduler evidence is inconsistent")
    offline_local_root = (
        paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve()
    )
    if mode in {
        "offline-local fixture", "offline-local operator confirmation",
    } and not offline_local_root:
        raise ValueError(
            "canonical cleared submission may not retain offline scheduler evidence"
        )
    if mode == "live scheduler absence query" and (
        not isinstance(value["squeue_command"], list)
        or not isinstance(value["sacct_command"], list)
        or not isinstance(value["squeue_output"], str)
        or not isinstance(value["sacct_output"], str)
    ):
        raise ValueError("cleared submission live scheduler evidence is invalid")
    if mode == "live scheduler absence query" and (
        scheduler_output_contains_job(value["squeue_output"], expected_name)
        or scheduler_output_contains_job(value["sacct_output"], expected_name)
    ):
        raise ValueError(
            "retained scheduler evidence still reports the prepared allocation"
        )
    if mode == "offline-local fixture" and (
        not isinstance(value["fixture"], dict)
        or value["fixture"].get("absent") is not True
    ):
        raise ValueError("cleared submission offline fixture is invalid")
    if mode == "break-glass after scheduler query failure" and (
        not isinstance(value["operator_evidence"], str)
        or not value["operator_evidence"].strip()
        or not isinstance(value["query_error"], str)
        or not value["query_error"].strip()
    ):
        raise ValueError("cleared submission break-glass evidence is invalid")


def transaction_expected_columns(kind: str) -> frozenset[str]:
    """Return the fixed journal schema for one transition kind."""

    if kind == "submit_pending":
        return TRANSACTION_COMMON_COLUMNS | frozenset({
            "prepared_manifest_sha256", "submission_audit",
        })
    columns = TRANSACTION_COMMON_COLUMNS | TRANSACTION_PAYLOAD_COLUMNS
    if kind == "submitted":
        return columns | frozenset({
            "prepared_manifest_sha256", "submission_audit",
            "job_id", "submitted_recorded_utc",
        })
    if kind == "submit_cleared":
        return columns | frozenset({
            "prepared_manifest_sha256", "submission_audit",
            "recovery_notes", "scheduler_absence_evidence",
        })
    return columns


def read_transaction(paths: dict[str, Path], path: Path) -> dict[str, object]:
    """Read and fully validate one durable Stage I transition journal."""

    resolved = path.expanduser().resolve()
    if resolved.parent != paths["transactions"].resolve() or resolved.suffix != ".json":
        raise ValueError(f"transaction journal is outside the E03 store: {resolved}")
    value = json.loads(resolved.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid Stage I transaction journal: {resolved}")
    kind = value.get("kind")
    if kind not in TRANSACTION_KINDS:
        raise ValueError(f"transaction journal has invalid kind: {resolved}")
    if frozenset(value) != transaction_expected_columns(str(kind)):
        raise ValueError(f"transaction journal has invalid schema for {kind}: {resolved}")
    if value.get("schema_version") != 1:
        raise ValueError(f"transaction journal has wrong schema version: {resolved}")
    if value.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError(f"transaction journal has wrong execution epoch: {resolved}")
    if value.get("transaction_id") != resolved.stem:
        raise ValueError(f"transaction journal ID differs from filename: {resolved}")
    parse_utc_timestamp(value.get("created_utc"), "transaction created_utc")
    prior_reservations_sha256 = value.get("prior_reservations_sha256")
    if (
        not isinstance(prior_reservations_sha256, str)
        or SHA256_PATTERN.fullmatch(prior_reservations_sha256) is None
    ):
        raise ValueError(f"transaction journal has invalid reservation baseline: {resolved}")
    prior_reservations = value.get("prior_reservations")
    if (
        not isinstance(prior_reservations, list)
        or stable_json_sha256(prior_reservations) != prior_reservations_sha256
    ):
        raise ValueError(f"transaction journal reservation baseline is invalid: {resolved}")
    reservation_records_by_manifest(
        paths, prior_reservations, "transaction prior reservation snapshot"
    )
    manifest_path = require_safe_manifest_path(paths, value["manifest_path"])
    if kind == "submit_pending":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("pending submission journal has invalid manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        authenticate_canonical_reservation_snapshot(paths, prior_reservations)
        return value
    manifest = value.get("manifest")
    reservations = value.get("reservations")
    if not isinstance(manifest, dict) or not isinstance(reservations, list):
        raise ValueError(f"transaction payload is incomplete: {resolved}")
    require_current_epoch(manifest, "transaction manifest")
    if Path(str(manifest.get("project_root", ""))).resolve() != paths["root"].resolve():
        raise ValueError("transaction manifest project root differs from E03 root")
    run = manifest.get("run")
    if (
        not isinstance(run, dict)
        or run.get("case_id") != manifest_path.parents[2].name
        or run.get("segment") != manifest_path.parents[1].name
    ):
        raise ValueError("transaction manifest identity differs from target path")
    validated = [
        validate_reservation_record(paths, reservation)
        for reservation in reservations
    ]
    validate_reservation_snapshot_transition(
        paths, str(kind), manifest_path, prior_reservations, validated
    )
    matches = [
        reservation for reservation in validated
        if Path(str(reservation["manifest"])).resolve() == manifest_path
    ]
    if len(matches) != 1 or matches[0].get("state") != manifest.get("state"):
        raise ValueError("transaction target reservation state differs from manifest")
    expected_state = {
        "prepared": "prepared",
        "submitted": "submitted",
        "submit_cleared": "prepared",
        "recorded": "recorded",
        "cancelled": "cancelled",
    }[str(kind)]
    if manifest.get("state") != expected_state:
        raise ValueError(f"transaction manifest state is invalid for {kind}")
    if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
        validate_prepared_resources(manifest, canonical_production=True)
        require_reservation_matches_manifest(matches[0], manifest)
        require_reserved_execution_intent(
            matches[0], manifest, allow_legacy_local=False
        )
        if kind == "prepared":
            require_prepare_case_policy(
                paths, str(run["case_id"]), int(matches[0]["nodes"]),
                offline_local_root=False,
            )
    row = value.get("ledger_row")
    if kind == "recorded":
        validate_transaction_ledger_row(row, manifest)
        if paths["root"].resolve() == DEFAULT_ROOT.expanduser().resolve():
            require_recorded_scheduler_evidence(paths, row, manifest)
        difference = abs(
            float(matches[0].get("actual_node_hours", -1.0))
            - float(row.get("actual_node_hours", -2.0))
        )
        if (
            matches[0].get("result") != row.get("result")
            or not math.isfinite(difference)
            or difference > 5.0e-7
        ):
            raise ValueError("recorded transaction reservation differs from ledger")
    elif row is not None:
        raise ValueError(f"{kind} transaction unexpectedly controls a ledger row")
    if kind == "submitted":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("submitted journal has invalid prepared manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        require_numeric_job_id(str(value.get("job_id", "")))
        if manifest.get("job_id") != value.get("job_id"):
            raise ValueError("submitted journal job ID differs from manifest")
        parse_utc_timestamp(
            value.get("submitted_recorded_utc"), "transaction submitted_recorded_utc"
        )
    if kind == "submit_cleared":
        digest = value.get("prepared_manifest_sha256")
        if not isinstance(digest, str) or SHA256_PATTERN.fullmatch(digest) is None:
            raise ValueError("cleared journal has invalid prepared manifest digest")
        validate_submission_audit(paths, value.get("submission_audit"))
        if not isinstance(value.get("recovery_notes"), str):
            raise ValueError("cleared submission transaction lacks recovery notes")
        validate_scheduler_absence_evidence(
            paths, value.get("scheduler_absence_evidence"), value
        )
    return value


def expected_pretransition_reservation(kind: str,
                                       reservation: dict[str, object],
                                       ) -> dict[str, object] | None:
    """Return the only permitted pre-transition form of one target reservation."""

    if kind == "prepared":
        return None
    result = dict(reservation)
    if kind == "submitted":
        result["state"] = "prepared"
        result.pop("job_id", None)
    elif kind == "recorded":
        result["state"] = "submitted"
        result.pop("actual_node_hours", None)
        result.pop("result", None)
    elif kind == "cancelled":
        result["state"] = "prepared"
        result.pop("notes", None)
    elif kind != "submit_cleared":
        raise ValueError(f"transaction kind has no reservation transition: {kind}")
    return result


def reservation_records_by_manifest(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    label: str,
) -> dict[Path, dict[str, object]]:
    """Index one validated reservation snapshot without duplicate targets."""

    indexed: dict[Path, dict[str, object]] = {}
    for reservation in reservations:
        validated = validate_reservation_record(paths, reservation)
        path = Path(str(validated["manifest"])).resolve()
        if path in indexed:
            raise ValueError(f"{label} contains duplicate manifest records")
        indexed[path] = validated
    return indexed


def authenticate_canonical_reservation_snapshot(
    paths: dict[str, Path],
    reservations: list[dict[str, object]],
    exempt_paths: frozenset[Path] = frozenset(),
) -> None:
    """Bind canonical baseline reservations to their retained launch intents."""

    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    indexed = reservation_records_by_manifest(
        paths, reservations, "canonical reservation snapshot"
    )
    for path, reservation in indexed.items():
        if path in exempt_paths:
            continue
        if not path.is_file():
            raise ValueError(
                f"canonical reservation lacks retained manifest during replay: {path}"
            )
        manifest = read_manifest(path)
        require_current_epoch(manifest, "canonical replay reservation manifest")
        validate_prepared_resources(manifest, canonical_production=True)
        require_reservation_matches_manifest(reservation, manifest)
        require_reserved_execution_intent(
            reservation, manifest, allow_legacy_local=False
        )


def validate_reservation_snapshot_transition(
    paths: dict[str, Path],
    kind: str,
    target: Path,
    prior: list[dict[str, object]],
    payload: list[dict[str, object]],
) -> None:
    """Require a journal snapshot to preserve every unrelated reservation."""

    prior_by_manifest = reservation_records_by_manifest(
        paths, prior, "transaction prior reservation snapshot"
    )
    payload_by_manifest = reservation_records_by_manifest(
        paths, payload, "transaction reservation snapshot"
    )
    authenticate_canonical_reservation_snapshot(
        paths, prior, exempt_paths=frozenset({target})
    )
    authenticate_canonical_reservation_snapshot(
        paths, payload, exempt_paths=frozenset({target})
    )
    if (
        {path: record for path, record in prior_by_manifest.items() if path != target}
        != {path: record for path, record in payload_by_manifest.items() if path != target}
    ):
        raise ValueError("transaction snapshot would alter unrelated reservations")
    payload_target = payload_by_manifest.get(target)
    if payload_target is None:
        raise ValueError("transaction snapshot lacks target reservation")
    if prior_by_manifest.get(target) != expected_pretransition_reservation(
        kind, payload_target
    ):
        raise ValueError("transaction snapshot target reservation transition is invalid")


def validate_transaction_reservation_baseline(
    paths: dict[str, Path],
    transaction: dict[str, object],
) -> None:
    """Authenticate reservation replay against its retained pre-transition store."""

    reservations = transaction.get("reservations")
    if not isinstance(reservations, list):
        raise ValueError("transaction reservation snapshot is invalid")
    current = read_reservations(paths)
    current_sha256 = sha256(paths["reservations"])
    payload_sha256 = stable_json_sha256(reservations)
    if current_sha256 not in {
        transaction.get("prior_reservations_sha256"), payload_sha256,
    }:
        raise ValueError("reservation store differs from transaction replay baseline")
    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    target = require_safe_manifest_path(paths, transaction["manifest_path"])
    reservation_records_by_manifest(paths, current, "reservation store")
    authenticate_canonical_reservation_snapshot(
        paths, current, exempt_paths=frozenset({target})
    )


def apply_transaction(paths: dict[str, Path], transaction_path: Path) -> None:
    """Idempotently complete one journaled metadata transition."""

    transaction = read_transaction(paths, transaction_path)
    kind = transaction.get("kind")
    if kind == "submit_pending":
        raise ValueError(
            "submission outcome is ambiguous; use recover-submit with the "
            f"scheduler job ID: {transaction_path}"
        )
    if kind not in {
        "prepared", "submitted", "submit_cleared", "recorded", "cancelled"
    }:
        raise ValueError(f"transaction journal has invalid kind: {transaction_path}")
    manifest_path = Path(str(transaction["manifest_path"])).resolve()
    manifest = transaction.get("manifest")
    reservations = transaction.get("reservations")
    if not isinstance(manifest, dict) or not isinstance(reservations, list):
        raise ValueError(f"transaction payload is incomplete: {transaction_path}")
    validate_transaction_reservation_baseline(paths, transaction)
    row = transaction.get("ledger_row")
    if row is not None:
        if not isinstance(row, dict):
            raise ValueError(f"transaction ledger row is invalid: {transaction_path}")
        ledger = read_ledger(paths)
        matches = [item for item in ledger if item.get("job_id") == row.get("job_id")]
        if len(matches) > 1:
            raise ValueError(f"transaction ledger job is duplicated: {transaction_path}")
        if matches and matches[0] != row:
            raise ValueError(f"transaction ledger row conflicts: {transaction_path}")
        if not matches:
            append_ledger_row(paths["ledger"], row, ledger)
    write_json(paths["reservations"], reservations)
    write_json(manifest_path, manifest)
    refresh_summary(paths)
    unlink_durable(transaction_path)


def durable_transition(paths: dict[str, Path], kind: str, manifest_path: Path,
                       manifest: dict[str, object],
                       reservations: list[dict[str, object]],
                       ledger_row: dict[str, object] | None = None) -> None:
    """Journal and apply one recoverable cross-file metadata transition."""

    require_no_pending_transactions(paths)
    require_safe_manifest_path(paths, manifest_path)
    prior_reservations = read_reservations(paths)
    transaction_path = (
        paths["transactions"] / f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}.json"
    )
    write_json(transaction_path, {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_path.stem,
        "kind": kind,
        "created_utc": utc_now(),
        "manifest_path": str(manifest_path),
        "prior_reservations": prior_reservations,
        "prior_reservations_sha256": stable_json_sha256(prior_reservations),
        "manifest": manifest,
        "reservations": reservations,
        "ledger_row": ledger_row,
    })
    apply_transaction(paths, transaction_path)


def write_submit_pending_transaction(paths: dict[str, Path],
                                     manifest_path: Path,
                                     submission_audit: dict[str, object],
                                     ) -> Path:
    """Persist an ambiguity barrier immediately before invoking sbatch."""

    require_no_pending_transactions(paths)
    require_safe_manifest_path(paths, manifest_path)
    prior_reservations = read_reservations(paths)
    transaction_path = (
        paths["transactions"] / f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}.json"
    )
    write_json(transaction_path, {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_path.stem,
        "kind": "submit_pending",
        "created_utc": utc_now(),
        "manifest_path": str(manifest_path),
        "prior_reservations": prior_reservations,
        "prior_reservations_sha256": stable_json_sha256(prior_reservations),
        "prepared_manifest_sha256": sha256(manifest_path),
        "submission_audit": submission_audit,
    })
    return transaction_path


def submit_pending_transaction(paths: dict[str, Path],
                               manifest_path: Path) -> Path:
    """Return the one ambiguous sbatch journal for a prepared manifest."""

    matches = []
    for path in pending_transaction_paths(paths):
        transaction = read_transaction(paths, path)
        if (
            transaction.get("kind") == "submit_pending"
            and Path(str(transaction.get("manifest_path", ""))).resolve()
            == manifest_path.resolve()
        ):
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(
            f"expected one pending submission journal for {manifest_path}, "
            f"found {len(matches)}"
        )
    return matches[0]


def finish_submit_transaction(paths: dict[str, Path], transaction_path: Path,
                              manifest_path: Path,
                              manifest: dict[str, object],
                              reservations: list[dict[str, object]],
                              job_id: str) -> None:
    """Attach a scheduler ID to an ambiguity journal and commit submission."""

    transaction = read_transaction(paths, transaction_path)
    if transaction.get("kind") != "submit_pending":
        raise ValueError(f"transaction is not an ambiguous submission: {transaction_path}")
    if transaction.get("prepared_manifest_sha256") != sha256(manifest_path):
        raise ValueError("prepared manifest changed after the sbatch boundary")
    if transaction.get("prior_reservations_sha256") != sha256(paths["reservations"]):
        raise ValueError("reservation store changed after the sbatch boundary")
    require_numeric_job_id(job_id)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if manifest.get("state") != "prepared" or reservation.get("state") != "prepared":
        raise ValueError("submission recovery requires matching prepared state")
    require_reserved_execution_intent(
        reservation, manifest,
        allow_legacy_local=paths["root"].resolve() != DEFAULT_ROOT.resolve(),
    )
    submitted_utc = utc_now()
    audit = transaction.get("submission_audit")
    if not isinstance(audit, dict):
        raise ValueError("submission transition lacks retained policy audit")
    audits = manifest.setdefault("submission_audits", [])
    if not isinstance(audits, list):
        raise ValueError("prepared manifest has invalid submission audits")
    audits.append(audit)
    manifest["state"] = "submitted"
    manifest["job_id"] = job_id
    manifest["submitted_recorded_utc"] = submitted_utc
    reservation["state"] = "submitted"
    reservation["job_id"] = job_id
    transaction.update({
        "kind": "submitted",
        "job_id": job_id,
        "submitted_recorded_utc": submitted_utc,
        "manifest": manifest,
        "reservations": reservations,
        "ledger_row": None,
    })
    write_json(transaction_path, transaction)
    apply_transaction(paths, transaction_path)


def active_reservations(reservations: list[dict[str, object]]
                        ) -> list[dict[str, object]]:
    """Return unaccounted prepared or submitted segment reservations."""

    return [
        item for item in reservations
        if item.get("state") in {"prepared", "submitted"}
    ]


def reservation_usage(paths: dict[str, Path]) -> tuple[float, float]:
    """Return actual and actively reserved current-epoch Stage I node-hours."""

    actual = sum(float(row["actual_node_hours"]) for row in read_ledger(paths))
    reserved = sum(
        float(item["reserved_node_hours"])
        for item in active_reservations(read_reservations(paths))
    )
    return actual, reserved


def refresh_summary(paths: dict[str, Path]) -> None:
    """Regenerate the human-readable current-epoch production summary."""

    with canonical_root_lock(paths["root"]):
        ledger = read_ledger(paths)
        reservations = read_reservations(paths)
        actual = sum(float(row["actual_node_hours"]) for row in ledger)
        active = active_reservations(reservations)
        reserved = sum(float(item["reserved_node_hours"]) for item in active)
        stage_remaining = max(
            0.0, CURRENT_STAGE_I_RESERVED_NODE_HOURS - actual - reserved
        )
        project_remaining = PROJECT_BUDGET_NODE_HOURS - actual - reserved
        qualification = qualification_approval_status(paths)
        if qualification["state"] == "approved":
            qualification_line = (
                "- E03 corrected-build Frontier qualification: approved by "
                f"`{qualification['approved_by']}` at "
                f"`{qualification['approved_utc']}` for executable "
                f"`{qualification['approved_executable_sha256']}`."
            )
        else:
            qualification_line = (
                "- E03 corrected-build Frontier qualification: pending until "
                f"a valid approval token exists at `{qualification['path']}` "
                f"({qualification['reason']})."
            )
        lines = [
            f"# MKS24 Stage I Frontier {EXECUTION_EPOCH} Budget",
            "",
            f"- Updated UTC: `{utc_now()}`",
            f"- Execution epoch: `{EXECUTION_EPOCH}`",
            qualification_line,
            f"- Fresh incremental project ceiling: "
            f"`{PROJECT_BUDGET_NODE_HOURS:.6f}` node-hours",
            f"- Historical debug qualification use, reported but not charged to E03: "
            f"`{HISTORICAL_DEBUG_NODE_HOURS:.6f}` node-hours",
            f"- Historical E01 Stage I use, reported but not charged to E03: "
            f"`{HISTORICAL_E01_STAGE_I_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 pipeline evidence, reported but not charged to E03: "
            f"`{HISTORICAL_E02_PIPELINE_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R16 use: `{COMPLETED_R16_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R02 standard-layout timing-pilot use: "
            f"`{COMPLETED_R02_STANDARD_LAYOUT_PILOT_NODE_HOURS:.6f}` node-hours",
            f"- Historical E02 R17 high-resolution timing-pilot use: "
            f"`{COMPLETED_R17_HIGH_RESOLUTION_PILOT_NODE_HOURS:.6f}` node-hours",
            f"- Current E03 mapped-matrix planning envelope: "
            f"`{MEASURED_STAGE_I_RESERVED_NODE_HOURS:.6f}` node-hours",
            f"- {EXECUTION_EPOCH} Stage I actual use: `{actual:.6f}` node-hours",
            f"- Active segment reservations: `{reserved:.6f}` node-hours",
            f"- Unreserved E03 mapped-matrix remainder: "
            f"`{stage_remaining:.6f}` node-hours",
            f"- Incremental project remainder after active E03 Stage I use: "
            f"`{project_remaining:.6f}` node-hours",
            "",
            "## Recorded Segments",
            "",
        ]
        if ledger:
            lines.extend([
                "| Job | Case/segment | State | Node-hours | Result |",
                "| --- | --- | --- | ---: | --- |",
            ])
            for row in ledger:
                lines.append(
                    "| `{job_id}` | `{case_id}/{segment}` | {state} | "
                    "`{actual_node_hours}` | {result} |".format(**row)
                )
        else:
            lines.append(
                f"No {EXECUTION_EPOCH} Stage I production allocation has been recorded."
            )
        lines.extend(["", "## Active Reservations", ""])
        if active:
            lines.extend([
                "| Case/segment | Nodes | Walltime | Node-hours | State |",
                "| --- | ---: | --- | ---: | --- |",
            ])
            for item in active:
                lines.append(
                    "| `{case_id}/{segment}` | `{nodes}` | `{requested_walltime}` | "
                    "`{reserved_node_hours:.6f}` | {state} |".format(**item)
                )
        else:
            lines.append(
                f"No prepared or submitted {EXECUTION_EPOCH} Stage I segment "
                "is reserved."
            )
        lines.extend([
            "",
            f"Only one {EXECUTION_EPOCH} Stage I segment may be prepared or "
            "submitted at a time. "
            "E02 is retained only as pipeline and cost evidence after the Phase A "
            "forcing-policy audit. "
            "Jobs use the `batch` partition with Frontier's default production "
            "`normal` QOS; the `debug` QOS is not used for paper production.",
            "",
        ])
        write_text(paths["summary"], "\n".join(lines))


def load_matrix(path: Path) -> dict[str, object]:
    """Read the mapped Stage I matrix."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid matrix manifest: {path}")
    return value


def parse_input_parameters(path: Path) -> dict[str, str]:
    """Read block-qualified input values while ignoring comments."""

    block = ""
    values: dict[str, str] = {}
    for original in path.read_text(encoding="utf-8").splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if "=" in line and block:
            parameter, value = line.split("=", 1)
            values[f"{block}/{parameter.strip()}"] = value.strip()
    return values


def time_tlim_override_target(overrides: list[str],
                              allow_missing: bool = False) -> float | None:
    """Return the single positive finite ``time/tlim`` execution target."""

    targets = [
        override.split("=", 1)[1]
        for override in overrides
        if "=" in override and override.split("=", 1)[0] == "time/tlim"
    ]
    if not targets and allow_missing:
        return None
    if len(targets) != 1:
        raise ValueError(
            "prepare requires exactly one --override time/tlim=<numeric-target>"
        )
    try:
        target = float(targets[0])
    except ValueError as error:
        raise ValueError("time/tlim override target must be numeric") from error
    if not math.isfinite(target) or target <= 0.0:
        raise ValueError("time/tlim override target must be positive and finite")
    return target


def validate_prepare_overrides(path: Path, overrides: list[str],
                               allow_missing_time_target: bool = False,
                               canonical_production: bool = False,
                               ) -> float | None:
    """Require declared input keys and one explicit segment time target."""

    keys = set(parse_input_parameters(path))
    if canonical_production and (
        len(overrides) != 1 or not overrides[0].startswith("time/tlim=")
    ):
        raise ValueError(
            "canonical production requires exactly one --override "
            "time/tlim=<positive-target>"
        )
    for override in overrides:
        if "=" not in override:
            raise ValueError(f"override must use block/name=value: {override}")
        key, _ = override.split("=", 1)
        if key not in keys:
            raise ValueError(
                f"override targets parameter absent from input deck: {key}"
            )
    return time_tlim_override_target(
        overrides, allow_missing=allow_missing_time_target
    )


def prepared_time_tlim_target(manifest: dict[str, object]) -> float:
    """Return and verify the prepared segment's retained execution target."""

    command = manifest.get("command")
    if not isinstance(command, dict):
        raise ValueError("prepared manifest lacks command metadata")
    overrides = command.get("overrides")
    if not isinstance(overrides, list) or not all(
        isinstance(item, str) for item in overrides
    ):
        raise ValueError("prepared manifest lacks command overrides")
    override_target = time_tlim_override_target(overrides)
    retained_target = command.get("time_tlim_target", override_target)
    try:
        target = float(retained_target)
    except (TypeError, ValueError) as error:
        raise ValueError("prepared manifest has invalid time/tlim target") from error
    if (
        not math.isfinite(target)
        or override_target is None
        or abs(target - override_target) > 1.0e-12
    ):
        raise ValueError("prepared manifest time/tlim target is inconsistent")
    return target


def validate_prepared_continuation_target(manifest: dict[str, object]) -> None:
    """Require a retained continuation target to advance beyond its parent."""

    command = manifest.get("command")
    if not isinstance(command, dict):
        raise ValueError("prepared manifest lacks command metadata")
    parent = command.get("parent_segment")
    if parent is None:
        return
    if not isinstance(parent, dict):
        raise ValueError("prepared continuation parent metadata is invalid")
    try:
        final_time = float(parent["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("prepared continuation lacks parent final time") from error
    target = prepared_time_tlim_target(manifest)
    if not math.isfinite(final_time) or target <= final_time + 1.0e-12:
        raise ValueError(
            "prepared continuation time/tlim target does not advance beyond "
            f"parent inspection time {final_time:.12g}"
        )


def validate_matrix(matrix_path: Path, source_dir: Path) -> dict[str, object]:
    """Validate unique case mapping and admitted file-level aliases."""

    matrix = load_matrix(matrix_path)
    cases = matrix.get("cases")
    if not isinstance(cases, list):
        raise ValueError("matrix cases must be a list")
    expected_ids = [f"R{number:02d}" for number in range(2, 18)]
    ids = [str(case.get("id")) for case in cases]
    if ids != expected_ids:
        raise ValueError(f"Stage I case identifiers must be {expected_ids}: {ids}")
    inputs = [str(case.get("input")) for case in cases]
    if len(set(inputs)) != len(inputs):
        raise ValueError("canonical Stage I inputs are not unique")
    if matrix.get("authorization", {}).get("mapped_unique_runs") != len(cases):
        raise ValueError("mapped run count disagrees with matrix cases")
    for path_text in inputs:
        if not (source_dir / path_text).is_file():
            raise ValueError(f"matrix input does not exist: {path_text}")
    excluded = {
        str(item.get("input"))
        for item in matrix.get("excluded_definitions", [])
    }
    if excluded.intersection(inputs):
        raise ValueError("an excluded input appears in the canonical run matrix")
    for alias in matrix.get("aliases", []):
        equivalent = alias.get("equivalent_input")
        if equivalent is None:
            continue
        canonical = source_dir / str(alias["canonical_input"])
        comparison = source_dir / str(equivalent)
        if not canonical.is_file() or not comparison.is_file():
            raise ValueError(f"alias inputs are missing for {alias['canonical_case']}")
        ignored = set(alias.get("comparison_ignore_parameters", []))
        canonical_values = parse_input_parameters(canonical)
        comparison_values = parse_input_parameters(comparison)
        for key in ignored:
            canonical_values.pop(key, None)
            comparison_values.pop(key, None)
        if canonical_values != comparison_values:
            raise ValueError(
                f"alias is not physically identical for {alias['canonical_case']}"
            )
    return matrix


def case_for_id(matrix: dict[str, object], case_id: str) -> dict[str, object]:
    """Select one unique canonical Stage I case."""

    matches = [
        case for case in matrix["cases"]
        if str(case.get("id")) == case_id
    ]
    if len(matches) != 1:
        raise ValueError(f"unknown or duplicate Stage I case ID: {case_id}")
    return matches[0]


def git_revision_for_input(source_dir: Path, input_path: Path,
                           matrix_path: Path) -> str:
    """Require submitted input and matrix to match committed source content."""

    revision = subprocess.run(
        ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    relative_paths = [
        str(input_path.relative_to(source_dir)),
        str(matrix_path.relative_to(source_dir)),
    ]
    for diff_args in (["diff", "--quiet", "--"], ["diff", "--cached", "--quiet", "--"]):
        result = subprocess.run(
            ["git", "-C", str(source_dir), *diff_args, *relative_paths],
            check=False,
        )
        if result.returncode != 0:
            raise ValueError(
                "submitted input and matrix must be committed before production"
            )
    return revision


def source_bundle_provenance(source_bundle_value: str | None,
                             revisions: list[str], root: Path,
                             allow_local_root: bool
                             ) -> dict[str, object] | None:
    """Require retained bundle provenance for real Stage I preparation."""

    if source_bundle_value is None:
        if allow_local_root:
            return None
        raise ValueError(
            "--source-bundle is required for retained Stage I source provenance"
        )
    source_bundle = Path(source_bundle_value).expanduser().resolve()
    require_beneath_root(
        source_bundle, root, "source bundle", allow_local_root
    )
    if not source_bundle.is_file():
        raise ValueError(f"source bundle is missing: {source_bundle}")
    with tempfile.TemporaryDirectory(prefix="cgl_lf_bundle_verify_") as directory:
        repository = Path(directory) / "source.git"
        try:
            subprocess.run(
                ["git", "clone", "--bare", "--quiet",
                 str(source_bundle), str(repository)],
                check=True, capture_output=True, text=True,
            )
        except subprocess.CalledProcessError as error:
            raise ValueError(
                f"source bundle cannot be cloned: {source_bundle}"
            ) from error
        for revision in revisions:
            present = subprocess.run(
                ["git", "-C", str(repository), "cat-file", "-e",
                 f"{revision}^{{commit}}"],
                check=False, capture_output=True, text=True,
            )
            if present.returncode != 0:
                raise ValueError(
                    f"source bundle does not contain revision {revision}: "
                    f"{source_bundle}"
                )
    return {
        "path": str(source_bundle),
        "sha256": sha256(source_bundle),
        "verified_revisions": revisions,
    }


def production_utility_provenance(allow_uncommitted: bool = False
                                  ) -> dict[str, object]:
    """Require the production-control script itself to be committed."""

    script_path = Path(__file__).resolve()
    relative = str(script_path.relative_to(ROOT_DIR))
    if Path(relative) != PRODUCTION_UTILITY_RELATIVE:
        raise ValueError("production utility path is inconsistent")
    revision = subprocess.run(
        ["git", "-C", str(ROOT_DIR), "rev-parse", "--verify", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if GIT_REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError("production utility revision is invalid")
    tracked = subprocess.run(
        ["git", "-C", str(ROOT_DIR), "ls-files", "--error-unmatch", "--", relative],
        check=False, capture_output=True, text=True,
    )
    if tracked.returncode != 0:
        raise ValueError("production utility is not tracked by Git")
    committed = True
    for diff_args in (["diff", "--quiet", "--"], ["diff", "--cached", "--quiet", "--"]):
        result = subprocess.run(
            ["git", "-C", str(ROOT_DIR), *diff_args, relative],
            check=False,
        )
        if result.returncode == 1:
            committed = False
        elif result.returncode != 0:
            raise ValueError("cannot determine production utility worktree status")
    if not committed and not allow_uncommitted:
        raise ValueError(
            "production utility must be committed before preparing a segment"
        )
    return {
        "path": str(script_path),
        "revision": revision,
        "sha256": sha256(script_path),
        "committed": committed,
    }


def authenticate_production_utility(
    record: object,
    *,
    source_bundle: object = None,
    allow_uncommitted: bool = False,
    allow_historical: bool = False,
) -> None:
    """Revalidate the retained production helper used to prepare a segment."""

    if not isinstance(record, dict):
        raise ValueError("prepared manifest lacks production utility provenance")
    path = Path(str(record.get("path", ""))).resolve()
    if path != Path(__file__).resolve():
        raise ValueError("prepared production utility path is inconsistent")
    try:
        relative = path.relative_to(ROOT_DIR)
    except ValueError as error:
        raise ValueError("prepared production utility path is inconsistent") from error
    if relative != PRODUCTION_UTILITY_RELATIVE:
        raise ValueError("prepared production utility path is inconsistent")
    revision = record.get("revision")
    if not isinstance(revision, str) or GIT_REVISION_PATTERN.fullmatch(revision) is None:
        raise ValueError("prepared production utility revision is invalid")
    if record.get("committed") is not True and not (
        allow_uncommitted and record.get("committed") is False
    ):
        raise ValueError("prepared production utility is not committed")
    expected = record.get("sha256")
    if not isinstance(expected, str) or SHA256_PATTERN.fullmatch(expected) is None:
        raise ValueError("prepared production utility checksum is invalid")
    production_utility_provenance(allow_uncommitted=allow_uncommitted)
    if not allow_historical:
        require_file_sha256(path, expected, "production utility")
        return
    if source_bundle is None and allow_uncommitted:
        require_file_sha256(path, expected, "production utility")
        return
    if not isinstance(source_bundle, dict):
        raise ValueError(
            "recorded manifest lacks historical production utility bundle provenance"
        )
    bundle_path = Path(str(source_bundle.get("path", ""))).resolve()
    require_file_sha256(bundle_path, source_bundle.get("sha256"), "source bundle")
    revisions = source_bundle.get("verified_revisions")
    if (
        not isinstance(revisions, list)
        or revision not in revisions
    ):
        raise ValueError(
            "recorded source bundle lacks the historical production utility revision"
        )
    with tempfile.TemporaryDirectory(prefix="cgl_lf_utility_verify_") as directory:
        repository = Path(directory) / "source.git"
        try:
            subprocess.run(
                ["git", "clone", "--bare", "--quiet",
                 str(bundle_path), str(repository)],
                check=True, capture_output=True, text=True,
            )
        except subprocess.CalledProcessError as error:
            raise ValueError(
                f"recorded source bundle cannot be cloned: {bundle_path}"
            ) from error
        historical = subprocess.run(
            ["git", "-C", str(repository), "show",
             f"{revision}:{PRODUCTION_UTILITY_RELATIVE}"],
            check=False, capture_output=True,
        )
    if (
        historical.returncode != 0
        or hashlib.sha256(historical.stdout).hexdigest() != expected
    ):
        raise ValueError(f"prepared production utility checksum has changed: {path}")


def read_build_provenance(executable: Path,
                          build_manifest: Path) -> dict[str, str]:
    """Verify an archived executable against its immutable build manifest."""

    digest_file = build_manifest / "athena.sha256"
    environment_file = build_manifest / "environment.txt"
    if not digest_file.is_file() or not environment_file.is_file():
        raise ValueError(f"incomplete build manifest directory: {build_manifest}")
    recorded_sha = digest_file.read_text(encoding="utf-8").split()[0]
    actual_sha = sha256(executable)
    if recorded_sha != actual_sha:
        raise ValueError("executable digest does not match build manifest")
    revision_match = re.search(
        r"(?m)^git_revision=([0-9a-f]{40})\s*$",
        environment_file.read_text(encoding="utf-8"),
    )
    if revision_match is None:
        raise ValueError("build manifest does not record a full git revision")
    return {
        "revision": revision_match.group(1),
        "sha256": actual_sha,
        "manifest_dir": str(build_manifest),
    }


@locked_root_action
def approve_qualification(args: argparse.Namespace) -> int:
    """Atomically retain reviewed corrected-build qualification approval."""

    if not args.confirm_corrected_build_frontier_qualified:
        raise ValueError("--confirm-corrected-build-frontier-qualified is required")
    if not args.approved_by.strip() or not args.review_notes.strip():
        raise ValueError("--approved-by and --review-notes must be nonempty")
    root = require_root(Path(args.root), args.allow_local_root)
    paths = initialize(root)
    require_reconciled_store_consistency(
        paths,
        allow_absent_qualification=True,
        allow_invalid_qualification=args.replace_existing_approval,
    )
    executable = Path(args.executable).expanduser().resolve()
    build_manifest = Path(args.build_manifest).expanduser().resolve()
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(
        build_manifest, root, "build manifest", args.allow_local_root
    )
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is unavailable: {executable}")
    provenance = read_build_provenance(executable, build_manifest)
    retained_manifests = sorted(
        paths["runs"].glob("*/*/manifest/prepared_run.json")
    )
    if retained_manifests:
        raise ValueError(
            "E03 qualification approval cannot change after segment preparation: "
            + ", ".join(str(path) for path in retained_manifests)
        )
    if paths["qualification"].exists() and not args.replace_existing_approval:
        raise ValueError(
            "E03 qualification approval token already exists; pass "
            "--replace-existing-approval only after reviewing the replacement build"
        )
    approval = {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "approval_scope": "corrected-build Frontier qualification for E03 prepare",
        "approved_utc": utc_now(),
        "approved_by": args.approved_by,
        "review_notes": args.review_notes,
        "approved_executable": str(executable),
        "approved_executable_sha256": provenance["sha256"],
        "approved_executable_revision": provenance["revision"],
        "build_manifest": provenance["manifest_dir"],
    }
    write_json(paths["qualification"], approval)
    refresh_summary(paths)
    print(f"Wrote E03 qualification approval token: {paths['qualification']}")
    print(f"Approved executable sha256: {provenance['sha256']}")
    print(f"Approved git revision: {provenance['revision']}")
    return 0


def quote(value: str | Path) -> str:
    """Quote a shell literal in the generated Slurm script."""

    return shlex.quote(str(value))


def expected_job_name(manifest: dict[str, object]) -> str:
    """Return the exact safe Slurm job name for one prepared segment."""

    run = manifest["run"]
    segment = require_safe_segment(str(run["segment"]))
    job_name = (
        f"cgl_mks24_{EXECUTION_EPOCH_SLUG}_{run['case_id']}_{segment}"
    )
    job_name = re.sub(r"[^A-Za-z0-9_]+", "_", job_name)[:60]
    if re.fullmatch(r"[A-Za-z0-9_]+", job_name) is None:
        raise ValueError(f"generated Stage I job name is unsafe: {job_name}")
    return job_name


def normalized_batch_script_text(value: str) -> str:
    """Normalize the embedded self-digest before hashing a batch script."""

    return BATCH_SCRIPT_DIGEST_PATTERN.sub(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}", value
    )


def normalized_batch_script_sha256(path: Path) -> str:
    """Return the digest authenticated by a generated batch script itself."""

    return hashlib.sha256(
        normalized_batch_script_text(path.read_text(encoding="utf-8")).encode("utf-8")
    ).hexdigest()


def finalize_batch_script(value: str) -> tuple[str, str]:
    """Embed a normalized self-digest in one generated batch script."""

    matches = BATCH_SCRIPT_DIGEST_PATTERN.findall(value)
    if matches != [BATCH_SCRIPT_DIGEST_PLACEHOLDER]:
        raise ValueError("generated batch script lacks one self-digest placeholder")
    digest = hashlib.sha256(normalized_batch_script_text(value).encode("utf-8")).hexdigest()
    return value.replace(
        f"BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}",
        f"BATCH_SCRIPT_SHA256={digest}",
        1,
    ), digest


def require_file_sha256(path: Path, expected: object, label: str) -> None:
    """Require one retained prepared artifact to preserve its digest."""

    if not path.is_file():
        raise ValueError(f"prepared {label} is missing: {path}")
    if not isinstance(expected, str) or sha256(path) != expected:
        raise ValueError(f"prepared {label} checksum has changed: {path}")


def prepared_restart_inventory(manifest_path: Path,
                               command: dict[str, object]) -> list[Path]:
    """Return the exact retained restart archive inventory."""

    records = command.get("restart_files")
    if not isinstance(records, list):
        raise ValueError("prepared manifest lacks restart sibling metadata")
    paths = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError("prepared restart sibling metadata is invalid")
        paths.append(Path(str(record.get("path", ""))).resolve())
    restart = command.get("restart_file")
    if restart is None:
        if paths:
            raise ValueError("prepared manifest retains restart siblings without a restart")
        return []
    restart_path = Path(str(restart)).resolve()
    if restart_path not in paths:
        raise ValueError("prepared primary restart is absent from retained siblings")
    archive_root = manifest_path.parent / "submitted_restart"
    if archive_root.is_dir():
        actual = sorted(path.resolve() for path in archive_root.rglob("*") if path.is_file())
    else:
        actual = [restart_path] if restart_path.is_file() else []
    if sorted(paths) != actual:
        raise ValueError("prepared restart archive inventory has changed")
    return paths


def restart_time_marker(path: Path) -> float:
    """Read the explicit physical-time marker from a restart parameter dump."""

    marker = b"<par_end>"
    header = b""
    with path.open("rb") as stream:
        while marker not in header and len(header) <= MAX_RESTART_PARAMETER_DUMP_BYTES:
            block = stream.read(65536)
            if not block:
                break
            header += block
    if len(header) > MAX_RESTART_PARAMETER_DUMP_BYTES:
        raise ValueError(f"restart parameter dump is implausibly large: {path}")
    end = header.find(marker)
    if end < 0:
        raise ValueError(f"restart parameter dump lacks <par_end>: {path}")
    try:
        text = header[:end].decode("utf-8")
    except UnicodeDecodeError as error:
        raise ValueError(f"restart parameter dump is not UTF-8 text: {path}") from error
    block = ""
    markers = []
    for original in text.splitlines():
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = line.split("=", 1)
            if key.strip() == "restart_time":
                markers.append(value.strip())
    if len(markers) != 1:
        raise ValueError(
            f"restart parameter dump must contain one time/restart_time marker: {path}"
        )
    try:
        result = float(markers[0])
    except ValueError as error:
        raise ValueError(f"restart time marker is not numeric: {path}") from error
    if not math.isfinite(result):
        raise ValueError(f"restart time marker is not finite: {path}")
    return result


def restart_product_time(paths: list[Path],
                         allow_missing_marker: bool = False) -> float | None:
    """Require every selected restart sibling to retain one physical time."""

    try:
        times = [restart_time_marker(path) for path in paths]
    except (OSError, ValueError):
        if allow_missing_marker:
            return None
        raise
    if not times:
        raise ValueError("restart product has no selected siblings")
    if any(abs(value - times[0]) > 1.0e-12 for value in times[1:]):
        raise ValueError("restart sibling physical-time markers disagree")
    return times[0]


def validate_prepared_resources(manifest: dict[str, object],
                                canonical_production: bool) -> None:
    """Require retained Slurm and Athena resource policy to remain valid."""

    allocation = manifest.get("allocation")
    command = manifest.get("command")
    if not isinstance(allocation, dict) or not isinstance(command, dict):
        raise ValueError("prepared manifest lacks resource metadata")
    try:
        nodes = int(allocation["nodes"])
        requested_seconds = parse_walltime(str(allocation["requested_walltime"]))
        athena_seconds = parse_walltime(str(command["athena_walltime"]))
        ranks_per_node = int(allocation["ranks_per_node"])
        cpus_per_task = int(allocation["cpus_per_task"])
        reserved_node_hours = float(allocation["reserved_node_hours"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("prepared manifest has invalid resource metadata") from error
    if (
        nodes < 1
        or requested_seconds <= 0
        or athena_seconds <= 0
        or not math.isfinite(reserved_node_hours)
        or reserved_node_hours <= 0.0
    ):
        raise ValueError("prepared manifest resource values must be positive")
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError("prepared Slurm walltime exceeds the Stage I limit")
    if athena_seconds >= requested_seconds:
        raise ValueError("prepared Athena walltime is not shorter than Slurm walltime")
    if athena_seconds > requested_seconds - 600:
        raise ValueError("prepared Athena walltime lacks the shutdown margin")
    if (
        allocation.get("requested_seconds") != requested_seconds
        or abs(reserved_node_hours - node_hours(nodes, requested_seconds)) > 5.0e-12
    ):
        raise ValueError("prepared resource accounting is inconsistent")
    if canonical_production and (
        ranks_per_node != EXPECTED_RANKS_PER_NODE
        or cpus_per_task != EXPECTED_CPUS_PER_TASK
    ):
        raise ValueError("canonical Frontier resource shape is inconsistent")
    if canonical_production:
        run = manifest.get("run")
        if not isinstance(run, dict):
            raise ValueError("prepared manifest lacks run metadata")
        require_case_node_count(str(run.get("case_id")), nodes)


def require_reservation_matches_manifest(reservation: dict[str, object],
                                         manifest: dict[str, object]) -> None:
    """Bind one canonical reservation snapshot to its manifest allocation."""

    allocation = manifest.get("allocation")
    run = manifest.get("run")
    if not isinstance(allocation, dict) or not isinstance(run, dict):
        raise ValueError("transaction manifest lacks reservation metadata")
    try:
        matches = (
            reservation.get("case_id") == run["case_id"]
            and reservation.get("case_name") == run["case_name"]
            and reservation.get("segment") == run["segment"]
            and int(reservation["nodes"]) == int(allocation["nodes"])
            and reservation.get("requested_walltime")
            == allocation["requested_walltime"]
            and abs(
                float(reservation["reserved_node_hours"])
                - float(allocation["reserved_node_hours"])
            ) <= 5.0e-12
        )
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "transaction manifest has invalid reservation metadata"
        ) from error
    if not matches:
        raise ValueError("transaction reservation allocation differs from manifest")


def authenticate_prepared_execution(manifest: dict[str, object],
                                    manifest_path: Path,
                                    allow_legacy_local: bool = False) -> None:
    """Authenticate every prepared artifact used by one production launch."""

    command = manifest.get("command")
    paths = manifest.get("paths")
    if allow_legacy_local and (
        not isinstance(command, dict) or not isinstance(paths, dict)
    ):
        return
    if not isinstance(command, dict) or not isinstance(paths, dict):
        raise ValueError("prepared manifest lacks execution metadata")
    required = (
        "batch_script_sha256",
        "input_file",
        "input_sha256",
        "matrix_file",
        "matrix_sha256",
        "executable",
        "executable_sha256",
        "restart_files",
        "production_utility",
    )
    if allow_legacy_local and any(key not in command for key in required):
        return
    overrides = command.get("overrides")
    if not allow_legacy_local and (
        not isinstance(overrides, list)
        or len(overrides) != 1
        or not isinstance(overrides[0], str)
        or not overrides[0].startswith("time/tlim=")
    ):
        raise ValueError(
            "canonical prepared manifest requires exactly one time/tlim override"
        )
    if not allow_legacy_local or isinstance(overrides, list):
        prepared_time_tlim_target(manifest)
    validate_prepared_continuation_target(manifest)
    validate_prepared_resources(manifest, canonical_production=not allow_legacy_local)
    bundle = command.get("source_bundle")
    recorded = manifest.get("state") == "recorded"
    authenticate_production_utility(
        command.get("production_utility"),
        source_bundle=bundle,
        allow_uncommitted=allow_legacy_local,
        allow_historical=recorded,
    )
    batch_script = Path(str(paths.get("batch_script", ""))).resolve()
    if batch_script != (manifest_path.parent / "cgl_lf_stage_i.sbatch").resolve():
        raise ValueError("prepared batch script path is inconsistent")
    if not batch_script.is_file():
        raise ValueError(f"prepared batch script is missing: {batch_script}")
    script_text = batch_script.read_text(encoding="utf-8")
    embedded = BATCH_SCRIPT_DIGEST_PATTERN.findall(script_text)
    if embedded != [command.get("batch_script_sha256")]:
        raise ValueError("prepared batch script self-digest is inconsistent")
    if normalized_batch_script_sha256(batch_script) != command.get(
        "batch_script_sha256"
    ):
        raise ValueError("prepared batch script normalized checksum has changed")
    if not recorded:
        expected_script = generated_batch_script(manifest, manifest_path)
        if normalized_batch_script_text(script_text) != normalized_batch_script_text(
            expected_script
        ):
            raise ValueError("prepared batch script differs from retained launch intent")
    require_file_sha256(
        Path(str(command.get("input_file", ""))).resolve(),
        command.get("input_sha256"),
        "input",
    )
    require_file_sha256(
        Path(str(command.get("matrix_file", ""))).resolve(),
        command.get("matrix_sha256"),
        "matrix",
    )
    executable = Path(str(command.get("executable", ""))).resolve()
    require_file_sha256(executable, command.get("executable_sha256"), "executable")
    if not os.access(executable, os.X_OK):
        raise ValueError(f"prepared executable is not executable: {executable}")
    qualification = command.get("qualification_approval")
    if qualification is None:
        if not allow_legacy_local:
            raise ValueError("prepared manifest lacks E03 qualification approval")
    elif not isinstance(qualification, dict):
        raise ValueError("prepared E03 qualification approval metadata is invalid")
    else:
        approval_path = Path(str(qualification.get("path", ""))).resolve()
        expected_path = layout(
            Path(str(manifest.get("project_root", ""))).resolve()
        )["qualification"].resolve()
        if approval_path != expected_path:
            raise ValueError("prepared E03 qualification token path is inconsistent")
        require_file_sha256(
            approval_path, qualification.get("sha256"), "E03 qualification token"
        )
        approval = read_qualification_approval(approval_path)
        if (
            qualification.get("execution_epoch") != EXECUTION_EPOCH
            or approval["approved_executable_sha256"]
            != command.get("executable_sha256")
            or approval["approved_executable_revision"]
            != command.get("executable_revision")
            or qualification.get("approved_executable_sha256")
            != approval["approved_executable_sha256"]
            or qualification.get("approved_executable_revision")
            != approval["approved_executable_revision"]
            or qualification.get("token") != approval
        ):
            raise ValueError(
                "prepared E03 qualification approval does not match the executable"
            )
    if bundle is None:
        if not allow_legacy_local:
            raise ValueError("prepared manifest lacks source bundle provenance")
    elif not isinstance(bundle, dict):
        raise ValueError("prepared source bundle metadata is invalid")
    else:
        revisions = bundle.get("verified_revisions")
        expected_revisions = {
            command.get("input_revision"),
            command.get("executable_revision"),
            command["production_utility"].get("revision"),
        }
        if (
            not isinstance(revisions, list)
            or not all(isinstance(item, str) for item in revisions)
            or not expected_revisions.issubset(set(revisions))
        ):
            raise ValueError(
                "prepared source bundle lacks launch provenance revisions"
            )
        require_file_sha256(
            Path(str(bundle.get("path", ""))).resolve(),
            bundle.get("sha256"),
            "source bundle",
        )
    records = command.get("restart_files")
    if not isinstance(records, list):
        raise ValueError("prepared manifest lacks restart sibling metadata")
    for record in records:
        revalidate_retained_file(record, label="prepared restart sibling")
    restart_paths = prepared_restart_inventory(manifest_path, command)
    if restart_paths:
        marker_bypass = (
            allow_legacy_local
            and command.get("allow_missing_restart_time_marker") is True
        )
        if (
            command.get("allow_missing_restart_time_marker") is True
            and not allow_legacy_local
        ):
            raise ValueError(
                "canonical prepared restart may not bypass physical-time markers"
            )
        parent = command.get("parent_segment")
        if not isinstance(parent, dict):
            if not allow_legacy_local:
                raise ValueError("prepared restart lacks inspected parent metadata")
        else:
            marker = restart_product_time(
                restart_paths, allow_missing_marker=marker_bypass
            )
            try:
                parent_final_time = float(parent["final_time"])
            except (KeyError, TypeError, ValueError) as error:
                raise ValueError("prepared restart parent lacks final time") from error
            if marker is not None and abs(marker - parent_final_time) > 1.0e-12:
                raise ValueError(
                    "prepared restart physical time differs from parent inspection"
                )


def generated_batch_script(manifest: dict[str, object],
                           manifest_path: Path) -> str:
    """Generate an authenticated normal-QOS production segment."""

    run = manifest["run"]
    allocation = manifest["allocation"]
    command = manifest["command"]
    paths = manifest["paths"]
    overrides = " ".join(quote(value) for value in command["overrides"])
    if overrides:
        overrides = " " + overrides
    restart = command.get("restart_file")
    restart_literal = quote(str(restart)) if restart else "''"
    job_name = expected_job_name(manifest)
    bundle = command.get("source_bundle")
    runtime_checks = [
        f"require_sha {quote(command['production_utility']['sha256'])} "
        f"{quote(command['production_utility']['path'])} production_utility",
        f"require_sha {quote(command['input_sha256'])} \"${{INPUT}}\" input",
        f"require_sha {quote(command['matrix_sha256'])} {quote(command['matrix_file'])} matrix",
        f"require_sha {quote(command['executable_sha256'])} \"${{ATHENA}}\" executable",
    ]
    if isinstance(bundle, dict):
        runtime_checks.append(
            f"require_sha {quote(bundle['sha256'])} {quote(bundle['path'])} source_bundle"
        )
    qualification = command.get("qualification_approval")
    if isinstance(qualification, dict):
        runtime_checks.append(
            f"require_sha {quote(qualification['sha256'])} "
            f"{quote(qualification['path'])} qualification_approval"
        )
    restart_records = command.get("restart_files", [])
    if not isinstance(restart_records, list):
        raise ValueError("prepared restart sibling metadata is invalid")
    for index, record in enumerate(restart_records):
        runtime_checks.append(
            f"require_sha {quote(record['sha256'])} {quote(record['path'])} "
            f"restart_{index:04d}"
        )
    archive_root = manifest_path.parent / "submitted_restart"
    if restart_records and archive_root.is_dir():
        runtime_checks.append(
            f'test "$(find {quote(archive_root)} -type f -print | wc -l)" '
            f'-eq {len(restart_records)} || {{ echo "restart inventory changed" >&2; '
            "exit 1; }"
        )
    runtime_checks_text = "\n".join(runtime_checks)
    return f"""#!/bin/bash
#SBATCH -J {job_name}
#SBATCH -A {ACCOUNT}
#SBATCH -o {paths["slurm_log"]}
#SBATCH -p {PARTITION}
# Frontier default normal QOS is required for production; do not add -q debug.
#SBATCH -t {allocation["requested_walltime"]}
#SBATCH -N {allocation["nodes"]}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail

RUN_MANIFEST={quote(manifest_path)}
BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}
ATHENA={quote(command["executable"])}
INPUT={quote(command["input_file"])}
RESTART={restart_literal}
OUT_DIR={quote(paths["output_dir"])}
ENV_LOG={quote(paths["environment_log"])}
RANKS_PER_NODE={allocation["ranks_per_node"]}
CPUS_PER_TASK={allocation["cpus_per_task"]}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * RANKS_PER_NODE))"

require_sha() {{
  local expected="$1"
  local path="$2"
  local label="$3"
  local actual
  test -f "${{path}}" || {{ echo "missing ${{label}}: ${{path}}" >&2; exit 1; }}
  actual="$(sha256sum "${{path}}" | awk '{{print $1}}')"
  test "${{actual}}" = "${{expected}}" || {{
    echo "checksum mismatch for ${{label}}: ${{path}}" >&2
    exit 1
  }}
}}

normalized_script_sha256() {{
  sed -E 's/^BATCH_SCRIPT_SHA256=[0-9a-f]{{64}}$/BATCH_SCRIPT_SHA256={BATCH_SCRIPT_DIGEST_PLACEHOLDER}/' "$0" |
    sha256sum | awk '{{print $1}}'
}}

test -f "${{RUN_MANIFEST}}"
test -x "${{ATHENA}}"
test -f "${{INPUT}}"
if [[ -n "${{RESTART}}" ]]; then
  test -f "${{RESTART}}"
fi
test "$(normalized_script_sha256)" = "${{BATCH_SCRIPT_SHA256}}" || {{
  echo "batch script checksum mismatch" >&2
  exit 1
}}
{runtime_checks_text}
mkdir -p "${{OUT_DIR}}"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_ENV_DISPLAY=1
export MPICH_VERSION_DISPLAY=1
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

{{
  date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ"
  echo "slurm_job_id=${{SLURM_JOB_ID:?Missing SLURM_JOB_ID}}"
  echo "prepared_manifest=${{RUN_MANIFEST}}"
  echo "nodes=${{NNODES}}"
  echo "ranks=${{NRANKS}}"
  module -t list 2>&1
  env | LC_ALL=C sort | grep -E '^(MPICH_|FI_|HSA_|OMP_|ROCR_|HIP_|CRAY_)' || true
}} > "${{ENV_LOG}}"

RUN_ARGS=(-i "${{INPUT}}")
if [[ -n "${{RESTART}}" ]]; then
  RUN_ARGS=(-r "${{RESTART}}")
fi

srun -N "${{NNODES}}" -n "${{NRANKS}}" --ntasks-per-node="${{RANKS_PER_NODE}}" \
  -c "${{CPUS_PER_TASK}}" --threads-per-core=1 --cpu-bind=threads \
  --gpus-per-task=1 --gpu-bind=closest \
  "${{ATHENA}}" "${{RUN_ARGS[@]}}" -d "${{OUT_DIR}}" \
  -t {quote(command["athena_walltime"])} \
  job/basename={quote(run["run_basename"])}{overrides}

date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "${{ENV_LOG}}"
"""


def execution_intent_sha256(manifest: dict[str, object]) -> str:
    """Digest the immutable prepared execution intent across lifecycle updates."""

    try:
        intent = {
            key: manifest[key]
            for key in (
                "schema_version", "execution_epoch", "project_root", "policy",
                "run", "allocation", "command", "paths",
            )
        }
    except KeyError as error:
        raise ValueError(
            f"prepared manifest lacks immutable execution intent field {error.args[0]}"
        ) from error
    return hashlib.sha256(
        json.dumps(intent, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def require_reserved_execution_intent(reservation: dict[str, object],
                                      manifest: dict[str, object],
                                      allow_legacy_local: bool) -> None:
    """Require reservation metadata to bind the immutable prepared launch."""

    expected = reservation.get("execution_intent_sha256")
    if expected is None and allow_legacy_local:
        return
    if (
        not isinstance(expected, str)
        or expected != execution_intent_sha256(manifest)
    ):
        raise ValueError("prepared execution intent differs from its reservation")


@locked_root_action
def prepare(args: argparse.Namespace) -> Path:
    """Create one retained, sequentially submitted production segment."""

    root = require_root(Path(args.root), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    if offline_local_root:
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    require_authorized_case(args.case_id)
    require_safe_segment(args.segment)
    if active_reservations(read_reservations(paths)):
        raise ValueError(
            "another Stage I segment is prepared or submitted; "
            "record or cancel it before preparing a new segment"
        )
    if args.nodes < 1:
        raise ValueError("--nodes must be positive")
    require_prepare_case_policy(
        paths, args.case_id, args.nodes, offline_local_root
    )
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    executable = Path(args.executable).expanduser().resolve()
    build_manifest = Path(args.build_manifest).expanduser().resolve()
    restart = (
        Path(args.restart_file).expanduser().resolve()
        if args.restart_file else None
    )
    require_beneath_root(executable, root, "executable", args.allow_local_root)
    require_beneath_root(
        build_manifest, root, "build manifest", args.allow_local_root
    )
    if restart is not None:
        require_beneath_root(restart, root, "restart file", args.allow_local_root)
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"executable is unavailable: {executable}")
    matrix = validate_matrix(matrix_path, source_dir)
    case = case_for_id(matrix, args.case_id)
    input_path = source_dir / str(case["input"])
    allow_missing_time_target = getattr(
        args, "allow_missing_time_target", False
    )
    if allow_missing_time_target and not offline_local_root:
        raise ValueError(
            "--allow-missing-time-target is restricted to offline validation"
        )
    time_tlim_target = validate_prepare_overrides(
        input_path, args.override,
        allow_missing_time_target=allow_missing_time_target,
        canonical_production=not offline_local_root,
    )
    input_revision = git_revision_for_input(source_dir, input_path, matrix_path)
    utility_provenance = production_utility_provenance(
        allow_uncommitted=offline_local_root
    )
    provenance = read_build_provenance(executable, build_manifest)
    qualification_approval = require_qualification_approval(
        paths, provenance["sha256"], provenance["revision"], offline_local_root
    )
    bundle_provenance = source_bundle_provenance(
        args.source_bundle,
        list(dict.fromkeys([
            input_revision,
            provenance["revision"],
            utility_provenance["revision"],
        ])),
        root,
        offline_local_root,
    )
    if restart is not None and not restart.is_file():
        raise ValueError(f"restart file is unavailable: {restart}")
    allow_missing_restart_time_marker = getattr(
        args, "allow_missing_restart_time_marker", False
    )
    if allow_missing_restart_time_marker and not offline_local_root:
        raise ValueError(
            "--allow-missing-restart-time-marker is restricted to offline validation"
        )
    parent_segment = (
        verify_continuation_restart(
            restart,
            allow_missing_restart_time_marker=allow_missing_restart_time_marker,
        )
        if restart else None
    )
    if parent_segment is not None:
        if parent_segment["case_id"] != args.case_id:
            raise ValueError("continuation restart belongs to a different case")
        if parent_segment["input_sha256"] != sha256(input_path):
            raise ValueError("continuation restart input differs from its parent")
        if parent_segment["executable_sha256"] != provenance["sha256"]:
            raise ValueError("continuation executable differs from its parent")
        if (
            time_tlim_target is not None
            and time_tlim_target <= float(parent_segment["final_time"]) + 1.0e-12
        ):
            raise ValueError(
                "continuation time/tlim target must advance beyond its parent "
                f"inspection time {float(parent_segment['final_time']):.12g}"
            )
    requested_seconds = parse_walltime(args.walltime)
    if requested_seconds <= 0:
        raise ValueError("--walltime must be positive")
    if requested_seconds > MAX_SEGMENT_SECONDS:
        raise ValueError(
            "a Stage I normal-QOS segment may not request more than two hours"
        )
    athena_seconds = parse_walltime(args.athena_walltime)
    if athena_seconds <= 0:
        raise ValueError("--athena-walltime must be positive")
    if athena_seconds >= requested_seconds:
        raise ValueError("Athena walltime must be shorter than the Slurm walltime")
    if athena_seconds > requested_seconds - 600:
        raise ValueError("Athena walltime must leave ten minutes for shutdown")
    if not offline_local_root and (
        args.ranks_per_node != EXPECTED_RANKS_PER_NODE
        or args.cpus_per_task != EXPECTED_CPUS_PER_TASK
    ):
        raise ValueError(
            "canonical Frontier production requires --ranks-per-node=8 "
            "and --cpus-per-task=7"
        )
    segment_hours = node_hours(args.nodes, requested_seconds)
    actual, reserved = reservation_usage(paths)
    if actual + reserved + segment_hours > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("proposed segment exceeds the Stage I reservation")
    if (
        actual + reserved + segment_hours > PROJECT_BUDGET_NODE_HOURS
    ):
        raise ValueError("proposed segment exceeds the project ceiling")
    run_dir = paths["runs"] / args.case_id / args.segment
    manifest_dir = run_dir / "manifest"
    manifest_path = manifest_dir / "prepared_run.json"
    if run_dir.exists():
        raise ValueError(f"segment run directory already exists: {run_dir}")
    output_dir = run_dir / "output"
    for path in (manifest_dir, output_dir):
        mkdir_durable(path)
    archived_input = manifest_dir / "submitted_input.athinput"
    archived_matrix = manifest_dir / "mks24_stage_i_manifest.json"
    copy_file(input_path, archived_input)
    copy_file(matrix_path, archived_matrix)
    archived_restart = None
    archived_restart_files: list[Path] = []
    if restart is not None:
        source_restart_files = (
            [Path(str(path)) for path in parent_segment["restart_files"]]
            if parent_segment is not None
            else [restart]
        )
        if restart.parent.name.startswith("rank_"):
            archive_root = manifest_dir / "submitted_restart"
            for source in source_restart_files:
                target = archive_root / source.parent.name / source.name
                mkdir_durable(target.parent)
                copy_file(source, target)
                archived_restart_files.append(target)
            archived_restart = (
                archive_root / restart.parent.name / restart.name
            )
        else:
            archived_restart = manifest_dir / "submitted_restart.rst"
            copy_file(restart, archived_restart)
            archived_restart_files.append(archived_restart)
    batch_script = manifest_dir / "cgl_lf_stage_i.sbatch"
    manifest: dict[str, object] = {
        "schema_version": 3,
        "execution_epoch": EXECUTION_EPOCH,
        "state": "prepared",
        "prepared_utc": utc_now(),
        "project_root": str(root),
        "policy": {
            "account": ACCOUNT,
            "partition": PARTITION,
            "qos": PRODUCTION_QOS,
            "project_budget_node_hours": PROJECT_BUDGET_NODE_HOURS,
            "historical_debug_node_hours": HISTORICAL_DEBUG_NODE_HOURS,
            "historical_e01_stage_i_node_hours": HISTORICAL_E01_STAGE_I_NODE_HOURS,
            "stage_i_reserved_node_hours": CURRENT_STAGE_I_RESERVED_NODE_HOURS,
            "stage_i_authorization": (
                "frozen mapped Stage I matrix R02-R17 under sequential inspection"
            ),
            "atomic_submission_required": True,
        },
        "run": {
            "case_id": args.case_id,
            "case_name": case["name"],
            "segment": args.segment,
            "run_basename": f"{EXECUTION_EPOCH_SLUG}_{case['name']}_{args.segment}",
            "resolution": case["resolution"],
            "figure_roles": case["figure_roles"],
            "acceptance_criterion": args.acceptance_criterion,
        },
        "allocation": {
            "nodes": args.nodes,
            "requested_walltime": args.walltime,
            "requested_seconds": requested_seconds,
            "reserved_node_hours": segment_hours,
            "ranks_per_node": args.ranks_per_node,
            "cpus_per_task": args.cpus_per_task,
        },
        "command": {
            "production_utility": utility_provenance,
            "qualification_approval": qualification_approval,
            "source_dir": str(source_dir),
            "source_bundle": bundle_provenance,
            "input_revision": input_revision,
            "source_input_file": str(input_path),
            "input_file": str(archived_input),
            "input_sha256": sha256(archived_input),
            "matrix_file": str(archived_matrix),
            "matrix_sha256": sha256(archived_matrix),
            "executable": str(executable),
            "executable_revision": provenance["revision"],
            "executable_sha256": provenance["sha256"],
            "build_manifest": provenance["manifest_dir"],
            "source_restart_file": str(restart) if restart else None,
            "restart_file": str(archived_restart) if archived_restart else None,
            "restart_sha256": sha256(archived_restart) if archived_restart else None,
            "restart_files": [
                retained_file(path) for path in archived_restart_files
            ],
            "parent_segment": parent_segment,
            "allow_missing_restart_time_marker": (
                allow_missing_restart_time_marker
            ),
            "athena_walltime": args.athena_walltime,
            "overrides": args.override,
            "time_tlim_target": time_tlim_target,
        },
        "paths": {
            "run_dir": str(run_dir),
            "output_dir": str(output_dir),
            "environment_log": str(manifest_dir / "run_environment.txt"),
            "batch_script": str(batch_script),
            "slurm_log": str(paths["logs_slurm"] / "%x.%j.log"),
        },
    }
    script, script_sha256 = finalize_batch_script(
        generated_batch_script(manifest, manifest_path)
    )
    manifest["command"]["batch_script_sha256"] = script_sha256
    write_text(batch_script, script, mode=0o750)
    authenticate_prepared_execution(
        manifest, manifest_path, allow_legacy_local=offline_local_root
    )
    reservations = read_reservations(paths)
    reservations.append({
        "execution_epoch": EXECUTION_EPOCH,
        "manifest": str(manifest_path),
        "case_id": args.case_id,
        "case_name": case["name"],
        "segment": args.segment,
        "nodes": args.nodes,
        "requested_walltime": args.walltime,
        "reserved_node_hours": segment_hours,
        "execution_intent_sha256": execution_intent_sha256(manifest),
        "state": "prepared",
        "prepared_utc": manifest["prepared_utc"],
    })
    durable_transition(paths, "prepared", manifest_path, manifest, reservations)
    print(f"Prepared Stage I segment: {manifest_path}")
    print(f"Reserved node-hours: {segment_hours:.6f}")
    return manifest_path


def read_manifest(path: Path) -> dict[str, object]:
    """Read one segment preparation manifest."""

    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"invalid run manifest: {path}")
    return value


def verify_continuation_restart(
    restart: Path,
    allow_missing_restart_time_marker: bool = False,
) -> dict[str, object]:
    """Require an inspected production parent for a continuation restart."""

    parent_manifest_path = None
    for ancestor in restart.parents:
        candidate = ancestor / "manifest" / "prepared_run.json"
        if candidate.is_file():
            parent_manifest_path = candidate
            break
    if parent_manifest_path is None:
        raise ValueError("restart does not belong to a retained Stage I segment")
    parent = read_manifest(parent_manifest_path)
    require_current_epoch(parent, "continuation parent")
    parent_root = Path(str(parent.get("project_root", ""))).resolve()
    if (
        allow_missing_restart_time_marker
        and parent_root == DEFAULT_ROOT.expanduser().resolve()
    ):
        raise ValueError(
            "canonical continuation may not bypass restart physical-time markers"
        )
    authenticate_prepared_execution(
        parent, parent_manifest_path,
        allow_legacy_local=parent_root != DEFAULT_ROOT.expanduser().resolve(),
    )
    if parent_root == DEFAULT_ROOT.expanduser().resolve():
        parent_paths = layout(parent_root)
        require_existing_layout(parent_paths)
        require_reserved_execution_intent(
            reservation_for_manifest(
                read_reservations(parent_paths), parent_manifest_path
            ),
            parent,
            allow_legacy_local=False,
        )
    accounting = parent.get("accounting", {})
    inspection = parent.get("scientific_inspection", {})
    if (
        parent.get("state") != "recorded"
        or not isinstance(accounting, dict)
        or accounting.get("result") not in {"accepted", "clean_partial"}
        or not isinstance(inspection, dict)
        or inspection.get("clean_for_continuation") is not True
    ):
        raise ValueError("restart parent has not passed continuation inspection")
    revalidate_inspection_files(inspection, parent)
    terminal = inspection.get("terminal_restart")
    if (
        not isinstance(terminal, dict)
        or Path(str(terminal.get("path", ""))).resolve() != restart
        or terminal.get("sha256") != sha256(restart)
    ):
        raise ValueError("continuation must use the inspected terminal restart")
    revalidate_retained_product(terminal)
    restart_files = retained_product_paths(terminal)
    restart_time = restart_product_time(
        restart_files,
        allow_missing_marker=allow_missing_restart_time_marker,
    )
    try:
        final_time = float(inspection["final_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("continuation parent inspection lacks final time") from error
    if (
        restart_time is not None
        and abs(restart_time - final_time) > 1.0e-12
    ):
        raise ValueError(
            "continuation terminal restart physical time differs from inspection"
        )
    return {
        "execution_epoch": EXECUTION_EPOCH,
        "manifest": str(parent_manifest_path),
        "case_id": parent["run"]["case_id"],
        "segment": parent["run"]["segment"],
        "result": accounting["result"],
        "restart_sha256": terminal["sha256"],
        "restart_files": [str(path) for path in restart_files],
        "final_time": final_time,
        "restart_time": restart_time,
        "input_sha256": parent["command"]["input_sha256"],
        "executable_sha256": parent["command"]["executable_sha256"],
    }


def reservation_for_manifest(reservations: list[dict[str, object]],
                             manifest_path: Path) -> dict[str, object]:
    """Return exactly one reservation corresponding to a run manifest."""

    matches = [
        item for item in reservations
        if Path(str(item.get("manifest", ""))).resolve() == manifest_path.resolve()
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one reservation for {manifest_path}, found {len(matches)}"
        )
    if matches[0].get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("reservation does not belong to the current execution epoch")
    return matches[0]


def scheduler_environment() -> dict[str, str]:
    """Return inherited process state without caller-selected Slurm routing."""

    return {
        key: value for key, value in os.environ.items()
        if not key.startswith("SLURM_")
    }


def scheduler_account_name() -> str:
    """Return the account name bound to this process effective UID."""

    try:
        return pwd.getpwuid(os.geteuid()).pw_name
    except KeyError as error:
        raise ValueError(
            "effective UID has no account name; cannot query scheduler"
        ) from error


def production_queue_output(args: argparse.Namespace,
                            offline_local_root: bool) -> str:
    """Query all of this user's queued jobs before a root-writing submission."""

    fixture = getattr(args, "squeue_file", None)
    if fixture:
        if not offline_local_root:
            raise ValueError("--squeue-file is restricted to offline local roots")
        return Path(fixture).read_text(encoding="utf-8")
    if offline_local_root:
        raise ValueError("offline local-root submission requires --squeue-file")
    user = scheduler_account_name()
    try:
        return subprocess.run(
            [str(SQUEUE), "-h", "-u", user, "-o", "%i|%P|%T|%j"],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (FileNotFoundError, subprocess.CalledProcessError) as error:
        raise ValueError(
            "squeue is unavailable; refusing Stage I submission"
        ) from error


def validate_submission_fixture_options(args: argparse.Namespace,
                                        offline_local_root: bool) -> None:
    """Keep scheduler fixture injection and test bypass out of production."""

    if getattr(args, "squeue_file", None) and not offline_local_root:
        raise ValueError("--squeue-file is restricted to offline local roots")
    if getattr(args, "skip_slurm_test", False) and not offline_local_root:
        raise ValueError("--skip-slurm-test is restricted to offline local roots")


def scheduler_test_only_output(script: Path) -> str:
    """Run the live Slurm submission validator with trusted routing."""

    return subprocess.run(
        [str(SBATCH), "--test-only", str(script)],
        check=True,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
    ).stdout


def shared_root_campaign_conflicts(root: Path,
                                   allowed: set[str]) -> list[str]:
    """Return top-level CGL campaigns requiring an explicit overlap review."""

    conflicts = []
    for manifest_path in sorted(
        (root / "runs").glob("*/manifest/prepared_run.json")
    ):
        manifest = read_manifest(manifest_path)
        state = str(manifest.get("state", "")).lower()
        if state not in SHARED_ROOT_ACTIVE_STATES:
            continue
        campaign_id = str(
            manifest.get("campaign_id", manifest_path.parents[1].name)
        )
        if campaign_id not in allowed:
            conflicts.append(f"{campaign_id}|{state}|{manifest_path}")
    return conflicts


def require_existing_layout(paths: dict[str, Path]) -> None:
    """Require initialized stores without modifying them."""

    for key in ("ledger", "reservations"):
        if not paths[key].is_file():
            raise ValueError(f"Stage I store is missing: {paths[key]}")
    for key in ("runs", "transactions"):
        if not paths[key].is_dir():
            raise ValueError(f"Stage I directory is missing: {paths[key]}")
    read_ledger(paths)
    read_reservations(paths)


def require_reconciled_store_consistency(
    paths: dict[str, Path],
    *,
    allow_absent_qualification: bool = False,
    allow_invalid_qualification: bool = False,
    ignored_issue_prefixes: tuple[str, ...] = (),
) -> None:
    """Fail closed before canonical metadata or bundle mutations."""

    if paths["root"].resolve() != DEFAULT_ROOT.expanduser().resolve():
        return
    report = reconcile_report(paths["root"])
    issues = []
    for issue in report["issues"]:
        if (
            (
                allow_absent_qualification
                and issue.startswith(
                    "E03 corrected-build Frontier qualification is pending:"
                )
            )
            or (
                allow_invalid_qualification
                and issue.startswith(
                    "E03 corrected-build Frontier qualification is invalid:"
                )
            )
        ):
            continue
        if any(issue.startswith(prefix) for prefix in ignored_issue_prefixes):
            continue
        issues.append(issue)
    if issues:
        raise ValueError(
            "canonical E03 store reconciliation failed before mutation: "
            + "; ".join(issues)
        )


def submission_preflight(args: argparse.Namespace, manifest_path: Path,
                         manifest: dict[str, object],
                         run_slurm_test: bool,
                         ) -> tuple[dict[str, Path], Path, dict[str, object]]:
    """Authenticate and check one prepared segment immediately before submission."""

    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    validate_submission_fixture_options(args, offline_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    if manifest.get("state") != "prepared":
        raise ValueError("only a prepared segment can be submission-checked")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is unavailable for submission")
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    active = active_reservations(reservations)
    if active != [reservation]:
        raise ValueError("submission requires exactly one matching active reservation")
    authenticate_prepared_execution(
        manifest, manifest_path, allow_legacy_local=offline_local_root
    )
    if (
        not offline_local_root
        or isinstance(manifest.get("command", {}).get("overrides"), list)
    ):
        prepared_time_tlim_target(manifest)
    actual, reserved = reservation_usage(paths)
    if actual + reserved > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("active Stage I reservation exceeds its ceiling")
    lines = [line for line in production_queue_output(
        args, offline_local_root
    ).splitlines()
             if line.strip()]
    if lines:
        raise ValueError(
            "another user job is queued; review shared-root concurrency before "
            f"submitting {EXECUTION_EPOCH}: " + "; ".join(lines)
        )
    conflicts = shared_root_campaign_conflicts(
        root, set(getattr(args, "allow_shared_root_campaign", []))
    )
    if conflicts:
        raise ValueError(
            "shared-root campaign records require explicit review; pass "
            "--allow-shared-root-campaign only after confirming isolation: "
            + "; ".join(conflicts)
        )
    script = Path(str(manifest["paths"]["batch_script"])).resolve()
    slurm_test_outcome = "not requested"
    if (
        run_slurm_test
        and not offline_local_root
        and not getattr(args, "skip_slurm_test", False)
    ):
        print(scheduler_test_only_output(script).strip())
        slurm_test_outcome = "passed"
    elif run_slurm_test and offline_local_root:
        slurm_test_outcome = "offline local-root fixture"
    elif run_slurm_test and getattr(args, "skip_slurm_test", False):
        slurm_test_outcome = "operator skipped with --skip-slurm-test"
    audit = {
        "created_utc": utc_now(),
        "offline_local_root": offline_local_root,
        "skip_slurm_test": bool(getattr(args, "skip_slurm_test", False)),
        "slurm_test_only": slurm_test_outcome,
        "acknowledged_shared_root_campaigns": sorted(
            set(getattr(args, "allow_shared_root_campaign", []))
        ),
    }
    return paths, script, audit


def check_submit(args: argparse.Namespace) -> int:
    """Print a read-only, explicitly non-authoritative submission preview."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    _, script, _ = submission_preflight(
        args, manifest_path, manifest, run_slurm_test=True
    )
    print("Non-authoritative preview passed. Re-run the checks atomically with:")
    acknowledgements = "".join(
        " " + shlex.quote(f"--allow-shared-root-campaign={campaign}")
        for campaign in sorted(set(getattr(args, "allow_shared_root_campaign", [])))
    )
    print(
        f"  python3 {shlex.quote(str(Path(__file__).resolve()))} submit "
        f"--manifest {shlex.quote(str(manifest_path))}{acknowledgements}"
    )
    print(f"Prepared script: {script}")
    return 0


def parse_sbatch_job_id(output: str) -> str:
    """Return the top-level numeric ID from ``sbatch --parsable`` output."""

    match = re.fullmatch(r"\s*([1-9][0-9]*)(?:;[^\s;]+)?\s*", output)
    if match is None:
        raise ValueError(f"sbatch did not return one numeric job ID: {output!r}")
    return require_numeric_job_id(match.group(1))


def scheduler_submit_time_evidence(value: object,
                                   transaction: dict[str, object]) -> str:
    """Require scheduler submit time within five minutes of the sbatch barrier."""

    submitted = parse_utc_timestamp(value, "scheduler submit time")
    barrier = parse_utc_timestamp(
        transaction.get("created_utc"), "submission ambiguity barrier time"
    )
    if abs((submitted - barrier).total_seconds()) > (
        SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS
    ):
        raise ValueError(
            "scheduler recovery submit time is outside the symmetric "
            f"{SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS}-second ambiguity "
            "barrier tolerance"
        )
    return submitted.isoformat()


def verify_recovered_scheduler_job(manifest: dict[str, object], job_id: str,
                                   offline_local_root: bool,
                                   transaction: dict[str, object],
                                   ) -> dict[str, object]:
    """Require an operator-supplied recovery ID to name the prepared job."""

    require_numeric_job_id(job_id)
    if offline_local_root:
        return {
            "mode": "offline-local fixture",
            "checked_utc": utc_now(),
            "job_id": job_id,
        }
    expected_name = expected_job_name(manifest)
    expected_script = Path(str(manifest["paths"]["batch_script"])).resolve()
    control = subprocess.run(
        [str(SCONTROL), "show", "job", "-o", job_id],
        check=False,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
    )
    if control.returncode == 0:
        fields = {
            key: value
            for key, value in (
                token.split("=", 1)
                for token in shlex.split(control.stdout)
                if "=" in token
            )
        }
        command = fields.get("Command")
        if (
            fields.get("JobId") == job_id
            and fields.get("JobName") == expected_name
            and str(fields.get("Account", "")).casefold() == ACCOUNT.casefold()
            and fields.get("Partition") == PARTITION
            and (
                command in {None, "", "(null)", "N/A"}
                or Path(command).expanduser().resolve() == expected_script
            )
        ):
            submitted = scheduler_submit_time_evidence(
                fields.get("SubmitTime"), transaction
            )
            return {
                "mode": "scontrol",
                "checked_utc": utc_now(),
                "job_id": job_id,
                "job_name": expected_name,
                "account": fields["Account"],
                "partition": fields["Partition"],
                "command": command,
                "submit_time": submitted,
            }
    accounting = subprocess.run(
        [
            str(SACCT), "-X", "-j", job_id,
            "--format=JobIDRaw,JobName,Account,Partition,Submit", "-n", "-P",
        ],
        check=True,
        capture_output=True,
        text=True,
        env=scheduler_environment(),
    ).stdout
    rows = [
        row for row in csv.reader(accounting.splitlines(), delimiter="|")
        if row and row[0] == job_id
    ]
    if (
        len(rows) != 1
        or len(rows[0]) < 5
        or rows[0][1] != expected_name
        or rows[0][2].casefold() != ACCOUNT.casefold()
        or rows[0][3] != PARTITION
    ):
        raise ValueError("scheduler recovery job does not match the prepared segment")
    submitted = scheduler_submit_time_evidence(rows[0][4], transaction)
    return {
        "mode": "sacct",
        "checked_utc": utc_now(),
        "job_id": job_id,
        "job_name": expected_name,
        "account": rows[0][2],
        "partition": rows[0][3],
        "submit_time": submitted,
    }


def scheduler_absence_evidence(args: argparse.Namespace,
                               manifest: dict[str, object],
                               transaction: dict[str, object],
                               offline_local_root: bool,
                               ) -> dict[str, object]:
    """Retain scheduler-side evidence before clearing an ambiguous submission."""

    expected_name = expected_job_name(manifest)
    checked_utc = utc_now()
    fixture = getattr(args, "scheduler_absence_evidence_file", None)
    if fixture:
        if not offline_local_root:
            raise ValueError(
                "--scheduler-absence-evidence-file is restricted to offline validation"
            )
        value = json.loads(Path(fixture).read_text(encoding="utf-8"))
        if not isinstance(value, dict) or value.get("absent") is not True:
            raise ValueError("offline scheduler absence fixture must assert absent=true")
        return {
            "mode": "offline-local fixture",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
            "fixture": value,
        }
    if offline_local_root:
        return {
            "mode": "offline-local operator confirmation",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
        }
    user = scheduler_account_name()
    barrier = parse_utc_timestamp(
        transaction.get("created_utc"), "submission ambiguity barrier time"
    )
    scheduler_start = (
        barrier - timedelta(seconds=SCHEDULER_SUBMIT_BARRIER_TOLERANCE_SECONDS)
    ).strftime(
        "%Y-%m-%dT%H:%M:%S"
    )
    squeue_command = [
        str(SQUEUE), "-h", "-u", user, "-o", "%i|%j|%a|%P|%V|%o",
    ]
    sacct_command = [
        str(SACCT), "-X", "-S", scheduler_start,
        "--format=JobIDRaw,JobName,Account,Partition,Submit", "-n", "-P",
    ]
    try:
        queued = subprocess.run(
            squeue_command,
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
        accounted = subprocess.run(
            sacct_command,
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        break_glass = getattr(args, "break_glass_clear_evidence", "")
        if (
            not getattr(args, "confirm_break_glass_clear", False)
            or not str(break_glass).strip()
        ):
            raise ValueError(
                "scheduler absence query failed; explicit break-glass evidence "
                "and --confirm-break-glass-clear are required"
            ) from error
        return {
            "mode": "break-glass after scheduler query failure",
            "checked_utc": checked_utc,
            "expected_job_name": expected_name,
            "ambiguity_created_utc": transaction["created_utc"],
            "operator_evidence": break_glass,
            "query_error": str(error),
        }
    if (
        scheduler_output_contains_job(queued, expected_name)
        or scheduler_output_contains_job(accounted, expected_name)
    ):
        raise ValueError("scheduler still reports a matching ambiguous submission")
    return {
        "mode": "live scheduler absence query",
        "checked_utc": checked_utc,
        "expected_job_name": expected_name,
        "ambiguity_created_utc": transaction["created_utc"],
        "squeue_command": squeue_command,
        "squeue_output": queued,
        "sacct_command": sacct_command,
        "sacct_output": accounted,
    }


@locked_manifest_action
def submit(args: argparse.Namespace) -> int:
    """Atomically authenticate, submit, and retain one scheduler job ID."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    output_file = getattr(args, "sbatch_output_file", None)
    if output_file and not offline_local_root:
        raise ValueError("--sbatch-output-file is restricted to offline validation")
    paths, script, audit = submission_preflight(
        args, manifest_path, manifest, run_slurm_test=True
    )
    transaction_path = write_submit_pending_transaction(
        paths, manifest_path, audit
    )
    if output_file:
        output = Path(output_file).read_text(encoding="utf-8")
    else:
        output = subprocess.run(
            [str(SBATCH), "--parsable", str(script)],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    job_id = parse_sbatch_job_id(output)
    finish_submit_transaction(
        paths, transaction_path, manifest_path, manifest,
        read_reservations(paths), job_id,
    )
    print(f"Submitted Stage I job {job_id}: {script}")
    return 0


@locked_manifest_action
def mark_submitted(args: argparse.Namespace) -> int:
    """Attach a fake scheduler ID only for legacy offline validation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    if not is_offline_local_root(root, args.allow_local_root):
        raise ValueError("production submissions must use the atomic submit action")
    paths = initialize(root)
    require_numeric_job_id(args.job_id)
    if manifest.get("state") != "prepared":
        raise ValueError("manifest is not in prepared state")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation is not in prepared state")
    transaction_path = write_submit_pending_transaction(
        paths, manifest_path, {
            "created_utc": utc_now(),
            "offline_local_root": True,
            "legacy_mark_submitted": True,
            "skip_slurm_test": True,
            "slurm_test_only": "offline local-root fixture",
            "acknowledged_shared_root_campaigns": [],
        }
    )
    finish_submit_transaction(
        paths, transaction_path, manifest_path, manifest, reservations, args.job_id
    )
    print(f"Marked Stage I job {args.job_id} submitted.")
    return 0


@locked_manifest_action
def recover_submit(args: argparse.Namespace) -> int:
    """Resolve an ambiguous sbatch boundary with its retained scheduler ID."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    transaction_path = submit_pending_transaction(paths, manifest_path)
    transaction = read_transaction(paths, transaction_path)
    evidence = verify_recovered_scheduler_job(
        manifest, args.job_id, offline_local_root, transaction
    )
    retained_evidence = manifest.setdefault("scheduler_recovery_evidence", [])
    if not isinstance(retained_evidence, list):
        raise ValueError("prepared manifest has invalid scheduler recovery evidence")
    retained_evidence.append(evidence)
    finish_submit_transaction(
        paths, transaction_path, manifest_path, manifest,
        read_reservations(paths), args.job_id,
    )
    print(f"Recovered submitted Stage I job {args.job_id}.")
    return 0


@locked_manifest_action
def clear_submit_pending(args: argparse.Namespace) -> int:
    """Clear an ambiguous submit barrier after confirming no job was launched."""

    if not args.confirm_no_job_submitted:
        raise ValueError("--confirm-no-job-submitted is required")
    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    transaction_path = submit_pending_transaction(paths, manifest_path)
    transaction = read_transaction(paths, transaction_path)
    if transaction.get("prepared_manifest_sha256") != sha256(manifest_path):
        raise ValueError("prepared manifest changed after the sbatch boundary")
    if transaction.get("prior_reservations_sha256") != sha256(paths["reservations"]):
        raise ValueError("reservation store changed after the sbatch boundary")
    absence_evidence = scheduler_absence_evidence(
        args, manifest, transaction, offline_local_root
    )
    notes = manifest.setdefault("submission_recovery_notes", [])
    if not isinstance(notes, list):
        raise ValueError("prepared manifest has invalid submission recovery notes")
    notes.append({
        "cleared_utc": utc_now(),
        "notes": args.notes,
        "outcome": "operator confirmed no scheduler job was submitted",
        "scheduler_absence_evidence": absence_evidence,
    })
    transaction.update({
        "kind": "submit_cleared",
        "manifest": manifest,
        "reservations": read_reservations(paths),
        "ledger_row": None,
        "recovery_notes": args.notes,
        "scheduler_absence_evidence": absence_evidence,
    })
    write_json(transaction_path, transaction)
    apply_transaction(paths, transaction_path)
    print(f"Cleared ambiguous submission barrier: {manifest_path}")
    return 0


@locked_root_action
def recover_transactions(args: argparse.Namespace) -> int:
    """Replay deterministic journals while preserving ambiguous submissions."""

    root = require_root(Path(args.root), args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    recovered = 0
    for transaction_path in pending_transaction_paths(paths):
        apply_transaction(paths, transaction_path)
        recovered += 1
    print(f"Recovered {recovered} deterministic Stage I transaction(s).")
    return 0


def parse_history(path: Path) -> dict[str, list[float]]:
    """Read one AthenaK history file keyed by its labeled columns."""

    labels: list[str] = []
    rows: list[list[float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            found = re.findall(r"\[\d+\]=([A-Za-z0-9_]+)", line)
            if found:
                labels = found
            continue
        if line.strip():
            rows.append([float(value) for value in line.split()])
    if not labels or not rows:
        raise ValueError(f"history file is missing labels or data: {path}")
    if any(len(row) != len(labels) for row in rows):
        raise ValueError(f"history row width does not match labels: {path}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def retained_file(path: Path) -> dict[str, object]:
    """Describe an output retained for acceptance review."""

    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def output_product_groups(directory: Path, pattern: str,
                          expected_ranks: int | None = None) -> list[list[Path]]:
    """Collect shared outputs or complete rank-local product sets."""

    groups = [[path] for path in sorted(directory.glob(pattern))]
    rank0_dir = directory / "rank_00000000"
    if rank0_dir.is_dir():
        for rank0 in sorted(rank0_dir.glob(pattern)):
            rank_files = sorted(directory.glob(f"rank_*/{rank0.name}"))
            expected_names = [
                f"rank_{rank:08d}" for rank in range(len(rank_files))
            ]
            if [path.parent.name for path in rank_files] != expected_names:
                raise ValueError(f"rank-local output set is not contiguous: {rank0}")
            if expected_ranks is not None and len(rank_files) != expected_ranks:
                raise ValueError(
                    f"rank-local output set has {len(rank_files)} files; "
                    f"expected {expected_ranks}: {rank0}"
                )
            groups.append(rank_files)
    return sorted(groups, key=lambda group: str(group[0]))


def retained_product(group: list[Path]) -> dict[str, object]:
    """Describe one shared or rank-local retained output product."""

    representative = retained_file(group[0])
    if group[0].parent.name.startswith("rank_"):
        representative["storage"] = "per_rank"
        representative["rank_files"] = [retained_file(path) for path in group]
    else:
        representative["storage"] = "shared_mpiio"
    return representative


def retained_product_paths(record: dict[str, object]) -> list[Path]:
    """Return every file belonging to an inspected retained product."""

    rank_files = record.get("rank_files")
    if isinstance(rank_files, list):
        if not rank_files or not all(isinstance(item, dict) for item in rank_files):
            raise ValueError("inspection-retained rank-local product is invalid")
        return [Path(str(item["path"])).resolve() for item in rank_files]
    return [Path(str(record["path"])).resolve()]


def revalidate_retained_file(record: object,
                             label: str = "inspection-retained file") -> None:
    """Require a retained file to preserve size and digest."""

    if not isinstance(record, dict):
        raise ValueError(f"{label} record is invalid")
    path = Path(str(record.get("path", ""))).resolve()
    if not path.is_file():
        raise ValueError(f"{label} is missing: {path}")
    if record.get("size_bytes") != path.stat().st_size:
        raise ValueError(f"{label} size has changed: {path}")
    if record.get("sha256") != sha256(path):
        raise ValueError(f"{label} checksum has changed: {path}")


def revalidate_retained_product(record: object) -> None:
    """Require every member of one inspected output product to be unchanged."""

    revalidate_retained_file(record)
    if not isinstance(record, dict):
        raise ValueError("inspection-retained product record is invalid")
    rank_files = record.get("rank_files")
    if rank_files is None:
        return
    if not isinstance(rank_files, list) or not rank_files:
        raise ValueError("inspection-retained rank-local product is invalid")
    for rank_file in rank_files:
        revalidate_retained_file(rank_file)


def recorded_product_groups(records: object, label: str) -> list[list[Path]]:
    """Return the exact inspected path groups for one retained output class."""

    if not isinstance(records, list):
        raise ValueError(f"segment inspection lacks retained {label}")
    groups = []
    for record in records:
        if not isinstance(record, dict):
            raise ValueError(f"segment inspection has invalid retained {label}")
        groups.append(retained_product_paths(record))
    return groups


def output_group_signature(groups: list[list[Path]]) -> list[tuple[str, ...]]:
    """Return a stable path-only signature for grouped output inventory."""

    return sorted(tuple(str(path.resolve()) for path in group) for group in groups)


def require_complete_output_inventory(directory: Path, pattern: str,
                                      groups: list[list[Path]],
                                      label: str) -> None:
    """Reject product files omitted by shared or rank-local grouping."""

    actual = sorted(
        path.resolve() for path in directory.rglob(pattern) if path.is_file()
    )
    grouped = sorted(path.resolve() for group in groups for path in group)
    if actual != grouped:
        raise ValueError(f"{label} output inventory contains ungrouped files")


def revalidate_inspection_inventory(manifest: dict[str, object],
                                    inspection: dict[str, object]) -> None:
    """Reject additions, removals, or regrouping after formal inspection."""

    output_dir = Path(str(manifest["paths"]["output_dir"])).resolve()
    expected_ranks = (
        int(manifest["allocation"]["nodes"])
        * int(manifest["allocation"].get("ranks_per_node", 1))
    )
    histories = {
        "mhd_history": sorted(path.resolve() for path in output_dir.glob("*.mhd.hst")),
        "user_history": sorted(path.resolve() for path in output_dir.glob("*.user.hst")),
    }
    for key, actual in histories.items():
        record = inspection.get(key)
        if not isinstance(record, dict):
            raise ValueError(f"segment inspection lacks retained {key}")
        expected = [Path(str(record.get("path", ""))).resolve()]
        if actual != expected:
            raise ValueError(f"inspection-retained {key} inventory has changed")
    for key, directory, pattern in (
        ("snapshots", output_dir / "bin", "*.bin"),
        ("restarts", output_dir / "rst", "*.rst"),
    ):
        expected = recorded_product_groups(inspection.get(key), key)
        actual = output_product_groups(directory, pattern, expected_ranks)
        require_complete_output_inventory(directory, pattern, actual, key)
        if output_group_signature(actual) != output_group_signature(expected):
            raise ValueError(f"inspection-retained {key} inventory has changed")


def revalidate_inspection_restart_times(inspection: dict[str, object],
                                        allow_legacy_local: bool) -> None:
    """Reparse retained restart markers and bind the terminal product to final time."""

    if allow_legacy_local and "restart_times" not in inspection:
        return
    groups = recorded_product_groups(inspection.get("restarts"), "restarts")
    bypass = (
        allow_legacy_local
        and inspection.get("restart_time_marker_bypass") is True
    )
    parsed = [
        restart_product_time(group, allow_missing_marker=bypass)
        for group in groups
    ]
    retained = inspection.get("restart_times")
    if not isinstance(retained, list) or len(retained) != len(parsed):
        raise ValueError("segment inspection restart-time evidence is incomplete")
    for expected, actual in zip(retained, parsed):
        if expected is None or actual is None:
            if not bypass or expected is not None or actual is not None:
                raise ValueError("segment inspection restart-time bypass is inconsistent")
        else:
            try:
                difference = abs(float(expected) - actual)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    "segment inspection restart-time evidence is invalid"
                ) from error
            if difference > 1.0e-12:
                raise ValueError("segment inspection restart-time evidence has changed")
    if bypass:
        return
    try:
        final_time = float(inspection["final_time"])
        terminal_time = float(inspection["terminal_restart_time"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError("segment inspection lacks terminal restart time") from error
    matches = [
        index for index, value in enumerate(parsed)
        if value is not None and abs(value - final_time) <= 1.0e-10
    ]
    if len(matches) != 1 or abs(terminal_time - final_time) > 1.0e-10:
        raise ValueError("segment inspection terminal restart time differs from final time")
    terminal = inspection.get("terminal_restart")
    if not isinstance(terminal, dict):
        raise ValueError("segment inspection lacks terminal restart product")
    if (
        output_group_signature([retained_product_paths(terminal)])
        != output_group_signature([groups[matches[0]]])
    ):
        raise ValueError("segment inspection terminal restart product is inconsistent")


def revalidate_inspection_files(inspection: dict[str, object],
                                manifest: dict[str, object] | None = None) -> None:
    """Recheck every file retained by an accepted or clean-partial inspection."""

    for key in ("mhd_history", "user_history"):
        revalidate_retained_file(inspection.get(key))
    for key in ("snapshots", "restarts"):
        records = inspection.get(key)
        if not isinstance(records, list):
            raise ValueError(f"segment inspection lacks retained {key}")
        for record in records:
            revalidate_retained_product(record)
    allow_legacy_local = (
        manifest is not None
        and Path(str(manifest.get("project_root", ""))).resolve()
        != DEFAULT_ROOT.expanduser().resolve()
    )
    revalidate_inspection_restart_times(inspection, allow_legacy_local)
    if manifest is not None:
        revalidate_inspection_inventory(manifest, inspection)


def merge_history_files(sources: list[Path], destination: Path) -> None:
    """Merge restart-segment histories while removing repeated boundary rows."""

    reference_labels: list[str] | None = None
    header: list[str] = []
    retained_rows: list[str] = []
    last_time = float("-inf")
    for source_index, source in enumerate(sources):
        source_header: list[str] = []
        source_labels: list[str] = []
        source_rows: list[str] = []
        for line in source.read_text(encoding="utf-8").splitlines():
            if line.startswith("#"):
                source_header.append(line)
                found = re.findall(r"\[\d+\]=([A-Za-z0-9_]+)", line)
                if found:
                    source_labels = found
            elif line.strip():
                source_rows.append(line)
        if not source_labels:
            raise ValueError(f"history file has no labeled header: {source}")
        if reference_labels is None:
            reference_labels = source_labels
            header = source_header
        elif source_labels != reference_labels:
            raise ValueError("restart-segment history columns do not match")
        for line in source_rows:
            time = float(line.split()[0])
            if time > last_time + 1.0e-12:
                retained_rows.append(line)
                last_time = time
        if source_index == 0 and not retained_rows:
            raise ValueError(f"history file has no retained rows: {source}")
    write_text(destination, "\n".join([*header, *retained_rows]) + "\n")


def binary_snapshot_time(path: Path) -> float:
    """Read a retained Athena binary snapshot time without loading field data."""

    with path.open("rb") as stream:
        code_header = stream.readline().split()
        if not code_header or code_header[0] != b"Athena":
            raise ValueError(f"invalid Athena binary snapshot: {path}")
        pheader_count = int(stream.readline().split(b"=")[-1])
        values: dict[str, str] = {}
        for _ in range(pheader_count - 1):
            key, value = [
                token.strip()
                for token in stream.readline().decode("utf-8").split("=", 1)
            ]
            values[key] = value
    if "time" not in values:
        raise ValueError(f"binary snapshot has no time field: {path}")
    return float(values["time"])


def binary_product_time(group: list[Path]) -> float:
    """Require every rank-local member of a snapshot to record one time."""

    times = [binary_snapshot_time(path) for path in group]
    if any(abs(time - times[0]) > 1.0e-12 for time in times[1:]):
        raise ValueError(f"rank-local snapshot times disagree: {group[0]}")
    return times[0]


def analysis_model_choices(input_path: Path) -> dict[str, str]:
    """Use the workflow's interpretation mapping for a production bundle."""

    workflow_path = ROOT_DIR / "scripts/cgl_lf_workflow.py"
    spec = importlib.util.spec_from_file_location(
        "_cgl_lf_workflow_stage_i", workflow_path
    )
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load workflow model mapping: {workflow_path}")
    workflow = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = workflow
    spec.loader.exec_module(workflow)
    return workflow.model_choices(input_path.read_text(encoding="utf-8"), [])


@locked_manifest_action
def inspect_segment(args: argparse.Namespace) -> int:
    """Inspect terminal-time output evidence before scientific acceptance."""

    manifest_path = Path(args.manifest).expanduser().resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "submitted segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    allow_missing_restart_time_marker = getattr(
        args, "allow_missing_restart_time_marker", False
    )
    if allow_missing_restart_time_marker and not offline_local_root:
        raise ValueError(
            "--allow-missing-restart-time-marker is restricted to offline validation"
        )
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    reservation = reservation_for_manifest(read_reservations(paths), manifest_path)
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    authenticate_prepared_execution(
        manifest, manifest_path,
        allow_legacy_local=offline_local_root,
    )
    if manifest.get("state") not in {"submitted", "recorded"}:
        raise ValueError("only a submitted or recorded segment may be inspected")
    required_time = float(args.required_time)
    prepared_target = prepared_time_tlim_target(manifest)
    if (
        not math.isfinite(required_time)
        or abs(required_time - prepared_target) > 1.0e-12
    ):
        raise ValueError(
            "--required-time must match the prepared time/tlim target "
            f"{prepared_target:.12g}"
        )
    output_dir = Path(str(manifest["paths"]["output_dir"])).resolve()
    mhd_histories = sorted(output_dir.glob("*.mhd.hst"))
    user_histories = sorted(output_dir.glob("*.user.hst"))
    expected_ranks = (
        int(manifest["allocation"]["nodes"])
        * int(manifest["allocation"].get("ranks_per_node", 1))
    )
    snapshots = output_product_groups(
        output_dir / "bin", "*.bin", expected_ranks
    )
    restarts = output_product_groups(
        output_dir / "rst", "*.rst", expected_ranks
    )
    require_complete_output_inventory(output_dir / "bin", "*.bin", snapshots, "snapshot")
    require_complete_output_inventory(output_dir / "rst", "*.rst", restarts, "restart")
    if len(mhd_histories) != 1 or len(user_histories) != 1:
        raise ValueError("segment must retain exactly one MHD and one user history")
    history = parse_history(mhd_histories[0])
    if "time" not in history:
        raise ValueError("MHD history does not contain a time column")
    missing = [
        label for label in STRICT_LF_FAILURE_COLUMNS
        if label not in history
    ]
    if missing:
        raise ValueError(f"MHD history lacks strict safety counters: {missing}")
    final_time = history["time"][-1]
    maximum_failure_counts = {
        label: max(history[label])
        for label in STRICT_LF_FAILURE_COLUMNS
    }
    snapshot_times = [binary_product_time(group) for group in snapshots]
    restart_times = [
        restart_product_time(
            group,
            allow_missing_marker=allow_missing_restart_time_marker,
        )
        for group in restarts
    ]
    restart_records = [retained_product(group) for group in restarts]
    restart_markers_verified = bool(restart_times) and all(
        value is not None for value in restart_times
    )
    terminal_restart_matches = [
        index for index, value in enumerate(restart_times)
        if value is not None and abs(value - final_time) <= 1.0e-10
    ]
    if restart_markers_verified:
        if len(terminal_restart_matches) != 1:
            raise ValueError(
                "segment must retain exactly one restart product whose explicit "
                "physical time matches the inspected final history time"
            )
        terminal_restart = restart_records[terminal_restart_matches[0]]
        terminal_restart_time = restart_times[terminal_restart_matches[0]]
    elif allow_missing_restart_time_marker and restart_records:
        terminal_restart = restart_records[-1]
        terminal_restart_time = None
    else:
        terminal_restart = None
        terminal_restart_time = None
    checks = {
        "required_time_reached": final_time >= required_time - 1.0e-10,
        "strict_lf_failure_counters_zero": all(
            value == 0.0 for value in maximum_failure_counts.values()
        ),
        "snapshots_retained": bool(snapshots),
        "terminal_snapshot_retained": bool(snapshot_times)
        and max(snapshot_times) >= final_time - 1.0e-10,
        "restart_retained": bool(restarts),
        "terminal_restart_physical_time_matches_final": (
            restart_markers_verified and len(terminal_restart_matches) == 1
        ) or (
            offline_local_root and allow_missing_restart_time_marker
            and bool(restart_records)
        ),
    }
    accepted = all(checks.values())
    clean_for_continuation = all(
        checks[key]
        for key in (
            "strict_lf_failure_counters_zero",
            "snapshots_retained",
            "terminal_snapshot_retained",
            "restart_retained",
            "terminal_restart_physical_time_matches_final",
        )
    )
    inspection: dict[str, object] = {
        "schema_version": 3,
        "execution_epoch": EXECUTION_EPOCH,
        "inspected_utc": utc_now(),
        "manifest": str(manifest_path),
        "job_id": manifest.get("job_id"),
        "case_id": manifest["run"]["case_id"],
        "segment": manifest["run"]["segment"],
        "required_time": required_time,
        "final_time": final_time,
        "maximum_strict_failure_counts": maximum_failure_counts,
        "checks": checks,
        "accepted": accepted,
        "clean_for_continuation": clean_for_continuation,
        "mhd_history": retained_file(mhd_histories[0]),
        "user_history": retained_file(user_histories[0]),
        "snapshots": [retained_product(group) for group in snapshots],
        "snapshot_times": snapshot_times,
        "restarts": restart_records,
        "restart_times": restart_times,
        "terminal_restart": terminal_restart,
        "terminal_restart_time": terminal_restart_time,
        "restart_time_marker_bypass": (
            offline_local_root
            and allow_missing_restart_time_marker
            and not restart_markers_verified
        ),
    }
    if "lf_hwproj" in history:
        inspection["final_hardwall_projection_count"] = history["lf_hwproj"][-1]
    inspection_path = manifest_path.parent / "segment_inspection.json"
    write_json(inspection_path, inspection)
    status = (
        "accepted"
        if accepted
        else "clean partial" if clean_for_continuation else "not accepted"
    )
    print(
        f"Segment inspection {status}: final_time={final_time:.12g}, "
        f"required_time={required_time:.12g}"
    )
    print(f"Wrote {inspection_path}")
    return 0 if accepted else 1


def sacct_output(args: argparse.Namespace, paths: dict[str, Path],
                 allow_fixture: bool) -> str:
    """Read or query a top-level Slurm allocation record."""

    if args.sacct_file:
        if not allow_fixture:
            raise ValueError("--sacct-file is allowed only for offline local roots")
        output = Path(args.sacct_file).read_text(encoding="utf-8")
    else:
        output = subprocess.run(
            [
                str(SACCT), "-X", "-j", args.job_id,
                "--format=JobIDRaw,JobName,State,ExitCode,AllocNodes,"
                "ElapsedRaw,Submit,End", "-n", "-P",
            ],
            check=True,
            capture_output=True,
            text=True,
            env=scheduler_environment(),
        ).stdout
    write_text(paths["accounting"] / f"{args.job_id}.stage_i.sacct.txt", output)
    return output


def parse_sacct(output: str, job_id: str) -> dict[str, str]:
    """Select a completed top-level allocation row."""

    require_numeric_job_id(job_id)
    records: list[list[str]] = []
    for row in csv.reader(output.splitlines(), delimiter="|"):
        if row and row[0] == job_id:
            while row and not row[-1]:
                row.pop()
            records.append(row)
    if len(records) != 1 or len(records[0]) != 8:
        raise ValueError(f"expected one eight-column allocation row for {job_id}")
    keys = (
        "job_id", "job_name", "state", "exit_code", "nodes",
        "elapsed_seconds", "submitted_utc", "completed_utc",
    )
    result = dict(zip(keys, records[0]))
    state = result["state"].split()[0].split("+")[0]
    if state in NONTERMINAL_STATES:
        raise ValueError(f"job {job_id} is not complete: {result['state']}")
    result["state"] = state
    return result


def require_recorded_scheduler_evidence(paths: dict[str, Path],
                                        row: dict[str, object],
                                        manifest: dict[str, object]) -> None:
    """Bind canonical recorded accounting to its archived Slurm allocation row."""

    job_id = str(row["job_id"])
    evidence = paths["accounting"] / f"{job_id}.stage_i.sacct.txt"
    if not evidence.is_file():
        raise ValueError(f"recorded transaction lacks scheduler evidence: {evidence}")
    sacct = parse_sacct(evidence.read_text(encoding="utf-8"), job_id)
    if sacct["job_name"] != expected_job_name(manifest):
        raise ValueError("recorded scheduler evidence job name differs from manifest")
    expected = {
        "job_id": row["job_id"],
        "state": row["state"],
        "exit_code": row["exit_code"],
        "nodes": row["nodes"],
        "elapsed_seconds": row["elapsed_seconds"],
        "submitted_utc": row["submitted_utc"],
        "completed_utc": row["completed_utc"],
    }
    for key, value in expected.items():
        if sacct[key] != str(value):
            raise ValueError(
                f"recorded scheduler evidence {key} differs from ledger"
            )


@locked_manifest_action
def record(args: argparse.Namespace) -> int:
    """Account a completed segment and release its reservation."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "submitted segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    paths = layout(root)
    require_existing_layout(paths)
    require_no_pending_transactions(paths)
    require_no_orphaned_segment_runs(paths)
    require_reconciled_store_consistency(paths)
    require_numeric_job_id(args.job_id)
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    require_reserved_execution_intent(
        reservation, manifest, allow_legacy_local=offline_local_root
    )
    authenticate_prepared_execution(
        manifest, manifest_path,
        allow_legacy_local=offline_local_root,
    )
    if manifest.get("state") != "submitted":
        raise ValueError("only a submitted segment can be accounted")
    if str(manifest.get("job_id")) != args.job_id:
        raise ValueError("job ID does not match the submitted segment")
    ledger = read_ledger(paths)
    if any(row["job_id"] == args.job_id for row in ledger):
        raise ValueError(f"job {args.job_id} is already accounted")
    sacct = parse_sacct(
        sacct_output(args, paths, allow_fixture=offline_local_root), args.job_id
    )
    if sacct["job_name"] != expected_job_name(manifest):
        raise ValueError(
            "sacct job name does not match the prepared segment: "
            f"{sacct['job_name']!r}"
        )
    nodes = int(sacct["nodes"])
    if nodes != int(manifest["allocation"]["nodes"]):
        raise ValueError("allocated nodes differ from the prepared reservation")
    inspection = None
    if args.result in {"accepted", "clean_partial"}:
        if sacct["state"] != "COMPLETED" or sacct["exit_code"] != "0:0":
            raise ValueError(
                "scientific continuation output requires a clean COMPLETED job"
            )
        inspection_path = manifest_path.parent / "segment_inspection.json"
        if not inspection_path.is_file():
            raise ValueError(
                "scientific continuation output requires inspect-segment evidence"
            )
        inspection = json.loads(inspection_path.read_text(encoding="utf-8"))
        if (
            not isinstance(inspection, dict)
            or inspection.get("execution_epoch") != EXECUTION_EPOCH
            or inspection.get("job_id") != args.job_id
            or Path(str(inspection.get("manifest", ""))).resolve()
            != manifest_path
        ):
            raise ValueError("segment inspection does not match this submitted job")
        try:
            inspected_target = float(inspection["required_time"])
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("segment inspection lacks its prepared target") from error
        if (
            not math.isfinite(inspected_target)
            or abs(inspected_target - prepared_time_tlim_target(manifest)) > 1.0e-12
        ):
            raise ValueError("segment inspection target differs from preparation")
        if args.result == "accepted" and inspection.get("accepted") is not True:
            raise ValueError("segment inspection does not accept this submitted job")
        if (
            args.result == "clean_partial"
            and (
                inspection.get("accepted") is True
                or inspection.get("clean_for_continuation") is not True
            )
        ):
            raise ValueError("clean_partial requires clean output short of its target")
        revalidate_inspection_files(inspection, manifest)
    actual = node_hours(nodes, int(sacct["elapsed_seconds"]))
    cumulative = sum(float(row["actual_node_hours"]) for row in ledger) + actual
    if cumulative > CURRENT_STAGE_I_RESERVED_NODE_HOURS:
        raise ValueError("actual use exceeds the Stage I reservation")
    if cumulative > PROJECT_BUDGET_NODE_HOURS:
        raise ValueError("actual use exceeds the incremental project ceiling")
    command = manifest["command"]
    run = manifest["run"]
    allocation = manifest["allocation"]
    row = {
        "execution_epoch": EXECUTION_EPOCH,
        "job_id": args.job_id,
        "submitted_utc": sacct["submitted_utc"],
        "completed_utc": sacct["completed_utc"],
        "case_id": run["case_id"],
        "case_name": run["case_name"],
        "segment": run["segment"],
        "state": sacct["state"],
        "exit_code": sacct["exit_code"],
        "nodes": str(nodes),
        "requested_walltime": allocation["requested_walltime"],
        "elapsed_seconds": sacct["elapsed_seconds"],
        "reserved_node_hours": f"{float(allocation['reserved_node_hours']):.6f}",
        "actual_node_hours": f"{actual:.6f}",
        "cumulative_stage_i_node_hours": f"{cumulative:.6f}",
        "executable_revision": command["executable_revision"],
        "executable_sha256": command["executable_sha256"],
        "input_revision": command["input_revision"],
        "input_file": command["input_file"],
        "output_dir": manifest["paths"]["output_dir"],
        "result": args.result,
        "notes": args.notes,
    }
    reservation["state"] = "recorded"
    reservation["actual_node_hours"] = actual
    reservation["result"] = args.result
    manifest["state"] = "recorded"
    manifest["accounting"] = row
    if inspection is not None:
        manifest["scientific_inspection"] = inspection
    durable_transition(
        paths, "recorded", manifest_path, manifest, reservations, ledger_row=row
    )
    print(
        f"Recorded {actual:.6f} node-hours for {run['case_id']}/{run['segment']}; "
        f"Stage I cumulative={cumulative:.6f}."
    )
    return 0


def accepted_case_segments(paths: dict[str, Path],
                           case_id: str) -> list[dict[str, object]]:
    """Load scientifically qualifying, accounted segments for one mapped case."""

    segments: list[dict[str, object]] = []
    manifests = sorted(
        (paths["runs"] / case_id).glob("*/manifest/prepared_run.json")
    )
    allow_legacy_local = paths["root"].resolve() != DEFAULT_ROOT.resolve()
    reservations = read_reservations(paths)
    for manifest_path in manifests:
        manifest = read_manifest(manifest_path)
        require_current_epoch(manifest, "retained segment")
        accounting = manifest.get("accounting", {})
        if manifest.get("state") != "recorded" or not isinstance(accounting, dict):
            continue
        result = accounting.get("result")
        if result not in {"accepted", "clean_partial"}:
            continue
        inspection = manifest.get("scientific_inspection")
        if (
            not isinstance(inspection, dict)
            or (
                result == "accepted"
                and inspection.get("accepted") is not True
            )
            or (
                result == "clean_partial"
                and inspection.get("clean_for_continuation") is not True
            )
        ):
            raise ValueError(
                f"retained segment lacks qualifying inspection: {manifest_path}"
            )
        authenticate_prepared_execution(
            manifest, manifest_path, allow_legacy_local=allow_legacy_local
        )
        require_reserved_execution_intent(
            reservation_for_manifest(reservations, manifest_path),
            manifest,
            allow_legacy_local=allow_legacy_local,
        )
        if not allow_legacy_local or "mhd_history" in inspection:
            revalidate_inspection_files(inspection, manifest)
        manifest["_manifest_path"] = str(manifest_path)
        segments.append(manifest)
    segments.sort(
        key=lambda item: float(item["scientific_inspection"]["final_time"])
    )
    return segments


def accepted_case_lineage(paths: dict[str, Path],
                          case_id: str) -> list[dict[str, object]]:
    """Select the recorded restart lineage ending at the latest accepted segment."""

    segments = accepted_case_segments(paths, case_id)
    if not segments:
        return []
    accepted = [
        segment for segment in segments
        if segment["accounting"]["result"] == "accepted"
    ]
    if not accepted:
        raise ValueError(f"{case_id} has no accepted terminal segment")
    final_time = max(
        float(segment["scientific_inspection"]["final_time"])
        for segment in accepted
    )
    terminals = [
        segment for segment in accepted
        if abs(
            float(segment["scientific_inspection"]["final_time"]) - final_time
        ) <= 1.0e-10
    ]
    if len(terminals) != 1:
        raise ValueError(
            f"{case_id} has multiple latest accepted terminal segments"
        )
    indexed = {
        Path(str(segment["_manifest_path"])).resolve(): segment
        for segment in segments
    }
    lineage: list[dict[str, object]] = []
    current = terminals[0]
    visited: set[Path] = set()
    while True:
        path = Path(str(current["_manifest_path"])).resolve()
        if path in visited:
            raise ValueError(f"{case_id} accepted restart lineage contains a cycle")
        visited.add(path)
        lineage.append(current)
        parent = current["command"].get("parent_segment")
        if parent is None:
            break
        if not isinstance(parent, dict) or "manifest" not in parent:
            raise ValueError(
                f"{case_id} accepted restart lineage lacks a parent manifest"
            )
        parent_path = Path(str(parent["manifest"])).resolve()
        if parent_path not in indexed:
            raise ValueError(
                f"{case_id} accepted restart lineage parent is not a "
                f"qualifying recorded segment: {parent_path}"
            )
        current = indexed[parent_path]
    lineage.reverse()
    return lineage


def one_segment_output(manifest: dict[str, object], pattern: str) -> Path:
    """Select the unique retained output matching one segment product."""

    output_dir = Path(str(manifest["paths"]["output_dir"]))
    matches = sorted(output_dir.glob(pattern))
    if len(matches) != 1:
        raise ValueError(
            f"segment output must have one {pattern} product: {output_dir}"
        )
    return matches[0]


def inspected_history_path(manifest: dict[str, object], key: str,
                           pattern: str, allow_legacy_local: bool) -> Path:
    """Return one authenticated inspection-retained history path."""

    inspection = manifest["scientific_inspection"]
    record = inspection.get(key)
    if isinstance(record, dict):
        revalidate_retained_file(record)
        return Path(str(record["path"])).resolve()
    if allow_legacy_local:
        return one_segment_output(manifest, pattern)
    raise ValueError(f"retained segment inspection lacks {key}")


def inspected_snapshot_groups(manifest: dict[str, object],
                              allow_legacy_local: bool) -> list[list[Path]]:
    """Return only snapshot groups retained by formal inspection."""

    inspection = manifest["scientific_inspection"]
    records = inspection.get("snapshots")
    if isinstance(records, list):
        groups = recorded_product_groups(records, "snapshots")
        for record in records:
            revalidate_retained_product(record)
        return groups
    if allow_legacy_local:
        output_dir = Path(str(manifest["paths"]["output_dir"]))
        expected_ranks = (
            int(manifest["allocation"]["nodes"])
            * int(manifest["allocation"].get("ranks_per_node", 1))
        )
        return output_product_groups(output_dir / "bin", "*.bin", expected_ranks)
    raise ValueError("retained segment inspection lacks snapshots")


def link_distinct_snapshots(sources: list[list[Path]], destination: Path,
                            case_name: str) -> list[Path]:
    """Link one retained binary per physical time into an analysis bundle."""

    timed = sorted((binary_product_time(group), group) for group in sources)
    linked: list[Path] = []
    last_time = float("-inf")
    for time, group in timed:
        if time <= last_time + 1.0e-12:
            continue
        name = f"{case_name}.{len(linked):05d}.bin"
        if group[0].parent.name.startswith("rank_"):
            rank0_target = None
            for source in group:
                rank_dir = destination / source.parent.name
                mkdir_durable(rank_dir)
                target = rank_dir / name
                target.symlink_to(source)
                fsync_directory(target.parent)
                if source.parent.name == "rank_00000000":
                    rank0_target = target
            if rank0_target is None:
                raise ValueError("rank-local snapshot set has no rank zero file")
            linked.append(rank0_target)
        else:
            target = destination / name
            target.symlink_to(group[0])
            fsync_directory(target.parent)
            linked.append(target)
        last_time = time
    if not linked:
        raise ValueError("accepted segments retain no distinct snapshots")
    return linked


@locked_root_action
def bundle_case(args: argparse.Namespace) -> int:
    """Assemble accepted restart segments as one analyzer-compatible bundle."""

    required_final_time = require_positive_finite_float(
        args.required_final_time, "--required-final-time"
    )
    root = require_root(Path(args.root), args.allow_local_root)
    if is_offline_local_root(root, args.allow_local_root):
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    matrix = validate_matrix(matrix_path, source_dir)
    case = case_for_id(matrix, args.case_id)
    segments = accepted_case_lineage(paths, args.case_id)
    if not segments:
        raise ValueError(f"no accepted segments are recorded for {args.case_id}")
    final_time = float(segments[-1]["scientific_inspection"]["final_time"])
    if final_time < required_final_time - 1.0e-10:
        raise ValueError(
            f"{args.case_id} reaches only t={final_time}; "
            f"required final time is {required_final_time}"
        )
    first_command = segments[0]["command"]
    expected_digests = (
        first_command["input_sha256"],
        first_command["executable_sha256"],
    )
    submitted_input = Path(str(first_command["input_file"]))
    allow_legacy_local = root.resolve() != DEFAULT_ROOT.resolve()
    if not allow_legacy_local:
        require_file_sha256(
            submitted_input, first_command["input_sha256"], "accepted input"
        )
    model_choices = analysis_model_choices(submitted_input)
    if model_choices.get("output1_file_type") != "hst":
        raise ValueError("accepted production input does not retain output1 history")
    try:
        history_interval = float(model_choices["output1_dt"])
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(
            "accepted production input has no numeric history cadence"
        ) from error
    if not math.isfinite(history_interval) or history_interval <= 0.0:
        raise ValueError("accepted production input has invalid history cadence")
    previous_final = None
    previous_history: dict[str, list[float]] | None = None
    for segment in segments:
        command = segment["command"]
        if (
            command["input_sha256"],
            command["executable_sha256"],
        ) != expected_digests:
            raise ValueError("accepted segments do not share input/executable digests")
        history = parse_history(
            inspected_history_path(
                segment, "mhd_history", "*.mhd.hst", allow_legacy_local
            )
        )
        first_time = history["time"][0]
        segment_final = history["time"][-1]
        if previous_final is None:
            if first_time > 1.0e-10:
                raise ValueError("accepted case history does not begin at t=0")
        else:
            if (
                previous_history is None
                or "dt" not in previous_history
                or "dt" not in history
            ):
                raise ValueError(
                    "accepted histories lack timesteps needed to validate restarts"
                )
            boundary_step = max(previous_history["dt"][-1], history["dt"][0])
            if not math.isfinite(boundary_step) or boundary_step < 0.0:
                raise ValueError("accepted histories have invalid boundary timestep")
            # A walltime terminal row can replace one scheduled output at restart.
            permitted_gap = 2.0 * history_interval + boundary_step + 1.0e-10
            if first_time > previous_final + permitted_gap:
                raise ValueError(
                    "accepted restart histories exceed configured sampling cadence"
                )
            if segment_final <= previous_final + 1.0e-10:
                raise ValueError("accepted restart segment does not advance time")
        previous_final = segment_final
        previous_history = history
    bundle = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else paths["runs"] / "bundles" / args.case_id
    )
    require_beneath_root(bundle, root, "analysis bundle", args.allow_local_root)
    if bundle.exists():
        if not args.replace:
            raise ValueError(f"analysis bundle already exists: {bundle}")
        shutil.rmtree(bundle)
    history_dir = bundle / "history"
    input_dir = bundle / "inputs"
    snapshot_dir = bundle / "cases" / str(case["name"]) / "bin"
    for directory in (history_dir, input_dir, snapshot_dir):
        mkdir_durable(directory)
    mhd_history = history_dir / f"{case['name']}.mhd.hst"
    user_history = history_dir / f"{case['name']}.user.hst"
    merge_history_files(
        [
            inspected_history_path(
                segment, "mhd_history", "*.mhd.hst", allow_legacy_local
            )
            for segment in segments
        ],
        mhd_history,
    )
    merge_history_files(
        [
            inspected_history_path(
                segment, "user_history", "*.user.hst", allow_legacy_local
            )
            for segment in segments
        ],
        user_history,
    )
    snapshot_sources: list[list[Path]] = []
    for segment in segments:
        snapshot_sources.extend(inspected_snapshot_groups(segment, allow_legacy_local))
    snapshots = link_distinct_snapshots(
        snapshot_sources, snapshot_dir, str(case["name"])
    )
    archived_input = input_dir / submitted_input.name
    copy_file(submitted_input, archived_input)
    case_entry = {
        "name": case["name"],
        "input": case["input"],
        "execution_input": str(archived_input.relative_to(bundle)),
        "command": ["archived production segments"],
        "overrides": [],
        "log": "retained in production segment archive",
        "status": "passed",
        "outputs": {
            "mhd_history": str(mhd_history.relative_to(bundle)),
            "user_history": str(user_history.relative_to(bundle)),
            "snapshot_paths": [
                str(path.relative_to(bundle)) for path in snapshots
            ],
        },
        "lf_active": True,
        "amr": False,
        "paper_smoke": False,
        "model_choices": model_choices,
    }
    manifest = {
        "workflow": "paper-mks24-stage-i-production",
        "execution_epoch": EXECUTION_EPOCH,
        "created_utc": utc_now(),
        "status": "accepted_for_analysis",
        "git_revision": first_command["input_revision"],
        "git_worktree_dirty": (
            "not recorded; submitted input and matrix verified committed"
        ),
        "executable": first_command["executable"],
        "production_case_id": args.case_id,
        "required_final_time": required_final_time,
        "accepted_final_time": final_time,
        "production_segment_manifests": [
            str(segment["_manifest_path"]) for segment in segments
        ],
        "cases": [case_entry],
    }
    write_json(bundle / "manifest.json", manifest)
    print(f"Wrote accepted Stage I analysis bundle: {bundle}")
    print(f"  final_time={final_time:.12g}, snapshots={len(snapshots)}")
    return 0


def prefix_bundle_case_paths(case: dict[str, object], prefix: Path
                             ) -> dict[str, object]:
    """Rewrite nested per-case paths for a campaign-level bundle manifest."""

    copied = json.loads(json.dumps(case))
    copied["execution_input"] = str(prefix / str(copied["execution_input"]))
    outputs = copied["outputs"]
    for key in ("mhd_history", "user_history"):
        outputs[key] = str(prefix / str(outputs[key]))
    outputs["snapshot_paths"] = [
        str(prefix / str(path)) for path in outputs["snapshot_paths"]
    ]
    return copied


@locked_root_action
def bundle_campaign(args: argparse.Namespace) -> int:
    """Assemble all accepted mapped cases into one paper-analysis bundle."""

    required_final_time = require_positive_finite_float(
        args.required_final_time, "--required-final-time"
    )
    root = require_root(Path(args.root), args.allow_local_root)
    if is_offline_local_root(root, args.allow_local_root):
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
        require_reconciled_store_consistency(paths)
    source_dir = Path(args.source_dir).expanduser().resolve()
    matrix_path = Path(args.matrix).expanduser().resolve()
    matrix = validate_matrix(matrix_path, source_dir)
    for case in matrix["cases"]:
        segments = accepted_case_lineage(paths, str(case["id"]))
        if not segments:
            raise ValueError(f"no accepted segments are recorded for {case['id']}")
        final_time = float(segments[-1]["scientific_inspection"]["final_time"])
        if final_time < required_final_time - 1.0e-10:
            raise ValueError(
                f"{case['id']} reaches only t={final_time}; "
                f"required final time is {required_final_time}"
            )
    bundle = (
        Path(args.output_dir).expanduser().resolve()
        if args.output_dir
        else paths["runs"] / "bundles" / "mks24-stage-i-campaign"
    )
    require_beneath_root(bundle, root, "analysis bundle", args.allow_local_root)
    if bundle.exists():
        if not args.replace:
            raise ValueError(f"analysis bundle already exists: {bundle}")
        shutil.rmtree(bundle)
    mkdir_durable(bundle)
    cases: list[dict[str, object]] = []
    segment_manifests: list[str] = []
    case_times: dict[str, float] = {}
    for case in matrix["cases"]:
        case_id = str(case["id"])
        nested = bundle / "case_bundles" / case_id
        case_args = argparse.Namespace(**vars(args))
        case_args.case_id = case_id
        case_args.output_dir = str(nested)
        case_args.replace = False
        bundle_case(case_args)
        nested_manifest = read_manifest(nested / "manifest.json")
        prefix = nested.relative_to(bundle)
        cases.append(prefix_bundle_case_paths(nested_manifest["cases"][0], prefix))
        segment_manifests.extend(nested_manifest["production_segment_manifests"])
        case_times[case_id] = float(nested_manifest["accepted_final_time"])
    campaign_manifest = {
        "workflow": "paper-mks24-stage-i-production",
        "execution_epoch": EXECUTION_EPOCH,
        "created_utc": utc_now(),
        "status": "accepted_for_analysis",
        "git_revision": subprocess.run(
            ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
            check=True, capture_output=True, text=True,
        ).stdout.strip(),
        "git_worktree_dirty": (
            "not recorded; submitted inputs and matrices verified committed"
        ),
        "executable": "recorded per accepted production segment",
        "required_final_time": required_final_time,
        "accepted_case_final_times": case_times,
        "production_segment_manifests": segment_manifests,
        "cases": cases,
    }
    write_json(bundle / "manifest.json", campaign_manifest)
    print(f"Wrote accepted Stage I campaign bundle: {bundle}")
    print(f"  accepted cases={len(cases)}")
    return 0


@locked_manifest_action
def cancel(args: argparse.Namespace) -> int:
    """Release a segment that was prepared but never submitted."""

    manifest_path = Path(args.manifest).resolve()
    manifest = read_manifest(manifest_path)
    require_current_epoch(manifest, "prepared segment")
    root = require_root(Path(str(manifest["project_root"])), args.allow_local_root)
    offline_local_root = is_offline_local_root(root, args.allow_local_root)
    if offline_local_root:
        paths = initialize(root)
    else:
        paths = layout(root)
        require_existing_layout(paths)
        require_no_pending_transactions(paths)
        require_no_orphaned_segment_runs(paths)
    if manifest.get("state") != "prepared":
        raise ValueError("only an unsubmitted segment may be cancelled")
    reservations = read_reservations(paths)
    reservation = reservation_for_manifest(reservations, manifest_path)
    if reservation.get("state") != "prepared":
        raise ValueError("reservation has progressed beyond preparation")
    break_glass = getattr(args, "break_glass_cancel_evidence", "")
    if getattr(args, "confirm_break_glass_cancel", False):
        if not str(break_glass).strip():
            raise ValueError("--break-glass-cancel-evidence must be nonempty")
        cancellation_mode = "break-glass"
    else:
        if break_glass:
            raise ValueError(
                "--confirm-break-glass-cancel is required with break-glass evidence"
            )
        require_reconciled_store_consistency(paths)
        require_reserved_execution_intent(
            reservation, manifest, allow_legacy_local=offline_local_root
        )
        authenticate_prepared_execution(
            manifest, manifest_path, allow_legacy_local=offline_local_root
        )
        cancellation_mode = "authenticated"
    reservation["state"] = "cancelled"
    reservation["notes"] = args.notes
    manifest["state"] = "cancelled"
    manifest["cancellation_notes"] = args.notes
    manifest["cancellation"] = {
        "cancelled_utc": utc_now(),
        "mode": cancellation_mode,
        "notes": args.notes,
        "break_glass_evidence": break_glass or None,
    }
    durable_transition(paths, "cancelled", manifest_path, manifest, reservations)
    print(f"Cancelled Stage I reservation: {manifest_path}")
    return 0


def reconcile_report(root: Path) -> dict[str, object]:
    """Read retained stores and report consistency without modifying them."""

    paths = layout(root)
    issues: list[str] = []
    offline_local_root = root.resolve() != DEFAULT_ROOT.expanduser().resolve()
    qualification = qualification_approval_status(paths)
    if (
        qualification["state"] == "invalid"
        or (
            not offline_local_root
            and qualification["state"] != "approved"
        )
    ):
        issues.append(
            "E03 corrected-build Frontier qualification is "
            f"{qualification['state']}: {qualification.get('reason', qualification['path'])}"
        )
    if not paths["transactions"].is_dir():
        issues.append(f"transaction store is missing: {paths['transactions']}")
    transactions = pending_transaction_paths(paths)
    for transaction_path in transactions:
        try:
            transaction = read_transaction(paths, transaction_path)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            issues.append(f"cannot read transaction {transaction_path}: {error}")
        else:
            issues.append(
                f"pending {transaction.get('kind')} transaction requires recovery: "
                f"{transaction_path}"
            )
    if paths["ledger"].is_file():
        try:
            ledger = read_ledger(paths)
        except (OSError, ValueError) as error:
            ledger = []
            issues.append(f"cannot read ledger: {error}")
    else:
        ledger = []
        issues.append(f"ledger is missing: {paths['ledger']}")
    if paths["reservations"].is_file():
        try:
            reservations = read_reservations(paths)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            reservations = []
            issues.append(f"cannot read reservations: {error}")
    else:
        reservations = []
        issues.append(f"reservation store is missing: {paths['reservations']}")

    manifests: dict[Path, dict[str, object]] = {}
    if paths["runs"].is_dir():
        for manifest_path in sorted(
            paths["runs"].glob("*/*/manifest/prepared_run.json")
        ):
            resolved = manifest_path.resolve()
            try:
                manifests[resolved] = read_manifest(resolved)
            except (OSError, ValueError, json.JSONDecodeError) as error:
                issues.append(f"cannot read manifest {resolved}: {error}")
    else:
        issues.append(f"run store is missing: {paths['runs']}")
    for run_dir in orphaned_segment_run_directories(paths):
        issues.append(f"orphaned segment run directory: {run_dir}")

    active = [
        reservation for reservation in reservations
        if isinstance(reservation, dict)
        and reservation.get("state") in {"prepared", "submitted"}
    ]
    if len(active) > 1:
        issues.append(f"reservation store has {len(active)} active segments")

    reservations_by_manifest: dict[Path, list[dict[str, object]]] = {}
    for reservation in reservations:
        if not isinstance(reservation, dict):
            issues.append(f"reservation record is invalid: {reservation!r}")
            continue
        try:
            validate_reservation_record(paths, reservation)
        except (KeyError, TypeError, ValueError) as error:
            issues.append(f"reservation record is invalid: {error}")
            continue
        manifest_path = Path(str(reservation.get("manifest", ""))).resolve()
        reservations_by_manifest.setdefault(manifest_path, []).append(reservation)
        if reservation.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"reservation has wrong execution epoch: {manifest_path}")
        manifest = manifests.get(manifest_path)
        if manifest is None:
            issues.append(f"reservation lacks retained manifest: {manifest_path}")
            continue
        run = manifest.get("run")
        allocation = manifest.get("allocation")
        if not isinstance(run, dict) or not isinstance(allocation, dict):
            issues.append(f"manifest lacks run or allocation metadata: {manifest_path}")
            continue
        for key in ("case_id", "case_name", "segment"):
            if reservation.get(key) != run.get(key):
                issues.append(f"reservation {key} differs from manifest: {manifest_path}")
        for key in ("nodes", "requested_walltime"):
            if reservation.get(key) != allocation.get(key):
                issues.append(
                    f"reservation {key} differs from allocation: {manifest_path}"
                )
        try:
            reserved_difference = abs(
                float(reservation["reserved_node_hours"])
                - float(allocation["reserved_node_hours"])
            )
        except (KeyError, TypeError, ValueError):
            issues.append(f"reservation node-hours are invalid: {manifest_path}")
        else:
            if (
                not math.isfinite(reserved_difference)
                or reserved_difference > 5.0e-12
            ):
                issues.append(
                    f"reservation node-hours differ from allocation: {manifest_path}"
                )
        if reservation.get("state") != manifest.get("state"):
            issues.append(f"reservation state differs from manifest: {manifest_path}")
        try:
            require_reserved_execution_intent(
                reservation, manifest, allow_legacy_local=offline_local_root
            )
        except (KeyError, TypeError, ValueError) as error:
            issues.append(
                f"reservation execution intent differs from manifest "
                f"{manifest_path}: {error}"
            )
        if (
            reservation.get("state") in {"submitted", "recorded"}
            and reservation.get("job_id") != manifest.get("job_id")
        ):
            issues.append(f"reservation job ID differs from manifest: {manifest_path}")

    ledger_by_job: dict[str, list[dict[str, str]]] = {}
    for row in ledger:
        job_id = row.get("job_id", "")
        ledger_by_job.setdefault(job_id, []).append(row)
        if row.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"ledger row has wrong execution epoch: {job_id}")
        try:
            require_numeric_job_id(job_id)
        except ValueError:
            issues.append(f"ledger row has invalid job ID: {job_id!r}")
    for job_id, rows in ledger_by_job.items():
        if len(rows) != 1:
            issues.append(f"ledger has {len(rows)} rows for job {job_id}")

    manifests_by_job: dict[str, list[Path]] = {}
    recorded_manifests_by_job: dict[str, list[Path]] = {}
    for manifest_path, manifest in manifests.items():
        if manifest.get("execution_epoch") != EXECUTION_EPOCH:
            issues.append(f"manifest has wrong execution epoch: {manifest_path}")
        matches = reservations_by_manifest.get(manifest_path, [])
        if len(matches) != 1:
            issues.append(
                f"manifest has {len(matches)} reservation records: {manifest_path}"
            )
        state = manifest.get("state")
        if state not in {"prepared", "submitted", "recorded", "cancelled"}:
            issues.append(f"manifest has invalid state {state!r}: {manifest_path}")
        if state in {"prepared", "submitted", "recorded"}:
            try:
                command = manifest.get("command", {})
                if (
                    not offline_local_root
                    or (
                        isinstance(command, dict)
                        and isinstance(command.get("overrides"), list)
                    )
                ):
                    prepared_time_tlim_target(manifest)
                authenticate_prepared_execution(
                    manifest, manifest_path,
                    allow_legacy_local=offline_local_root,
                )
            except (OSError, ValueError, KeyError, TypeError) as error:
                issues.append(f"prepared artifact drift for {manifest_path}: {error}")
        if state == "prepared":
            if manifest.get("job_id") is not None:
                issues.append(f"prepared manifest unexpectedly has a job ID: {manifest_path}")
            continue
        if state not in {"submitted", "recorded"}:
            continue
        job_id = str(manifest.get("job_id", ""))
        try:
            require_numeric_job_id(job_id)
        except ValueError:
            issues.append(f"manifest has invalid job ID: {manifest_path}")
        manifests_by_job.setdefault(job_id, []).append(manifest_path)
        if state != "recorded":
            continue
        recorded_manifests_by_job.setdefault(job_id, []).append(manifest_path)
        accounting = manifest.get("accounting")
        if not isinstance(accounting, dict) or accounting.get("job_id") != job_id:
            issues.append(f"recorded manifest lacks matching accounting: {manifest_path}")
        rows = ledger_by_job.get(job_id, [])
        if len(rows) != 1:
            issues.append(f"recorded manifest lacks one ledger row: {manifest_path}")
            continue
        row = rows[0]
        if accounting != row:
            issues.append(f"manifest accounting differs from ledger: {manifest_path}")
        try:
            validate_transaction_ledger_row(row, manifest)
        except (KeyError, TypeError, ValueError) as error:
            issues.append(
                f"recorded ledger provenance differs for {manifest_path}: {error}"
            )
        if not offline_local_root:
            try:
                require_recorded_scheduler_evidence(paths, row, manifest)
            except (OSError, KeyError, TypeError, ValueError) as error:
                issues.append(
                    f"recorded scheduler evidence differs for {manifest_path}: {error}"
                )
        run = manifest.get("run")
        if not isinstance(run, dict):
            issues.append(f"recorded manifest lacks run metadata: {manifest_path}")
        else:
            for key in ("case_id", "case_name", "segment"):
                if row.get(key) != run.get(key):
                    issues.append(f"ledger {key} differs from manifest: {manifest_path}")
        if len(matches) == 1:
            reservation = matches[0]
            if reservation.get("result") != row.get("result"):
                issues.append(f"reservation result differs from ledger: {manifest_path}")
            try:
                difference = abs(
                    float(reservation["actual_node_hours"])
                    - float(row["actual_node_hours"])
                )
            except (KeyError, TypeError, ValueError):
                issues.append(
                    f"reservation actual node-hours are invalid: {manifest_path}"
                )
            else:
                if not math.isfinite(difference) or difference > 5.0e-7:
                    issues.append(
                        f"reservation actual node-hours differ from ledger: "
                        f"{manifest_path}"
                    )
        if isinstance(accounting, dict) and accounting.get("result") in {
            "accepted", "clean_partial"
        }:
            inspection = manifest.get("scientific_inspection")
            if not isinstance(inspection, dict):
                issues.append(f"recorded manifest lacks inspection: {manifest_path}")
            else:
                try:
                    inspected_target = float(inspection["required_time"])
                    if (
                        not math.isfinite(inspected_target)
                        or abs(inspected_target - prepared_time_tlim_target(manifest))
                        > 1.0e-12
                    ):
                        raise ValueError("inspection target differs from preparation")
                    revalidate_inspection_files(inspection, manifest)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    issues.append(f"inspection drift for {manifest_path}: {error}")

    for job_id, manifest_paths in manifests_by_job.items():
        if len(manifest_paths) != 1:
            issues.append(
                f"job {job_id} is attached to {len(manifest_paths)} manifests"
            )
    for job_id in ledger_by_job:
        matches = recorded_manifests_by_job.get(job_id, [])
        if len(matches) != 1:
            issues.append(f"ledger job {job_id} has {len(matches)} recorded manifests")

    return {
        "execution_epoch": EXECUTION_EPOCH,
        "root": str(root),
        "qualification": qualification,
        "consistent": not issues,
        "counts": {
            "transactions": len(transactions),
            "reservations": len(reservations),
            "active_reservations": len(active),
            "ledger_rows": len(ledger),
            "manifests": len(manifests),
        },
        "issues": issues,
    }


def reconcile(args: argparse.Namespace) -> int:
    """Print a read-only reservations, ledger, and manifest consistency report."""

    root = require_root(Path(args.root), args.allow_local_root)
    report = reconcile_report(root)
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["consistent"] else 1


def parser() -> argparse.ArgumentParser:
    """Build command-line parsing."""

    command = argparse.ArgumentParser(
        description=__doc__,
        epilog=(
            "E03 production qualification remains pending until the canonical "
            "corrected-build approval token exists. Create it only after review "
            "with approve-qualification."
        ),
    )
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument(
        "--allow-local-root", action="store_true",
        help="Permit a non-project root only for offline validation.",
    )
    actions = command.add_subparsers(dest="action", required=True)
    validate = actions.add_parser("validate-matrix")
    validate.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    validate.add_argument("--source-dir", default=str(ROOT_DIR))
    actions.add_parser("init")
    approved = actions.add_parser(
        "approve-qualification",
        help=(
            "Atomically create the E03 corrected-build qualification token "
            "after Frontier review."
        ),
    )
    approved.add_argument("--executable", required=True)
    approved.add_argument("--build-manifest", required=True)
    approved.add_argument("--approved-by", required=True)
    approved.add_argument("--review-notes", required=True)
    approved.add_argument(
        "--confirm-corrected-build-frontier-qualified", action="store_true"
    )
    approved.add_argument("--replace-existing-approval", action="store_true")
    prepare_parser = actions.add_parser(
        "prepare",
        help=(
            "Prepare one E03 segment; canonical production rejects until the "
            "qualification token exists."
        ),
    )
    prepare_parser.add_argument("--case-id", required=True)
    prepare_parser.add_argument("--segment", type=require_safe_segment, required=True)
    prepare_parser.add_argument("--acceptance-criterion", required=True)
    prepare_parser.add_argument("--executable", required=True)
    prepare_parser.add_argument("--build-manifest", required=True)
    prepare_parser.add_argument("--source-dir", default=str(ROOT_DIR))
    prepare_parser.add_argument(
        "--source-bundle",
        help=(
            "Retained Git bundle beneath the CGL root containing the input "
            "and executable revisions."
        ),
    )
    prepare_parser.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    prepare_parser.add_argument("--restart-file")
    prepare_parser.add_argument("--nodes", type=int, required=True)
    prepare_parser.add_argument("--walltime", required=True)
    prepare_parser.add_argument("--athena-walltime", required=True)
    prepare_parser.add_argument("--ranks-per-node", type=int, default=8)
    prepare_parser.add_argument("--cpus-per-task", type=int, default=7)
    prepare_parser.add_argument("--override", action="append", default=[])
    prepare_parser.add_argument(
        "--allow-missing-time-target", action="store_true",
        help=(
            "Permit offline local-root preparation without a time/tlim "
            "override. This is never allowed for real production."
        ),
    )
    prepare_parser.add_argument(
        "--allow-missing-restart-time-marker", action="store_true",
        help=(
            "Permit an offline legacy continuation fixture without explicit "
            "time/restart_time markers. This is never allowed for production."
        ),
    )
    checked = actions.add_parser("check-submit")
    submitted_atomically = actions.add_parser("submit")
    for submit_parser in (checked, submitted_atomically):
        submit_parser.add_argument("--manifest", required=True)
        submit_parser.add_argument("--squeue-file")
        submit_parser.add_argument("--skip-slurm-test", action="store_true")
        submit_parser.add_argument(
            "--allow-shared-root-campaign", action="append", default=[],
            help=(
                "Acknowledge one reviewed top-level CGL-root campaign record. "
                "Queued user jobs still fail closed."
            ),
        )
    submitted_atomically.add_argument(
        "--sbatch-output-file",
        help="Use retained sbatch --parsable output only for offline validation.",
    )
    submitted = actions.add_parser("mark-submitted")
    submitted.add_argument("--manifest", required=True)
    submitted.add_argument("--job-id", type=require_numeric_job_id, required=True)
    recovered_submit = actions.add_parser("recover-submit")
    recovered_submit.add_argument("--manifest", required=True)
    recovered_submit.add_argument(
        "--job-id", type=require_numeric_job_id, required=True
    )
    cleared_submit = actions.add_parser("clear-submit-pending")
    cleared_submit.add_argument("--manifest", required=True)
    cleared_submit.add_argument("--notes", required=True)
    cleared_submit.add_argument("--confirm-no-job-submitted", action="store_true")
    cleared_submit.add_argument(
        "--scheduler-absence-evidence-file",
        help="Use machine-readable absence evidence only for offline local fixtures.",
    )
    cleared_submit.add_argument("--break-glass-clear-evidence")
    cleared_submit.add_argument("--confirm-break-glass-clear", action="store_true")
    inspected = actions.add_parser("inspect-segment")
    inspected.add_argument("--manifest", required=True)
    inspected.add_argument("--required-time", type=float, required=True)
    inspected.add_argument(
        "--allow-missing-restart-time-marker", action="store_true",
        help=(
            "Permit offline inspection of legacy restart fixtures without "
            "explicit time/restart_time markers."
        ),
    )
    recorded = actions.add_parser("record")
    recorded.add_argument("--manifest", required=True)
    recorded.add_argument("--job-id", type=require_numeric_job_id, required=True)
    recorded.add_argument(
        "--result",
        choices=("accepted", "clean_partial", "rejected", "failed", "aborted"),
        required=True,
    )
    recorded.add_argument("--notes", default="")
    recorded.add_argument("--sacct-file")
    bundled = actions.add_parser("bundle-case")
    bundled.add_argument("--case-id", required=True)
    bundled.add_argument(
        "--required-final-time", type=positive_finite_float_arg,
        default=REQUIRED_CASE_FINAL_TIME
    )
    bundled.add_argument("--source-dir", default=str(ROOT_DIR))
    bundled.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    bundled.add_argument("--output-dir")
    bundled.add_argument("--replace", action="store_true")
    campaign = actions.add_parser("bundle-campaign")
    campaign.add_argument(
        "--required-final-time", type=positive_finite_float_arg,
        default=REQUIRED_CASE_FINAL_TIME
    )
    campaign.add_argument("--source-dir", default=str(ROOT_DIR))
    campaign.add_argument("--matrix", default=str(DEFAULT_MATRIX))
    campaign.add_argument("--output-dir")
    campaign.add_argument("--replace", action="store_true")
    cancelled = actions.add_parser("cancel")
    cancelled.add_argument("--manifest", required=True)
    cancelled.add_argument("--notes", required=True)
    cancelled.add_argument("--break-glass-cancel-evidence")
    cancelled.add_argument("--confirm-break-glass-cancel", action="store_true")
    actions.add_parser("summary")
    actions.add_parser("reconcile")
    actions.add_parser("recover-transactions")
    return command


def main() -> int:
    """Command-line entry point."""

    args = parser().parse_args()
    try:
        if args.action == "validate-matrix":
            matrix = validate_matrix(
                Path(args.matrix).expanduser().resolve(),
                Path(args.source_dir).expanduser().resolve(),
            )
            print(f"Validated {len(matrix['cases'])} mapped Stage I cases.")
            return 0
        root = require_root(Path(args.root), args.allow_local_root)
        if args.action == "init":
            initialize(root)
            print(f"Initialized Stage I production accounting beneath {root}.")
            return 0
        if args.action == "approve-qualification":
            return approve_qualification(args)
        if args.action == "prepare":
            prepare(args)
            return 0
        if args.action == "check-submit":
            return check_submit(args)
        if args.action == "submit":
            return submit(args)
        if args.action == "mark-submitted":
            return mark_submitted(args)
        if args.action == "recover-submit":
            return recover_submit(args)
        if args.action == "clear-submit-pending":
            return clear_submit_pending(args)
        if args.action == "inspect-segment":
            return inspect_segment(args)
        if args.action == "record":
            return record(args)
        if args.action == "bundle-case":
            return bundle_case(args)
        if args.action == "bundle-campaign":
            return bundle_campaign(args)
        if args.action == "cancel":
            return cancel(args)
        if args.action == "summary":
            with canonical_root_lock(root):
                if is_offline_local_root(root, args.allow_local_root):
                    paths = initialize(root)
                else:
                    paths = layout(root)
                    require_existing_layout(paths)
                    require_no_pending_transactions(paths)
                    require_no_orphaned_segment_runs(paths)
                    require_reconciled_store_consistency(
                        paths, allow_absent_qualification=True
                    )
                refresh_summary(paths)
            print(f"Wrote {paths['summary']}")
            return 0
        if args.action == "reconcile":
            return reconcile(args)
        if args.action == "recover-transactions":
            return recover_transactions(args)
        raise ValueError(f"unsupported action: {args.action}")
    except (KeyError, OSError, TypeError, ValueError,
            subprocess.CalledProcessError) as error:
        print(f"Stage I production utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
