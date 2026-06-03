#!/usr/bin/env python3
"""Verify and publish retained Stage I recost checkpoints.

The companion is intentionally generic: callers supply the retained artifact
name, authenticated file digests, the sole authorized next-segment profile,
expected reconcile counts, and JSON pointers locating those bindings in the
artifact.  Canonical publication uses a same-directory link/fsync/unlink/fsync
sequence under the Stage I lock.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import datetime, timedelta, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import stat
import subprocess
import sys
import uuid


DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
CANONICAL_REPOSITORY_ROOT = Path("/autofs/nccs-svm1_home2/dfielding/athenak-df")
SQUEUE = Path("/usr/bin/squeue")
GIT = Path("/usr/bin/git")
OFFLINE_ALLOWED_PREFIXES = (Path("/tmp"),)
EXECUTION_EPOCH = "E03-forcing-policy"
EXECUTION_EPOCH_SLUG = "E03_forcing_policy"
STAGE_I_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i.py")
UTILITY_RELATIVE = Path("scripts/frontier/cgl_lf_stage_i_checkpoint.py")
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
ARTIFACT_NAME_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,199}\.json")
SELF_DESCRIPTOR_ENV = "_CGL_LF_RECOST_UTILITY_DESCRIPTOR"
SELF_SOURCE_ENV = "_CGL_LF_RECOST_UTILITY_SOURCE"
ROOT_DIR_ENV = "_CGL_LF_RECOST_REPOSITORY_ROOT"
COUNT_KEYS = (
    "transactions",
    "reservations",
    "active_reservations",
    "ledger_rows",
    "manifests",
)
CGL_JOB_NAME_PREFIX = "cgl_"


class LinkAttemptError(RuntimeError):
    """Report an os.link attempt whose post-state requires recovery review."""


def utc_now() -> str:
    """Return one stable UTC timestamp."""

    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def require_utc_timestamp(value: object, label: str) -> None:
    """Require one UTC ISO-8601 timestamp."""

    if not isinstance(value, str):
        raise ValueError(f"{label} must be an ISO-8601 timestamp")
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise ValueError(f"{label} must be an ISO-8601 timestamp") from error
    if parsed.tzinfo is None:
        raise ValueError(f"{label} must include a timezone")
    if parsed.utcoffset() != timedelta(0):
        raise ValueError(f"{label} must use UTC")


def sha256_bytes(value: bytes) -> str:
    """Return the SHA-256 digest for retained bytes."""

    return hashlib.sha256(value).hexdigest()


def read_descriptor_bytes(descriptor: int) -> bytes:
    """Read one descriptor without changing its caller-visible offset."""

    offset = os.lseek(descriptor, 0, os.SEEK_CUR)
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
        blocks = []
        while True:
            block = os.read(descriptor, 1024 * 1024)
            if not block:
                return b"".join(blocks)
            blocks.append(block)
    finally:
        os.lseek(descriptor, offset, os.SEEK_SET)


def sha256_descriptor(descriptor: int) -> str:
    """Return the digest of one open descriptor."""

    return sha256_bytes(read_descriptor_bytes(descriptor))


def sha256(path: Path) -> str:
    """Return the digest of one regular retained file."""

    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        return sha256_descriptor(descriptor)
    finally:
        os.close(descriptor)


def require_sha256(value: object, label: str) -> str:
    """Require a lowercase SHA-256 digest."""

    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise ValueError(f"{label} must be a lowercase SHA-256 digest")
    return value


def sha256_arg(value: str) -> str:
    """Parse a SHA-256 command argument."""

    try:
        return require_sha256(value, "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError(str(error)) from error


def mode_arg(value: str) -> int:
    """Parse a four-digit retained-file mode."""

    if re.fullmatch(r"0?[0-7]{3}", value) is None:
        raise argparse.ArgumentTypeError("mode must be an octal value such as 0644")
    return int(value, 8)


def capped_mode_arg(value: str, maximum: int, label: str) -> int:
    """Parse a retained-file mode with no permissions beyond one ceiling."""

    mode = mode_arg(value)
    if mode & ~maximum:
        raise argparse.ArgumentTypeError(
            f"{label} mode must not exceed {maximum:04o}"
        )
    return mode


def artifact_mode_arg(value: str) -> int:
    """Parse a retained recost-artifact mode."""

    return capped_mode_arg(value, 0o644, "artifact")


def generator_mode_arg(value: str) -> int:
    """Parse a retained recost-generator mode."""

    return capped_mode_arg(value, 0o755, "generator")


def scheduler_mode_arg(value: str) -> int:
    """Parse a retained scheduler-evidence mode."""

    return capped_mode_arg(value, 0o644, "scheduler evidence")


def json_object_arg(value: str) -> dict[str, object]:
    """Parse one explicit JSON object argument."""

    try:
        parsed = json.loads(value)
    except json.JSONDecodeError as error:
        raise argparse.ArgumentTypeError("value must be a JSON object") from error
    if not isinstance(parsed, dict) or not parsed:
        raise argparse.ArgumentTypeError("value must be a nonempty JSON object")
    return parsed


def nonnegative_int_arg(value: str) -> int:
    """Parse one non-negative retained count."""

    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("count must be a non-negative integer") from error
    if parsed < 0:
        raise argparse.ArgumentTypeError("count must be a non-negative integer")
    return parsed


def expected_self_sha256(argv: list[str]) -> str:
    """Read the externally supplied self digest before general CLI parsing."""

    option = "--expected-utility-sha256"
    matches = [
        argv[index + 1]
        for index, value in enumerate(argv[:-1])
        if value == option
    ]
    if len(matches) != 1:
        raise ValueError(f"{option} must be supplied exactly once")
    return require_sha256(matches[0], option)


def initial_source_path() -> Path:
    """Return the retained source path across descriptor re-execution."""

    source = os.environ.get(SELF_SOURCE_ENV)
    if source is not None:
        return Path(source).expanduser().resolve()
    return Path(__file__).expanduser().resolve()


def repository_root(source: Path) -> Path:
    """Return the repository root retained across descriptor re-execution."""

    root = os.environ.get(ROOT_DIR_ENV)
    if root is not None:
        return Path(root).expanduser().resolve()
    try:
        return source.parents[2]
    except IndexError as error:
        raise ValueError(f"utility source path is invalid: {source}") from error


def authenticate_retained_self_source(source: Path, expected: str) -> None:
    """Require the retained source profile even after descriptor re-execution."""

    require_file_sha256(
        source,
        expected,
        "retained utility",
        expected_mode=0o755,
        expected_links=1,
    )


def authenticate_self(expected: str) -> tuple[Path, Path]:
    """Authenticate this utility and re-execute immutable descriptor bytes."""

    expected = require_sha256(expected, "utility SHA-256")
    inherited = os.environ.get(SELF_DESCRIPTOR_ENV)
    if inherited is not None:
        try:
            descriptor = int(inherited)
        except ValueError as error:
            raise ValueError("utility descriptor marker is invalid") from error
        if __file__ != f"/proc/self/fd/{descriptor}":
            raise ValueError("utility descriptor marker is not attached to this execution")
        require_regular_profile(
            os.fstat(descriptor), Path(__file__), "authenticated utility descriptor"
        )
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("authenticated utility descriptor checksum has changed")
        source = initial_source_path()
        root = repository_root(source)
        if source != root / UTILITY_RELATIVE:
            raise ValueError(f"utility retained path is invalid: {source}")
        authenticate_retained_self_source(source, expected)
        return source, root

    source = initial_source_path()
    root = repository_root(source)
    if source != root / UTILITY_RELATIVE:
        raise ValueError(f"utility retained path is invalid: {source}")
    authenticate_retained_self_source(source, expected)
    value = require_file_sha256(
        source,
        expected,
        "retained utility",
        expected_mode=0o755,
        expected_links=1,
    )
    if not hasattr(os, "memfd_create"):
        raise ValueError("descriptor re-execution requires os.memfd_create")
    retained = os.memfd_create("cgl_lf_stage_i_checkpoint")
    os.write(retained, value)
    os.lseek(retained, 0, os.SEEK_SET)
    os.set_inheritable(retained, True)
    environment = dict(os.environ)
    environment[SELF_DESCRIPTOR_ENV] = str(retained)
    environment[SELF_SOURCE_ENV] = str(source)
    environment[ROOT_DIR_ENV] = str(root)
    os.execve(
        sys.executable,
        [sys.executable, f"/proc/self/fd/{retained}", *sys.argv[1:]],
        environment,
    )
    raise AssertionError("descriptor re-execution returned unexpectedly")


def fsync_directory(path: Path) -> None:
    """Persist directory-entry changes."""

    descriptor = os.open(path, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def create_exclusive(path: Path, flags: int, mode: int) -> int:
    """Create one exact-mode entry without allowing the caller umask to narrow it."""

    previous = os.umask(0)
    try:
        return os.open(path, flags | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW, mode)
    finally:
        os.umask(previous)


def write_json(path: Path, value: object, mode: int = 0o644, *,
               simulate_interruption_before_directory_fsync: bool = False) -> None:
    """Write one JSON file atomically and durably."""

    temporary = path.with_name(f".{path.name}.{os.getpid()}.{uuid.uuid4().hex}.tmp")
    descriptor = create_exclusive(
        temporary,
        os.O_WRONLY,
        mode,
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
            os.fchmod(stream.fileno(), mode)
            stream.write(json.dumps(value, indent=2, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        if simulate_interruption_before_directory_fsync:
            raise ValueError("simulated interruption before JSON directory fsync")
        fsync_directory(path.parent)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def json_temporary_entries(path: Path) -> list[Path]:
    """Return narrowly matched atomic-write temporaries for one JSON path."""

    pattern = re.compile(
        rf"\.{re.escape(path.name)}\.[0-9]+\.[0-9a-f]{{32}}\.tmp"
    )
    return [
        entry for entry in directory_entries(path.parent)
        if pattern.fullmatch(entry.name) is not None
    ]


def recost_transaction_temporary_entries(paths: dict[str, Path]) -> list[Path]:
    """Return narrowly matched orphan journal temporaries."""

    directory = paths["recost_transactions"]
    require_managed_directory(directory, "recost transaction directory")
    pattern = re.compile(r"\..+\.json\.[0-9]+\.[0-9a-f]{32}\.tmp")
    return [
        entry for entry in directory_entries(directory)
        if pattern.fullmatch(entry.name) is not None
    ]


def remove_json_temporaries(entries: list[Path], label: str) -> None:
    """Remove authenticated utility-owned atomic-write temporaries."""

    for entry in entries:
        with regular_descriptor(
            entry,
            label,
            expected_links=1,
        ) as descriptor:
            require_regular_mode_subset(
                os.fstat(descriptor),
                entry,
                label,
                maximum_mode=0o644,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
        unlink_durable(entry)


def mkdir_durable(path: Path) -> None:
    """Create one directory tree and persist each added entry."""

    missing = []
    current = path
    while not current.exists():
        missing.append(current)
        current = current.parent
    for directory in reversed(missing):
        previous = os.umask(0)
        try:
            directory.mkdir(mode=0o755)
        finally:
            os.umask(previous)
        fsync_directory(directory)
        fsync_directory(directory.parent)


def unlink_durable(path: Path) -> None:
    """Remove one entry and persist its parent directory."""

    path.unlink()
    fsync_directory(path.parent)


def entry_exists(path: Path) -> bool:
    """Return whether one directory entry exists without following symlinks."""

    try:
        path.lstat()
    except FileNotFoundError:
        return False
    return True


def require_regular_profile(value: os.stat_result, path: Path, label: str, *,
                            expected_mode: int | None = None,
                            expected_links: int | None = None,
                            expected_uid: int | None = None) -> None:
    """Require one descriptor-backed regular-file profile."""

    if not stat.S_ISREG(value.st_mode):
        raise ValueError(f"{label} is not a regular file: {path}")
    mode = stat.S_IMODE(value.st_mode)
    if expected_mode is not None and mode != expected_mode:
        raise ValueError(
            f"{label} mode is {mode:04o}, expected {expected_mode:04o}: {path}"
        )
    if expected_links is not None and value.st_nlink != expected_links:
        raise ValueError(
            f"{label} has {value.st_nlink} links, expected {expected_links}: {path}"
        )
    if expected_uid is not None and value.st_uid != expected_uid:
        raise ValueError(
            f"{label} owner is UID {value.st_uid}, expected {expected_uid}: {path}"
        )


def require_regular_mode_subset(value: os.stat_result, path: Path, label: str, *,
                                maximum_mode: int,
                                expected_links: int,
                                expected_uid: int | None = None) -> None:
    """Require a retained file with no permissions beyond a selected cap."""

    require_regular_profile(
        value,
        path,
        label,
        expected_links=expected_links,
        expected_uid=expected_uid,
    )
    mode = stat.S_IMODE(value.st_mode)
    if mode & ~maximum_mode:
        raise ValueError(
            f"{label} mode is {mode:04o}, exceeds {maximum_mode:04o}: {path}"
        )


@contextmanager
def regular_descriptor(path: Path, label: str, *,
                       flags: int = os.O_RDONLY,
                       expected_mode: int | None = None,
                       expected_links: int | None = None):
    """Yield one O_NOFOLLOW descriptor after validating its regular profile."""

    try:
        descriptor = os.open(path, flags | os.O_NOFOLLOW)
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    try:
        require_regular_profile(
            os.fstat(descriptor),
            path,
            label,
            expected_mode=expected_mode,
            expected_links=expected_links,
        )
        yield descriptor
    finally:
        os.close(descriptor)


def require_file_sha256(path: Path, expected: str, label: str, *,
                        expected_mode: int | None = None,
                        expected_links: int | None = None) -> bytes:
    """Return descriptor-read bytes after authenticating one retained file."""

    with regular_descriptor(
        path,
        label,
        expected_mode=expected_mode,
        expected_links=expected_links,
    ) as descriptor:
        value = read_descriptor_bytes(descriptor)
    if sha256_bytes(value) != require_sha256(expected, f"{label} SHA-256"):
        raise ValueError(f"{label} checksum has changed: {path}")
    return value


def safe_relative_path(value: str, label: str) -> Path:
    """Require a relative path without traversal."""

    path = Path(value)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise ValueError(f"{label} must be a relative path without '..': {value}")
    return path


def absolute_path(value: Path) -> Path:
    """Return one lexical absolute path without resolving symlinks."""

    return Path(os.path.abspath(os.path.expanduser(str(value))))


def require_offline_path(path: Path, label: str) -> Path:
    """Require offline fixture paths beneath an explicit local prefix."""

    path = absolute_path(path)
    for allowed in OFFLINE_ALLOWED_PREFIXES:
        allowed = absolute_path(allowed)
        if path == allowed or allowed in path.parents:
            return path
    allowed = ", ".join(str(absolute_path(path)) for path in OFFLINE_ALLOWED_PREFIXES)
    raise ValueError(f"{label} must be under local fixture storage {allowed}: {path}")


@contextmanager
def absolute_descriptor(path: Path, label: str, *, flags: int):
    """Open an absolute path without accepting symlinks in any component."""

    path = absolute_path(path)
    if not path.is_absolute() or len(path.parts) < 2:
        raise ValueError(f"{label} path is invalid: {path}")
    parent = os.open("/", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [parent]
    try:
        for component in path.parts[1:-1]:
            parent = os.open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            opened.append(parent)
        descriptor = os.open(path.parts[-1], flags | os.O_NOFOLLOW, dir_fd=parent)
        opened.append(descriptor)
        yield descriptor
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)


def require_confined_file_sha256(root: Path, relative_value: str, expected: str,
                                 label: str, *, expected_mode: int) -> Path:
    """Authenticate one in-root file without accepting symlink traversal."""

    relative = safe_relative_path(relative_value, f"{label} relative path")
    parent = os.open(root, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    opened = [parent]
    try:
        for component in relative.parts[:-1]:
            parent = os.open(
                component,
                os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                dir_fd=parent,
            )
            opened.append(parent)
        descriptor = os.open(
            relative.parts[-1],
            os.O_RDONLY | os.O_NOFOLLOW,
            dir_fd=parent,
        )
        opened.append(descriptor)
        require_regular_profile(
            os.fstat(descriptor),
            root / relative,
            label,
            expected_mode=expected_mode,
            expected_links=1,
        )
        if sha256_descriptor(descriptor) != require_sha256(expected, f"{label} SHA-256"):
            raise ValueError(f"{label} checksum has changed: {root / relative}")
    finally:
        for descriptor in reversed(opened):
            os.close(descriptor)
    return root / relative


def require_artifact_name(value: str) -> str:
    """Require a basename suitable for exact namespace checks."""

    if ARTIFACT_NAME_PATTERN.fullmatch(value) is None:
        raise ValueError(
            "artifact name must be a JSON basename using ASCII letters, digits, "
            "dots, underscores, or hyphens"
        )
    return value


def require_root(root: Path, allow_local_root: bool) -> tuple[Path, bool]:
    """Select the canonical root or an explicitly enabled offline fixture."""

    requested = absolute_path(root)
    canonical = absolute_path(DEFAULT_ROOT)
    if requested == canonical:
        if allow_local_root:
            raise ValueError("--allow-local-root may not relax the canonical root")
        if root.expanduser().resolve() != canonical:
            raise ValueError(f"canonical Stage I root resolves unexpectedly: {root}")
        return canonical, False
    if not allow_local_root:
        raise ValueError(
            f"Stage I root must be {canonical}; use --allow-local-root only "
            "for offline fixtures"
        )
    requested = require_offline_path(requested, "offline fixture root")
    with absolute_descriptor(
        requested, "offline fixture root", flags=os.O_RDONLY | os.O_DIRECTORY
    ):
        pass
    return requested, True


def layout(root: Path, artifact_name: str) -> dict[str, Path]:
    """Return companion and Stage I retained paths."""

    accounting = root / "accounting"
    canonical = accounting / require_artifact_name(artifact_name)
    return {
        "root": root,
        "accounting": accounting,
        "stage_i_transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_transactions"
        ),
        "recost_transactions": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_transactions"
        ),
        "recost_forensics_root": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_forensics"
        ),
        "recost_forensics": (
            accounting / f"mks24_stage_i_{EXECUTION_EPOCH_SLUG}_recost_forensics"
            / canonical.name
        ),
        "lock": root / f".mks24_stage_i_{EXECUTION_EPOCH_SLUG}.lock",
        "canonical": canonical,
        "staged": accounting / f"{artifact_name}.staged",
        "audit": accounting / f"{artifact_name}.publication_audit.json",
    }


def require_directory(path: Path, label: str) -> None:
    """Require one retained directory without accepting a symlink."""

    try:
        value = path.lstat()
    except FileNotFoundError as error:
        raise ValueError(f"{label} is missing: {path}") from error
    if not stat.S_ISDIR(value.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")


def require_trusted_directory(path: Path, label: str) -> None:
    """Require one owned directory without unsafe replacement permissions."""

    require_directory(path, label)
    require_trusted_directory_profile(path.lstat(), path, label)


def require_trusted_directory_profile(value: os.stat_result, path: Path,
                                      label: str) -> None:
    """Require one descriptor- or pathname-backed trusted directory profile."""

    if not stat.S_ISDIR(value.st_mode):
        raise ValueError(f"{label} is not a directory: {path}")
    if value.st_uid != os.geteuid():
        raise ValueError(
            f"{label} owner is UID {value.st_uid}, expected {os.geteuid()}: {path}"
        )
    mode = stat.S_IMODE(value.st_mode)
    allowed = 0o755 | stat.S_ISGID
    if mode & ~allowed:
        raise ValueError(
            f"{label} mode is {mode:04o}, exceeds trusted profile 0755: {path}"
        )


def require_managed_directory(path: Path, label: str) -> None:
    """Require one companion-managed trusted directory."""

    require_trusted_directory(path, label)


def directory_entries(path: Path) -> list[Path]:
    """Return sorted direct children of one directory."""

    return sorted(path.iterdir(), key=lambda item: item.name)


def require_empty_directory(path: Path, label: str, *,
                            allow_absent: bool = False,
                            managed: bool = False,
                            trusted: bool = False) -> None:
    """Require an empty retained directory."""

    if not entry_exists(path) and allow_absent:
        return
    if managed or trusted:
        require_managed_directory(path, label)
    else:
        require_directory(path, label)
    entries = directory_entries(path)
    if entries:
        raise ValueError(
            f"{label} is not empty: " + ", ".join(str(item) for item in entries)
        )


def require_exact_directory_entries(path: Path, label: str,
                                    expected: set[Path], *,
                                    managed: bool = False) -> None:
    """Require one retained directory to contain exactly the selected entries."""

    if managed:
        require_managed_directory(path, label)
    else:
        require_directory(path, label)
    retained = set(directory_entries(path))
    if retained != expected:
        raise ValueError(
            f"{label} entries differ; found: "
            + ", ".join(str(item) for item in sorted(retained))
        )


def artifact_namespace_entries(paths: dict[str, Path]) -> list[Path]:
    """Return direct accounting entries related to one artifact basename."""

    accounting = paths["accounting"]
    require_trusted_directory(accounting, "accounting directory")
    name = paths["canonical"].name
    return [
        path for path in directory_entries(accounting)
        if (
            path.name == name
            or path.name.startswith(f"{name}.")
            or path.name == f".{name}"
            or path.name.startswith(f".{name}.")
        )
    ]


def require_artifact_namespace(paths: dict[str, Path], state: str) -> None:
    """Require an exact staged-only or canonical-only artifact namespace."""

    if state == "staged":
        expected = {paths["staged"]}
    elif state == "promoted":
        expected = {paths["canonical"], paths["audit"]}
    elif state == "linked-pair":
        expected = {paths["staged"], paths["canonical"]}
    elif state == "canonical-pending-audit":
        expected = {paths["canonical"]}
    else:
        raise ValueError(f"unsupported namespace state: {state}")
    retained = set(artifact_namespace_entries(paths))
    if retained != expected:
        raise ValueError(
            f"recost namespace is not {state}-only; found: "
            + ", ".join(str(path) for path in sorted(retained))
        )


def json_pointer(value: object, pointer: str, label: str) -> object:
    """Resolve one RFC 6901 JSON pointer."""

    if pointer == "":
        return value
    if not pointer.startswith("/"):
        raise ValueError(f"{label} JSON pointer must start with '/': {pointer}")
    current = value
    for raw in pointer[1:].split("/"):
        token = raw.replace("~1", "/").replace("~0", "~")
        if isinstance(current, dict) and token in current:
            current = current[token]
            continue
        if isinstance(current, list) and re.fullmatch(r"0|[1-9][0-9]*", token):
            index = int(token)
            if index < len(current):
                current = current[index]
                continue
        raise ValueError(f"{label} JSON pointer is absent: {pointer}")
    return current


def expected_counts(args: argparse.Namespace) -> dict[str, int]:
    """Return the explicitly authorized reconcile counts."""

    counts = {
        key: int(getattr(args, f"expected_{key}"))
        for key in COUNT_KEYS
    }
    if counts["transactions"] != 0:
        raise ValueError("expected Stage I transaction count must be zero")
    if counts["active_reservations"] != 0:
        raise ValueError("expected active reservation count must be zero")
    return counts


def require_selected_counts(value: object, expected: dict[str, int],
                            label: str) -> None:
    """Require the selected artifact or reconcile report counts."""

    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    selected = {key: value.get(key) for key in COUNT_KEYS}
    if selected != expected:
        raise ValueError(f"{label} differ: expected {expected}, found {selected}")


def validate_artifact_payload(retained: bytes, path: Path,
                              args: argparse.Namespace) -> dict[str, object]:
    """Validate authenticated recost-artifact bytes."""

    try:
        value = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"recost artifact is not valid JSON: {path}") from error
    if not isinstance(value, dict):
        raise ValueError(f"recost artifact must be a JSON object: {path}")
    bindings = (
        ("execution epoch", args.artifact_epoch_pointer, EXECUTION_EPOCH),
        (
            "sole next-segment profile",
            args.artifact_profile_pointer,
            args.authorized_next_segment_profile_json,
        ),
        (
            "Stage I helper SHA-256",
            args.artifact_stage_i_sha256_pointer,
            args.expected_stage_i_sha256,
        ),
        (
            "generator SHA-256",
            args.artifact_generator_sha256_pointer,
            args.expected_generator_sha256,
        ),
        (
            "scheduler SHA-256",
            args.artifact_scheduler_sha256_pointer,
            args.expected_scheduler_sha256,
        ),
    )
    for label, pointer, expected in bindings:
        retained = json_pointer(value, pointer, label)
        if retained != expected:
            raise ValueError(
                f"recost artifact {label} differs: expected {expected!r}, "
                f"found {retained!r}"
            )
    require_selected_counts(
        json_pointer(value, args.artifact_counts_pointer, "counts"),
        expected_counts(args),
        "recost artifact counts",
    )
    return value


def read_artifact(path: Path, args: argparse.Namespace, *,
                  expected_links: int = 1) -> dict[str, object]:
    """Authenticate and validate one staged or promoted recost artifact."""

    retained = require_file_sha256(
        path,
        args.expected_artifact_sha256,
        "recost artifact",
        expected_mode=args.expected_artifact_mode,
        expected_links=expected_links,
    )
    return validate_artifact_payload(retained, path, args)


def git_run(root: Path, arguments: list[str], *,
            capture_output: bool = False) -> subprocess.CompletedProcess:
    """Run one Git query against the source repository."""

    environment = {
        key: value for key, value in os.environ.items()
        if not key.startswith("GIT_")
    }
    return subprocess.run(
        [str(GIT), "--no-replace-objects", "-C", str(root), *arguments],
        check=False,
        capture_output=capture_output,
        env=environment,
    )


def authenticate_committed_stage_i(root_dir: Path, expected: str) -> Path:
    """Authenticate the committed live Stage I helper bytes."""

    expected = require_sha256(expected, "Stage I helper SHA-256")
    path = root_dir / STAGE_I_RELATIVE
    require_file_sha256(
        path, expected, "Stage I helper", expected_mode=0o644, expected_links=1
    )
    relative = STAGE_I_RELATIVE.as_posix()
    if git_run(root_dir, ["ls-files", "--error-unmatch", "--", relative]).returncode:
        raise ValueError("Stage I helper is not tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(root_dir, arguments).returncode:
            raise ValueError("Stage I helper must be committed before recost publication")
    retained = git_run(
        root_dir, ["show", f"HEAD:{relative}"], capture_output=True
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError("committed Stage I helper checksum differs from live bytes")
    return path


def authenticate_committed_utility(root_dir: Path, expected: str) -> Path:
    """Authenticate the tracked, clean companion before canonical use."""

    expected = require_sha256(expected, "utility SHA-256")
    path = root_dir / UTILITY_RELATIVE
    require_file_sha256(
        path, expected, "retained utility", expected_mode=0o755, expected_links=1
    )
    relative = UTILITY_RELATIVE.as_posix()
    if git_run(root_dir, ["ls-files", "--error-unmatch", "--", relative]).returncode:
        raise ValueError("retained utility is not tracked by Git")
    for arguments in (
        ["diff", "--quiet", "--", relative],
        ["diff", "--cached", "--quiet", "--", relative],
    ):
        if git_run(root_dir, arguments).returncode:
            raise ValueError("retained utility must be committed before canonical use")
    retained = git_run(
        root_dir, ["show", f"HEAD:{relative}"], capture_output=True
    )
    if retained.returncode or sha256_bytes(retained.stdout) != expected:
        raise ValueError("committed utility checksum differs from live bytes")
    return path


@contextmanager
def stage_i_descriptor(root_dir: Path, expected: str):
    """Yield an authenticated descriptor for the committed Stage I helper."""

    path = authenticate_committed_stage_i(root_dir, expected)
    with regular_descriptor(
        path, "Stage I helper", expected_mode=0o644, expected_links=1
    ) as descriptor:
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("Stage I helper descriptor checksum has changed")
        os.set_inheritable(descriptor, True)
        yield descriptor
        if sha256_descriptor(descriptor) != expected:
            raise ValueError("Stage I helper descriptor checksum changed during reconcile")


def run_reconcile(root_dir: Path, root: Path, offline: bool,
                  args: argparse.Namespace) -> dict[str, object]:
    """Execute reconcile from the authenticated Stage I helper descriptor."""

    with stage_i_descriptor(root_dir, args.expected_stage_i_sha256) as descriptor:
        command = [
            sys.executable,
            f"/proc/self/fd/{descriptor}",
            "--root",
            str(root),
        ]
        if offline:
            command.append("--allow-local-root")
        command.append("reconcile")
        completed = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            pass_fds=(descriptor,),
        )
    if completed.returncode:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise ValueError(f"authenticated Stage I reconcile failed: {detail}")
    try:
        report = json.loads(completed.stdout)
    except json.JSONDecodeError as error:
        raise ValueError("authenticated Stage I reconcile did not return JSON") from error
    if not isinstance(report, dict) or report.get("consistent") is not True:
        raise ValueError(f"authenticated Stage I reconcile is inconsistent: {report!r}")
    if report.get("execution_epoch") != EXECUTION_EPOCH:
        raise ValueError("authenticated Stage I reconcile returned the wrong epoch")
    require_selected_counts(
        report.get("counts"), expected_counts(args), "Stage I reconcile counts"
    )
    return report


def scheduler_environment() -> dict[str, str]:
    """Return inherited process state without caller-selected Slurm routing."""

    return {
        key: value for key, value in os.environ.items()
        if not key.startswith("SLURM_")
    }


def require_empty_queue(root: Path, offline: bool, queue_file: str | None) -> None:
    """Require no queued CGL workflow job, with fixture injection only offline."""

    if offline and queue_file is None:
        raise ValueError("offline fixture root requires --squeue-file")
    if queue_file is not None:
        if not offline:
            raise ValueError("--squeue-file is permitted only with a local fixture root")
        path = require_offline_path(Path(queue_file), "queue fixture")
        with absolute_descriptor(path, "queue fixture", flags=os.O_RDONLY) as descriptor:
            require_regular_profile(os.fstat(descriptor), path, "queue fixture")
            output = read_descriptor_bytes(descriptor).decode("utf-8")
    else:
        try:
            user = pwd.getpwuid(os.geteuid()).pw_name
        except KeyError as error:
            raise ValueError("effective UID has no account name; cannot check the queue") from error
        try:
            output = subprocess.run(
                [str(SQUEUE), "-h", "-u", user, "-o", "%i|%j|%T"],
                check=True,
                capture_output=True,
                text=True,
                env=scheduler_environment(),
            ).stdout
        except (FileNotFoundError, subprocess.CalledProcessError) as error:
            raise ValueError("squeue is unavailable; refusing recost publication") from error
    queued = []
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split("|")
        if len(fields) != 3:
            raise ValueError("squeue output row is not exactly three fields: " + line)
        if any(not field or field.strip() != field for field in fields):
            raise ValueError("squeue output row has malformed fields: " + line)
        if fields[1].startswith(CGL_JOB_NAME_PREFIX):
            queued.append(line)
    if queued:
        raise ValueError("another CGL job is queued: " + "; ".join(queued))


def validate_fixture_options(args: argparse.Namespace, offline: bool) -> None:
    """Keep mutation simulation and queue injection out of canonical use."""

    fixture_values = (
        getattr(args, "squeue_file", None),
        getattr(args, "pre_link_squeue_file", None),
        getattr(args, "pre_unlink_squeue_file", None),
        getattr(args, "simulate_post_link_failure", False),
        getattr(args, "simulate_link_return_failure", False),
        getattr(args, "simulate_adoption_interruption_after_journal", False),
        getattr(
            args,
            "simulate_adoption_interruption_before_forensic_directory_fsync",
            False,
        ),
        getattr(args, "simulate_adoption_interruption_after_forensic", False),
        getattr(
            args,
            "simulate_adoption_interruption_before_audit_directory_fsync",
            False,
        ),
        getattr(args, "simulate_adoption_interruption_after_audit", False),
    )
    if any(value not in (None, False) for value in fixture_values) and not offline:
        raise ValueError("fixture injection options are permitted only for a local root")


def authenticate_external_files(root: Path, args: argparse.Namespace
                                ) -> tuple[Path, Path]:
    """Authenticate retained generator and scheduler evidence files."""

    generator = require_confined_file_sha256(
        root,
        args.generator_relative_path,
        args.expected_generator_sha256,
        "recost generator",
        expected_mode=args.expected_generator_mode,
    )
    scheduler = require_confined_file_sha256(
        root,
        args.scheduler_relative_path,
        args.expected_scheduler_sha256,
        "scheduler evidence",
        expected_mode=args.expected_scheduler_mode,
    )
    return generator, scheduler


@contextmanager
def promotion_lock(paths: dict[str, Path]):
    """Take the nonblocking cooperative Stage I root lock."""

    root = paths["root"]
    require_trusted_directory(root, "Stage I root")
    created = False
    try:
        descriptor = create_exclusive(
            paths["lock"],
            os.O_RDWR,
            0o644,
        )
        created = True
    except FileExistsError:
        descriptor = os.open(paths["lock"], os.O_RDWR | os.O_NOFOLLOW)
    try:
        require_regular_mode_subset(
            os.fstat(descriptor),
            paths["lock"],
            "Stage I lock",
            maximum_mode=0o644,
            expected_links=1,
            expected_uid=os.geteuid(),
        )
    except BaseException:
        os.close(descriptor)
        raise
    if created:
        fsync_directory(root)
    stream = os.fdopen(descriptor, "a+", encoding="utf-8")
    try:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as error:
            if error.errno not in (errno.EACCES, errno.EAGAIN):
                raise
            raise ValueError(f"another Stage I mutation holds {paths['lock']}") from error
        if stat.S_IMODE(os.fstat(stream.fileno()).st_mode) != 0o644:
            os.fchmod(stream.fileno(), 0o644)
            os.fsync(stream.fileno())
            fsync_directory(root)
        yield
    finally:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
        finally:
            stream.close()


def copy_forensic(source: Path, destination: Path, expected: str, *,
                  expected_mode: int,
                  source_label: str = "staged recost artifact",
                  simulate_interruption_before_directory_fsync: bool = False) -> None:
    """Retain one independent authenticated forensic copy."""

    source_descriptor = os.open(source, os.O_RDONLY | os.O_NOFOLLOW)
    try:
        require_regular_profile(
            os.fstat(source_descriptor),
            source,
            source_label,
            expected_mode=expected_mode,
            expected_links=1,
        )
        destination_descriptor = create_exclusive(
            destination,
            os.O_WRONLY,
            0o444,
        )
        try:
            os.fchmod(destination_descriptor, 0o444)
            while True:
                block = os.read(source_descriptor, 1024 * 1024)
                if not block:
                    break
                view = memoryview(block)
                while view:
                    written = os.write(destination_descriptor, view)
                    view = view[written:]
            os.fsync(destination_descriptor)
        finally:
            os.close(destination_descriptor)
    finally:
        os.close(source_descriptor)
    if simulate_interruption_before_directory_fsync:
        raise ValueError("simulated interruption before forensic directory fsync")
    fsync_directory(destination.parent)
    require_file_sha256(
        destination,
        expected,
        "recost forensic copy",
        expected_mode=0o444,
        expected_links=1,
    )


def recover_or_copy_forensic(source: Path, destination: Path, expected: str, *,
                             expected_mode: int,
                             simulate_interruption_before_directory_fsync: bool
                             ) -> None:
    """Repair one interrupted provisional forensic copy, then persist it."""

    if entry_exists(destination):
        with regular_descriptor(
            destination,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        ) as descriptor:
            profile = os.fstat(descriptor)
            require_regular_profile(
                profile,
                destination,
                "recost forensic copy",
                expected_mode=0o444,
                expected_links=1,
                expected_uid=os.geteuid(),
            )
            retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != require_sha256(expected, "recost forensic copy SHA-256"):
            unlink_durable(destination)
    if not entry_exists(destination):
        copy_forensic(
            source,
            destination,
            expected,
            expected_mode=expected_mode,
            source_label="canonical recost artifact",
            simulate_interruption_before_directory_fsync=(
                simulate_interruption_before_directory_fsync
            ),
        )
    require_file_sha256(
        destination,
        expected,
        "recost forensic copy",
        expected_mode=0o444,
        expected_links=1,
    )
    fsync_directory(destination.parent)


def descriptor_identity(descriptor: int) -> tuple[int, int]:
    """Return the stable device/inode identity for one open descriptor."""

    value = os.fstat(descriptor)
    return value.st_dev, value.st_ino


def require_named_artifact_identity(path: Path, args: argparse.Namespace, *,
                                    expected_links: int,
                                    expected_identity: tuple[int, int]) -> None:
    """Require one artifact name still references the selected inode."""

    read_artifact(path, args, expected_links=expected_links)
    with regular_descriptor(
        path,
        "recost artifact",
        expected_mode=args.expected_artifact_mode,
        expected_links=expected_links,
    ) as descriptor:
        if descriptor_identity(descriptor) != expected_identity:
            raise ValueError(f"recost artifact inode identity has changed: {path}")


def link_after_empty_queue(root: Path, offline: bool, queue_file: str | None,
                           directory_descriptor: int,
                           paths: dict[str, Path],
                           args: argparse.Namespace,
                           expected_identity: tuple[int, int]) -> None:
    """Check the live queue and immediately attempt same-directory publication."""

    require_artifact_namespace_in_directory(directory_descriptor, paths, "staged")
    require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["staged"].name,
        paths["staged"],
        args,
        expected_links=1,
        expected_identity=expected_identity,
    )
    require_empty_queue(root, offline, queue_file)
    try:
        os.link(
            paths["staged"].name,
            paths["canonical"].name,
            src_dir_fd=directory_descriptor,
            dst_dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
    except OSError as error:
        raise LinkAttemptError(f"recost publication link failed: {error}") from error
    if args.simulate_link_return_failure:
        raise LinkAttemptError("simulated failure immediately after recost publication link")


def unlink_linked_pair_after_empty_queue(root: Path, offline: bool,
                                         queue_file: str | None,
                                         paths: dict[str, Path],
                                         args: argparse.Namespace,
                                         before_unlink=None) -> None:
    """Revalidate one linked pair and immediately remove its staged entry."""

    if before_unlink is not None:
        before_unlink()
    directory_descriptor = os.open(
        paths["accounting"],
        os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
    )
    try:
        require_trusted_directory_profile(
            os.fstat(directory_descriptor),
            paths["accounting"],
            "accounting directory",
        )
        require_linked_pair_in_directory(directory_descriptor, paths, args)
        require_empty_queue(root, offline, queue_file)
        os.unlink(paths["staged"].name, dir_fd=directory_descriptor)
        os.fsync(directory_descriptor)
    finally:
        os.close(directory_descriptor)


def audit_record(paths: dict[str, Path], args: argparse.Namespace,
                 source: Path, generator: Path, scheduler: Path,
                 forensic: Path, transaction_id: str) -> dict[str, object]:
    """Build the finalized durable publication audit."""

    return {
        "schema_version": 1,
        "record_type": "observed-publication",
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "published_utc": utc_now(),
        "artifact": {
            "path": str(paths["canonical"]),
            "sha256": args.expected_artifact_sha256,
            "mode": f"{args.expected_artifact_mode:04o}",
            "links": 1,
        },
        "authorized_sole_next_segment_profile": (
            args.authorized_next_segment_profile_json
        ),
        "counts": expected_counts(args),
        "generator": {
            "path": str(generator),
            "sha256": args.expected_generator_sha256,
            "mode": f"{args.expected_generator_mode:04o}",
        },
        "scheduler_evidence": {
            "path": str(scheduler),
            "sha256": args.expected_scheduler_sha256,
            "mode": f"{args.expected_scheduler_mode:04o}",
        },
        "stage_i_helper": {
            "path": str(source.parent / "cgl_lf_stage_i.py"),
            "sha256": args.expected_stage_i_sha256,
            "committed": True,
            "reconcile_execution": "descriptor",
        },
        "utility": {
            "path": str(source),
            "sha256": args.expected_utility_sha256,
            "execution": "authenticated-descriptor",
            "committed": not args.allow_local_root,
        },
        "forensic_copy": {
            "path": str(forensic),
            "sha256": args.expected_artifact_sha256,
            "mode": "0444",
            "links": 1,
        },
        "publication": "same-directory-link-fsync-unlink-fsync",
    }


def adoption_audit_record(paths: dict[str, Path], args: argparse.Namespace,
                          source: Path, generator: Path, scheduler: Path,
                          forensic: Path, transaction_id: str) -> dict[str, object]:
    """Build a present-time audit for a legacy canonical-only artifact."""

    return {
        "schema_version": 1,
        "record_type": "legacy-canonical-adoption",
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "adopted_utc": utc_now(),
        "artifact": {
            "path": str(paths["canonical"]),
            "sha256": args.expected_artifact_sha256,
            "mode": f"{args.expected_artifact_mode:04o}",
            "links": 1,
        },
        "authorized_sole_next_segment_profile": (
            args.authorized_next_segment_profile_json
        ),
        "counts": expected_counts(args),
        "generator": {
            "path": str(generator),
            "sha256": args.expected_generator_sha256,
            "mode": f"{args.expected_generator_mode:04o}",
        },
        "scheduler_evidence": {
            "path": str(scheduler),
            "sha256": args.expected_scheduler_sha256,
            "mode": f"{args.expected_scheduler_mode:04o}",
        },
        "stage_i_helper": {
            "path": str(source.parent / "cgl_lf_stage_i.py"),
            "sha256": args.expected_stage_i_sha256,
            "committed": True,
            "reconcile_execution": "descriptor",
        },
        "utility": {
            "path": str(source),
            "sha256": args.expected_utility_sha256,
            "execution": "authenticated-descriptor",
            "committed": not args.allow_local_root,
        },
        "forensic_copy": {
            "path": str(forensic),
            "sha256": args.expected_artifact_sha256,
            "mode": "0444",
            "links": 1,
        },
        "original_publication_transition_observed": False,
        "original_publication_method": "unknown",
        "original_publisher": "unknown",
        "present_authentication": {
            "canonical_artifact": "authenticated-under-stage-i-lock",
            "empty_user_queue": "required-under-stage-i-lock",
            "reconciliation": "clean-required-under-stage-i-lock",
        },
    }


def validate_audit(paths: dict[str, Path], args: argparse.Namespace,
                   source: Path, generator: Path, scheduler: Path) -> None:
    """Validate the finalized audit and its retained forensic copy."""

    with regular_descriptor(
        paths["audit"],
        "recost publication audit",
        expected_mode=0o644,
        expected_links=1,
    ) as descriptor:
        retained = read_descriptor_bytes(descriptor)
    try:
        audit = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("recost publication audit is not valid JSON") from error
    if not isinstance(audit, dict):
        raise ValueError("recost publication audit must be an object")
    forensic_record = audit.get("forensic_copy")
    if not isinstance(forensic_record, dict):
        raise ValueError("recost publication audit forensic binding is invalid")
    transaction_id = audit.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("recost publication audit transaction ID is invalid")
    forensic = Path(str(forensic_record.get("path", "")))
    expected_forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    if forensic != expected_forensic:
        raise ValueError("recost forensic copy path does not match its transaction")
    record_type = audit.get("record_type")
    if record_type == "observed-publication":
        expected = audit_record(
            paths, args, source, generator, scheduler, forensic, transaction_id
        )
        timestamp = "published_utc"
        timestamp_label = "publication audit timestamp"
    elif record_type == "legacy-canonical-adoption":
        expected = adoption_audit_record(
            paths, args, source, generator, scheduler, forensic, transaction_id
        )
        timestamp = "adopted_utc"
        timestamp_label = "adoption audit timestamp"
    else:
        raise ValueError("recost publication audit record type is invalid")
    if frozenset(audit) != frozenset(expected):
        raise ValueError("recost publication audit schema differs")
    require_utc_timestamp(audit.get(timestamp), timestamp_label)
    for key in expected:
        if key == timestamp:
            continue
        if audit.get(key) != expected[key]:
            raise ValueError(f"recost publication audit {key} binding differs")
    require_exact_directory_entries(
        paths["recost_forensics"],
        "recost forensic directory",
        {forensic},
        managed=True,
    )
    require_file_sha256(
        forensic,
        args.expected_artifact_sha256,
        "recost forensic copy",
        expected_mode=0o444,
        expected_links=1,
    )


def verify_staged_state(paths: dict[str, Path], root_dir: Path, root: Path,
                        offline: bool, args: argparse.Namespace) -> None:
    """Validate a staged-only recost boundary without mutation."""

    validate_fixture_options(args, offline)
    require_empty_queue(root, offline, args.squeue_file)
    require_empty_directory(
        paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
    )
    require_empty_directory(
        paths["recost_transactions"],
        "recost transaction directory",
        allow_absent=True,
        managed=True,
    )
    require_artifact_namespace(paths, "staged")
    authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
    authenticate_external_files(root, args)
    read_artifact(paths["staged"], args)
    run_reconcile(root_dir, root, offline, args)


def verify_promoted_state(paths: dict[str, Path], root_dir: Path, root: Path,
                          offline: bool, args: argparse.Namespace) -> None:
    """Validate the canonical-only post-publication boundary."""

    validate_fixture_options(args, offline)
    require_empty_queue(root, offline, args.squeue_file)
    require_empty_directory(
        paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
    )
    require_empty_directory(
        paths["recost_transactions"], "recost transaction directory", managed=True
    )
    require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
    require_managed_directory(paths["recost_forensics"], "recost forensic directory")
    require_artifact_namespace(paths, "promoted")
    authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
    generator, scheduler = authenticate_external_files(root, args)
    read_artifact(paths["canonical"], args)
    validate_audit(paths, args, initial_source_path(), generator, scheduler)
    run_reconcile(root_dir, root, offline, args)


def remove_prelink_records(journal: Path, forensic: Path) -> None:
    """Remove retry-safe records when publication was never attempted."""

    if entry_exists(forensic):
        unlink_durable(forensic)
    if entry_exists(journal):
        unlink_durable(journal)


def recost_journal(paths: dict[str, Path]) -> tuple[Path, dict[str, object]]:
    """Read the sole durable companion journal for linked-pair recovery."""

    require_managed_directory(paths["recost_transactions"], "recost transaction directory")
    entries = directory_entries(paths["recost_transactions"])
    if len(entries) != 1 or entries[0].suffix != ".json":
        raise ValueError(
            "linked-pair recovery requires exactly one recost journal: "
            + ", ".join(str(path) for path in entries)
        )
    journal = entries[0]
    with regular_descriptor(
        journal, "recost transaction journal", expected_mode=0o644, expected_links=1
    ) as descriptor:
        retained = read_descriptor_bytes(descriptor)
    try:
        record = json.loads(retained)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError(f"recost transaction journal is not valid JSON: {journal}") from error
    if not isinstance(record, dict):
        raise ValueError(f"recost transaction journal must be an object: {journal}")
    if record.get("transaction_id") != journal.stem:
        raise ValueError("recost transaction journal ID differs from its filename")
    return journal, record


def recovery_bindings(paths: dict[str, Path], args: argparse.Namespace,
                      transaction_id: str, forensic: Path) -> dict[str, object]:
    """Return every authorization-relevant immutable journal binding."""

    return {
        "schema_version": 1,
        "execution_epoch": EXECUTION_EPOCH,
        "transaction_id": transaction_id,
        "staged_path": str(paths["staged"]),
        "canonical_path": str(paths["canonical"]),
        "artifact_sha256": args.expected_artifact_sha256,
        "artifact_mode": f"{args.expected_artifact_mode:04o}",
        "artifact_epoch_pointer": args.artifact_epoch_pointer,
        "artifact_profile_pointer": args.artifact_profile_pointer,
        "artifact_counts_pointer": args.artifact_counts_pointer,
        "artifact_stage_i_sha256_pointer": args.artifact_stage_i_sha256_pointer,
        "artifact_generator_sha256_pointer": args.artifact_generator_sha256_pointer,
        "artifact_scheduler_sha256_pointer": args.artifact_scheduler_sha256_pointer,
        "authorized_sole_next_segment_profile": (
            args.authorized_next_segment_profile_json
        ),
        "counts": expected_counts(args),
        "forensic_path": str(forensic),
        "utility_sha256": args.expected_utility_sha256,
        "stage_i_helper_sha256": args.expected_stage_i_sha256,
        "generator_relative_path": args.generator_relative_path,
        "generator_sha256": args.expected_generator_sha256,
        "generator_mode": f"{args.expected_generator_mode:04o}",
        "scheduler_relative_path": args.scheduler_relative_path,
        "scheduler_sha256": args.expected_scheduler_sha256,
        "scheduler_mode": f"{args.expected_scheduler_mode:04o}",
    }


def adoption_bindings(paths: dict[str, Path], args: argparse.Namespace,
                      transaction_id: str, forensic: Path) -> dict[str, object]:
    """Return immutable bindings for one legacy canonical adoption."""

    return {
        **recovery_bindings(paths, args, transaction_id, forensic),
        "operation": "legacy-canonical-adoption",
    }


def validate_adoption_journal(paths: dict[str, Path], args: argparse.Namespace,
                              record: dict[str, object]) -> tuple[str, Path, str]:
    """Validate one resumable legacy canonical adoption journal."""

    transaction_id = record.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("legacy adoption journal transaction ID is invalid")
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    expected = adoption_bindings(paths, args, transaction_id, forensic)
    state = record.get("state")
    if state not in {"preparing", "forensic-copied", "audit-written"}:
        raise ValueError("legacy adoption journal is not in an allowed recovery state")
    if frozenset(record) != frozenset({*expected, "state", "created_utc"}):
        raise ValueError("legacy adoption journal schema differs")
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"legacy adoption journal {key} binding differs")
    require_utc_timestamp(record.get("created_utc"), "legacy adoption journal timestamp")
    return transaction_id, forensic, state


def validate_recovery_journal(paths: dict[str, Path], args: argparse.Namespace,
                              record: dict[str, object], *,
                              allowed_states: set[str] | None = None,
                              require_forensic_sha256: bool = True
                              ) -> tuple[str, Path]:
    """Validate immutable recovery bindings and its forensic copy."""

    transaction_id = record.get("transaction_id")
    if not isinstance(transaction_id, str) or not transaction_id:
        raise ValueError("recost recovery journal transaction ID is invalid")
    forensic = (
        paths["recost_forensics"]
        / f"{transaction_id}.{paths['canonical'].name}.forensic"
    )
    expected = recovery_bindings(paths, args, transaction_id, forensic)
    state = record.get("state")
    if allowed_states is None:
        allowed_states = {"link-pending", "ambiguous-after-link-attempt"}
    if state not in allowed_states:
        raise ValueError("recost recovery journal is not in an allowed recovery state")
    keys = {*expected, "state", "created_utc"}
    if state == "ambiguous-after-link-attempt":
        keys.add("ambiguity_recorded_utc")
    if frozenset(record) != frozenset(keys):
        raise ValueError("recost recovery journal schema differs")
    for key, value in expected.items():
        if record.get(key) != value:
            raise ValueError(f"recost recovery journal {key} binding differs")
    require_utc_timestamp(record.get("created_utc"), "recost journal creation timestamp")
    if state == "ambiguous-after-link-attempt":
        require_utc_timestamp(
            record.get("ambiguity_recorded_utc"), "recost journal ambiguity timestamp"
        )
    require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
    require_exact_directory_entries(
        paths["recost_forensics"],
        "recost forensic directory",
        {forensic},
        managed=True,
    )
    if require_forensic_sha256:
        require_file_sha256(
            forensic,
            args.expected_artifact_sha256,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        )
    else:
        with regular_descriptor(
            forensic,
            "recost forensic copy",
            expected_mode=0o444,
            expected_links=1,
        ):
            pass
    return transaction_id, forensic


def require_linked_pair(paths: dict[str, Path], args: argparse.Namespace, *,
                        expected_identity: tuple[int, int] | None = None) -> None:
    """Require staged and canonical names to be one authenticated two-link inode."""

    require_artifact_namespace(paths, "linked-pair")
    read_artifact(paths["staged"], args, expected_links=2)
    read_artifact(paths["canonical"], args, expected_links=2)
    with regular_descriptor(
        paths["staged"], "staged recost artifact", expected_links=2
    ) as staged_descriptor:
        with regular_descriptor(
            paths["canonical"], "canonical recost artifact", expected_links=2
        ) as canonical_descriptor:
            staged = os.fstat(staged_descriptor)
            canonical = os.fstat(canonical_descriptor)
            if (staged.st_dev, staged.st_ino) != (canonical.st_dev, canonical.st_ino):
                raise ValueError(
                    "staged and canonical recost entries are not one linked pair"
                )
            if (
                expected_identity is not None
                and (staged.st_dev, staged.st_ino) != expected_identity
            ):
                raise ValueError("linked recost pair differs from authenticated staged inode")


def require_artifact_namespace_in_directory(directory_descriptor: int,
                                            paths: dict[str, Path],
                                            state: str) -> None:
    """Require one artifact namespace through a stable accounting descriptor."""

    if state == "staged":
        expected = {paths["staged"].name}
    elif state == "linked-pair":
        expected = {paths["staged"].name, paths["canonical"].name}
    else:
        raise ValueError(f"unsupported descriptor-relative namespace state: {state}")
    name = paths["canonical"].name
    retained = {
        entry for entry in os.listdir(directory_descriptor)
        if (
            entry == name
            or entry.startswith(f"{name}.")
            or entry == f".{name}"
            or entry.startswith(f".{name}.")
        )
    }
    if retained != expected:
        raise ValueError(
            f"descriptor-relative recost namespace is not {state}-only; found: "
            + ", ".join(sorted(retained))
        )


def require_named_artifact_identity_in_directory(
    directory_descriptor: int,
    name: str,
    path: Path,
    args: argparse.Namespace,
    *,
    expected_links: int,
    expected_identity: tuple[int, int] | None = None,
) -> tuple[int, int]:
    """Authenticate one artifact name through a stable accounting descriptor."""

    descriptor = os.open(name, os.O_RDONLY | os.O_NOFOLLOW, dir_fd=directory_descriptor)
    try:
        value = os.fstat(descriptor)
        require_regular_profile(
            value,
            path,
            "recost artifact",
            expected_mode=args.expected_artifact_mode,
            expected_links=expected_links,
        )
        retained = read_descriptor_bytes(descriptor)
        if sha256_bytes(retained) != args.expected_artifact_sha256:
            raise ValueError(f"recost artifact checksum has changed: {path}")
        validate_artifact_payload(retained, path, args)
        identity = descriptor_identity(descriptor)
        if expected_identity is not None and identity != expected_identity:
            raise ValueError(f"recost artifact inode identity has changed: {path}")
        return identity
    finally:
        os.close(descriptor)


def require_linked_pair_in_directory(directory_descriptor: int,
                                     paths: dict[str, Path],
                                     args: argparse.Namespace, *,
                                     expected_identity: tuple[int, int] | None = None) -> None:
    """Require one authenticated linked pair through a stable directory."""

    require_artifact_namespace_in_directory(directory_descriptor, paths, "linked-pair")
    staged_identity = require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["staged"].name,
        paths["staged"],
        args,
        expected_links=2,
        expected_identity=expected_identity,
    )
    canonical_identity = require_named_artifact_identity_in_directory(
        directory_descriptor,
        paths["canonical"].name,
        paths["canonical"],
        args,
        expected_links=2,
        expected_identity=expected_identity,
    )
    if staged_identity != canonical_identity:
        raise ValueError("staged and canonical recost entries are not one linked pair")


def link_failure_is_retry_safe(paths: dict[str, Path], args: argparse.Namespace,
                               expected_identity: tuple[int, int]) -> bool:
    """Return whether a failed link attempt provably left staged-only state."""

    try:
        require_artifact_namespace(paths, "staged")
        require_named_artifact_identity(
            paths["staged"],
            args,
            expected_links=1,
            expected_identity=expected_identity,
        )
    except (OSError, TypeError, ValueError):
        return False
    return True


def finalize_linked_pair(paths: dict[str, Path], root_dir: Path, root: Path,
                         offline: bool, args: argparse.Namespace) -> None:
    """Finish an authenticated publication left ambiguous after os.link."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        require_empty_queue(root, offline, args.squeue_file)
        require_empty_directory(
            paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
        )
        authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
        generator, scheduler = authenticate_external_files(root, args)
        run_reconcile(root_dir, root, offline, args)
        remove_json_temporaries(
            recost_transaction_temporary_entries(paths),
            "recost transaction temporary",
        )
        journal, record = recost_journal(paths)
        transaction_id, forensic = validate_recovery_journal(paths, args, record)
        remove_json_temporaries(
            json_temporary_entries(paths["audit"]),
            "recost publication audit temporary",
        )
        retained = set(artifact_namespace_entries(paths))
        if retained == {paths["staged"], paths["canonical"]}:
            require_linked_pair(paths, args)
            validate_recovery_journal(paths, args, record)
            unlink_linked_pair_after_empty_queue(
                root, offline, args.pre_unlink_squeue_file or args.squeue_file,
                paths,
                args,
                before_unlink=lambda: validate_recovery_journal(paths, args, record),
            )
        elif retained == {paths["canonical"]}:
            read_artifact(paths["canonical"], args)
            fsync_directory(paths["accounting"])
        elif retained == {paths["canonical"], paths["audit"]}:
            read_artifact(paths["canonical"], args)
            validate_audit(
                paths,
                args,
                initial_source_path(),
                generator,
                scheduler,
            )
        else:
            raise ValueError(
                "recost recovery namespace is not a linked pair or canonical-only: "
                + ", ".join(str(path) for path in sorted(retained))
            )
        read_artifact(paths["canonical"], args)
        if not paths["audit"].exists():
            require_artifact_namespace_after_unlink(paths)
            write_json(
                paths["audit"],
                audit_record(
                    paths,
                    args,
                    initial_source_path(),
                    generator,
                    scheduler,
                    forensic,
                    transaction_id,
                ),
            )
        validate_audit(paths, args, initial_source_path(), generator, scheduler)
        unlink_durable(journal)
        require_empty_directory(
            paths["recost_transactions"], "recost transaction directory", managed=True
        )
    verify_promoted_state(paths, root_dir, root, offline, args)


def retire_preparing(paths: dict[str, Path], root_dir: Path, root: Path,
                     offline: bool, args: argparse.Namespace) -> None:
    """Retire a staged-only transaction interrupted before publication."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        require_empty_queue(root, offline, args.squeue_file)
        require_empty_directory(
            paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
        )
        authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
        authenticate_external_files(root, args)
        run_reconcile(root_dir, root, offline, args)
        require_artifact_namespace(paths, "staged")
        read_artifact(paths["staged"], args)
        require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
        require_managed_directory(paths["recost_forensics"], "recost forensic directory")
        remove_json_temporaries(
            recost_transaction_temporary_entries(paths),
            "recost transaction temporary",
        )
        entries = directory_entries(paths["recost_transactions"])
        if not entries:
            require_empty_directory(
                paths["recost_forensics"], "recost forensic directory", managed=True
            )
            require_empty_queue(root, offline, args.squeue_file)
            verify_staged_state(paths, root_dir, root, offline, args)
            return
        journal, record = recost_journal(paths)
        transaction_id = record.get("transaction_id")
        if not isinstance(transaction_id, str) or not transaction_id:
            raise ValueError("recost recovery journal transaction ID is invalid")
        forensic = (
            paths["recost_forensics"]
            / f"{transaction_id}.{paths['canonical'].name}.forensic"
        )
        expected = recovery_bindings(paths, args, transaction_id, forensic)
        if frozenset(record) != frozenset({*expected, "state", "created_utc"}):
            raise ValueError("recost prepublication journal schema differs")
        for key, value in expected.items():
            if record.get(key) != value:
                raise ValueError(f"recost prepublication journal {key} binding differs")
        state = record.get("state")
        if state not in {"preparing", "link-pending"}:
            raise ValueError("recost transaction is not an interrupted prepublication state")
        require_utc_timestamp(record.get("created_utc"), "recost journal creation timestamp")
        retained = set(directory_entries(paths["recost_forensics"]))
        if retained not in (set(), {forensic}):
            raise ValueError(
                "recost forensic directory entries differ; found: "
                + ", ".join(str(item) for item in sorted(retained))
            )
        if state == "link-pending":
            if entry_exists(forensic):
                require_file_sha256(
                    forensic,
                    args.expected_artifact_sha256,
                    "recost forensic copy",
                    expected_mode=0o444,
                    expected_links=1,
                )
        elif entry_exists(forensic):
            with regular_descriptor(
                forensic,
                "recost forensic copy",
                expected_links=1,
            ) as descriptor:
                require_regular_mode_subset(
                    os.fstat(descriptor),
                    forensic,
                    "recost forensic copy",
                    maximum_mode=0o444,
                    expected_links=1,
                    expected_uid=os.geteuid(),
                )
        require_empty_queue(root, offline, args.squeue_file)
        remove_prelink_records(journal, forensic)
        require_empty_directory(
            paths["recost_transactions"], "recost transaction directory", managed=True
        )
    verify_staged_state(paths, root_dir, root, offline, args)


def promote(paths: dict[str, Path], root_dir: Path, root: Path, offline: bool,
            args: argparse.Namespace) -> None:
    """Publish one staged recost artifact exactly once."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        verify_staged_state(paths, root_dir, root, offline, args)
        mkdir_durable(paths["recost_transactions"])
        require_empty_directory(
            paths["recost_transactions"], "recost transaction directory", managed=True
        )
        if not paths["recost_forensics_root"].exists():
            mkdir_durable(paths["recost_forensics_root"])
        require_managed_directory(paths["recost_forensics_root"], "recost forensic root directory")
        if not paths["recost_forensics"].exists():
            mkdir_durable(paths["recost_forensics"])
        require_empty_directory(
            paths["recost_forensics"], "recost forensic directory", managed=True
        )
        transaction_id = f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}"
        journal = paths["recost_transactions"] / f"{transaction_id}.json"
        forensic = (
            paths["recost_forensics"]
            / f"{transaction_id}.{paths['canonical'].name}.forensic"
        )
        record = {
            **recovery_bindings(paths, args, transaction_id, forensic),
            "state": "preparing",
            "created_utc": utc_now(),
        }
        write_json(journal, record)
        staged_identity = None
        try:
            copy_forensic(
                paths["staged"],
                forensic,
                args.expected_artifact_sha256,
                expected_mode=args.expected_artifact_mode,
            )
            record["state"] = "link-pending"
            write_json(journal, record)
            queue_file = args.pre_link_squeue_file or args.squeue_file
            with regular_descriptor(
                paths["staged"],
                "staged recost artifact",
                expected_mode=args.expected_artifact_mode,
                expected_links=1,
            ) as staged_descriptor:
                if sha256_descriptor(staged_descriptor) != args.expected_artifact_sha256:
                    raise ValueError("staged recost artifact checksum changed before link")
                staged_identity = descriptor_identity(staged_descriptor)
                directory_descriptor = os.open(
                    paths["accounting"],
                    os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW,
                )
                try:
                    require_trusted_directory_profile(
                        os.fstat(directory_descriptor),
                        paths["accounting"],
                        "accounting directory",
                    )
                    link_after_empty_queue(
                        root,
                        offline,
                        queue_file,
                        directory_descriptor,
                        paths,
                        args,
                        staged_identity,
                    )
                    require_linked_pair_in_directory(
                        directory_descriptor,
                        paths,
                        args,
                        expected_identity=staged_identity,
                    )
                    os.fsync(directory_descriptor)
                    if args.simulate_post_link_failure:
                        raise ValueError("simulated post-link publication failure")
                    os.unlink(paths["staged"].name, dir_fd=directory_descriptor)
                    os.fsync(directory_descriptor)
                finally:
                    os.close(directory_descriptor)
                require_named_artifact_identity(
                    paths["canonical"],
                    args,
                    expected_links=1,
                    expected_identity=staged_identity,
                )
            require_artifact_namespace_after_unlink(paths)
            generator, scheduler = authenticate_external_files(root, args)
            write_json(
                paths["audit"],
                audit_record(
                    paths,
                    args,
                    initial_source_path(),
                    generator,
                    scheduler,
                    forensic,
                    transaction_id,
                ),
            )
            validate_audit(
                paths,
                args,
                initial_source_path(),
                generator,
                scheduler,
            )
            unlink_durable(journal)
            require_empty_directory(
                paths["recost_transactions"], "recost transaction directory", managed=True
            )
        except BaseException as error:
            if (
                staged_identity is not None
                and link_failure_is_retry_safe(paths, args, staged_identity)
            ):
                remove_prelink_records(journal, forensic)
            elif record.get("state") == "link-pending":
                record["state"] = "ambiguous-after-link-attempt"
                record["ambiguity_recorded_utc"] = utc_now()
                write_json(journal, record)
            else:
                remove_prelink_records(journal, forensic)
            raise
    verify_promoted_state(paths, root_dir, root, offline, args)


def adopt_legacy_canonical(paths: dict[str, Path], root_dir: Path, root: Path,
                           offline: bool, args: argparse.Namespace) -> None:
    """Attest an authenticated canonical-only artifact without publication claims."""

    validate_fixture_options(args, offline)
    with promotion_lock(paths):
        require_empty_queue(root, offline, args.squeue_file)
        require_empty_directory(
            paths["stage_i_transactions"], "Stage I transaction directory", trusted=True
        )
        authenticate_committed_stage_i(root_dir, args.expected_stage_i_sha256)
        generator, scheduler = authenticate_external_files(root, args)
        read_artifact(paths["canonical"], args)
        run_reconcile(root_dir, root, offline, args)
        transactions_existed = entry_exists(paths["recost_transactions"])
        if not transactions_existed:
            require_artifact_namespace(paths, "canonical-pending-audit")
            if entry_exists(paths["recost_forensics_root"]):
                raise ValueError(
                    "legacy adoption requires an absent recost forensic root "
                    "before initial attestation"
                )
        else:
            require_managed_directory(
                paths["recost_transactions"], "recost transaction directory"
            )
            if not directory_entries(paths["recost_transactions"]):
                raise ValueError(
                    "legacy adoption requires an absent transaction directory "
                    "or one resumable journal"
                )
        mkdir_durable(paths["recost_transactions"])
        require_managed_directory(
            paths["recost_transactions"], "recost transaction directory"
        )
        if not paths["recost_forensics_root"].exists():
            mkdir_durable(paths["recost_forensics_root"])
        require_managed_directory(
            paths["recost_forensics_root"], "recost forensic root directory"
        )
        if not paths["recost_forensics"].exists():
            mkdir_durable(paths["recost_forensics"])
        require_managed_directory(
            paths["recost_forensics"], "recost forensic directory"
        )
        remove_json_temporaries(
            recost_transaction_temporary_entries(paths),
            "recost transaction temporary",
        )
        remove_json_temporaries(
            json_temporary_entries(paths["audit"]),
            "recost publication audit temporary",
        )
        if directory_entries(paths["recost_transactions"]):
            journal, record = recost_journal(paths)
            transaction_id, forensic, state = validate_adoption_journal(
                paths, args, record
            )
        else:
            if transactions_existed:
                raise ValueError("legacy adoption transaction journal disappeared")
            require_artifact_namespace(paths, "canonical-pending-audit")
            require_empty_directory(
                paths["recost_forensics"], "recost forensic directory", managed=True
            )
            transaction_id = f"{utc_now().replace(':', '')}-{uuid.uuid4().hex}"
            journal = paths["recost_transactions"] / f"{transaction_id}.json"
            forensic = (
                paths["recost_forensics"]
                / f"{transaction_id}.{paths['canonical'].name}.forensic"
            )
            record = {
                **adoption_bindings(paths, args, transaction_id, forensic),
                "state": "preparing",
                "created_utc": utc_now(),
            }
            write_json(journal, record)
            state = "preparing"
            if args.simulate_adoption_interruption_after_journal:
                raise ValueError("simulated interruption after legacy adoption journal")

        if state == "preparing":
            require_artifact_namespace(paths, "canonical-pending-audit")
            retained = set(directory_entries(paths["recost_forensics"]))
            if retained == set():
                recover_or_copy_forensic(
                    paths["canonical"],
                    forensic,
                    args.expected_artifact_sha256,
                    expected_mode=args.expected_artifact_mode,
                    simulate_interruption_before_directory_fsync=(
                        args.simulate_adoption_interruption_before_forensic_directory_fsync
                    ),
                )
            elif retained != {forensic}:
                raise ValueError(
                    "recost forensic directory entries differ; found: "
                    + ", ".join(str(item) for item in sorted(retained))
                )
            else:
                recover_or_copy_forensic(
                    paths["canonical"],
                    forensic,
                    args.expected_artifact_sha256,
                    expected_mode=args.expected_artifact_mode,
                    simulate_interruption_before_directory_fsync=(
                        args.simulate_adoption_interruption_before_forensic_directory_fsync
                    ),
                )
            record["state"] = "forensic-copied"
            write_json(journal, record)
            state = "forensic-copied"
            if args.simulate_adoption_interruption_after_forensic:
                raise ValueError("simulated interruption after legacy adoption forensic copy")

        if state == "forensic-copied":
            require_exact_directory_entries(
                paths["recost_forensics"],
                "recost forensic directory",
                {forensic},
                managed=True,
            )
            require_file_sha256(
                forensic,
                args.expected_artifact_sha256,
                "recost forensic copy",
                expected_mode=0o444,
                expected_links=1,
            )
            if not entry_exists(paths["audit"]):
                require_artifact_namespace(paths, "canonical-pending-audit")
                read_artifact(paths["canonical"], args)
                require_empty_queue(root, offline, args.squeue_file)
                write_json(
                    paths["audit"],
                    adoption_audit_record(
                        paths,
                        args,
                        initial_source_path(),
                        generator,
                        scheduler,
                        forensic,
                        transaction_id,
                    ),
                    simulate_interruption_before_directory_fsync=(
                        args.simulate_adoption_interruption_before_audit_directory_fsync
                    ),
                )
                if args.simulate_adoption_interruption_after_audit:
                    raise ValueError("simulated interruption after legacy adoption audit")
            fsync_directory(paths["accounting"])
            require_artifact_namespace(paths, "promoted")
            validate_audit(
                paths, args, initial_source_path(), generator, scheduler
            )
            record["state"] = "audit-written"
            write_json(journal, record)
            state = "audit-written"

        if state == "audit-written":
            require_artifact_namespace(paths, "promoted")
            validate_audit(
                paths, args, initial_source_path(), generator, scheduler
            )
        require_empty_queue(root, offline, args.squeue_file)
        unlink_durable(journal)
        require_empty_directory(
            paths["recost_transactions"], "recost transaction directory", managed=True
        )
    verify_promoted_state(paths, root_dir, root, offline, args)


def require_artifact_namespace_after_unlink(paths: dict[str, Path]) -> None:
    """Require only the canonical artifact before its audit is written."""

    retained = set(artifact_namespace_entries(paths))
    expected = {paths["canonical"]}
    if retained != expected:
        raise ValueError(
            "recost namespace is not canonical-only after publication; found: "
            + ", ".join(str(path) for path in sorted(retained))
        )


def add_artifact_arguments(command: argparse.ArgumentParser) -> None:
    """Add explicit generic recost binding arguments."""

    command.add_argument("--artifact-name", required=True)
    command.add_argument("--expected-artifact-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-artifact-mode", type=artifact_mode_arg, default=0o644)
    command.add_argument("--generator-relative-path", required=True)
    command.add_argument("--expected-generator-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-generator-mode", type=generator_mode_arg, default=0o755)
    command.add_argument("--scheduler-relative-path", required=True)
    command.add_argument("--expected-scheduler-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-scheduler-mode", type=scheduler_mode_arg, default=0o644)
    command.add_argument(
        "--authorized-next-segment-profile-json",
        type=json_object_arg,
        required=True,
    )
    for key in COUNT_KEYS:
        command.add_argument(
            f"--expected-{key.replace('_', '-')}",
            dest=f"expected_{key}",
            type=nonnegative_int_arg,
            required=True,
        )
    command.add_argument("--artifact-epoch-pointer", default="/execution_epoch")
    command.add_argument(
        "--artifact-profile-pointer",
        default="/authorization/sole_next_segment_profile",
    )
    command.add_argument("--artifact-counts-pointer", default="/reconcile/counts")
    command.add_argument(
        "--artifact-stage-i-sha256-pointer",
        default="/provenance/stage_i_helper_sha256",
    )
    command.add_argument(
        "--artifact-generator-sha256-pointer",
        default="/provenance/generator_sha256",
    )
    command.add_argument(
        "--artifact-scheduler-sha256-pointer",
        default="/provenance/scheduler_sha256",
    )


def parser() -> argparse.ArgumentParser:
    """Build the retained companion CLI."""

    command = argparse.ArgumentParser(prog=UTILITY_RELATIVE.name, description=__doc__)
    command.add_argument("--root", default=str(DEFAULT_ROOT))
    command.add_argument("--allow-local-root", action="store_true")
    command.add_argument("--expected-utility-sha256", type=sha256_arg, required=True)
    command.add_argument("--expected-stage-i-sha256", type=sha256_arg, required=True)
    command.add_argument(
        "--squeue-file",
        help="Offline local-root fixture input in place of querying squeue.",
    )
    actions = command.add_subparsers(dest="action", required=True)
    for name in (
        "audit-recost",
        "verify-staged-recost",
        "promote-recost",
        "verify-promoted-recost",
        "finalize-linked-pair",
        "retire-preparing",
        "adopt-legacy-canonical",
    ):
        action = actions.add_parser(name)
        add_artifact_arguments(action)
        if name == "audit-recost":
            action.add_argument(
                "--publication-state",
                choices=("staged", "promoted"),
                required=True,
            )
        if name == "promote-recost":
            action.add_argument(
                "--pre-link-squeue-file",
                help="Offline queue fixture read immediately before publication.",
            )
            action.add_argument(
                "--simulate-post-link-failure",
                action="store_true",
                help="Offline fixture hook retaining an ambiguous link attempt.",
            )
            action.add_argument(
                "--simulate-link-return-failure",
                action="store_true",
                help="Offline fixture hook interrupting immediately after os.link().",
            )
        if name == "finalize-linked-pair":
            action.add_argument(
                "--pre-unlink-squeue-file",
                help="Offline queue fixture read immediately before staged unlink.",
            )
        if name == "adopt-legacy-canonical":
            action.add_argument(
                "--simulate-adoption-interruption-after-journal",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption journal.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-before-forensic-directory-fsync",
                action="store_true",
                help="Offline fixture hook interrupting before forensic-directory fsync.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-after-forensic",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption forensic copy.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-before-audit-directory-fsync",
                action="store_true",
                help="Offline fixture hook interrupting before audit-directory fsync.",
            )
            action.add_argument(
                "--simulate-adoption-interruption-after-audit",
                action="store_true",
                help="Offline fixture hook interrupting after the adoption audit.",
            )
    return command


def main() -> int:
    """Run one authenticated retained-companion operation."""

    try:
        utility_sha = expected_self_sha256(sys.argv[1:])
        _, root_dir = authenticate_self(utility_sha)
        args = parser().parse_args()
        root, offline = require_root(Path(args.root), args.allow_local_root)
        require_trusted_directory(root, "Stage I root")
        if not offline:
            if root_dir != CANONICAL_REPOSITORY_ROOT:
                raise ValueError(
                    "canonical use requires the retained repository root "
                    f"{CANONICAL_REPOSITORY_ROOT}: {root_dir}"
                )
            authenticate_committed_utility(root_dir, args.expected_utility_sha256)
        paths = layout(root, args.artifact_name)
        if args.action == "audit-recost":
            if args.publication_state == "staged":
                verify_staged_state(paths, root_dir, root, offline, args)
                print(f"Audited staged-only recost artifact: {paths['staged']}")
            else:
                verify_promoted_state(paths, root_dir, root, offline, args)
                print(f"Audited canonical-only recost artifact: {paths['canonical']}")
            return 0
        if args.action == "verify-staged-recost":
            verify_staged_state(paths, root_dir, root, offline, args)
            print(f"Verified staged-only recost artifact: {paths['staged']}")
            return 0
        if args.action == "promote-recost":
            promote(paths, root_dir, root, offline, args)
            print(f"Promoted and verified recost artifact: {paths['canonical']}")
            return 0
        if args.action == "verify-promoted-recost":
            verify_promoted_state(paths, root_dir, root, offline, args)
            print(f"Verified canonical-only recost artifact: {paths['canonical']}")
            return 0
        if args.action == "finalize-linked-pair":
            finalize_linked_pair(paths, root_dir, root, offline, args)
            print(f"Finalized and verified recost artifact: {paths['canonical']}")
            return 0
        if args.action == "retire-preparing":
            retire_preparing(paths, root_dir, root, offline, args)
            print(f"Retired interrupted preparing transaction: {paths['staged']}")
            return 0
        if args.action == "adopt-legacy-canonical":
            adopt_legacy_canonical(paths, root_dir, root, offline, args)
            print(f"Adopted and verified legacy canonical artifact: {paths['canonical']}")
            return 0
        raise ValueError(f"unsupported action: {args.action}")
    except (KeyError, OSError, TypeError, ValueError,
            subprocess.CalledProcessError) as error:
        print(f"Stage I recost checkpoint utility failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
