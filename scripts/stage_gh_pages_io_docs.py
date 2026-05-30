#!/usr/bin/env python3
"""Stage deferred IO documentation into a detached gh-pages worktree.

The helper intentionally never fetches, commits, pushes, or switches branches.
Strict staging accepts only a clean detached worktree at the frozen baseline.
"""

from __future__ import annotations

import argparse
import difflib
import hashlib
import json
import os
import re
import stat
import subprocess
import sys
import tempfile
from pathlib import Path, PurePosixPath
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
BUNDLE_DIR = (
    REPO_ROOT / "deferred_docs" / "gh-pages" / "io-output-formats-and-sharding"
)
MANIFEST_PATH = BUNDLE_DIR / "manifest.json"
MARKER_PREFIX = "RCP-07 io-output-formats-and-sharding"
GIT_TIMEOUT_SECONDS = 30
DOCS_BUILD_TIMEOUT_SECONDS = 300


class StageError(RuntimeError):
    """Raised when staging cannot proceed without an explicit reconciliation."""


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _git(
    repo: Path, *args: str, check: bool = True, text: bool = True
) -> subprocess.CompletedProcess[Any]:
    proc = subprocess.run(
        ["git", "-C", str(repo), *args],
        check=False,
        capture_output=True,
        text=text,
        timeout=GIT_TIMEOUT_SECONDS,
    )
    if check and proc.returncode != 0:
        stderr = proc.stderr.strip() if text else proc.stderr.decode(errors="replace")
        raise StageError(f"git {' '.join(args)} failed in {repo}: {stderr}")
    return proc


def _safe_relative_path(value: str, *, label: str) -> PurePosixPath:
    path = PurePosixPath(value)
    if path.is_absolute() or not path.parts or ".." in path.parts:
        raise StageError(f"{label} must be a safe relative path: {value!r}")
    return path


def _resolve_git_path(repo: Path, value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo / path
    return path.resolve()


def _registered_worktrees(repo: Path) -> set[Path]:
    worktrees = _git(repo, "worktree", "list", "--porcelain").stdout
    return {
        Path(line.removeprefix("worktree ")).resolve()
        for line in worktrees.splitlines()
        if line.startswith("worktree ")
    }


def _is_within(path: Path, directory: Path) -> bool:
    try:
        path.relative_to(directory)
    except ValueError:
        return False
    return True


def _git_blob(repo: Path, revision: str, path: str) -> str | None:
    proc = _git(repo, "rev-parse", f"{revision}:{path}", check=False)
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def _git_show(repo: Path, revision: str, path: str) -> bytes | None:
    proc = _git(repo, "show", f"{revision}:{path}", check=False, text=False)
    if proc.returncode != 0:
        return None
    return proc.stdout


def _hash_worktree_file(repo: Path, path: str) -> str | None:
    target = repo / path
    if not target.is_file():
        return None
    return _git(repo, "hash-object", "--", path).stdout.strip()


def _status_entries(repo: Path) -> list[tuple[str, str]]:
    proc = _git(
        repo,
        "status",
        "--porcelain=v1",
        "-z",
        "--untracked-files=all",
        text=False,
    )
    records = proc.stdout.split(b"\0")
    entries: list[tuple[str, str]] = []
    index = 0
    while index < len(records):
        record = records[index]
        index += 1
        if not record:
            continue
        if len(record) < 4 or record[2:3] != b" ":
            raise StageError(f"unexpected git status record: {record!r}")
        state = record[:2].decode("ascii", errors="replace")
        path = record[3:].decode(errors="surrogateescape")
        if "R" in state or "C" in state:
            if index < len(records):
                index += 1
            raise StageError(f"renamed or copied target path is not supported: {path}")
        entries.append((state, path))
    return entries


def _status_paths(entries: list[tuple[str, str]]) -> set[str]:
    return {path for _, path in entries}


def _index_entries(repo: Path) -> list[tuple[str, str, str, str]]:
    proc = _git(repo, "ls-files", "--stage", "-z", text=False)
    entries: list[tuple[str, str, str, str]] = []
    for record in proc.stdout.split(b"\0"):
        if not record:
            continue
        try:
            metadata, raw_path = record.split(b"\t", 1)
            mode, object_id, stage = metadata.decode("ascii").split()
        except ValueError as exc:
            raise StageError(f"unexpected git index record: {record!r}") from exc
        entries.append(
            (mode, object_id, stage, raw_path.decode(errors="surrogateescape"))
        )
    return entries


def _verify_canonical_index(repo: Path) -> None:
    index_entries = _index_entries(repo)
    for mode, object_id, stage, path in index_entries:
        if stage != "0":
            raise StageError(f"target index has non-stage-zero entry: {path}")
        if set(object_id) == {"0"}:
            raise StageError(f"target index has intent-to-add entry: {path}")
        if mode not in {"100644", "100755", "120000", "160000"}:
            raise StageError(f"target index has unsupported mode {mode}: {path}")
    head_entries = _git(repo, "ls-tree", "-r", "-z", "--full-tree", "HEAD", text=False)
    expected: list[tuple[str, str, str]] = []
    for record in head_entries.stdout.split(b"\0"):
        if not record:
            continue
        try:
            metadata, raw_path = record.split(b"\t", 1)
            mode, _, object_id = metadata.decode("ascii").split()
        except ValueError as exc:
            raise StageError(f"unexpected HEAD tree record: {record!r}") from exc
        expected.append((mode, object_id, raw_path.decode(errors="surrogateescape")))
    actual = [(mode, object_id, path) for mode, object_id, _, path in index_entries]
    if actual != expected:
        raise StageError(
            "target index differs from HEAD, including possible intent-to-add entries"
        )
    proc = _git(repo, "ls-files", "-v", "-z", text=False)
    for record in proc.stdout.split(b"\0"):
        if not record:
            continue
        tag = chr(record[0])
        path = record[2:].decode(errors="surrogateescape")
        if tag != "H":
            raise StageError(
                f"target index has noncanonical flag tag {tag!r}: {path}"
            )


def _git_tree_mode(repo: Path, revision: str, path: str) -> int | None:
    proc = _git(repo, "ls-tree", revision, "--", path, check=False)
    if proc.returncode != 0 or not proc.stdout:
        return None
    mode = proc.stdout.split(maxsplit=1)[0]
    if mode not in {"100644", "100755"}:
        raise StageError(f"{revision} target is not a regular file: {path}")
    return int(mode[-3:], 8)


def _verify_regular_file_mode(path: Path, expected_mode: int, *, label: str) -> None:
    try:
        info = path.lstat()
    except FileNotFoundError as exc:
        raise StageError(f"{label} is missing: {path}") from exc
    if not stat.S_ISREG(info.st_mode):
        raise StageError(f"{label} must be a regular file: {path}")
    actual_mode = stat.S_IMODE(info.st_mode)
    if actual_mode != expected_mode:
        raise StageError(
            f"{label} mode drift for {path}: expected {expected_mode:o}, "
            f"found {actual_mode:o}"
        )


def _verify_target_file_modes(
    target: Path, manifest: dict[str, Any], *, require_additions: bool
) -> None:
    for item in manifest["targets"]:
        path = item["path"]
        target_path = target / path
        if item.get("expected_absent") and not target_path.exists():
            if require_additions:
                raise StageError(f"expected staged add target is missing: {path}")
            continue
        expected_mode = _git_tree_mode(target, "HEAD", path)
        if expected_mode is None:
            expected_mode = 0o644
        _verify_regular_file_mode(
            target_path, expected_mode, label="target worktree file"
        )


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        manifest = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise StageError(f"cannot load manifest {path}: {exc}") from exc
    if manifest.get("schema_version") != 1:
        raise StageError("unsupported manifest schema")
    return manifest


def _validate_manifest(manifest: dict[str, Any], bundle_dir: Path) -> None:
    allowlist = manifest.get("allowlist")
    targets = manifest.get("targets")
    if not isinstance(allowlist, list) or not isinstance(targets, list):
        raise StageError("manifest allowlist and targets must be lists")
    if len(allowlist) != len(set(allowlist)):
        raise StageError("manifest allowlist contains duplicates")
    target_paths = [item.get("path") for item in targets]
    if set(allowlist) != set(target_paths) or len(target_paths) != len(set(target_paths)):
        raise StageError("manifest target paths must equal the exact allowlist")
    for path in allowlist:
        _safe_relative_path(path, label="target path")
    for path in manifest.get("protected_blobs", {}):
        _safe_relative_path(path, label="protected path")
        if path in allowlist:
            raise StageError(f"protected path must not be staged: {path}")
    seen_markers: set[str] = set()
    for target in targets:
        if bool(target.get("expected_absent")) == bool(target.get("baseline_blob")):
            raise StageError(
                f"{target['path']} needs exactly one of baseline_blob or expected_absent"
            )
        operations = target.get("operations")
        if not isinstance(operations, list) or not operations:
            raise StageError(f"{target['path']} has no operations")
        for operation in operations:
            kind = operation.get("kind")
            if kind not in {
                "replace_section",
                "insert_before",
                "replace_file",
                "add_file",
            }:
                raise StageError(f"unsupported operation {kind!r} for {target['path']}")
            source = _safe_relative_path(operation.get("source", ""), label="source")
            source_path = (bundle_dir / source).resolve()
            try:
                source_path.relative_to(bundle_dir.resolve())
            except ValueError as exc:
                raise StageError(f"source escapes bundle: {source}") from exc
            if kind in {"replace_section", "insert_before"}:
                marker = operation.get("marker")
                if not marker or marker in seen_markers:
                    raise StageError(f"missing or duplicate marker ID: {marker!r}")
                seen_markers.add(marker)


def _load_sources(manifest: dict[str, Any], bundle_dir: Path) -> dict[str, bytes]:
    sources: dict[str, bytes] = {}
    for target in manifest["targets"]:
        for operation in target["operations"]:
            source = operation["source"]
            if source in sources:
                data = sources[source]
            else:
                path = bundle_dir / source
                try:
                    data = path.read_bytes()
                except OSError as exc:
                    raise StageError(
                        f"cannot read source fragment {path}: {exc}"
                    ) from exc
                sources[source] = data
            actual = _sha256(data)
            if actual != operation["source_sha256"]:
                raise StageError(
                    f"source SHA256 drift for {source}: expected "
                    f"{operation['source_sha256']}, found {actual}"
                )
            if b"Merge this reviewed section" in data:
                raise StageError(
                    f"internal merge prose leaked into public payload: {source}"
                )
    return sources


def _validate_target(
    source_repo: Path,
    target_arg: Path,
    manifest: dict[str, Any],
    *,
    strict_baseline: bool,
) -> Path:
    target_input = target_arg.expanduser()
    target = target_input.resolve()
    if not target_input.is_absolute():
        raise StageError("target must be an absolute worktree path")
    if not target.is_dir():
        raise StageError(f"target worktree does not exist: {target}")
    top = Path(_git(target, "rev-parse", "--show-toplevel").stdout.strip()).resolve()
    if top != target:
        raise StageError(f"target must be the canonical worktree root: {target}")
    source_common = _resolve_git_path(
        source_repo, _git(source_repo, "rev-parse", "--git-common-dir").stdout.strip()
    )
    target_common = _resolve_git_path(
        target, _git(target, "rev-parse", "--git-common-dir").stdout.strip()
    )
    if source_common != target_common:
        raise StageError("target is not a linked worktree of the source repository")
    symbolic = _git(target, "symbolic-ref", "--quiet", "HEAD", check=False)
    if symbolic.returncode == 0:
        raise StageError(
            f"target must be detached, not checked out at {symbolic.stdout.strip()}"
        )
    listed = _registered_worktrees(source_repo)
    if target not in listed:
        raise StageError(f"target is not registered as a linked worktree: {target}")
    if strict_baseline:
        expected = manifest["expected_baseline"]
        head = _git(target, "rev-parse", "HEAD").stdout.strip()
        origin = _git(
            source_repo, "rev-parse", "refs/remotes/origin/gh-pages"
        ).stdout.strip()
        if head != expected:
            raise StageError(f"target HEAD drift: expected {expected}, found {head}")
        if origin != expected:
            raise StageError(
                f"local origin/gh-pages drift: expected {expected}, found {origin}"
            )
    return target


def _verify_inventory(
    target: Path,
    manifest: dict[str, Any],
    *,
    require_clean_payloads: bool,
) -> None:
    for path, expected in manifest.get("protected_blobs", {}).items():
        head_blob = _git_blob(target, "HEAD", path)
        worktree_blob = _hash_worktree_file(target, path)
        if head_blob != expected:
            raise StageError(
                f"protected HEAD blob drift for {path}: expected {expected}, "
                f"found {head_blob}"
            )
        if worktree_blob != expected:
            raise StageError(
                f"protected worktree blob drift for {path}: expected {expected}, "
                f"found {worktree_blob}"
            )
    for item in manifest["targets"]:
        path = item["path"]
        head_blob = _git_blob(target, "HEAD", path)
        if item.get("expected_absent"):
            if head_blob is not None:
                raise StageError(f"expected add target already exists in HEAD: {path}")
            if require_clean_payloads and (target / path).exists():
                raise StageError(
                    f"expected add target already exists in worktree: {path}"
                )
            continue
        expected = item["baseline_blob"]
        if head_blob != expected:
            raise StageError(
                f"target HEAD blob drift for {path}: expected {expected}, "
                f"found {head_blob}"
            )
        if require_clean_payloads:
            worktree_blob = _hash_worktree_file(target, path)
            if worktree_blob != expected:
                raise StageError(
                    f"target worktree blob drift for {path}: expected {expected}, "
                    f"found {worktree_blob}"
                )


def _unique_index(text: str, anchor: str, *, label: str) -> int:
    count = text.count(anchor)
    if count != 1:
        raise StageError(f"{label} must occur exactly once; found {count}")
    return text.index(anchor)


def _marker_block(marker: str, payload: str) -> str:
    begin = f"<!-- BEGIN {MARKER_PREFIX}: {marker} -->"
    end = f"<!-- END {MARKER_PREFIX}: {marker} -->"
    return f"{begin}\n{payload.rstrip()}\n{end}"


def _join_marker_suffix(block: str, suffix: str) -> str:
    """Keep retained Markdown from sharing the closing-marker line."""
    if suffix and not suffix.startswith("\n"):
        return block + "\n" + suffix
    return block + suffix


def _transform_bounded(
    text: str,
    operation: dict[str, Any],
    payload: str,
    *,
    allow_section_drift: bool,
) -> str:
    marker = operation["marker"]
    begin = f"<!-- BEGIN {MARKER_PREFIX}: {marker} -->"
    end = f"<!-- END {MARKER_PREFIX}: {marker} -->"
    begin_count = text.count(begin)
    end_count = text.count(end)
    if begin_count != end_count or begin_count > 1:
        raise StageError(
            f"marker {marker!r} must be absent or appear once as a matched pair; "
            f"found begin={begin_count}, end={end_count}"
        )
    block = _marker_block(marker, payload)
    if begin_count == 1:
        start = text.index(begin)
        finish = text.index(end, start) + len(end)
        return text[:start] + _join_marker_suffix(block, text[finish:])
    kind = operation["kind"]
    if kind == "replace_section":
        start = _unique_index(
            text, operation["start_anchor"], label=f"{marker} start anchor"
        )
        finish = _unique_index(
            text, operation["end_anchor"], label=f"{marker} end anchor"
        )
        if finish <= start:
            raise StageError(f"{marker} end anchor precedes its start anchor")
        section = text[start:finish].encode()
        actual = _sha256(section)
        expected = operation["baseline_section_sha256"]
        if actual != expected and not allow_section_drift:
            raise StageError(
                f"baseline section SHA256 drift for {marker}: expected {expected}, "
                f"found {actual}"
            )
        return text[:start] + _join_marker_suffix(block, text[finish:])
    if kind == "insert_before":
        anchor = operation["anchor"]
        index = _unique_index(text, anchor, label=f"{marker} insertion anchor")
        actual = _sha256(anchor.encode())
        if actual != operation["anchor_sha256"]:
            raise StageError(
                f"anchor SHA256 drift for {marker}: expected "
                f"{operation['anchor_sha256']}, "
                f"found {actual}"
            )
        return text[:index] + _join_marker_suffix(block, text[index:])
    raise StageError(f"bounded transform received unsupported kind: {kind}")


def _transform_target(
    item: dict[str, Any],
    current: bytes | None,
    sources: dict[str, bytes],
    *,
    allow_section_drift: bool,
) -> bytes:
    result = current
    for operation in item["operations"]:
        kind = operation["kind"]
        payload = sources[operation["source"]]
        if kind in {"replace_file", "add_file"}:
            if len(item["operations"]) != 1:
                raise StageError(
                    f"{kind} must be the sole operation for {item['path']}"
                )
            result = payload
            continue
        if result is None:
            raise StageError(f"bounded-edit target is missing: {item['path']}")
        try:
            text = result.decode()
            fragment = payload.decode()
        except UnicodeDecodeError as exc:
            raise StageError(
                f"documentation payload is not UTF-8: {item['path']}"
            ) from exc
        result = _transform_bounded(
            text,
            operation,
            fragment,
            allow_section_drift=allow_section_drift,
        ).encode()
    if result is None:
        raise StageError(f"no transformed content generated for {item['path']}")
    return result


def _plan_transformations(
    target: Path,
    manifest: dict[str, Any],
    sources: dict[str, bytes],
    *,
    allow_section_drift: bool,
    base_revision: str | None = None,
) -> dict[str, bytes]:
    plan: dict[str, bytes] = {}
    for item in manifest["targets"]:
        path = item["path"]
        if base_revision is None:
            worktree_path = target / path
            current = worktree_path.read_bytes() if worktree_path.is_file() else None
        else:
            current = _git_show(target, base_revision, path)
        once = _transform_target(
            item, current, sources, allow_section_drift=allow_section_drift
        )
        twice = _transform_target(
            item, once, sources, allow_section_drift=allow_section_drift
        )
        if once != twice:
            raise StageError(f"pure second transformation is not idempotent: {path}")
        plan[path] = once
    if set(plan) != set(manifest["allowlist"]):
        raise StageError("planned file set differs from manifest allowlist")
    return plan


def _virtual_markdown_documents(target: Path, plan: dict[str, bytes]) -> dict[str, str]:
    documents: dict[str, str] = {}
    source_root = target / "docs" / "source"
    if source_root.is_dir():
        for path in source_root.rglob("*.md"):
            relative = path.relative_to(target).as_posix()
            documents[relative] = path.read_text()
    for path, data in plan.items():
        if path.endswith(".md"):
            documents[path] = data.decode()
    return documents


def _contradiction_report(
    target: Path, plan: dict[str, bytes], manifest: dict[str, Any]
) -> list[dict[str, Any]]:
    documents = _virtual_markdown_documents(target, plan)
    report: list[dict[str, Any]] = []
    failures: list[str] = []
    for check in manifest.get("contradiction_checks", []):
        regex = re.compile(check["pattern"])
        matches: list[str] = []
        for path, text in sorted(documents.items()):
            for match in regex.finditer(text):
                line = text.count("\n", 0, match.start()) + 1
                matches.append(f"{path}:{line}:{match.group(0)}")
        row = {
            "pattern": check["pattern"],
            "description": check["description"],
            "count": len(matches),
            "matches": matches,
        }
        report.append(row)
        minimum = check.get("minimum", 0)
        maximum = check.get("maximum")
        if len(matches) < minimum:
            failures.append(
                f"{check['description']}: expected at least {minimum} matches, "
                f"found {len(matches)}"
            )
        if maximum is not None and len(matches) > maximum:
            failures.append(
                f"{check['description']}: expected at most {maximum} matches, "
                f"found {len(matches)}"
            )
    print("Contradiction search report:")
    for row in report:
        print(f"  {row['description']}: {row['count']} match(es)")
        for match in row["matches"][:10]:
            print(f"    {match}")
        if len(row["matches"]) > 10:
            print(f"    ... {len(row['matches']) - 10} more")
    if failures:
        raise StageError("contradiction search failed: " + "; ".join(failures))
    return report


def _atomic_write(path: Path, data: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    mode = path.stat().st_mode & 0o777 if path.exists() else 0o644
    with tempfile.NamedTemporaryFile(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", delete=False
    ) as handle:
        temporary = Path(handle.name)
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    try:
        os.chmod(temporary, mode)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _verify_exact_status(target: Path, manifest: dict[str, Any]) -> None:
    entries = _status_entries(target)
    paths = _status_paths(entries)
    expected = set(manifest["allowlist"])
    if paths != expected:
        raise StageError(
            "staged status paths differ from exact allowlist: "
            f"expected={sorted(expected)}, found={sorted(paths)}"
        )


def _verify_clean_index(target: Path) -> None:
    _verify_canonical_index(target)
    proc = _git(target, "diff", "--cached", "--quiet", "HEAD", "--", check=False)
    if proc.returncode == 1:
        raise StageError("target index differs from HEAD")
    if proc.returncode != 0:
        raise StageError(f"cannot compare target index to HEAD: {proc.stderr.strip()}")
    index_entries = [
        (state, path)
        for state, path in _status_entries(target)
        if state[0] not in {" ", "?"}
    ]
    if index_entries:
        raise StageError(f"target index has non-HEAD entries: {index_entries}")


def _print_review_diff(target: Path, manifest: dict[str, Any]) -> None:
    print("Review diff stat:")
    for path in manifest["allowlist"]:
        before = _git_show(target, "HEAD", path) or b""
        after_path = target / path
        after = after_path.read_bytes() if after_path.is_file() else b""
        before_lines = before.decode(errors="replace").splitlines()
        after_lines = after.decode(errors="replace").splitlines()
        additions = sum(
            1
            for line in difflib.ndiff(before_lines, after_lines)
            if line.startswith("+ ")
        )
        deletions = sum(
            1
            for line in difflib.ndiff(before_lines, after_lines)
            if line.startswith("- ")
        )
        print(f"  {path} | +{additions} -{deletions}")
    print("Review diff:")
    for path in manifest["allowlist"]:
        before = (_git_show(target, "HEAD", path) or b"").decode(errors="replace")
        after = (target / path).read_text()
        diff = difflib.unified_diff(
            before.splitlines(),
            after.splitlines(),
            fromfile=f"a/{path}",
            tofile=f"b/{path}",
            lineterm="",
        )
        for line in diff:
            print(line)


def _print_build_commands(target: Path) -> None:
    docs = target / "docs"
    print("Strict documentation commands:")
    print(f"  cd {docs}")
    print("  make clean html SPHINXOPTS='-W --keep-going'")
    print("  make linkcheck SPHINXOPTS='-W --keep-going'")


def _run_builds(target: Path) -> None:
    docs = target / "docs"
    for command in (
        ["make", "clean", "html", "SPHINXOPTS=-W --keep-going"],
        ["make", "linkcheck", "SPHINXOPTS=-W --keep-going"],
    ):
        print(f"Running: {' '.join(command)}")
        subprocess.run(
            command, cwd=docs, check=True, timeout=DOCS_BUILD_TIMEOUT_SECONDS
        )


def _write_packet_file(packet: Path, relative: str, data: bytes | None) -> None:
    path = packet / relative
    if data is None:
        path = path.with_name(path.name + ".ABSENT")
        data = b"ABSENT\n"
    _atomic_write(path, data)


def _fingerprint_entry(path: Path) -> dict[str, Any]:
    info = path.lstat()
    entry: dict[str, Any] = {
        "mode": info.st_mode,
        "size": info.st_size,
        "mtime_ns": info.st_mtime_ns,
        "ctime_ns": info.st_ctime_ns,
    }
    if stat.S_ISREG(info.st_mode):
        entry["sha256"] = _sha256(path.read_bytes())
    elif stat.S_ISLNK(info.st_mode):
        entry["symlink"] = os.readlink(path)
    return entry


def _target_fingerprint(target: Path) -> dict[str, Any]:
    # Git status hides ignored files and collapses repeated edits to a dirty file.
    # Snapshot the whole worktree namespace so packet generation is provably write-free.
    tree: dict[str, dict[str, Any]] = {}
    for root, directories, files in os.walk(target, followlinks=False):
        root_path = Path(root)
        for name in sorted([*directories, *files]):
            path = root_path / name
            tree[path.relative_to(target).as_posix()] = _fingerprint_entry(path)
        directories[:] = [
            name for name in directories if not (root_path / name).is_symlink()
        ]
    return {
        "root": _fingerprint_entry(target),
        "status": _status_entries(target),
        "tree": tree,
    }


def _write_reviewed_drift_packet(
    source_repo: Path,
    target: Path,
    packet_arg: Path,
    manifest: dict[str, Any],
    sources: dict[str, bytes],
) -> None:
    packet_input = packet_arg.expanduser()
    packet = packet_input.resolve()
    if not packet_input.is_absolute():
        raise StageError("reviewed-drift packet must use an absolute path")
    excluded = _registered_worktrees(source_repo)
    excluded.add(
        _resolve_git_path(
            source_repo,
            _git(source_repo, "rev-parse", "--git-common-dir").stdout.strip(),
        )
    )
    if any(_is_within(packet, directory) for directory in excluded):
        raise StageError(
            "reviewed-drift packet must be outside every registered worktree "
            "and repository administrative directory"
        )
    if packet.exists():
        raise StageError(f"reviewed-drift packet path already exists: {packet}")
    before = _target_fingerprint(target)
    packet.mkdir(parents=True)
    expected = manifest["expected_baseline"]
    report: dict[str, Any] = {
        "mode": "reviewed-drift",
        "expected_baseline": expected,
        "target_head": _git(target, "rev-parse", "HEAD").stdout.strip(),
        "local_origin_gh_pages": _git(
            source_repo, "rev-parse", "refs/remotes/origin/gh-pages", check=False
        ).stdout.strip(),
        "target_status": before["status"],
        "files": [],
    }
    proposal: dict[str, bytes] = {}
    for item in manifest["targets"]:
        path = item["path"]
        current_path = target / path
        current = current_path.read_bytes() if current_path.is_file() else None
        base = _git_show(source_repo, expected, path)
        error = None
        proposed = None
        try:
            proposed = _transform_target(
                item, current, sources, allow_section_drift=True
            )
            proposal[path] = proposed
        except StageError as exc:
            error = str(exc)
        _write_packet_file(packet, f"base/{path}", base)
        _write_packet_file(packet, f"current/{path}", current)
        _write_packet_file(packet, f"proposed/{path}", proposed)
        if proposed is not None:
            before_text = (current or b"").decode(errors="replace").splitlines()
            after_text = proposed.decode(errors="replace").splitlines()
            diff = "\n".join(
                difflib.unified_diff(
                    before_text,
                    after_text,
                    fromfile=f"current/{path}",
                    tofile=f"proposed/{path}",
                    lineterm="",
                )
            )
            _write_packet_file(
                packet, f"diffs/{path}.diff", (diff + "\n").encode()
            )
        report["files"].append({"path": path, "error": error})
    if len(proposal) == len(manifest["allowlist"]):
        try:
            report["contradiction_report"] = _contradiction_report(
                target, proposal, manifest
            )
        except StageError as exc:
            report["contradiction_error"] = str(exc)
    for source, data in sources.items():
        _write_packet_file(packet, f"sources/{source}", data)
    _atomic_write(
        packet / "report.json", (json.dumps(report, indent=2) + "\n").encode()
    )
    _atomic_write(
        packet / "README.md",
        (
            "# RCP-07 Reviewed Drift Reconciliation Packet\n\n"
            "This packet is write-free with respect to the target worktree. Review "
            "`base/`, `current/`, `proposed/`, `diffs/`, and `report.json` before "
            "approving any restaged Pages payload.\n"
        ).encode(),
    )
    after = _target_fingerprint(target)
    if before != after:
        raise StageError("reviewed-drift packet generation modified the target worktree")
    print(f"Reviewed-drift packet written outside target: {packet}")


def execute(
    target_arg: Path,
    *,
    source_repo: Path = REPO_ROOT,
    manifest_path: Path = MANIFEST_PATH,
    verify_staged: bool = False,
    reviewed_drift: Path | None = None,
    run_builds: bool = False,
) -> None:
    source_repo = source_repo.resolve()
    manifest_path = manifest_path.resolve()
    bundle_dir = manifest_path.parent
    manifest = _load_manifest(manifest_path)
    _validate_manifest(manifest, bundle_dir)
    sources = _load_sources(manifest, bundle_dir)
    target = _validate_target(
        source_repo,
        target_arg,
        manifest,
        strict_baseline=reviewed_drift is None,
    )
    _verify_canonical_index(target)
    if reviewed_drift is not None:
        _write_reviewed_drift_packet(
            source_repo, target, reviewed_drift, manifest, sources
        )
        _print_build_commands(target)
        return
    _verify_inventory(
        target, manifest, require_clean_payloads=not verify_staged
    )
    if verify_staged:
        _verify_clean_index(target)
        _verify_exact_status(target, manifest)
        _verify_target_file_modes(target, manifest, require_additions=True)
        plan = _plan_transformations(
            target,
            manifest,
            sources,
            allow_section_drift=False,
            base_revision="HEAD",
        )
        for path, expected in plan.items():
            actual = (target / path).read_bytes()
            if actual != expected:
                raise StageError(f"staged payload differs from frozen HEAD plan: {path}")
        _contradiction_report(target, plan, manifest)
        print(
            "Verified staged payload: exact allowlist and idempotence checks passed."
        )
    else:
        dirty = _status_entries(target)
        if dirty:
            raise StageError(f"strict target worktree must be clean: {dirty}")
        _verify_target_file_modes(target, manifest, require_additions=False)
        plan = _plan_transformations(
            target, manifest, sources, allow_section_drift=False
        )
        for path, data in plan.items():
            current_path = target / path
            current = current_path.read_bytes() if current_path.is_file() else None
            if current == data:
                raise StageError(f"strict transformation made no change: {path}")
        _contradiction_report(target, plan, manifest)
        for path, data in plan.items():
            _atomic_write(target / path, data)
        _verify_exact_status(target, manifest)
        for path, expected in plan.items():
            if (target / path).read_bytes() != expected:
                raise StageError(f"atomic write verification failed: {path}")
        print(
            "Staged payload: strict baseline, allowlist, and atomic-write checks "
            "passed."
        )
    _print_review_diff(target, manifest)
    _print_build_commands(target)
    if run_builds:
        _run_builds(target)


def _parse_args(argv: list[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", type=Path, help="absolute detached gh-pages worktree")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--verify-staged",
        action="store_true",
        help="allow intended allowlist dirtiness and verify the staged payload",
    )
    mode.add_argument(
        "--reviewed-drift",
        type=Path,
        metavar="PACKET_DIR",
        help="write a three-way reconciliation packet outside the target; never stage",
    )
    parser.add_argument(
        "--run-builds",
        action="store_true",
        help="run strict Sphinx HTML and link-check commands after verification",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)
    try:
        execute(
            args.target,
            verify_staged=args.verify_staged,
            reviewed_drift=args.reviewed_drift,
            run_builds=args.run_builds,
        )
    except (OSError, StageError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
