"""Diff-driven disposition gate for direct MPI calls touched by the IO branch."""

from collections import Counter
from pathlib import Path
import re
import subprocess


ROOT = Path(__file__).resolve().parents[3]
DISPOSITIONS = ROOT / "tst" / "test_suite" / "io" / "rcp02_mpi_dispositions.tsv"
MPI_CALL = re.compile(r"\b(MPI_[A-Za-z0-9_]+)\s*\(")
ALLOWED_DISPOSITIONS = {
    "checked",
    "intentional-fatal-shutdown",
    "wrapper-owned-status",
    "inherited-startup-lifecycle",
    "inherited-outside-rcp02",
}


def _subprocess_run(*args, **kwargs):
    kwargs.setdefault("timeout", 30)
    return subprocess.run(*args, **kwargs)


def _diff_touched_source_files():
    tracked_output = _subprocess_run(
        ["git", "-C", str(ROOT), "diff", "--name-only", "origin/main", "--", "src"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    untracked_output = _subprocess_run(
        [
            "git",
            "-C",
            str(ROOT),
            "ls-files",
            "--others",
            "--exclude-standard",
            "--",
            "src",
        ],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    paths = set(tracked_output.splitlines()) | set(untracked_output.splitlines())
    return [ROOT / line for line in sorted(paths) if line]


def _actual_inventory():
    inventory = Counter()
    for path in _diff_touched_source_files():
        if not path.is_file():
            continue
        relative_path = str(path.relative_to(ROOT))
        for symbol in MPI_CALL.findall(path.read_text()):
            inventory[(relative_path, symbol)] += 1
    return inventory


def _expected_inventory():
    inventory = Counter()
    for number, line in enumerate(DISPOSITIONS.read_text().splitlines(), start=1):
        if not line or line.startswith("#"):
            continue
        path, symbol, count, disposition, note = line.split("\t")
        assert disposition in ALLOWED_DISPOSITIONS, (
            f"{DISPOSITIONS}:{number} has unknown disposition {disposition!r}"
        )
        assert note, f"{DISPOSITIONS}:{number} requires an audit note"
        key = (path, symbol)
        assert key not in inventory, f"{DISPOSITIONS}:{number} duplicates {key}"
        inventory[key] = int(count)
    return inventory


def test_diff_touched_mpi_calls_have_explicit_dispositions():
    assert _actual_inventory() == _expected_inventory()
