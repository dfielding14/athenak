"""Local coverage for immutable external IO evidence-packet finalization."""

from pathlib import Path
import importlib.util
import os
import signal
import stat
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[3]
FINALIZER = ROOT / "scripts" / "finalize_external_io_qualification_packet.sh"
VALIDATOR = ROOT / "scripts" / "validate_external_io_packet_index.py"


def _write_executable(path: Path, payload: str):
    path.write_text(payload)
    path.chmod(0o755)


def _packet(tmp_path: Path):
    packet = tmp_path / "packet"
    (packet / "logs").mkdir(parents=True)
    (packet / "packet-index.tsv").write_text(
        "row\taction\tattempt_id\tsample\tmeasured\texit_code\tdisposition\t"
        "monotonic_elapsed_ns\toutput_dir\tfs_accounting_hook\n"
        "ED-1\tgenerate\tattempt-1\trank-map\t0\t0\tpassed-rank-map\t1\t"
        "N/A\tN/A\n"
    )
    (packet / "logs" / "retained file.txt").write_text("retained\n")
    (packet / "logs" / "artifacts.sha256").write_text("nested retained artifact\n")
    return packet


def _run_finalizer(packet: Path, archive_record: Path, check=True, env=None):
    return subprocess.run(
        ["/bin/bash", str(FINALIZER), str(packet), str(archive_record)],
        check=check,
        capture_output=True,
        text=True,
        env=env,
    )


def _run_validator(*arguments: str, check=True):
    return subprocess.run(
        [str(VALIDATOR), *arguments],
        check=check,
        capture_output=True,
        text=True,
    )


def _validator_module():
    spec = importlib.util.spec_from_file_location("external_io_packet_index", VALIDATOR)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _identity(path: Path):
    metadata = path.stat()
    return f"{metadata.st_dev}:{metadata.st_ino}"


def _open_descriptor_count():
    return len(os.listdir("/dev/fd"))


def _restore_packet_permissions(packet: Path):
    for root, directories, files in os.walk(packet):
        Path(root).chmod(0o755)
        for directory in directories:
            (Path(root) / directory).chmod(0o755)
        for filename in files:
            (Path(root) / filename).chmod(0o644)


def test_packet_finalizer_writes_acyclic_digests_and_removes_write_permissions(
    tmp_path,
):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    completed = _run_finalizer(packet, archive_record)

    try:
        manifest = packet / "artifacts.sha256"
        index = packet / "packet-index.md"
        manifest_text = manifest.read_text()
        archive_text = archive_record.read_text()
        assert "./packet-index.tsv" in manifest_text
        assert "./logs/retained file.txt" in manifest_text
        assert "./logs/artifacts.sha256" in manifest_text
        assert "./artifacts.sha256" not in manifest_text
        assert "./packet-index.md" not in manifest_text
        assert "artifacts_sha256=" in completed.stdout
        assert "packet_index_sha256=" in completed.stdout
        assert "artifacts.sha256" in archive_text
        assert "packet-index.md" in archive_text
        assert not (packet.stat().st_mode & stat.S_IWUSR)
        assert not (manifest.stat().st_mode & stat.S_IWUSR)
        assert not (index.stat().st_mode & stat.S_IWUSR)
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_repeats_idempotently(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    _run_finalizer(packet, archive_record)
    _restore_packet_permissions(packet)

    completed = _run_finalizer(packet, archive_record, check=False)

    try:
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
        assert not (packet.stat().st_mode & stat.S_IWUSR)
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_archive_record_inside_packet(tmp_path):
    packet = _packet(tmp_path)

    completed = _run_finalizer(packet, packet / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "outside the immutable packet" in completed.stderr


def test_packet_finalizer_rejects_packet_symlinks(tmp_path):
    packet = _packet(tmp_path)
    (packet / "logs" / "linked").symlink_to(packet / "packet-index.tsv")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "packet must not contain symlinks" in completed.stderr


def test_packet_finalizer_rejects_archive_record_symlinks(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_record.symlink_to(packet / "packet-index.tsv")

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "archive record must not be a symlink" in completed.stderr


def test_packet_finalizer_rejects_dangling_archive_record_symlinks(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    external_target = tmp_path / "missing-external.tsv"
    archive_record.symlink_to(external_target)

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "archive record must not be a symlink" in completed.stderr
    assert not external_target.exists()


def test_packet_finalizer_rejects_archive_record_hard_links(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    os.link(packet / "packet-index.tsv", archive_record)

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "packet must not contain hard-linked files" in completed.stderr


def test_packet_finalizer_rejects_external_archive_hard_link_alias(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    external_alias = tmp_path / "external-alias.tsv"
    external_alias.write_text("existing\n")
    os.link(external_alias, archive_record)

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "archive record must not have hard-link aliases" in completed.stderr


def test_packet_finalizer_rejects_malformed_existing_archive_record(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_record.write_text("malformed\n")

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


def test_packet_finalizer_rejects_conflicting_existing_archive_packet_claims(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_record.write_text(
        f"{tmp_path / 'claimed-packet'}\tartifacts.sha256\t{'a' * 64}\t"
        f"packet-index.md\t{'b' * 64}\n"
        f"{tmp_path / 'claimed-packet'}\tartifacts.sha256\t{'c' * 64}\t"
        f"packet-index.md\t{'d' * 64}\n"
    )

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode != 0
    assert "conflicts with an earlier packet record" in completed.stderr


def test_packet_finalizer_rejects_noncanonical_existing_archive_packet_claim(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_record.write_text(
        f"{tmp_path / 'packet' / '..' / 'packet'}\tartifacts.sha256\t{'a' * 64}\t"
        f"packet-index.md\t{'b' * 64}\n"
    )

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode != 0
    assert "noncanonical packet path" in completed.stderr


def test_packet_finalizer_rejects_control_characters_in_paths(tmp_path):
    packet = _packet(tmp_path / "control\tpath")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "control characters" in completed.stderr


def test_packet_finalizer_rejects_control_characters_in_artifact_paths(tmp_path):
    packet = _packet(tmp_path)
    (packet / "logs" / "bad\tartifact").write_text("bad\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "artifact paths must not contain control characters" in completed.stderr


def test_packet_finalizer_rejects_markdown_delimiters_in_packet_path(tmp_path):
    packet = _packet(tmp_path / "markdown|path")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "Markdown table delimiters" in completed.stderr


def test_packet_finalizer_rebuilds_writable_partial_metadata(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    (packet / "artifacts.sha256").write_text("partial\n")

    completed = _run_finalizer(packet, archive_record)

    try:
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
        assert not (packet.stat().st_mode & stat.S_IWUSR)
    finally:
        _restore_packet_permissions(packet)


@pytest.mark.parametrize(
    "relative_path",
    (
        ".artifacts.sha256.tmp.operator",
        ".packet-index.tsv.tmp.operator",
        "logs/.packet-index.md.tmp.operator",
    ),
)
def test_packet_finalizer_rejects_stale_helper_temporary_artifacts(
    tmp_path, relative_path
):
    packet = _packet(tmp_path)
    (packet / relative_path).write_text("stale\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode == 2
    assert "unsafe alias or reserved helper path" in completed.stderr


def test_packet_finalizer_rejects_ambiguous_owned_shape_helper_temporary(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    stale = packet / ".artifacts.sha256.tmp.0123456789abcdef"
    stale.write_text("interrupted publication\n")

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "unsafe alias or reserved helper path" in completed.stderr
    assert stale.exists()
    assert not archive_record.exists() or not archive_record.read_text()


def test_packet_finalizer_rebuilds_writable_inconsistent_manifest_pair(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    (packet / "artifacts.sha256").write_text("")
    (packet / "packet-index.md").write_text("incomplete\n")

    completed = _run_finalizer(packet, archive_record)

    try:
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
        assert "./packet-index.tsv" in (packet / "artifacts.sha256").read_text()
        assert "Inner artifact manifest" in (packet / "packet-index.md").read_text()
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_recovers_after_chmod_failure(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    bin_dir = tmp_path / "bin"
    marker = tmp_path / "chmod-failed"
    bin_dir.mkdir()
    (bin_dir / "chmod").write_text(
        f"""#!/usr/bin/env bash
if [[ ! -e "{marker}" ]]; then
  : > "{marker}"
  exit 9
fi
exec /bin/chmod "$@"
"""
    )
    (bin_dir / "chmod").chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    failed = _run_finalizer(packet, archive_record, check=False, env=env)
    assert not archive_record.exists() or not archive_record.read_text()
    completed = _run_finalizer(packet, archive_record)

    try:
        assert failed.returncode == 9
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
        assert not (packet.stat().st_mode & stat.S_IWUSR)
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_recovers_late_artifact_after_chmod_failure(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    bin_dir = tmp_path / "bin"
    marker = tmp_path / "chmod-failed"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "chmod",
        f"""#!/usr/bin/env bash
if [[ ! -e "{marker}" ]]; then
  : > "{marker}"
  exit 9
fi
exec /bin/chmod "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    failed = _run_finalizer(packet, archive_record, check=False, env=env)
    (packet / "late-artifact").write_text("late\n")
    completed = _run_finalizer(packet, archive_record)

    try:
        assert failed.returncode == 9
        assert completed.returncode == 0
        assert "./late-artifact" in (packet / "artifacts.sha256").read_text()
        assert len(archive_record.read_text().splitlines()) == 1
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_read_only_artifact_created_during_chmod(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "chmod",
        f"""#!/usr/bin/env bash
if [[ "$*" == *"{packet}"* && ! -e "{packet / 'late-during-chmod.txt'}" ]]; then
  printf 'late\\n' > "{packet / 'late-during-chmod.txt'}"
fi
exec /bin/chmod "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, check=False, env=env)

    try:
        assert completed.returncode == 2
        assert "inventory changed after write-permission removal" in completed.stderr
        assert not archive_record.exists() or not archive_record.read_text()
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_symlink_created_during_chmod(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    outside = tmp_path / "outside.txt"
    outside.write_text("outside\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "chmod",
        f"""#!/usr/bin/env bash
if [[ "$*" == *"{packet}"* && ! -L "{packet / 'late-link'}" ]]; then
  ln -s "{outside}" "{packet / 'late-link'}"
fi
exec /bin/chmod "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, check=False, env=env)

    try:
        assert completed.returncode == 2
        assert "unsafe alias or reserved helper path" in completed.stderr
        assert not archive_record.exists() or not archive_record.read_text()
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_symlink_created_during_manifest_scan(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    outside = tmp_path / "outside.txt"
    outside.write_text("outside\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "find",
        f"""#!/usr/bin/env bash
if [[ "$*" == *"-type f ! -path"* && ! -L "{packet / 'late-link'}" ]]; then
  ln -s "{outside}" "{packet / 'late-link'}"
fi
exec /usr/bin/find "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, check=False, env=env)

    try:
        assert completed.returncode == 2
        assert "unsafe alias or reserved helper path" in completed.stderr
        assert not archive_record.read_text()
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_permission_scan_failure(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "chmod",
        """#!/usr/bin/env bash
exit 0
""",
    )
    _write_executable(
        bin_dir / "find",
        """#!/usr/bin/env bash
if [[ " $* " == *" -perm "* ]]; then
  exit 9
fi
exec /usr/bin/find "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, check=False, env=env)

    assert completed.returncode != 0
    assert "write-permission scan failed" in completed.stderr
    assert packet.stat().st_mode & stat.S_IWUSR
    assert not archive_record.exists() or not archive_record.read_text()


def test_packet_finalizer_rejects_malformed_runner_index(tmp_path):
    packet = _packet(tmp_path)
    (packet / "packet-index.tsv").write_text("wrong\theader\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


def test_packet_finalizer_rejects_header_only_runner_index(tmp_path):
    packet = _packet(tmp_path)
    header = (packet / "packet-index.tsv").read_text().splitlines()[0]
    (packet / "packet-index.tsv").write_text(f"{header}\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode != 0
    assert "must contain at least one data row" in completed.stderr


@pytest.mark.parametrize(
    "row",
    (
        "ED-1\tgenerate\tattempt-1\trank-map\t0\t9\tpassed-rank-map\t1\tN/A\tN/A",
        "ED-1\tgenerate\tattempt-1\trank-map\t0\t0\tfailed-rank-map\t1\tN/A\tN/A",
        "ED-1\tgenerate\tattempt-1\trank-map\t0\t-1\tfailed-rank-map\t1\tN/A\tN/A",
    ),
)
def test_packet_finalizer_rejects_contradictory_runner_index_exit_status(tmp_path, row):
    packet = _packet(tmp_path)
    header = (packet / "packet-index.tsv").read_text().splitlines()[0]
    (packet / "packet-index.tsv").write_text(f"{header}\n{row}\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


@pytest.mark.parametrize(
    "rows",
    (
        (
            "ED-1\tgenerate\tattempt-1\trun\t0\t0\tpassed\t1\t"
            "/tmp/output\tINCOMPLETE-unset"
        ),
        (
            "MR-1\tresume\tattempt-1\trank-map\t0\t0\tpassed-rank-map\t1\tN/A\tN/A\n"
            "MR-1\tresume\tattempt-1\tsample-5\t1\t0\tpassed\t1\t"
            "/tmp/output\t/tmp/hook"
        ),
        (
            "ED-1\tgenerate\tattempt-1\trank-map\t0\t2\tfailed-rank-map\t1\tN/A\tN/A\n"
            "ED-1\tgenerate\tattempt-1\trun\t0\t0\tpassed\t1\t"
            "/tmp/output\tINCOMPLETE-unset"
        ),
    ),
)
def test_packet_finalizer_rejects_invalid_runner_lifecycle_prefix(tmp_path, rows):
    packet = _packet(tmp_path)
    header = (packet / "packet-index.tsv").read_text().splitlines()[0]
    (packet / "packet-index.tsv").write_text(f"{header}\n{rows}\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


@pytest.mark.parametrize(
    "rows",
    (
        (
            "ED-1\tgenerate\tattempt-1\trank-map\t0\t0\tpassed-rank-map\t1\tN/A\tN/A\n"
            "ED-1\tgenerate\tattempt-1\trun\t0\t3\tincomplete-accounting-before\t1\t"
            "/tmp/output\tINCOMPLETE-unset"
        ),
        (
            "MR-1\tresume\tattempt-1\trank-map\t0\t0\tpassed-rank-map\t1\tN/A\tN/A\n"
            "MR-1\tresume\tattempt-1\twarmup\t0\t9\tincomplete-accounting-after\t1\t"
            "/tmp/output\t/tmp/hook"
        ),
        (
            "MR-1\tresume\tattempt-1\trank-map\t0\t0\tpassed-rank-map\t1\tN/A\tN/A\n"
            "MR-1\tresume\tattempt-1\twarmup\t0\t0\tpassed\t1\t"
            "/tmp/output\tINCOMPLETE-unset\n"
            "MR-1\tresume\tattempt-1\tsample-1\t1\t0\tpassed\t1\t"
            "/tmp/output\tINCOMPLETE-unset"
        ),
    ),
)
def test_packet_finalizer_rejects_impossible_accounting_history(tmp_path, rows):
    packet = _packet(tmp_path)
    header = (packet / "packet-index.tsv").read_text().splitlines()[0]
    (packet / "packet-index.tsv").write_text(f"{header}\n{rows}\n")

    completed = _run_finalizer(packet, tmp_path / "archive.tsv", check=False)

    assert completed.returncode != 0
    assert "invalid packet index" in completed.stderr


def test_packet_index_validator_rejects_candidate_row_without_mutating_index(tmp_path):
    packet = _packet(tmp_path)
    index = packet / "packet-index.tsv"
    before = index.read_text()

    completed = _run_validator(
        "--append-row",
        str(index),
        f"{index.parent.stat().st_dev}:{index.parent.stat().st_ino}",
        f"{index.stat().st_dev}:{index.stat().st_ino}",
        "ED-1",
        "generate",
        "attempt-1",
        "run",
        "0",
        "9",
        "passed",
        "1",
        "/tmp/output",
        "INCOMPLETE-unset",
        check=False,
    )

    assert completed.returncode != 0
    assert "passing disposition with a nonzero exit code" in completed.stderr
    assert index.read_text() == before


def test_packet_finalizer_recovers_writable_retained_markdown_index_mutation(tmp_path):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    _run_finalizer(packet, archive_record)
    _restore_packet_permissions(packet)
    with (packet / "packet-index.md").open("a") as handle:
        handle.write("\nmutated\n")

    completed = _run_finalizer(packet, archive_record)

    try:
        assert completed.returncode == 0
        assert "mutated" not in (packet / "packet-index.md").read_text()
        assert len(archive_record.read_text().splitlines()) == 1
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_recovers_after_archive_append_failure(tmp_path):
    packet = _packet(tmp_path)
    archive_dir = tmp_path / "archive"
    archive_record = archive_dir / "archive.tsv"
    bin_dir = tmp_path / "bin"
    archive_dir.mkdir()
    bin_dir.mkdir()
    (bin_dir / "chmod").write_text(
        f"""#!/usr/bin/env bash
if [[ "$*" == *"{packet}"* ]]; then
  /bin/chmod a-w "{archive_record}"
fi
exec /bin/chmod "$@"
"""
    )
    (bin_dir / "chmod").chmod(0o755)
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    try:
        failed = _run_finalizer(packet, archive_record, check=False, env=env)
    finally:
        archive_record.chmod(0o644)
    completed = _run_finalizer(packet, archive_record)

    try:
        assert failed.returncode != 0
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
        assert not (packet.stat().st_mode & stat.S_IWUSR)
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_has_no_postpublication_grep_mutation_window(
    tmp_path,
):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_backup = tmp_path / "archive-backup.tsv"
    outside = tmp_path / "outside.tsv"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    archive_record.write_text(
        f"{tmp_path / 'existing-packet'}\tartifacts.sha256\t{'a' * 64}\t"
        f"packet-index.md\t{'b' * 64}\n"
    )
    _write_executable(
        bin_dir / "grep",
        f"""#!/usr/bin/env bash
if [[ ! -e "{archive_backup}" ]]; then
  mv "{archive_record}" "{archive_backup}"
  ln -s "{outside}" "{archive_record}"
fi
exec /usr/bin/grep "$@"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, env=env)

    try:
        assert completed.returncode == 0
        assert archive_record.is_file()
        assert not archive_record.is_symlink()
        assert not archive_backup.exists()
        assert not outside.exists()
    finally:
        _restore_packet_permissions(packet)


def test_packet_finalizer_rejects_clean_archive_record_replacement_before_append(
    tmp_path,
):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_backup = tmp_path / "archive-backup.tsv"
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    _write_executable(
        bin_dir / "chmod",
        f"""#!/usr/bin/env bash
/bin/chmod "$@"
status=$?
if [[ "$*" == *"{packet}"* && ! -e "{archive_backup}" ]]; then
  mv "{archive_record}" "{archive_backup}"
  cp "{archive_backup}" "{archive_record}"
fi
exit "$status"
""",
    )
    env = os.environ.copy()
    env["PATH"] = f"{bin_dir}:{env['PATH']}"

    completed = _run_finalizer(packet, archive_record, check=False, env=env)

    try:
        assert completed.returncode == 2
        assert "archive record identity changed" in completed.stderr
        assert archive_backup.is_file()
        assert not archive_record.read_text()
    finally:
        _restore_packet_permissions(packet)


def test_packet_index_atomic_append_preserves_original_after_failed_short_write(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    packet = _packet(tmp_path)
    index = packet / "packet-index.tsv"
    before = index.read_text()
    original_write = validator.os.write
    calls = 0

    def partial_then_fail(descriptor, payload):
        nonlocal calls
        calls += 1
        if calls == 1:
            return original_write(descriptor, payload[:3])
        raise OSError("injected packet-index write failure")

    monkeypatch.setattr(validator.os, "write", partial_then_fail)
    fields = [
        "ED-1",
        "generate",
        "attempt-1",
        "run",
        "0",
        "0",
        "passed",
        "2",
        "/tmp/output",
        "INCOMPLETE-unset",
    ]

    with pytest.raises(OSError, match="injected packet-index write failure"):
        validator.append_packet_index(
            index, _identity(packet), _identity(index), fields
        )

    assert index.read_text() == before
    monkeypatch.setattr(validator.os, "write", original_write)
    validator.append_packet_index(index, _identity(packet), _identity(index), fields)
    assert index.read_text().splitlines()[-1] == "\t".join(fields)


def test_archive_atomic_append_preserves_original_after_failed_short_write(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    prior = (
        f"{tmp_path / 'existing-packet'}\tartifacts.sha256\t{'c' * 64}\t"
        f"packet-index.md\t{'d' * 64}\n"
    )
    archive_record.write_text(prior)
    original_write = validator.os.write
    calls = 0

    def partial_then_fail(descriptor, payload):
        nonlocal calls
        calls += 1
        if calls == 1:
            return original_write(descriptor, payload[:3])
        raise OSError("injected archive-record write failure")

    monkeypatch.setattr(validator.os, "write", partial_then_fail)
    fields = [
        str(packet.resolve()),
        "artifacts.sha256",
        "a" * 64,
        "packet-index.md",
        "b" * 64,
    ]
    arguments = (
        archive_record,
        _identity(tmp_path),
        _identity(archive_record),
        _identity(packet),
        fields,
    )

    with pytest.raises(OSError, match="injected archive-record write failure"):
        validator.append_archive_record(*arguments)

    assert archive_record.read_text() == prior
    monkeypatch.setattr(validator.os, "write", original_write)
    validator.append_archive_record(*arguments)
    assert archive_record.read_text().splitlines() == [
        prior.rstrip("\n"),
        "\t".join(fields),
    ]


def test_packet_local_metadata_publication_recovers_partial_writes(tmp_path, monkeypatch):
    validator = _validator_module()
    source = tmp_path / "source"
    packet = tmp_path / "packet"
    destination = packet / "packet-index.md"
    source.write_text("metadata publication payload\n")
    packet.mkdir()
    original_write = validator.os.write

    def short_write(descriptor, payload):
        return original_write(descriptor, payload[: max(1, len(payload) // 2)])

    monkeypatch.setattr(validator.os, "write", short_write)
    validator.publish_file(source, destination, _identity(packet))

    assert destination.read_text() == source.read_text()
    assert list(packet.glob(".packet-index.md.tmp.*")) == []


def test_packet_local_metadata_can_resync_after_rename_directory_sync_failure(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    source = tmp_path / "source"
    packet = tmp_path / "packet"
    destination = packet / "packet-index.md"
    source.write_text("metadata publication payload\n")
    packet.mkdir()
    original_fsync = validator.os.fsync

    def fail_directory_sync(descriptor):
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise OSError("injected packet-directory sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "fsync", fail_directory_sync)
    with pytest.raises(OSError, match="injected packet-directory sync failure"):
        validator.publish_file(source, destination, _identity(packet))

    assert destination.read_text() == source.read_text()
    monkeypatch.setattr(validator.os, "fsync", original_fsync)
    validator.sync_directory(packet, _identity(packet))


def test_archive_sink_creation_resyncs_existing_sink_after_parent_sync_failure(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    archive_record = tmp_path / "archive.tsv"
    original_fsync = validator.os.fsync

    def fail_directory_sync(descriptor):
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise OSError("injected archive-directory sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "fsync", fail_directory_sync)
    with pytest.raises(OSError, match="injected archive-directory sync failure"):
        validator.ensure_archive_record(archive_record, _identity(tmp_path))

    assert archive_record.is_file()
    monkeypatch.setattr(validator.os, "fsync", original_fsync)
    validator.ensure_archive_record(archive_record, _identity(tmp_path))


def test_packet_index_atomic_append_is_complete_after_parent_sync_failure(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    packet = _packet(tmp_path)
    index = packet / "packet-index.tsv"
    original_fsync = validator.os.fsync

    def fail_directory_sync(descriptor):
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise OSError("injected index-directory sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "fsync", fail_directory_sync)
    fields = [
        "ED-1",
        "generate",
        "attempt-1",
        "run",
        "0",
        "0",
        "passed",
        "2",
        "/tmp/output",
        "INCOMPLETE-unset",
    ]

    with pytest.raises(OSError, match="injected index-directory sync failure"):
        validator.append_packet_index(
            index, _identity(packet), _identity(index), fields
        )

    monkeypatch.setattr(validator.os, "fsync", original_fsync)
    validator.validate_packet_index(index)
    assert index.read_text().splitlines()[-1] == "\t".join(fields)
    validator.sync_directory(packet, _identity(packet))


def test_archive_atomic_append_is_complete_after_parent_sync_failure(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    prior = (
        f"{tmp_path / 'existing-packet'}\tartifacts.sha256\t{'c' * 64}\t"
        f"packet-index.md\t{'d' * 64}\n"
    )
    archive_record.write_text(prior)
    original_fsync = validator.os.fsync
    directory_syncs = 0

    def fail_replacement_directory_sync(descriptor):
        nonlocal directory_syncs
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            directory_syncs += 1
            if directory_syncs == 2:
                raise OSError("injected archive-replacement sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "fsync", fail_replacement_directory_sync)
    fields = [
        str(packet.resolve()),
        "artifacts.sha256",
        "a" * 64,
        "packet-index.md",
        "b" * 64,
    ]

    with pytest.raises(OSError, match="injected archive-replacement sync failure"):
        validator.append_archive_record(
            archive_record,
            _identity(tmp_path),
            _identity(archive_record),
            _identity(packet),
            fields,
        )

    monkeypatch.setattr(validator.os, "fsync", original_fsync)
    validator.append_archive_record(
        archive_record,
        _identity(tmp_path),
        _identity(archive_record),
        _identity(packet),
        fields,
    )
    assert archive_record.read_text().splitlines() == [
        prior.rstrip("\n"),
        "\t".join(fields),
    ]


@pytest.mark.parametrize(
    "temporary_name",
    (".archive.tsv.tmp.0123456789abcdef", ".archive.tsv.tmp.operator"),
)
def test_packet_finalizer_rejects_archive_sibling_replacement_temporaries(
    tmp_path, temporary_name
):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    temporary = tmp_path / temporary_name
    temporary.write_text("ambiguous archive replacement\n")

    completed = _run_finalizer(packet, archive_record, check=False)

    assert completed.returncode == 2
    assert "reserved replacement temporary requires operator adjudication" in (
        completed.stderr
    )
    assert temporary.read_text() == "ambiguous archive replacement\n"
    assert not archive_record.exists()


def test_packet_finalizer_succeeds_after_operator_removes_archive_sibling_temporary(
    tmp_path,
):
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    temporary = tmp_path / ".archive.tsv.tmp.0123456789abcdef"
    temporary.write_text("ambiguous archive replacement\n")

    rejected = _run_finalizer(packet, archive_record, check=False)
    temporary.unlink()
    completed = _run_finalizer(packet, archive_record)

    try:
        assert rejected.returncode == 2
        assert completed.returncode == 0
        assert len(archive_record.read_text().splitlines()) == 1
    finally:
        _restore_packet_permissions(packet)


def test_archive_atomic_append_kill_before_rename_requires_operator_adjudication(
    tmp_path,
):
    validator = _validator_module()
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    prior = (
        f"{tmp_path / 'existing-packet'}\tartifacts.sha256\t{'c' * 64}\t"
        f"packet-index.md\t{'d' * 64}\n"
    )
    archive_record.write_text(prior)
    fields = [
        str(packet.resolve()),
        "artifacts.sha256",
        "a" * 64,
        "packet-index.md",
        "b" * 64,
    ]
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            f"""
from pathlib import Path
import importlib.util
import os
import signal
spec = importlib.util.spec_from_file_location(
    "external_io_packet_index", {str(VALIDATOR)!r}
)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)
def kill_before_rename(*args, **kwargs):
    os.kill(os.getpid(), signal.SIGKILL)
validator.os.rename = kill_before_rename
validator.append_archive_record(
    Path({str(archive_record)!r}),
    {_identity(tmp_path)!r},
    {_identity(archive_record)!r},
    {_identity(packet)!r},
    {fields!r},
)
""",
        ],
        check=False,
    )

    assert child.returncode == -signal.SIGKILL
    assert archive_record.read_text() == prior
    temporaries = list(tmp_path.glob(".archive.tsv.tmp.*"))
    assert len(temporaries) == 1
    with pytest.raises(SystemExit, match="operator adjudication"):
        validator.append_archive_record(
            archive_record,
            _identity(tmp_path),
            _identity(archive_record),
            _identity(packet),
            fields,
        )
    temporaries[0].unlink()
    validator.append_archive_record(
        archive_record,
        _identity(tmp_path),
        _identity(archive_record),
        _identity(packet),
        fields,
    )
    assert archive_record.read_text().splitlines() == [
        prior.rstrip("\n"),
        "\t".join(fields),
    ]


def test_archive_atomic_append_kill_after_rename_retries_idempotently(tmp_path):
    validator = _validator_module()
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    prior = (
        f"{tmp_path / 'existing-packet'}\tartifacts.sha256\t{'c' * 64}\t"
        f"packet-index.md\t{'d' * 64}\n"
    )
    archive_record.write_text(prior)
    fields = [
        str(packet.resolve()),
        "artifacts.sha256",
        "a" * 64,
        "packet-index.md",
        "b" * 64,
    ]
    child = subprocess.run(
        [
            sys.executable,
            "-c",
            f"""
from pathlib import Path
import importlib.util
import os
import signal
import stat
spec = importlib.util.spec_from_file_location(
    "external_io_packet_index", {str(VALIDATOR)!r}
)
validator = importlib.util.module_from_spec(spec)
spec.loader.exec_module(validator)
original_fsync = validator.os.fsync
directory_syncs = 0
def kill_after_rename(descriptor):
    global directory_syncs
    if stat.S_ISDIR(os.fstat(descriptor).st_mode):
        directory_syncs += 1
        if directory_syncs == 2:
            os.kill(os.getpid(), signal.SIGKILL)
    return original_fsync(descriptor)
validator.os.fsync = kill_after_rename
validator.append_archive_record(
    Path({str(archive_record)!r}),
    {_identity(tmp_path)!r},
    {_identity(archive_record)!r},
    {_identity(packet)!r},
    {fields!r},
)
""",
        ],
        check=False,
    )

    assert child.returncode == -signal.SIGKILL
    expected = [prior.rstrip("\n"), "\t".join(fields)]
    assert archive_record.read_text().splitlines() == expected
    assert list(tmp_path.glob(".archive.tsv.tmp.*")) == []
    validator.append_archive_record(
        archive_record,
        _identity(tmp_path),
        _identity(archive_record),
        _identity(packet),
        fields,
    )
    assert archive_record.read_text().splitlines() == expected


def test_archive_validation_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    archive_record = tmp_path / "archive.tsv"
    archive_record.write_text("malformed\n")
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="archive line 1"):
            validator.ensure_archive_record(archive_record, _identity(tmp_path))

    assert _open_descriptor_count() == before


def test_archive_identity_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    packet = _packet(tmp_path)
    archive_record = tmp_path / "archive.tsv"
    archive_record.write_text("")
    fields = [
        str(packet.resolve()),
        "artifacts.sha256",
        "a" * 64,
        "packet-index.md",
        "b" * 64,
    ]
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="archive record identity changed"):
            validator.append_archive_record(
                archive_record,
                _identity(tmp_path),
                "0:0",
                _identity(packet),
                fields,
            )

    assert _open_descriptor_count() == before


def test_directory_identity_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="parent directory identity changed"):
            validator.sync_directory(tmp_path, "0:0")

    assert _open_descriptor_count() == before


def test_directory_sync_rejections_do_not_leak_descriptors(tmp_path, monkeypatch):
    validator = _validator_module()
    original_fsync = validator.os.fsync

    def fail_directory_sync(descriptor):
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise OSError("injected directory sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "fsync", fail_directory_sync)
    before = _open_descriptor_count()
    for _ in range(16):
        with pytest.raises(OSError, match="injected directory sync failure"):
            validator.sync_directory(tmp_path, _identity(tmp_path))

    assert _open_descriptor_count() == before


def test_publication_target_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.write_text("source\n")
    destination.write_text("destination\n")
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="publication target already exists"):
            validator.publish_file(source, destination, _identity(tmp_path))

    assert _open_descriptor_count() == before


def test_publication_rollback_sync_rejections_do_not_leak_descriptors(
    tmp_path, monkeypatch
):
    validator = _validator_module()
    source = tmp_path / "source"
    destination = tmp_path / "destination"
    source.write_text("source\n")
    original_fsync = validator.os.fsync

    def fail_write(descriptor, payload):
        raise OSError("injected publication write failure")

    def fail_cleanup_sync(descriptor):
        if stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise OSError("injected publication cleanup sync failure")
        return original_fsync(descriptor)

    monkeypatch.setattr(validator.os, "write", fail_write)
    monkeypatch.setattr(validator.os, "fsync", fail_cleanup_sync)
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(OSError, match="injected publication cleanup sync failure"):
            validator.publish_file(source, destination, _identity(tmp_path))

    assert list(tmp_path.glob(".destination.tmp.*")) == []
    assert _open_descriptor_count() == before


def test_metadata_reset_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    manifest = tmp_path / "artifacts.sha256"
    index = tmp_path / "packet-index.md"
    alias = tmp_path / "manifest-alias"
    manifest.write_text("manifest\n")
    os.link(manifest, alias)
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="metadata changed or has aliases"):
            validator.reset_metadata_pair(manifest, index, _identity(tmp_path))

    assert _open_descriptor_count() == before


def test_packet_index_identity_rejections_do_not_leak_descriptors(tmp_path):
    validator = _validator_module()
    packet = _packet(tmp_path)
    index = packet / "packet-index.tsv"
    fields = [
        "ED-1",
        "generate",
        "attempt-1",
        "run",
        "0",
        "0",
        "passed",
        "2",
        "/tmp/output",
        "INCOMPLETE-unset",
    ]
    before = _open_descriptor_count()

    for _ in range(16):
        with pytest.raises(SystemExit, match="packet index changed or has aliases"):
            validator.append_packet_index(index, _identity(packet), "0:0", fields)

    assert _open_descriptor_count() == before
