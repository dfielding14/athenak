#!/usr/bin/env python3
"""Validate the strict TSV grammar for external IO qualification packet rows."""

from pathlib import Path
import os
import re
import secrets
import stat
import sys


HEADER = (
    "row\taction\tattempt_id\tsample\tmeasured\texit_code\tdisposition\t"
    "monotonic_elapsed_ns\toutput_dir\tfs_accounting_hook"
)
TOKEN = re.compile(r"[A-Za-z0-9._-]+\Z")
INTEGER = re.compile(r"[0-9]+\Z")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
ALLOWED_ROWS = {"ED-1", "MR-1", "MR-2", "scaling"}
ALLOWED_ACTIONS = {"generate", "resume"}
ALLOWED_MEASURED = {"0", "1"}
ALLOWED_ROW_ACTIONS = {
    ("ED-1", "generate"),
    ("MR-1", "generate"),
    ("MR-1", "resume"),
    ("MR-2", "generate"),
    ("MR-2", "resume"),
    ("scaling", "resume"),
}
RANK_MAP_DISPOSITIONS = {
    "failed-invalid-rank-map",
    "failed-invalid-rank-map-timing",
    "failed-rank-map",
    "passed-rank-map",
    "timeout-rank-map",
}
LAUNCH_DISPOSITIONS = {
    "failed",
    "failed-and-incomplete-accounting-after",
    "failed-forbidden-assembled",
    "failed-inventory-scan",
    "failed-invalid-launch-timing",
    "failed-staging-scan",
    "incomplete-accounting-after",
    "incomplete-accounting-before",
    "passed",
    "timeout-incomplete",
}
PASSING_DISPOSITIONS = {"passed", "passed-rank-map"}
TIMEOUT_DISPOSITIONS = {"timeout-incomplete", "timeout-rank-map"}
ACCOUNTING_DISPOSITIONS = {
    "failed-and-incomplete-accounting-after",
    "incomplete-accounting-after",
    "incomplete-accounting-before",
}


def fail(message: str):
    raise SystemExit(f"invalid packet index: {message}")


def has_control(value: str):
    return any(ord(character) < 32 or ord(character) == 127 for character in value)


def validate_field(value: str, label: str):
    if not value or has_control(value):
        fail(f"{label} contains an empty value or control character")


def validate_row(fields, line_number: int):
    if len(fields) != 10:
        fail(f"line {line_number} has {len(fields)} columns, expected 10")
    row, action, attempt, sample, measured = fields[:5]
    exit_code, disposition, elapsed, output, hook = fields[5:]
    if row not in ALLOWED_ROWS:
        fail(f"line {line_number} has unsupported row {row!r}")
    if action not in ALLOWED_ACTIONS:
        fail(f"line {line_number} has unsupported action {action!r}")
    if (row, action) not in ALLOWED_ROW_ACTIONS:
        fail(f"line {line_number} has unsupported row/action pair")
    for value, label in (
        (attempt, "attempt_id"),
        (sample, "sample"),
        (disposition, "disposition"),
    ):
        if TOKEN.fullmatch(value) is None:
            fail(f"line {line_number} has invalid {label} {value!r}")
    if measured not in ALLOWED_MEASURED:
        fail(f"line {line_number} has invalid measured flag {measured!r}")
    if INTEGER.fullmatch(exit_code) is None:
        fail(f"line {line_number} has invalid exit code {exit_code!r}")
    if INTEGER.fullmatch(elapsed) is None:
        fail(f"line {line_number} has invalid elapsed interval {elapsed!r}")
    validate_field(output, "output_dir")
    validate_field(hook, "fs_accounting_hook")
    if sample == "rank-map":
        if measured != "0" or disposition not in RANK_MAP_DISPOSITIONS:
            fail(f"line {line_number} has invalid rank-map semantics")
        if output != "N/A" or hook != "N/A":
            fail(f"line {line_number} has invalid rank-map metadata")
    else:
        if disposition not in LAUNCH_DISPOSITIONS:
            fail(f"line {line_number} has invalid launcher disposition")
        if action == "generate":
            if sample != "run" or measured != "0":
                fail(f"line {line_number} has invalid generate semantics")
        elif sample == "warmup":
            if measured != "0":
                fail(f"line {line_number} has invalid warmup semantics")
        elif sample not in {f"sample-{index}" for index in range(1, 6)}:
            fail(f"line {line_number} has invalid resume sample")
        elif measured != "1":
            fail(f"line {line_number} has invalid measured-resume semantics")
    if disposition in PASSING_DISPOSITIONS and exit_code != "0":
        fail(f"line {line_number} has a passing disposition with a nonzero exit code")
    if disposition in TIMEOUT_DISPOSITIONS and exit_code != "124":
        fail(f"line {line_number} has a timeout disposition without exit code 124")
    if (
        disposition not in PASSING_DISPOSITIONS | TIMEOUT_DISPOSITIONS
        and exit_code == "0"
    ):
        fail(f"line {line_number} has a non-passing disposition with exit code 0")
    if disposition in ACCOUNTING_DISPOSITIONS and measured != "1":
        fail(f"line {line_number} has accounting disposition for an unmeasured launch")
    if (
        measured == "1"
        and hook == "INCOMPLETE-unset"
        and disposition != "incomplete-accounting-before"
    ):
        fail(f"line {line_number} has measured launch without an accounting hook")


def validate_lifecycle(groups):
    measured_samples = [f"sample-{index}" for index in range(1, 6)]
    for key, rows in groups.items():
        row, action, attempt = key
        expected = ["rank-map", "run"] if action == "generate" else [
            "rank-map",
            "warmup",
            *measured_samples,
        ]
        samples = [fields[3] for fields in rows]
        if samples != expected[: len(samples)]:
            fail(f"{row} {action} {attempt} has invalid lifecycle prefix {samples!r}")
        for fields in rows[:-1]:
            if fields[6] not in PASSING_DISPOSITIONS:
                fail(
                    f"{row} {action} {attempt} continues after terminal "
                    f"disposition {fields[6]!r}"
                )


def validate_packet_index_lines(lines, require_rows=False):
    if not lines or lines[0] != HEADER:
        fail("header does not match the version-1 packet index grammar")
    if require_rows and len(lines) == 1:
        fail("finalized packet index must contain at least one data row")
    seen = set()
    groups = {}
    for line_number, line in enumerate(lines[1:], start=2):
        fields = line.split("\t")
        validate_row(fields, line_number)
        launch_key = tuple(fields[:4])
        if launch_key in seen:
            fail(f"line {line_number} duplicates launch key {launch_key!r}")
        seen.add(launch_key)
        groups.setdefault(tuple(fields[:3]), []).append(fields)
    validate_lifecycle(groups)


def validate_packet_index(index: Path, require_rows=False):
    if not index.is_file():
        fail(f"not a regular file: {index}")
    validate_packet_index_lines(index.read_text().splitlines(), require_rows)


def validate_archive_lines(lines):
    seen = {}
    for line_number, line in enumerate(lines, start=1):
        fields = line.split("\t")
        if len(fields) != 5:
            fail(f"archive line {line_number} has {len(fields)} columns, expected 5")
        packet, manifest, manifest_digest, index, index_digest = fields
        validate_field(packet, "archive packet path")
        if not Path(packet).is_absolute():
            fail(f"archive line {line_number} has a non-absolute packet path")
        canonical_packet = os.path.realpath(packet)
        if packet != canonical_packet:
            fail(f"archive line {line_number} has a noncanonical packet path")
        if manifest != "artifacts.sha256" or index != "packet-index.md":
            fail(f"archive line {line_number} has unsupported artifact names")
        if (
            SHA256.fullmatch(manifest_digest) is None
            or SHA256.fullmatch(index_digest) is None
        ):
            fail(f"archive line {line_number} has invalid SHA-256 digests")
        if canonical_packet in seen and seen[canonical_packet] != line:
            fail(f"archive line {line_number} conflicts with an earlier packet record")
        seen[canonical_packet] = line


def validate_archive_record(record: Path):
    if not record.is_file():
        fail(f"not a regular file: {record}")
    validate_archive_lines(record.read_text().splitlines())


def read_descriptor_text(descriptor):
    os.lseek(descriptor, 0, os.SEEK_SET)
    chunks = []
    while True:
        chunk = os.read(descriptor, 65536)
        if not chunk:
            return b"".join(chunks).decode()
        chunks.append(chunk)


def descriptor_identity(metadata):
    return f"{metadata.st_dev}:{metadata.st_ino}"


def open_verified_directory(directory: Path, expected_identity, label):
    directory_flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0)
    directory_flags |= getattr(os, "O_NOFOLLOW", 0)
    directory_fd = os.open(directory, directory_flags)
    try:
        directory_metadata = os.fstat(directory_fd)
        if descriptor_identity(directory_metadata) != expected_identity:
            fail(f"{label} parent directory identity changed")
        current_directory_metadata = directory.lstat()
        if (
            not stat.S_ISDIR(current_directory_metadata.st_mode)
            or not os.path.samestat(directory_metadata, current_directory_metadata)
        ):
            fail(f"{label} parent directory changed")
        return directory_fd
    except BaseException:
        os.close(directory_fd)
        raise


def write_all(descriptor, payload):
    remaining = memoryview(payload)
    while remaining:
        written = os.write(descriptor, remaining)
        if written <= 0:
            raise OSError("write returned no progress")
        remaining = remaining[written:]


def sync_directory(directory: Path, expected_directory_identity):
    directory_fd = open_verified_directory(
        directory, expected_directory_identity, "directory sync"
    )
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def atomic_replace_payload(directory_fd, destination: Path, payload, mode=0o644):
    temporary_name = f".{destination.name}.tmp.{secrets.token_hex(8)}"
    file_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    file_flags |= getattr(os, "O_NOFOLLOW", 0)
    temporary_fd = os.open(temporary_name, file_flags, mode, dir_fd=directory_fd)
    try:
        write_all(temporary_fd, payload)
        os.fsync(temporary_fd)
        os.close(temporary_fd)
        temporary_fd = None
        os.rename(
            temporary_name,
            destination.name,
            src_dir_fd=directory_fd,
            dst_dir_fd=directory_fd,
        )
        os.fsync(directory_fd)
    finally:
        if temporary_fd is not None:
            os.close(temporary_fd)
        try:
            os.unlink(temporary_name, dir_fd=directory_fd)
            os.fsync(directory_fd)
        except FileNotFoundError:
            pass


def reject_destination_temporaries(directory_fd, destination: Path):
    temporary_prefix = f".{destination.name}.tmp."
    for name in os.listdir(directory_fd):
        if name.startswith(temporary_prefix):
            fail(
                "reserved replacement temporary requires operator adjudication: "
                f"{destination.parent / name}"
            )


def publish_file(source: Path, destination: Path, expected_directory_identity):
    directory_fd = open_verified_directory(
        destination.parent, expected_directory_identity, "publication"
    )
    temporary_name = f".{destination.name}.tmp.{secrets.token_hex(8)}"
    file_flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    file_flags |= getattr(os, "O_NOFOLLOW", 0)
    temporary_fd = None
    try:
        if destination.exists() or destination.is_symlink():
            fail(f"publication target already exists: {destination}")
        temporary_fd = os.open(
            temporary_name, file_flags, 0o600, dir_fd=directory_fd
        )
        write_all(temporary_fd, source.read_bytes())
        os.fsync(temporary_fd)
        os.close(temporary_fd)
        temporary_fd = None
        os.rename(
            temporary_name,
            destination.name,
            src_dir_fd=directory_fd,
            dst_dir_fd=directory_fd,
        )
        os.fsync(directory_fd)
        published_metadata = destination.lstat()
        if (
            not stat.S_ISREG(published_metadata.st_mode)
            or published_metadata.st_nlink != 1
        ):
            fail(f"published metadata changed or has aliases: {destination}")
    finally:
        try:
            if temporary_fd is not None:
                os.close(temporary_fd)
            try:
                os.unlink(temporary_name, dir_fd=directory_fd)
                os.fsync(directory_fd)
            except FileNotFoundError:
                pass
        finally:
            os.close(directory_fd)


def reset_metadata_pair(
    manifest: Path, index: Path, expected_directory_identity
):
    if manifest.parent != index.parent:
        fail("finalization metadata must share one parent directory")
    directory_fd = open_verified_directory(
        manifest.parent, expected_directory_identity, "metadata reset"
    )
    try:
        for path in (manifest, index):
            try:
                metadata = os.stat(
                    path.name, dir_fd=directory_fd, follow_symlinks=False
                )
            except FileNotFoundError:
                continue
            if not stat.S_ISREG(metadata.st_mode) or metadata.st_nlink != 1:
                fail(f"finalization metadata changed or has aliases: {path}")
            os.unlink(path.name, dir_fd=directory_fd)
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def append_packet_index(
    index: Path, expected_directory_identity, expected_index_identity, fields
):
    validate_packet_index(index)
    validate_packet_index_lines(index.read_text().splitlines() + ["\t".join(fields)])
    directory_fd = open_verified_directory(
        index.parent, expected_directory_identity, "packet-index"
    )
    file_flags = os.O_RDWR | os.O_APPEND | getattr(os, "O_NOFOLLOW", 0)
    index_fd = None
    try:
        index_fd = os.open(index.name, file_flags, dir_fd=directory_fd)
        index_metadata = os.fstat(index_fd)
        current_index_metadata = index.lstat()
        if (
            not stat.S_ISREG(index_metadata.st_mode)
            or index_metadata.st_nlink != 1
            or descriptor_identity(index_metadata) != expected_index_identity
            or not os.path.samestat(index_metadata, current_index_metadata)
        ):
            fail("packet index changed or has aliases before append")
        lines = read_descriptor_text(index_fd).splitlines() + ["\t".join(fields)]
        validate_packet_index_lines(lines)
        os.close(index_fd)
        index_fd = None
        atomic_replace_payload(
            directory_fd, index, ("\n".join(lines) + "\n").encode()
        )
        replacement_metadata = index.lstat()
        if (
            not stat.S_ISREG(replacement_metadata.st_mode)
            or replacement_metadata.st_nlink != 1
        ):
            fail("packet index replacement changed or has aliases")
        replacement_identity = descriptor_identity(replacement_metadata)
        validate_packet_index(index)
    finally:
        if index_fd is not None:
            os.close(index_fd)
        os.close(directory_fd)
    print(replacement_identity)


def open_archive_record(record: Path, expected_directory_identity):
    directory_fd = open_verified_directory(
        record.parent, expected_directory_identity, "archive-record"
    )
    record_fd = None
    try:
        reject_destination_temporaries(directory_fd, record)
        file_flags = os.O_RDWR | os.O_APPEND | os.O_CREAT
        file_flags |= getattr(os, "O_NOFOLLOW", 0)
        record_fd = os.open(record.name, file_flags, 0o644, dir_fd=directory_fd)
        record_metadata = os.fstat(record_fd)
        current_record_metadata = record.lstat()
        if (
            not stat.S_ISREG(record_metadata.st_mode)
            or record_metadata.st_nlink != 1
            or not os.path.samestat(record_metadata, current_record_metadata)
        ):
            fail("archive record changed or has aliases before append")
        os.fsync(directory_fd)
        return directory_fd, record_fd, record_metadata
    except BaseException:
        if record_fd is not None:
            os.close(record_fd)
        os.close(directory_fd)
        raise


def ensure_archive_record(record: Path, expected_directory_identity):
    directory_fd, record_fd, record_metadata = open_archive_record(
        record, expected_directory_identity
    )
    try:
        validate_archive_lines(read_descriptor_text(record_fd).splitlines())
        record_identity = descriptor_identity(record_metadata)
    finally:
        os.close(record_fd)
        os.close(directory_fd)
    print(record_identity)


def append_archive_record(
    record: Path,
    expected_directory_identity,
    expected_record_identity,
    expected_packet_identity,
    fields,
):
    directory_fd, record_fd, record_metadata = open_archive_record(
        record, expected_directory_identity
    )
    try:
        if descriptor_identity(record_metadata) != expected_record_identity:
            fail("archive record identity changed before append")
        existing_lines = read_descriptor_text(record_fd).splitlines()
        archive_line = "\t".join(fields)
        validate_archive_lines(existing_lines + [archive_line])
        packet_metadata = Path(fields[0]).lstat()
        if (
            not stat.S_ISDIR(packet_metadata.st_mode)
            or descriptor_identity(packet_metadata) != expected_packet_identity
        ):
            fail("packet directory identity changed before archive append")
        if archive_line not in existing_lines:
            lines = existing_lines + [archive_line]
            validate_archive_lines(lines)
            os.close(record_fd)
            record_fd = None
            atomic_replace_payload(
                directory_fd, record, ("\n".join(lines) + "\n").encode()
            )
            replacement_metadata = record.lstat()
            if (
                not stat.S_ISREG(replacement_metadata.st_mode)
                or replacement_metadata.st_nlink != 1
            ):
                fail("archive record replacement changed or has aliases")
            record_fd = os.open(
                record.name,
                os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=directory_fd,
            )
            validate_archive_lines(read_descriptor_text(record_fd).splitlines())
    finally:
        if record_fd is not None:
            os.close(record_fd)
        os.close(directory_fd)


def main():
    if len(sys.argv) == 2:
        validate_packet_index(Path(sys.argv[1]))
    elif len(sys.argv) == 3 and sys.argv[1] == "--require-rows":
        validate_packet_index(Path(sys.argv[2]), require_rows=True)
    elif len(sys.argv) == 3 and sys.argv[1] == "--archive-record":
        validate_archive_record(Path(sys.argv[2]))
    elif len(sys.argv) == 4 and sys.argv[1] == "--ensure-archive-record":
        ensure_archive_record(Path(sys.argv[2]), sys.argv[3])
    elif len(sys.argv) == 5 and sys.argv[1] == "--publish-file":
        publish_file(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4])
    elif len(sys.argv) == 5 and sys.argv[1] == "--reset-metadata-pair":
        reset_metadata_pair(Path(sys.argv[2]), Path(sys.argv[3]), sys.argv[4])
    elif len(sys.argv) == 4 and sys.argv[1] == "--sync-directory":
        sync_directory(Path(sys.argv[2]), sys.argv[3])
    elif len(sys.argv) == 15 and sys.argv[1] == "--append-row":
        append_packet_index(Path(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5:])
    elif len(sys.argv) == 11 and sys.argv[1] == "--append-archive-row":
        append_archive_record(
            Path(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:]
        )
    else:
        fail(
            "usage: validate_external_io_packet_index.py "
            "[--archive-record|--require-rows] INDEX_OR_ARCHIVE; or "
            "--ensure-archive-record RECORD EXPECTED_DIR_ID; or "
            "--publish-file SOURCE DESTINATION EXPECTED_DIR_ID; or "
            "--reset-metadata-pair MANIFEST INDEX EXPECTED_DIR_ID; or "
            "--sync-directory DIRECTORY EXPECTED_DIR_ID; or "
            "--append-row INDEX EXPECTED_DIR_ID EXPECTED_INDEX_ID FIELD...; or "
            "--append-archive-row RECORD EXPECTED_DIR_ID EXPECTED_RECORD_ID "
            "EXPECTED_PACKET_ID FIELD..."
        )


if __name__ == "__main__":
    main()
