"""Private shared limits and preflight helpers for AthenaK Python readers."""

from dataclasses import dataclass
from numbers import Integral
import argparse
import os
import sys


MIB = 1024 * 1024
DEFAULT_MAX_LIVE_BYTES = 512 * MIB
DEFAULT_MAX_HEADER_READ_BYTES = 16 * MIB
DEFAULT_MAX_PAYLOAD_READ_BYTES = 512 * MIB


def _positive_integer(value, label):
    if (
        isinstance(value, bool)
        or not isinstance(value, Integral)
        or value <= 0
    ):
        raise ValueError(f"{label} must be a positive integer")
    return int(value)


@dataclass(frozen=True)
class ReaderLimits:
    """Practical reader budgets. Defaults preserve the established limits."""

    max_live_bytes: int = DEFAULT_MAX_LIVE_BYTES
    max_header_read_bytes: int = DEFAULT_MAX_HEADER_READ_BYTES
    max_payload_read_bytes: int = DEFAULT_MAX_PAYLOAD_READ_BYTES

    def __post_init__(self):
        for name in (
            "max_live_bytes",
            "max_header_read_bytes",
            "max_payload_read_bytes",
        ):
            value = _positive_integer(getattr(self, name), name)
            object.__setattr__(self, name, value)


def normalize_reader_limits(
    limits,
    *,
    max_live_bytes=DEFAULT_MAX_LIVE_BYTES,
    max_header_read_bytes=DEFAULT_MAX_HEADER_READ_BYTES,
    max_payload_read_bytes=DEFAULT_MAX_PAYLOAD_READ_BYTES,
):
    """Return validated limits while allowing legacy defaults."""
    if limits is None:
        return ReaderLimits(
            max_live_bytes=max_live_bytes,
            max_header_read_bytes=max_header_read_bytes,
            max_payload_read_bytes=max_payload_read_bytes,
        )
    if not isinstance(limits, ReaderLimits):
        raise TypeError("limits must be a ReaderLimits instance or None")
    return limits


def checked_sum(values, label):
    """Add non-negative integer values without fixed-width overflow."""
    total = 0
    for value in values:
        if (
            isinstance(value, bool)
            or not isinstance(value, Integral)
            or value < 0
        ):
            raise ValueError(f"{label} has invalid extent {value!r}")
        total += int(value)
    return total


def checked_product(values, label):
    """Multiply non-negative integer values without fixed-width overflow."""
    product = 1
    for value in values:
        if (
            isinstance(value, bool)
            or not isinstance(value, Integral)
            or value < 0
        ):
            raise ValueError(f"{label} has invalid extent {value!r}")
        product *= int(value)
    return product


def require_live_bytes(nbytes, label, limits):
    """Reject an aggregate retained peak above the configured budget."""
    nbytes = checked_sum((nbytes,), label)
    if nbytes > limits.max_live_bytes:
        raise ValueError(
            f"{label} requires {nbytes} bytes, exceeding the practical "
            "allocation "
            f"limit of {limits.max_live_bytes} bytes"
        )
    return nbytes


def conservative_text_metadata_bytes(nbytes):
    """Conservatively charge decoded text and retained Python string storage."""
    nbytes = checked_sum((nbytes,), "text metadata")
    return checked_product((nbytes, 4), "expanded text metadata")


def preflight_allocation(values, itemsize, label, limits, retained_bytes=0):
    """Reject an allocation and its retained peak before materialization."""
    count = checked_product(values, label)
    nbytes = checked_product((count, itemsize), label + " byte count")
    peak_bytes = checked_sum((retained_bytes, nbytes), label)
    require_live_bytes(peak_bytes, label, limits)
    return count


def preflight_file_read(path, max_read_bytes, label):
    """Reject a file read above its configured practical budget."""
    max_read_bytes = _positive_integer(max_read_bytes, "max_read_bytes")
    size = os.path.getsize(path)
    if size > max_read_bytes:
        raise ValueError(
            f"{label} {path!r} requires reading {size} bytes, exceeding the "
            f"practical file-read limit of {max_read_bytes} bytes"
        )
    return size


def require_read_bytes(nbytes, max_read_bytes, label):
    """Reject one declared bulk read before passing it to read()."""
    nbytes = checked_sum((nbytes,), label)
    max_read_bytes = _positive_integer(max_read_bytes, "max_read_bytes")
    if nbytes > max_read_bytes:
        raise ValueError(
            f"{label} requires reading {nbytes} bytes, exceeding the "
            "practical "
            f"file-read limit of {max_read_bytes} bytes"
        )
    return nbytes


def bounded_partition_files(
    path,
    kind,
    label,
    limits,
    partition_info,
    externally_retained_bytes=0,
):
    """Discover sibling shards incrementally under a conservative live budget."""
    path = os.path.abspath(path)
    if kind == "shared":
        files = [path]
        retained_bytes = sys.getsizeof(files) + sys.getsizeof(path)
        require_live_bytes(
            externally_retained_bytes + retained_bytes,
            f"{label} sibling inventory",
            limits,
        )
        return files, retained_bytes
    shard_dir = os.path.dirname(path)
    parent = os.path.dirname(shard_dir)
    basename = os.path.basename(path)
    files = []
    retained_member_bytes = 0
    retained_bytes = (
        sys.getsizeof(files)
        + sys.getsizeof(parent)
        + sys.getsizeof(basename)
    )
    require_live_bytes(
        externally_retained_bytes + retained_bytes,
        f"{label} sibling inventory",
        limits,
    )
    with os.scandir(parent) as entries:
        for entry in entries:
            if not entry.name.startswith(kind + "_"):
                continue
            candidate = os.path.join(parent, entry.name, basename)
            if not os.path.isfile(candidate):
                continue
            candidate_kind, _ = partition_info(candidate)
            if candidate_kind != kind:
                raise ValueError(
                    f"{label} shard {candidate!r} does not match {kind!r} inventory"
                )
            item_bytes = sys.getsizeof(candidate)
            projected_bytes = retained_bytes + item_bytes + 2 * sys.getsizeof(None)
            require_live_bytes(
                externally_retained_bytes + projected_bytes,
                f"{label} sibling inventory",
                limits,
            )
            files.append(candidate)
            retained_member_bytes += item_bytes
            retained_bytes = (
                retained_member_bytes
                + sys.getsizeof(files)
                + sys.getsizeof(parent)
                + sys.getsizeof(basename)
            )
            require_live_bytes(
                externally_retained_bytes + retained_bytes,
                f"{label} sibling inventory",
                limits,
            )
    if not files:
        raise FileNotFoundError(f"no {label} {kind} shards found for {path!r}")
    require_live_bytes(
        externally_retained_bytes
        + retained_bytes
        + len(files) * sys.getsizeof(None),
        f"{label} sibling inventory sort peak",
        limits,
    )
    files.sort()
    for expected, candidate in enumerate(files):
        _, actual = partition_info(candidate)
        if actual != expected:
            raise ValueError(
                f"{label} {kind} shard inventory is incomplete: "
                f"expected ID {expected}, found {actual}"
            )
    return files, retained_bytes


def shallow_mapping_bytes(mapping):
    """Conservatively charge one retained metadata summary mapping."""
    return sys.getsizeof(mapping) + sum(
        sys.getsizeof(key) + sys.getsizeof(value)
        for key, value in mapping.items()
    )


def metadata_record_bytes(key, value=None, *, charge_value=True):
    """Conservatively charge one retained dictionary or set metadata record."""
    return (
        sys.getsizeof({})
        + 2 * sys.getsizeof(None)
        + sys.getsizeof(key)
        + (sys.getsizeof(value) if charge_value else 0)
    )


def count_ascii_tokens(text):
    """Count whitespace-delimited tokens without first materializing a list."""
    count = 0
    in_token = False
    for character in text:
        if character.isspace():
            in_token = False
        elif not in_token:
            count += 1
            in_token = True
    return count


def parse_ascii_floats(text, label):
    """Return strictly validated ASCII floating-point tokens."""
    if not text.isascii():
        raise ValueError(f"{label} must contain ASCII numeric tokens")
    try:
        return [float(token) for token in text.split()]
    except ValueError as exc:
        raise ValueError(f"{label} contains an invalid numeric token") from exc


def _positive_integer_argument(value):
    try:
        return _positive_integer(int(value), "reader limit")
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def add_reader_limit_arguments(parser):
    """Add consistent reader-budget flags to one ArgumentParser."""
    parser.add_argument(
        "--max-live-bytes",
        type=_positive_integer_argument,
        default=DEFAULT_MAX_LIVE_BYTES,
        help="maximum aggregate live/materialized reader bytes",
    )
    parser.add_argument(
        "--max-header-read-bytes",
        type=_positive_integer_argument,
        default=DEFAULT_MAX_HEADER_READ_BYTES,
        help="maximum bytes read for one metadata/header region",
    )
    parser.add_argument(
        "--max-payload-read-bytes",
        type=_positive_integer_argument,
        default=DEFAULT_MAX_PAYLOAD_READ_BYTES,
        help="maximum bytes read for one bulk payload",
    )


def reader_limits_from_args(args):
    """Build immutable reader limits from shared CLI arguments."""
    return ReaderLimits(
        max_live_bytes=args.max_live_bytes,
        max_header_read_bytes=args.max_header_read_bytes,
        max_payload_read_bytes=args.max_payload_read_bytes,
    )
