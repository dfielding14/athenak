#!/usr/bin/env python

"""Read AthenaK spherical-slice binary output.

Shared files contain a complete ``(variable, theta, phi)`` payload.  Rank-
and node-sharded files contain sparse angular ownership records and are
reassembled when any sibling shard path is supplied.
"""

import glob
import os
import re
import sys

import numpy as np


_FIRST_LINE_RE = re.compile(r"^Athena spherical slice version=(.+)$")
_SHARD_DIRECTORY_RE = re.compile(r"^(rank|node)_([0-9]{8})$")
_SUPPORTED_VERSION = "1.0"
_MAX_DENSE_ALLOCATION_BYTES = 512 * 1024 * 1024
_MAX_FILE_READ_BYTES = 512 * 1024 * 1024
_MAX_HEADER_READ_BYTES = 16 * 1024 * 1024
_INT_KEYS = frozenset(
    (
        "cycle",
        "ntheta",
        "nphi",
        "number_of_variables",
        "size_of_variable",
        "npoints",
        "rank",
        "node",
        "number_of_nodes",
        "number_of_ranks",
        "single_file_per_rank",
        "header_offset",
    )
)
_FLOAT_KEYS = frozenset(("time", "radius"))
_REQUIRED_KEYS = frozenset(
    (
        "time",
        "cycle",
        "radius",
        "ntheta",
        "nphi",
        "size_of_variable",
        "number_of_variables",
        "npoints",
        "header_offset",
        "variables",
    )
)


def _checked_product(values, label):
    product = 1
    for value in values:
        if not isinstance(value, (int, np.integer)) or value < 0:
            raise ValueError(f"{label} has invalid extent {value!r}")
        product *= int(value)
    return product


def _require_allocation(values, itemsize, label):
    count = _checked_product(values, label)
    nbytes = _checked_product((count, itemsize), label + " byte count")
    if nbytes > _MAX_DENSE_ALLOCATION_BYTES:
        raise ValueError(
            f"{label} requires {nbytes} bytes, exceeding the practical allocation "
            f"limit of {_MAX_DENSE_ALLOCATION_BYTES} bytes"
        )
    return count


def _require_retained_bytes(nbytes, label):
    if nbytes > _MAX_DENSE_ALLOCATION_BYTES:
        raise ValueError(
            f"{label} requires {nbytes} bytes, exceeding the practical allocation "
            f"limit of {_MAX_DENSE_ALLOCATION_BYTES} bytes"
        )


def _require_file_size(path):
    size = os.path.getsize(path)
    if size > _MAX_FILE_READ_BYTES:
        raise ValueError(
            f"spherical-slice file {path!r} requires reading {size} bytes, "
            f"exceeding the practical file-read limit of {_MAX_FILE_READ_BYTES} bytes"
        )


def _read_limited_header_line(handle, budget):
    """Read one metadata line without accepting an unbounded header."""
    raw_line = handle.readline(_MAX_HEADER_READ_BYTES + 1)
    budget[0] += len(raw_line)
    if len(raw_line) > _MAX_HEADER_READ_BYTES or budget[0] > _MAX_HEADER_READ_BYTES:
        raise ValueError(
            f"spherical-slice header in {handle.name!r} exceeds the practical "
            f"metadata limit of {_MAX_HEADER_READ_BYTES} bytes"
        )
    return raw_line


def _count_ascii_tokens(text):
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


def _variable_token_peak_bytes(text, token_count):
    """Conservatively bound split variable strings and their pointer list."""
    pointer_bytes = np.dtype(np.intp).itemsize
    return (
        sys.getsizeof([])
        + 2 * token_count * pointer_bytes
        + token_count * sys.getsizeof("")
        + 4 * len(text)
    )


def _partition_info(path):
    directory = os.path.basename(os.path.dirname(os.path.abspath(path)))
    match = _SHARD_DIRECTORY_RE.fullmatch(directory)
    if match is not None:
        return match.group(1), int(match.group(2))
    if directory.startswith(("rank_", "node_")):
        raise ValueError(
            f"invalid spherical-slice shard directory {directory!r} for {path!r}"
        )
    return "shared", None


def _shard_kind(path):
    return _partition_info(path)[0]


def _glob_partition_files(path):
    kind = _shard_kind(path)
    if kind == "shared":
        return [os.path.abspath(path)]
    directory = os.path.dirname(os.path.abspath(path))
    pattern = os.path.join(
        os.path.dirname(directory), kind + "_*", os.path.basename(path)
    )
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(
            f"no spherical-slice {kind} shards found for pattern {pattern!r}"
        )
    shard_ids = []
    for candidate in files:
        candidate_kind, shard_id = _partition_info(candidate)
        if candidate_kind != kind:
            raise ValueError(
                f"spherical-slice shard {candidate!r} does not match {kind!r} inventory"
            )
        shard_ids.append(shard_id)
    if len(set(shard_ids)) != len(shard_ids):
        raise ValueError(
            f"spherical-slice {kind} shard inventory contains duplicate IDs"
        )
    expected_ids = set(range(len(shard_ids)))
    actual_ids = set(shard_ids)
    if actual_ids != expected_ids:
        raise ValueError(
            f"spherical-slice {kind} shard inventory is incomplete: "
            f"expected IDs {sorted(expected_ids)!r}, found {sorted(actual_ids)!r}"
        )
    return files


def _read_header(handle, externally_retained_bytes=0):
    header_budget = [0]
    first = _read_limited_header_line(handle, header_budget)
    if not first:
        raise ValueError(f"empty spherical-slice file {handle.name!r}")
    try:
        first_line = first.decode("ascii").rstrip()
    except UnicodeDecodeError as exc:
        raise ValueError(f"invalid spherical-slice header in {handle.name!r}") from exc
    match = _FIRST_LINE_RE.match(first_line)
    if match is None:
        raise ValueError(f"not a spherical-slice file: {handle.name!r}")
    if match.group(1) != _SUPPORTED_VERSION:
        raise ValueError(
            f"unsupported spherical-slice version {match.group(1)!r} in {handle.name!r}"
        )

    header = {"version": match.group(1)}
    while "header_offset" not in header:
        raw_line = _read_limited_header_line(handle, header_budget)
        if not raw_line:
            raise ValueError(f"truncated spherical-slice header in {handle.name!r}")
        try:
            line = raw_line.decode("ascii").strip()
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"invalid spherical-slice header text in {handle.name!r}"
            ) from exc
        if line.startswith("variables:"):
            if "variables" in header:
                raise ValueError(
                    f"duplicate spherical-slice variables metadata in {handle.name!r}"
                )
            variables_text = line[len("variables:"):]
            token_count = _count_ascii_tokens(variables_text)
            retained_variable_bytes = _variable_token_peak_bytes(
                variables_text, token_count
            )
            _require_retained_bytes(
                externally_retained_bytes + retained_variable_bytes,
                "spherical-slice variable tokenization peak",
            )
            header["variables"] = variables_text.split()
            header["_retained_metadata_bytes"] = retained_variable_bytes
            continue
        if "=" not in line:
            continue
        key, value = (part.strip() for part in line.partition("=")[::2])
        key = key.replace(" ", "_")
        try:
            if key in _INT_KEYS:
                header[key] = int(value)
            elif key in _FLOAT_KEYS:
                header[key] = float(value)
            else:
                header[key] = value
        except ValueError as exc:
            raise ValueError(
                f"invalid spherical-slice header value for {key!r} in {handle.name!r}"
            ) from exc

    missing = sorted(_REQUIRED_KEYS.difference(header))
    if missing:
        raise ValueError(
            f"spherical-slice header {handle.name!r} is missing {', '.join(missing)}"
        )
    header["nvars"] = header["number_of_variables"]
    if header["ntheta"] <= 0 or header["nphi"] <= 0 or header["nvars"] <= 0:
        raise ValueError(f"spherical-slice header {handle.name!r} has invalid dimensions")
    header["surface_points"] = _require_allocation(
        (header["ntheta"], header["nphi"]),
        np.dtype(bool).itemsize,
        "spherical-slice surface coverage",
    )
    dense_values = _require_allocation(
        (header["nvars"], header["surface_points"]),
        np.dtype(np.float32).itemsize,
        "spherical-slice dense values",
    )
    _require_allocation(
        (header["ntheta"],),
        np.dtype(np.float64).itemsize,
        "spherical-slice theta coordinates",
    )
    _require_allocation(
        (header["nphi"],),
        np.dtype(np.float64).itemsize,
        "spherical-slice phi coordinates",
    )
    header["_retained_dense_bytes"] = (
        header["surface_points"] * np.dtype(bool).itemsize
        + dense_values * np.dtype(np.float32).itemsize
        + (header["ntheta"] + header["nphi"]) * np.dtype(np.float64).itemsize
    )
    _require_retained_bytes(
        header["_retained_dense_bytes"], "spherical-slice retained arrays"
    )
    if header["size_of_variable"] != np.dtype(np.float32).itemsize:
        raise ValueError(
            f"spherical-slice file {handle.name!r} has unsupported variable size "
            f"{header['size_of_variable']}"
        )
    if len(header["variables"]) != header["nvars"]:
        raise ValueError(
            f"spherical-slice file {handle.name!r} variable count does not match metadata"
        )
    if header["header_offset"] < 0 or header["npoints"] < 0:
        raise ValueError(
            f"spherical-slice file {handle.name!r} has a negative size field"
        )
    if not np.isfinite(header["time"]) or not np.isfinite(header["radius"]):
        raise ValueError(f"spherical-slice file {handle.name!r} has non-finite metadata")
    return header


def _distribution_for(header, path):
    distribution = header.get("distribution")
    if distribution is None:
        distribution = "rank" if header.get("single_file_per_rank", 0) else "shared"
    if distribution not in ("shared", "rank", "node"):
        raise ValueError(
            f"spherical-slice file {path!r} has invalid distribution {distribution!r}"
        )
    kind, shard_id = _partition_info(path)
    if kind != distribution:
        raise ValueError(
            f"spherical-slice file {path!r} declares distribution={distribution!r} "
            f"but resides in a {kind!r} layout"
        )
    expected_layout = "dense" if distribution == "shared" else "sparse_angles"
    layout = header.get("layout", expected_layout)
    if layout != expected_layout:
        raise ValueError(
            f"spherical-slice file {path!r} declares layout={layout!r}, "
            f"expected {expected_layout!r} for distribution={distribution!r}"
        )
    header["layout"] = layout
    if distribution == "rank" and "rank" in header and header["rank"] != shard_id:
        raise ValueError(
            f"spherical-slice file {path!r} declares rank={header['rank']}, "
            f"but resides in rank_{shard_id:08d}"
        )
    if distribution == "node" and "node" in header and header["node"] != shard_id:
        raise ValueError(
            f"spherical-slice file {path!r} declares node={header['node']}, "
            f"but resides in node_{shard_id:08d}"
        )
    for count_key, expected_distribution in (
        ("number_of_ranks", "rank"),
        ("number_of_nodes", "node"),
    ):
        if count_key in header:
            if header[count_key] <= 0:
                raise ValueError(
                    f"spherical-slice file {path!r} has invalid {count_key}="
                    f"{header[count_key]}"
                )
            if distribution != expected_distribution:
                raise ValueError(
                    f"spherical-slice file {path!r} defines {count_key} for "
                    f"distribution={distribution!r}"
                )
    return distribution


def _validate_header_offset(handle, header):
    """Reject impossible embedded-input offsets before passing them to read()."""
    remaining = os.fstat(handle.fileno()).st_size - handle.tell()
    if header["header_offset"] > remaining:
        raise ValueError(
            f"truncated spherical-slice input-header block in {handle.name!r}"
        )


def _read_payload(handle, header, distribution, externally_retained_bytes=0):
    _validate_header_offset(handle, header)
    dump_bytes = header["header_offset"]
    npoints = header["npoints"]
    surface_points = header["surface_points"]
    nvars = header["nvars"]
    retained_header_bytes = (
        externally_retained_bytes
        + header["_retained_dense_bytes"]
        + header["_retained_metadata_bytes"]
    )
    if distribution == "shared":
        if npoints != surface_points:
            raise ValueError(
                f"shared spherical-slice file {handle.name!r} has npoints={npoints}, "
                f"expected {surface_points}"
            )
        expected = _checked_product(
            (nvars, surface_points, np.dtype(np.float32).itemsize),
            "spherical-slice dense payload",
        )
        _require_retained_bytes(
            retained_header_bytes + dump_bytes + 2 * expected,
            "spherical-slice dense reconstruction peak",
        )
        dump = handle.read(dump_bytes)
        if len(dump) != dump_bytes:
            raise ValueError(
                f"truncated spherical-slice input-header block in {handle.name!r}"
            )
        payload = handle.read()
        if len(payload) != expected:
            raise ValueError(
                f"truncated or overlong spherical-slice payload in {handle.name!r}: "
                f"expected {expected} bytes, found {len(payload)}"
            )
        values = np.frombuffer(payload, dtype=np.float32).copy()
        return None, values

    if npoints > surface_points:
        raise ValueError(
            f"spherical-slice shard {handle.name!r} has npoints={npoints}, "
            f"larger than surface size {surface_points}"
        )
    record_bytes = np.dtype(np.int32).itemsize + _checked_product(
        (nvars, np.dtype(np.float32).itemsize),
        "spherical-slice sparse value record",
    )
    expected = _checked_product(
        (npoints, record_bytes),
        "spherical-slice sparse payload",
    )
    index_bytes = _checked_product(
        (npoints, np.dtype(np.int32).itemsize),
        "spherical-slice sparse index copy",
    )
    _require_retained_bytes(
        retained_header_bytes + dump_bytes + 2 * expected + index_bytes,
        "spherical-slice sparse reconstruction peak",
    )
    value_bytes = _checked_product(
        (nvars, npoints, np.dtype(np.float32).itemsize),
        "spherical-slice sparse value copy",
    )
    unique_bytes = 2 * index_bytes + npoints * np.dtype(bool).itemsize
    _require_retained_bytes(
        retained_header_bytes
        + dump_bytes
        + expected
        + index_bytes
        + value_bytes
        + unique_bytes,
        "spherical-slice sparse duplicate-validation peak",
    )
    dump = handle.read(dump_bytes)
    if len(dump) != dump_bytes:
        raise ValueError(
            f"truncated spherical-slice input-header block in {handle.name!r}"
        )
    payload = handle.read()
    if len(payload) != expected:
        raise ValueError(
            f"truncated or overlong spherical-slice shard {handle.name!r}: "
            f"expected {expected} bytes, found {len(payload)}"
        )
    indices = np.frombuffer(payload, dtype=np.int32, count=npoints).copy()
    values = np.frombuffer(
        payload, dtype=np.float32, count=nvars * npoints, offset=4 * npoints
    ).copy()
    if npoints and (np.min(indices) < 0 or np.max(indices) >= surface_points):
        bad = int(np.min(indices) if np.min(indices) < 0 else np.max(indices))
        raise ValueError(
            f"spherical-slice shard {handle.name!r} has out-of-range angle index {bad}"
        )
    if npoints and np.unique(indices).size != npoints:
        raise ValueError(
            f"spherical-slice shard {handle.name!r} contains duplicate angle ownership"
        )
    return indices, values


def _read_one(path, externally_retained_bytes=0):
    _require_file_size(path)
    with open(path, "rb") as handle:
        header = _read_header(handle, externally_retained_bytes)
        _validate_header_offset(handle, header)
        distribution = _distribution_for(header, path)
        header["distribution"] = distribution
        indices, values = _read_payload(
            handle, header, distribution, externally_retained_bytes
        )
    return header, indices, values


def _compare_headers(reference, candidate, path):
    for key in (
        "version",
        "layout",
        "distribution",
        "time",
        "cycle",
        "radius",
        "ntheta",
        "nphi",
        "size_of_variable",
        "nvars",
        "variables",
    ):
        if candidate[key] != reference[key]:
            raise ValueError(
                f"spherical-slice shard metadata mismatch for {key!r} in {path!r}: "
                f"{candidate[key]!r} != {reference[key]!r}"
            )


def _validate_sibling_inventory(files, headers):
    distribution = headers[0]["distribution"]
    if distribution == "shared":
        return
    count_key = "number_of_ranks" if distribution == "rank" else "number_of_nodes"
    id_key = "rank" if distribution == "rank" else "node"
    count_values = [header.get(count_key) for header in headers]
    if all(value is None for value in count_values):
        return
    if any(value is None for value in count_values) or len(set(count_values)) != 1:
        raise ValueError(
            f"spherical-slice {distribution} shard inventory has inconsistent "
            f"{count_key} metadata"
        )
    expected_count = count_values[0]
    sibling_ids = []
    for path, header in zip(files, headers):
        kind, shard_id = _partition_info(path)
        if kind != distribution:
            raise ValueError(
                f"spherical-slice shard {path!r} does not match "
                f"{distribution!r} inventory"
            )
        if id_key not in header:
            raise ValueError(
                f"spherical-slice shard {path!r} is missing required {id_key} metadata"
            )
        if header[id_key] != shard_id:
            raise ValueError(
                f"spherical-slice shard {path!r} declares {id_key}={header[id_key]}, "
                f"but its directory identifies {shard_id}"
            )
        sibling_ids.append(shard_id)
    if len(set(sibling_ids)) != len(sibling_ids):
        raise ValueError(
            f"spherical-slice {distribution} shard inventory contains duplicate IDs"
        )
    if expected_count != len(sibling_ids):
        raise ValueError(
            f"spherical-slice {distribution} shard inventory is incomplete: "
            f"expected {expected_count} shards, found {len(sibling_ids)}"
        )
    expected_ids = set(range(len(sibling_ids)))
    actual_ids = set(sibling_ids)
    if actual_ids != expected_ids:
        raise ValueError(
            f"spherical-slice {distribution} shard inventory is incomplete: "
            f"expected IDs {sorted(expected_ids)!r}, found {sorted(actual_ids)!r}"
        )


def read_sphslice_header(path):
    """Read and validate only a spherical-slice preheader."""
    _require_file_size(path)
    with open(path, "rb") as handle:
        header = _read_header(handle)
        _validate_header_offset(handle, header)
    header["distribution"] = _distribution_for(header, path)
    header.pop("_retained_dense_bytes", None)
    header.pop("_retained_metadata_bytes", None)
    return header


def read_sphslice(path):
    """Read a shared or rank/node-sharded spherical slice as a dense surface.

    The returned ``data`` array is ordered ``(theta, phi, variable)``.
    Partitioned input must provide exactly one owner for every angular point;
    zero-point shard files are valid and are retained during validation.
    """
    files = _glob_partition_files(path)
    header = None
    full = None
    covered = None
    inventory_headers = []
    for shard in files:
        externally_retained_bytes = 0
        if header is not None:
            externally_retained_bytes = (
                header["_retained_dense_bytes"] + header["_retained_metadata_bytes"]
            )
        candidate, indices, values = _read_one(shard, externally_retained_bytes)
        inventory_headers.append(
            {
                key: candidate[key]
                for key in (
                    "distribution",
                    "number_of_ranks",
                    "number_of_nodes",
                    "rank",
                    "node",
                )
                if key in candidate
            }
        )
        if header is None:
            header = candidate
            surface_points = header["surface_points"]
            full = np.zeros((header["nvars"], surface_points), dtype=np.float32)
            covered = np.zeros(surface_points, dtype=bool)
        else:
            _compare_headers(header, candidate, shard)
        if indices is None:
            full[...] = values.reshape(header["nvars"], -1)
            covered[...] = True
            del indices, values
            if candidate is not header:
                del candidate
            continue
        if indices.size:
            retained_metadata_bytes = header["_retained_metadata_bytes"]
            if candidate is not header:
                retained_metadata_bytes += candidate["_retained_metadata_bytes"]
            _require_retained_bytes(
                header["_retained_dense_bytes"]
                + retained_metadata_bytes
                + indices.nbytes
                + values.nbytes
                + indices.size * np.dtype(bool).itemsize
                + indices.nbytes,
                "spherical-slice duplicate ownership check peak",
            )
            duplicate = indices[covered[indices]]
            if duplicate.size:
                raise ValueError(
                    f"spherical-slice duplicate ownership for angle index "
                    f"{int(duplicate[0])} in {shard!r}"
                )
            full[:, indices] = values.reshape(header["nvars"], indices.size)
            covered[indices] = True
        del indices, values
        if candidate is not header:
            del candidate

    if header is None:
        raise ValueError(f"no spherical-slice data found for {path!r}")
    _validate_sibling_inventory(files, inventory_headers)
    missing_count = covered.size - int(np.count_nonzero(covered))
    if missing_count:
        raise ValueError(
            f"spherical-slice reassembly from {path!r} is missing "
            f"{missing_count} of {covered.size} angular points "
            f"(first missing index {int(np.argmin(covered))})"
        )

    ntheta, nphi = header["ntheta"], header["nphi"]
    retained_dense_bytes = header["_retained_dense_bytes"]
    retained_metadata_bytes = header["_retained_metadata_bytes"]
    header.pop("_retained_dense_bytes", None)
    header.pop("rank", None)
    header.pop("node", None)
    header["npoints"] = header["surface_points"]
    data = full.reshape(header["nvars"], ntheta, nphi).transpose(1, 2, 0)
    coordinate_temporary_bytes = (
        4 * (ntheta + nphi) * np.dtype(np.float64).itemsize
    )
    _require_retained_bytes(
        retained_dense_bytes
        + retained_metadata_bytes
        + coordinate_temporary_bytes,
        "spherical-slice coordinate generation peak",
    )
    theta = np.arccos(-1.0 + 2.0 * (np.arange(ntheta) + 0.5) / ntheta)
    phi = 2.0 * np.pi * (np.arange(nphi) + 0.5) / nphi
    header.pop("_retained_metadata_bytes", None)
    return {
        "header": header,
        "data": data,
        "theta": theta,
        "phi": phi,
        "time": header["time"],
        "cycle": header["cycle"],
        "radius": header["radius"],
        "variables": header["variables"],
    }


__all__ = ["read_sphslice", "read_sphslice_header"]
