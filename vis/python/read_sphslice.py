#!/usr/bin/env python

"""Read AthenaK spherical-slice binary output.

Shared files contain a complete ``(variable, theta, phi)`` payload.  Rank-
and node-sharded files contain sparse angular ownership records and are
reassembled when any sibling shard path is supplied.
"""

import glob
import os
import re

import numpy as np


_FIRST_LINE_RE = re.compile(r"^Athena spherical slice version=(.+)$")
_SHARD_DIRECTORY_RE = re.compile(r"^(rank|node)_([0-9]+)$")
_SUPPORTED_VERSION = "1.0"
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
    for candidate in files:
        candidate_kind, _ = _partition_info(candidate)
        if candidate_kind != kind:
            raise ValueError(
                f"spherical-slice shard {candidate!r} does not match {kind!r} inventory"
            )
    return files


def _read_header(handle):
    first = handle.readline()
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
        raw_line = handle.readline()
        if not raw_line:
            raise ValueError(f"truncated spherical-slice header in {handle.name!r}")
        try:
            line = raw_line.decode("ascii").strip()
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"invalid spherical-slice header text in {handle.name!r}"
            ) from exc
        if line.startswith("variables:"):
            header["variables"] = line[len("variables:"):].split()
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


def _read_payload(handle, header, distribution):
    dump = handle.read(header["header_offset"])
    if len(dump) != header["header_offset"]:
        raise ValueError(
            f"truncated spherical-slice input-header block in {handle.name!r}"
        )
    npoints = header["npoints"]
    surface_points = header["ntheta"] * header["nphi"]
    nvars = header["nvars"]
    if distribution == "shared":
        if npoints != surface_points:
            raise ValueError(
                f"shared spherical-slice file {handle.name!r} has npoints={npoints}, "
                f"expected {surface_points}"
            )
        expected = nvars * surface_points * np.dtype(np.float32).itemsize
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
    expected = npoints * (
        np.dtype(np.int32).itemsize + nvars * np.dtype(np.float32).itemsize
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
    if npoints and (np.any(indices < 0) or np.any(indices >= surface_points)):
        bad = int(indices[(indices < 0) | (indices >= surface_points)][0])
        raise ValueError(
            f"spherical-slice shard {handle.name!r} has out-of-range angle index {bad}"
        )
    if npoints and np.unique(indices).size != npoints:
        raise ValueError(
            f"spherical-slice shard {handle.name!r} contains duplicate angle ownership"
        )
    return indices, values


def _read_one(path):
    with open(path, "rb") as handle:
        header = _read_header(handle)
        distribution = _distribution_for(header, path)
        header["distribution"] = distribution
        indices, values = _read_payload(handle, header, distribution)
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
    with open(path, "rb") as handle:
        header = _read_header(handle)
    header["distribution"] = _distribution_for(header, path)
    return header


def read_sphslice(path):
    """Read a shared or rank/node-sharded spherical slice as a dense surface.

    The returned ``data`` array is ordered ``(theta, phi, variable)``.
    Partitioned input must provide exactly one owner for every angular point;
    zero-point shard files are valid and are retained during validation.
    """
    files = _glob_partition_files(path)
    records = []
    for shard in files:
        candidate, indices, values = _read_one(shard)
        records.append((shard, candidate, indices, values))

    header = None
    full = None
    covered = None

    for shard, candidate, indices, values in records:
        if header is None:
            header = candidate
            surface_points = header["ntheta"] * header["nphi"]
            full = np.zeros((header["nvars"], surface_points), dtype=np.float32)
            covered = np.zeros(surface_points, dtype=bool)
        else:
            _compare_headers(header, candidate, shard)
        if indices is None:
            full[...] = values.reshape(header["nvars"], -1)
            covered[...] = True
            continue
        if indices.size:
            duplicate = indices[covered[indices]]
            if duplicate.size:
                raise ValueError(
                    f"spherical-slice duplicate ownership for angle index "
                    f"{int(duplicate[0])} in {shard!r}"
                )
            full[:, indices] = values.reshape(header["nvars"], indices.size)
            covered[indices] = True

    if header is None:
        raise ValueError(f"no spherical-slice data found for {path!r}")
    _validate_sibling_inventory(files, [record[1] for record in records])
    missing = np.flatnonzero(~covered)
    if missing.size:
        raise ValueError(
            f"spherical-slice reassembly from {path!r} is missing "
            f"{missing.size} of {covered.size} angular points "
            f"(first missing index {int(missing[0])})"
        )

    ntheta, nphi = header["ntheta"], header["nphi"]
    data = full.reshape(header["nvars"], ntheta, nphi).transpose(1, 2, 0)
    theta = np.arccos(-1.0 + 2.0 * (np.arange(ntheta) + 0.5) / ntheta)
    phi = 2.0 * np.pi * (np.arange(nphi) + 0.5) / nphi
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
