#!/usr/bin/env python

"""Read AthenaK probability-distribution-function output.

Supported inputs are the legacy shared ASCII ``.bins.pdf``/``.pdf`` pair and
the N-dimensional binary representation with a ``.header.pdf`` metadata file.
The binary format stores shared data densely and rank/node shards as sparse
COO records.  Passing any sparse shard discovers and sums its sibling shards.
"""

import glob
import os
import re
import struct

import numpy as np


_HEADER_KEY_RE = re.compile(r"^([A-Za-z0-9_]+)\s*=\s*(.*)$")
_BIN_EDGE_RE = re.compile(r"^bin_edges_(\d+)$")
_PDF_DATA_RE = re.compile(r"^(.*)\.(\d{5})\.pdf$")
_LEGACY_VARIABLE_RE = re.compile(r"^#\s*\[(\d+)\]\s*=\s*(\S+)")
_LEGACY_TIME_RE = re.compile(r"^#\s*time\s*=\s*(\S+)")
_SCALES = frozenset(("linear", "log", "symlog"))
_V2_MAGIC = b"AKPDFV2\0"
_V2_PREAMBLE = struct.Struct("=8sIIIIQdq")


def _parse_bool(value):
    normalized = value.strip().lower()
    if normalized in ("true", "t", "1", "yes", "y"):
        return True
    if normalized in ("false", "f", "0", "no", "n"):
        return False
    raise ValueError(f"PDF header has invalid boolean value {value!r}")


def _shard_kind(path):
    leaf = os.path.basename(os.path.dirname(os.path.abspath(path)))
    if leaf.startswith("rank_"):
        return "rank"
    if leaf.startswith("node_"):
        return "node"
    return "shared"


def _glob_partition_files(path):
    kind = _shard_kind(path)
    if kind == "shared":
        return [os.path.abspath(path)]
    shard_dir = os.path.dirname(os.path.abspath(path))
    parent = os.path.dirname(shard_dir)
    pattern = os.path.join(parent, kind + "_*", os.path.basename(path))
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"no PDF {kind} shards found for pattern {pattern!r}")
    return files


def _header_candidates(data_path):
    match = _PDF_DATA_RE.match(os.path.basename(data_path))
    if match is None:
        raise ValueError(f"unable to infer PDF header path from {data_path!r}")
    base = match.group(1)
    directory = os.path.dirname(os.path.abspath(data_path))
    directories = [directory]
    if _shard_kind(data_path) != "shared":
        directories.append(os.path.dirname(directory))
    for extension in (".header.pdf", ".bins.pdf"):
        for candidate_dir in directories:
            yield os.path.join(candidate_dir, base + extension)


def _infer_header_path(data_path):
    for candidate in _header_candidates(data_path):
        if os.path.exists(candidate):
            return candidate
    candidates = ", ".join(_header_candidates(data_path))
    raise FileNotFoundError(
        f"no PDF metadata file found for {data_path!r}; checked {candidates}"
    )


def _symlog_forward(values, linthresh):
    values = np.asarray(values, dtype=np.float64)
    absolute = np.abs(values)
    transformed = np.empty_like(absolute)
    linear = absolute <= linthresh
    transformed[linear] = absolute[linear] / linthresh
    transformed[~linear] = 1.0 + np.log10(absolute[~linear] / linthresh)
    return np.sign(values) * transformed


def _symlog_inverse(values, linthresh):
    values = np.asarray(values, dtype=np.float64)
    absolute = np.abs(values)
    physical = np.empty_like(absolute)
    linear = absolute <= 1.0
    physical[linear] = absolute[linear] * linthresh
    physical[~linear] = linthresh * np.power(10.0, absolute[~linear] - 1.0)
    return np.sign(values) * physical


def _generated_edges(info):
    transformed = np.linspace(info["bin_min"], info["bin_max"], info["nbin"] + 1)
    if info["scale"] == "log":
        transformed = np.linspace(
            np.log10(info["bin_min"]), np.log10(info["bin_max"]), info["nbin"] + 1
        )
        return np.power(10.0, transformed)
    if info["scale"] == "symlog":
        transformed = np.linspace(
            _symlog_forward(info["bin_min"], info["linthresh"]),
            _symlog_forward(info["bin_max"], info["linthresh"]),
            info["nbin"] + 1,
        )
        return _symlog_inverse(transformed, info["linthresh"])
    return transformed


def _bin_centers(info):
    edges = info["bin_edges"]
    if info["scale"] == "log":
        return np.sqrt(edges[:-1] * edges[1:])
    if info["scale"] == "symlog":
        transformed = _symlog_forward(edges, info["linthresh"])
        return _symlog_inverse(
            0.5 * (transformed[:-1] + transformed[1:]), info["linthresh"]
        )
    return 0.5 * (edges[:-1] + edges[1:])


def _read_legacy_header(header_path):
    dimensions = []
    variable_names = {}
    with open(header_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = _LEGACY_VARIABLE_RE.match(line)
            if match is not None:
                variable_names[int(match.group(1))] = match.group(2)
                continue
            if not line or line.startswith("#"):
                continue
            edges = np.fromstring(line, dtype=np.float64, sep=" ")
            if edges.size < 2:
                raise ValueError(f"legacy PDF header {header_path!r} has invalid bin row")
            if not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
                raise ValueError(
                    f"legacy PDF header {header_path!r} has non-increasing bin edges"
                )
            dimensions.append(
                {
                    "variable": variable_names.get(len(dimensions) + 1),
                    "nbin": int(edges.size - 1),
                    "nbin_with_overflow": int(edges.size + 1),
                    "bin_edges": edges,
                    "bin_centers": 0.5 * (edges[:-1] + edges[1:]),
                    "scale": None,
                }
            )
    if not dimensions or len(dimensions) > 2:
        raise ValueError(
            f"legacy PDF header {header_path!r} must contain one or two bin rows"
        )
    shape = tuple(item["nbin_with_overflow"] for item in dimensions)
    return {
        "format": "legacy_dense",
        "distribution": "shared",
        "ndim": len(dimensions),
        "dimensions": dimensions,
        "shape": shape,
        "total_bins": int(np.prod(shape)),
        "header_path": os.path.abspath(header_path),
    }


def read_pdf_header(header_path):
    """Parse a modern ``.header.pdf`` or legacy ``.bins.pdf`` metadata file."""
    if header_path.endswith(".bins.pdf"):
        return _read_legacy_header(header_path)

    header = {}
    dims = {}
    with open(header_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            match = _HEADER_KEY_RE.match(line)
            if match is None:
                raise ValueError(f"malformed PDF header line in {header_path!r}: {line!r}")
            key, value = match.group(1), match.group(2).strip()
            if key in ("format", "distribution", "weight", "weight_variable",
                       "binary_magic"):
                header[key] = value
            elif key == "layout":
                header["format"] = value
            elif key in ("ndim", "total_bins", "cycle"):
                header[key] = int(value)
            else:
                edge_match = _BIN_EDGE_RE.match(key)
                if edge_match is not None:
                    dimension = int(edge_match.group(1))
                    dims.setdefault(dimension, {})["bin_edges"] = np.array(
                        [float(entry) for entry in value.split()], dtype=np.float64
                    )
                    continue
                bound_match = re.match(r"^bin(\d+)_(min|max)$", key)
                if bound_match is not None:
                    dimension = int(bound_match.group(1))
                    dims.setdefault(dimension, {})[
                        "bin_" + bound_match.group(2)
                    ] = float(value)
                    continue
                suffix = re.search(r"(\d+)$", key)
                if suffix is None:
                    continue
                dimension = int(suffix.group(1))
                info = dims.setdefault(dimension, {})
                if key.startswith("variable_"):
                    info["variable"] = value
                elif key.startswith("nbin"):
                    info["nbin"] = int(value)
                elif key.startswith("scale"):
                    info["scale"] = value
                elif key.startswith("linthresh"):
                    info["linthresh"] = float(value)
                elif key.startswith("logscale"):
                    info["logscale"] = _parse_bool(value)
                elif key.startswith("stride"):
                    info["stride"] = int(value)

    if header.get("format", "dense") not in ("dense", "sparse_coo"):
        raise ValueError(f"unsupported PDF format in {header_path!r}: {header.get('format')!r}")
    header.setdefault("format", "dense")
    if "binary_magic" in header and header["binary_magic"] != "AKPDFV2":
        raise ValueError(
            f"unsupported PDF binary magic in {header_path!r}: "
            f"{header['binary_magic']!r}"
        )
    if "ndim" not in header or not 1 <= header["ndim"] <= 4:
        raise ValueError(f"PDF header {header_path!r} must declare ndim between 1 and 4")

    dimensions = []
    for dimension in range(1, header["ndim"] + 1):
        if dimension not in dims:
            raise ValueError(f"PDF header {header_path!r} is missing dimension {dimension}")
        info = dims[dimension]
        if not isinstance(info.get("nbin"), int) or info["nbin"] <= 0:
            raise ValueError(f"PDF header {header_path!r} has invalid nbin{dimension}")
        if "variable" not in info:
            raise ValueError(f"PDF header {header_path!r} is missing variable_{dimension}")
        scale = info.get("scale", "log" if info.get("logscale", False) else "linear")
        if scale not in _SCALES:
            raise ValueError(
                f"PDF header {header_path!r} has invalid scale{dimension}={scale!r}"
            )
        info["scale"] = scale
        if scale == "symlog":
            if info.get("linthresh", 0.0) <= 0.0:
                raise ValueError(
                    f"PDF header {header_path!r} has invalid linthresh{dimension}"
                )
        else:
            info.setdefault("linthresh", 1.0)
        if "bin_min" not in info or "bin_max" not in info:
            raise ValueError(f"PDF header {header_path!r} is missing bin bounds")
        if info["bin_min"] >= info["bin_max"]:
            raise ValueError(f"PDF header {header_path!r} has reversed bin bounds")
        if scale == "log" and info["bin_min"] <= 0.0:
            raise ValueError(f"PDF header {header_path!r} has non-positive log bounds")
        info.setdefault("bin_edges", _generated_edges(info))
        edges = info["bin_edges"]
        if edges.size != info["nbin"] + 1:
            raise ValueError(
                f"PDF header {header_path!r} has wrong number of bin_edges_{dimension}"
            )
        if not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
            raise ValueError(
                f"PDF header {header_path!r} has invalid bin_edges_{dimension}"
            )
        info["nbin_with_overflow"] = info["nbin"] + 2
        info["bin_centers"] = _bin_centers(info)
        dimensions.append(info)

    shape = tuple(item["nbin_with_overflow"] for item in dimensions)
    total_bins = int(np.prod(shape))
    if "total_bins" in header and header["total_bins"] != total_bins:
        raise ValueError(
            f"PDF header {header_path!r} total_bins={header['total_bins']} "
            f"does not match dimensions ({total_bins})"
        )
    header["total_bins"] = total_bins
    for dimension, info in enumerate(dimensions):
        expected_stride = int(np.prod(shape[dimension + 1 :], dtype=int))
        if "stride" in info and info["stride"] != expected_stride:
            raise ValueError(
                f"PDF header {header_path!r} has invalid stride{dimension + 1}"
            )
        info["stride"] = expected_stride
    header["dimensions"] = dimensions
    header["shape"] = shape
    header["header_path"] = os.path.abspath(header_path)
    return header


def _compare_headers(reference, candidate, path):
    for key in ("format", "distribution", "ndim", "weight", "weight_variable", "total_bins", "cycle"):
        if reference.get(key) != candidate.get(key):
            raise ValueError(
                f"PDF shard metadata mismatch for {key!r} in {path!r}: "
                f"{candidate.get(key)!r} != {reference.get(key)!r}"
            )
    for number, (left, right) in enumerate(
        zip(reference["dimensions"], candidate["dimensions"]), start=1
    ):
        for key in ("variable", "nbin", "scale", "linthresh", "stride"):
            if left.get(key) != right.get(key):
                raise ValueError(
                    f"PDF shard metadata mismatch for dimension {number} {key!r} "
                    f"in {path!r}"
                )
        if not np.array_equal(left["bin_edges"], right["bin_edges"]):
            raise ValueError(
                f"PDF shard metadata mismatch for dimension {number} bin edges in {path!r}"
            )


def _header_for_shard(shard_path, requested_header, reference_path):
    local_headers = list(_header_candidates(shard_path))
    for candidate in local_headers:
        if os.path.exists(candidate):
            return candidate
    if requested_header is not None:
        return reference_path
    raise FileNotFoundError(f"PDF shard {shard_path!r} has no readable metadata header")


def _read_v2_preamble(payload, path, header, expected_layout):
    if len(payload) < _V2_PREAMBLE.size:
        raise ValueError(f"truncated PDF V2 payload {path!r}: missing preamble")
    magic, version, layout, ndim, rank, count, time_value, cycle = (
        _V2_PREAMBLE.unpack_from(payload)
    )
    if magic != _V2_MAGIC or version != 2:
        raise ValueError(f"invalid PDF V2 magic or version in {path!r}")
    if layout != expected_layout:
        raise ValueError(
            f"PDF V2 layout mismatch in {path!r}: {layout} != {expected_layout}"
        )
    if ndim != header["ndim"]:
        raise ValueError(
            f"PDF V2 dimension mismatch in {path!r}: {ndim} != {header['ndim']}"
        )
    return int(count), float(time_value), int(cycle), _V2_PREAMBLE.size


def _read_sparse_file(path, header):
    with open(path, "rb") as handle:
        payload = handle.read()
    if payload.startswith(_V2_MAGIC):
        nnz, time_value, cycle, offset = _read_v2_preamble(payload, path, header, 1)
        expected_size = offset + nnz * (
            np.dtype(np.uint64).itemsize + np.dtype(np.float64).itemsize
        )
        if len(payload) != expected_size:
            qualifier = "truncated" if len(payload) < expected_size else "has trailing bytes"
            raise ValueError(
                f"sparse PDF shard {path!r} {qualifier}: nnz={nnz}, "
                f"expected {expected_size} bytes, found {len(payload)}"
            )
        records = np.frombuffer(
            payload,
            dtype=np.dtype([("index", np.uint64), ("value", np.float64)]),
            count=nnz,
            offset=offset,
        )
        indices = records["index"].copy()
        values = records["value"].copy()
        total_bins = header["total_bins"]
        if nnz and np.any(indices >= total_bins):
            bad = int(indices[indices >= total_bins][0])
            raise ValueError(
                f"sparse PDF shard {path!r} has out-of-range index {bad}; "
                f"total_bins={total_bins}"
            )
        if nnz and np.unique(indices).size != nnz:
            raise ValueError(f"sparse PDF shard {path!r} contains duplicate COO indices")
        return time_value, cycle, indices, values

    total_bins = header["total_bins"]
    if len(payload) < 12:
        raise ValueError(f"truncated sparse PDF shard {path!r}: missing time or nnz")
    time_value = float(np.frombuffer(payload, dtype=np.float64, count=1)[0])
    nnz = int(np.frombuffer(payload, dtype=np.uint32, count=1, offset=8)[0])
    expected_size = 12 + nnz * (np.dtype(np.uint32).itemsize + np.dtype(np.float64).itemsize)
    if len(payload) != expected_size:
        qualifier = "truncated" if len(payload) < expected_size else "has trailing bytes"
        raise ValueError(
            f"sparse PDF shard {path!r} {qualifier}: nnz={nnz}, "
            f"expected {expected_size} bytes, found {len(payload)}"
        )
    indices = np.frombuffer(payload, dtype=np.uint32, count=nnz, offset=12).copy()
    values = np.frombuffer(payload, dtype=np.float64, count=nnz, offset=12 + 4 * nnz).copy()
    if nnz and np.any(indices >= total_bins):
        bad = int(indices[indices >= total_bins][0])
        raise ValueError(
            f"sparse PDF shard {path!r} has out-of-range index {bad}; "
            f"total_bins={total_bins}"
        )
    if nnz and np.unique(indices).size != nnz:
        raise ValueError(f"sparse PDF shard {path!r} contains duplicate COO indices")
    return time_value, None, indices, values


def _read_dense_binary(path, header):
    with open(path, "rb") as handle:
        payload = handle.read()
    total_bins = header["total_bins"]
    if payload.startswith(_V2_MAGIC):
        count, time_value, cycle, offset = _read_v2_preamble(payload, path, header, 0)
        if count != total_bins:
            raise ValueError(
                f"dense PDF data {path!r} record count={count}, expected {total_bins}"
            )
        expected_size = offset + total_bins * np.dtype(np.float64).itemsize
        if len(payload) != expected_size:
            raise ValueError(
                f"dense PDF data {path!r} expected {expected_size} bytes, "
                f"found {len(payload)}"
            )
        return time_value, cycle, np.frombuffer(
            payload, dtype=np.float64, count=total_bins, offset=offset
        ).copy()
    raw = np.frombuffer(payload, dtype=np.float64)
    if raw.size != total_bins + 1:
        raise ValueError(
            f"dense PDF data {path!r} expected {total_bins + 1} float64 values, "
            f"found {raw.size}"
        )
    return float(raw[0]), None, raw[1:].copy()


def _read_legacy_data(data_path, header):
    time_value = None
    rows = []
    with open(data_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            time_match = _LEGACY_TIME_RE.match(line)
            if time_match is not None:
                if time_value is not None:
                    raise ValueError(f"legacy PDF data {data_path!r} contains multiple times")
                time_value = float(time_match.group(1))
                continue
            if line.startswith("#"):
                continue
            rows.append(np.fromstring(line, dtype=np.float64, sep=" "))
    if time_value is None:
        raise ValueError(f"legacy PDF data {data_path!r} is missing its time header")
    if header["ndim"] == 1:
        if len(rows) != 1 or rows[0].size != header["shape"][0]:
            raise ValueError(f"legacy PDF data {data_path!r} has an invalid 1-D shape")
        return time_value, rows[0]
    if len(rows) != header["shape"][1] or any(
        row.size != header["shape"][0] for row in rows
    ):
        raise ValueError(f"legacy PDF data {data_path!r} has an invalid 2-D shape")
    return time_value, np.vstack(rows).T


def read_pdf(data_path, header_path=None, reshape=True):
    """Read a PDF output and return ``time``, dense ``pdf``, and ``header``.

    Modern arrays are returned in declared dimension order.  Legacy two-axis
    ASCII data are transposed from their stored row order into the same
    ``(variable_1, variable_2)`` ordering.  Sparse repeated indices are valid
    across shards because shard histograms are summed; repeated indices within
    one shard are rejected as malformed.
    """
    explicit_header = header_path
    if header_path is None:
        header_path = _infer_header_path(data_path)
    header = read_pdf_header(header_path)
    fmt = header["format"]
    if fmt == "legacy_dense":
        if _shard_kind(data_path) != "shared":
            raise ValueError("legacy PDF data do not define rank/node shard reconstruction")
        time_value, data = _read_legacy_data(data_path, header)
        cycle = None
    elif fmt == "dense":
        if _shard_kind(data_path) != "shared":
            raise ValueError(f"dense PDF data {data_path!r} unexpectedly resides in a shard directory")
        time_value, cycle, data = _read_dense_binary(data_path, header)
    else:
        kind = _shard_kind(data_path)
        distribution = header.get("distribution", kind if kind != "shared" else None)
        if kind == "shared" or distribution != kind:
            raise ValueError(
                f"sparse PDF data {data_path!r} cannot be resolved as a {distribution!r} shard family"
            )
        header["distribution"] = distribution
        data = np.zeros(header["total_bins"], dtype=np.float64)
        time_value = None
        cycle = None
        for shard in _glob_partition_files(data_path):
            shard_header_path = _header_for_shard(shard, explicit_header, header_path)
            shard_header = read_pdf_header(shard_header_path)
            shard_header["distribution"] = shard_header.get("distribution", kind)
            _compare_headers(header, shard_header, shard_header_path)
            shard_time, shard_cycle, indices, values = _read_sparse_file(shard, header)
            if time_value is None:
                time_value = shard_time
                cycle = shard_cycle
            elif shard_time != time_value:
                raise ValueError(
                    f"PDF shard time mismatch in {shard!r}: {shard_time!r} != {time_value!r}"
                )
            if shard_cycle != cycle:
                raise ValueError(
                    f"PDF shard cycle mismatch in {shard!r}: {shard_cycle!r} != {cycle!r}"
                )
            np.add.at(data, indices, values)
        if time_value is None:
            raise ValueError(f"no PDF shards found for {data_path!r}")

    if not reshape:
        data = np.asarray(data).reshape(-1, order="C")
    else:
        data = np.asarray(data).reshape(header["shape"], order="C")
    result = {"time": time_value, "pdf": data, "header": header}
    if cycle is not None:
        result["cycle"] = cycle
    elif "cycle" in header:
        result["cycle"] = header["cycle"]
    return result


__all__ = ["read_pdf", "read_pdf_header"]
