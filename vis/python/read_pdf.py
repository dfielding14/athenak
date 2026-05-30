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
import sys

import numpy as np


_HEADER_KEY_RE = re.compile(r"^([A-Za-z0-9_]+)\s*=\s*(.*)$")
_BIN_EDGE_RE = re.compile(r"^bin_edges_(\d+)$")
_PDF_DATA_RE = re.compile(r"^(.*)\.(\d{5})\.pdf$")
_SHARD_DIRECTORY_RE = re.compile(r"^(rank|node)_([0-9]{8})$")
_LEGACY_VARIABLE_RE = re.compile(r"^#\s*\[(\d+)\]\s*=\s*(\S+)")
_LEGACY_TIME_RE = re.compile(r"^#\s*time\s*=\s*(\S+)\s*$")
_SCALES = frozenset(("linear", "log", "symlog"))
_V2_MAGIC = b"AKPDFV2\0"
_V2_PREAMBLE = struct.Struct("=8sIIIIQdq")
_MAX_DENSE_ALLOCATION_BYTES = 512 * 1024 * 1024
_MAX_HEADER_READ_BYTES = 16 * 1024 * 1024
_MAX_PAYLOAD_READ_BYTES = 512 * 1024 * 1024


def _parse_bool(value):
    normalized = value.strip().lower()
    if normalized in ("true", "t", "1", "yes", "y"):
        return True
    if normalized in ("false", "f", "0", "no", "n"):
        return False
    raise ValueError(f"PDF header has invalid boolean value {value!r}")


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


def _header_retained_bytes(header):
    return sum(
        dimension["bin_edges"].nbytes + dimension["bin_centers"].nbytes
        for dimension in header["dimensions"]
    )


def _require_finite(value, label):
    if not np.isfinite(value):
        raise ValueError(f"{label} must be finite")


def _require_file_size(path, limit, label):
    size = os.path.getsize(path)
    if size > limit:
        raise ValueError(
            f"{label} {path!r} requires reading {size} bytes, exceeding the "
            f"practical file-read limit of {limit} bytes"
        )
    return size


def _partition_info(path):
    leaf = os.path.basename(os.path.dirname(os.path.abspath(path)))
    match = _SHARD_DIRECTORY_RE.fullmatch(leaf)
    if match is not None:
        return match.group(1), int(match.group(2))
    if leaf.startswith(("rank_", "node_")):
        raise ValueError(f"invalid PDF shard directory {leaf!r} for {path!r}")
    return "shared", None


def _shard_kind(path):
    return _partition_info(path)[0]


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
    shard_ids = []
    for candidate in files:
        candidate_kind, shard_id = _partition_info(candidate)
        if candidate_kind != kind:
            raise ValueError(
                f"PDF shard {candidate!r} does not match {kind!r} inventory"
            )
        shard_ids.append(shard_id)
    if len(set(shard_ids)) != len(shard_ids):
        raise ValueError(f"PDF {kind} shard inventory contains duplicate IDs")
    expected_ids = set(range(len(shard_ids)))
    actual_ids = set(shard_ids)
    if actual_ids != expected_ids:
        raise ValueError(
            f"PDF {kind} shard inventory is incomplete: "
            f"expected IDs {sorted(expected_ids)!r}, found {sorted(actual_ids)!r}"
        )
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
    return np.linspace(info["bin_min"], info["bin_max"], info["nbin"] + 1)


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


def _symlog_metadata_peak_bytes(edge_count):
    """Conservatively bound vectorized symlog edge and center temporaries."""
    center_count = max(edge_count - 1, 0)
    float_bytes = np.dtype(np.float64).itemsize
    bool_bytes = np.dtype(bool).itemsize
    return (
        edge_count * (6 * float_bytes + bool_bytes)
        + center_count * (8 * float_bytes + bool_bytes)
    )


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


def _ascii_token_list_peak_bytes(text, token_count):
    """Conservatively bound split token strings and their pointer list."""
    pointer_bytes = np.dtype(np.intp).itemsize
    return (
        sys.getsizeof([])
        + 2 * token_count * pointer_bytes
        + token_count * sys.getsizeof("")
        + 4 * len(text)
    )


def _ascii_float_parse_peak_bytes(text, token_count):
    """Bound token strings, parsed floats, lists, and the NumPy array."""
    pointer_bytes = np.dtype(np.intp).itemsize
    float_list = sys.getsizeof([]) + token_count * (
        sys.getsizeof(0.0) + pointer_bytes
    )
    return (
        _ascii_token_list_peak_bytes(text, token_count)
        + float_list
        + token_count * np.dtype(np.float64).itemsize
    )


def _parse_ascii_floats(text, label, externally_retained_bytes=0):
    """Strictly parse one bounded ASCII float row."""
    if not text.isascii():
        raise ValueError(f"{label} must contain ASCII numeric tokens")
    token_count = _count_ascii_tokens(text)
    _require_allocation((token_count,), np.dtype(np.float64).itemsize, label)
    _require_retained_bytes(
        externally_retained_bytes + _ascii_float_parse_peak_bytes(text, token_count),
        label + " parsing peak",
    )
    try:
        return np.array([float(token) for token in text.split()], dtype=np.float64)
    except ValueError as exc:
        raise ValueError(f"{label} contains an invalid numeric token") from exc


def _edge_validation_peak_bytes(edge_count):
    """Bound isfinite and monotonicity temporaries while edges remain retained."""
    interval_count = max(edge_count - 1, 0)
    return (
        edge_count * np.dtype(bool).itemsize
        + interval_count * (np.dtype(np.float64).itemsize + np.dtype(bool).itemsize)
    )


def _generated_edge_peak_bytes(info):
    """Bound generated edges and scale-specific source temporaries."""
    edge_bytes = (info["nbin"] + 1) * np.dtype(np.float64).itemsize
    if info["scale"] == "symlog":
        return _symlog_metadata_peak_bytes(info["nbin"] + 1)
    if info["scale"] == "log":
        return 2 * edge_bytes
    return edge_bytes


def _bin_center_peak_bytes(info):
    """Bound bin-center output and vectorized source temporaries."""
    if info["scale"] == "symlog":
        return _symlog_metadata_peak_bytes(info["nbin"] + 1)
    return 2 * info["nbin"] * np.dtype(np.float64).itemsize


def _read_legacy_header(header_path, externally_retained_bytes=0):
    dimensions = []
    variable_names = {}
    retained_bytes = externally_retained_bytes
    _require_file_size(header_path, _MAX_HEADER_READ_BYTES, "legacy PDF header")
    with open(header_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            match = _LEGACY_VARIABLE_RE.match(line)
            if match is not None:
                variable_names[int(match.group(1))] = match.group(2)
                continue
            if not line or line.startswith("#"):
                continue
            edges = _parse_ascii_floats(
                line,
                "legacy PDF bin edges",
                externally_retained_bytes=retained_bytes,
            )
            if edges.size < 2:
                raise ValueError(f"legacy PDF header {header_path!r} has invalid bin row")
            _require_retained_bytes(
                retained_bytes
                + edges.nbytes
                + _edge_validation_peak_bytes(edges.size),
                "legacy PDF bin-edge validation peak",
            )
            if not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
                raise ValueError(
                    f"legacy PDF header {header_path!r} has non-increasing bin edges"
                )
            _require_retained_bytes(
                retained_bytes
                + edges.nbytes
                + 2 * (edges.size - 1) * np.dtype(np.float64).itemsize,
                "legacy PDF retained bin metadata",
            )
            centers = 0.5 * (edges[:-1] + edges[1:])
            retained_bytes += edges.nbytes + centers.nbytes
            dimensions.append(
                {
                    "variable": variable_names.get(len(dimensions) + 1),
                    "nbin": int(edges.size - 1),
                    "nbin_with_overflow": int(edges.size + 1),
                    "bin_edges": edges,
                    "bin_centers": centers,
                    "scale": None,
                }
            )
    if not dimensions or len(dimensions) > 2:
        raise ValueError(
            f"legacy PDF header {header_path!r} must contain one or two bin rows"
        )
    shape = tuple(item["nbin_with_overflow"] for item in dimensions)
    total_bins = _require_allocation(
        shape, np.dtype(np.float64).itemsize, "legacy PDF dense histogram"
    )
    _require_retained_bytes(
        retained_bytes + total_bins * np.dtype(np.float64).itemsize,
        "legacy PDF retained histogram and bin metadata",
    )
    return {
        "format": "legacy_dense",
        "distribution": "shared",
        "ndim": len(dimensions),
        "dimensions": dimensions,
        "shape": shape,
        "total_bins": total_bins,
        "header_path": os.path.abspath(header_path),
    }


def _read_pdf_header(header_path, externally_retained_bytes=0):
    """Parse a modern ``.header.pdf`` or legacy ``.bins.pdf`` metadata file."""
    if header_path.endswith(".bins.pdf"):
        return _read_legacy_header(header_path, externally_retained_bytes)

    header = {}
    dims = {}
    parsed_edge_bytes = externally_retained_bytes
    _require_file_size(header_path, _MAX_HEADER_READ_BYTES, "PDF header")
    with open(header_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            match = _HEADER_KEY_RE.match(line)
            if match is None:
                raise ValueError(
                    f"malformed PDF header line in {header_path!r}: {line!r}"
                )
            key, value = match.group(1), match.group(2).strip()
            if key in ("format", "distribution", "weight", "weight_variable",
                       "binary_magic"):
                header[key] = value
            elif key == "layout":
                header["format"] = value
            elif key in (
                "ndim",
                "total_bins",
                "cycle",
                "rank",
                "node",
                "number_of_ranks",
                "number_of_nodes",
            ):
                header[key] = int(value)
            else:
                edge_match = _BIN_EDGE_RE.match(key)
                if edge_match is not None:
                    dimension = int(edge_match.group(1))
                    edges = _parse_ascii_floats(
                        value,
                        f"PDF bin_edges_{dimension}",
                        externally_retained_bytes=parsed_edge_bytes,
                    )
                    dims.setdefault(dimension, {})["bin_edges"] = edges
                    parsed_edge_bytes += edges.nbytes
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
        raise ValueError(
            f"unsupported PDF format in {header_path!r}: {header.get('format')!r}"
        )
    header.setdefault("format", "dense")
    if header["format"] == "dense":
        distribution = header.get("distribution", "shared")
        if distribution != "shared":
            raise ValueError(
                f"dense PDF header {header_path!r} has invalid "
                f"distribution={distribution!r}"
            )
        header["distribution"] = distribution
    if "binary_magic" in header and header["binary_magic"] != "AKPDFV2":
        raise ValueError(
            f"unsupported PDF binary magic in {header_path!r}: "
            f"{header['binary_magic']!r}"
        )
    if "ndim" not in header or not 1 <= header["ndim"] <= 4:
        raise ValueError(f"PDF header {header_path!r} must declare ndim between 1 and 4")

    dimensions = []
    retained_bytes = externally_retained_bytes
    pending_explicit_edge_bytes = sum(
        info["bin_edges"].nbytes for info in dims.values() if "bin_edges" in info
    )
    for dimension in range(1, header["ndim"] + 1):
        if dimension not in dims:
            raise ValueError(
                f"PDF header {header_path!r} is missing dimension {dimension}"
            )
        info = dims[dimension]
        if not isinstance(info.get("nbin"), int) or info["nbin"] <= 0:
            raise ValueError(f"PDF header {header_path!r} has invalid nbin{dimension}")
        if "variable" not in info:
            raise ValueError(
                f"PDF header {header_path!r} is missing variable_{dimension}"
            )
        scale = info.get("scale", "log" if info.get("logscale", False) else "linear")
        if scale not in _SCALES:
            raise ValueError(
                f"PDF header {header_path!r} has invalid scale{dimension}={scale!r}"
            )
        info["scale"] = scale
        if scale == "symlog":
            if (
                not np.isfinite(info.get("linthresh", 0.0))
                or info.get("linthresh", 0.0) <= 0.0
            ):
                raise ValueError(
                    f"PDF header {header_path!r} has invalid linthresh{dimension}"
                )
        else:
            info.setdefault("linthresh", 1.0)
        if "bin_min" not in info or "bin_max" not in info:
            raise ValueError(f"PDF header {header_path!r} is missing bin bounds")
        if not np.isfinite(info["bin_min"]) or not np.isfinite(info["bin_max"]):
            raise ValueError(f"PDF header {header_path!r} has non-finite bin bounds")
        if info["bin_min"] >= info["bin_max"]:
            raise ValueError(f"PDF header {header_path!r} has reversed bin bounds")
        if scale == "log" and info["bin_min"] <= 0.0:
            raise ValueError(f"PDF header {header_path!r} has non-positive log bounds")
        _require_allocation(
            (info["nbin"] + 1,),
            np.dtype(np.float64).itemsize,
            f"PDF dimension {dimension} bin edges",
        )
        current_edges_are_explicit = "bin_edges" in info
        if not current_edges_are_explicit:
            _require_retained_bytes(
                retained_bytes
                + pending_explicit_edge_bytes
                + _generated_edge_peak_bytes(info),
                "PDF generated bin-edge peak",
            )
            info["bin_edges"] = _generated_edges(info)
        edges = info["bin_edges"]
        live_edge_bytes = retained_bytes + pending_explicit_edge_bytes
        if not current_edges_are_explicit:
            live_edge_bytes += edges.nbytes
        if edges.size != info["nbin"] + 1:
            raise ValueError(
                f"PDF header {header_path!r} has wrong number of bin_edges_{dimension}"
            )
        _require_retained_bytes(
            live_edge_bytes + _edge_validation_peak_bytes(edges.size),
            "PDF bin-edge validation peak",
        )
        if not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0.0):
            raise ValueError(
                f"PDF header {header_path!r} has invalid bin_edges_{dimension}"
            )
        info["nbin_with_overflow"] = info["nbin"] + 2
        _require_retained_bytes(
            live_edge_bytes + _bin_center_peak_bytes(info),
            "PDF retained bin metadata",
        )
        info["bin_centers"] = _bin_centers(info)
        if current_edges_are_explicit:
            pending_explicit_edge_bytes -= edges.nbytes
        retained_bytes += edges.nbytes + info["bin_centers"].nbytes
        dimensions.append(info)

    shape = tuple(item["nbin_with_overflow"] for item in dimensions)
    total_bins = _require_allocation(
        shape, np.dtype(np.float64).itemsize, "PDF dense histogram"
    )
    _require_retained_bytes(
        retained_bytes + total_bins * np.dtype(np.float64).itemsize,
        "PDF retained histogram and bin metadata",
    )
    if "total_bins" in header and header["total_bins"] != total_bins:
        raise ValueError(
            f"PDF header {header_path!r} total_bins={header['total_bins']} "
            f"does not match dimensions ({total_bins})"
        )
    header["total_bins"] = total_bins
    for dimension, info in enumerate(dimensions):
        expected_stride = _checked_product(
            shape[dimension + 1:], f"PDF stride{dimension + 1}"
        )
        if "stride" in info and info["stride"] != expected_stride:
            raise ValueError(
                f"PDF header {header_path!r} has invalid stride{dimension + 1}"
            )
        info["stride"] = expected_stride
    header["dimensions"] = dimensions
    header["shape"] = shape
    header["header_path"] = os.path.abspath(header_path)
    return header


def read_pdf_header(header_path):
    """Read and validate PDF metadata without exposing internal peak accounting."""
    return _read_pdf_header(header_path)


def _compare_headers(reference, candidate, path):
    for key in (
        "format",
        "distribution",
        "ndim",
        "weight",
        "weight_variable",
        "total_bins",
        "cycle",
        "binary_magic",
    ):
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
                f"PDF shard metadata mismatch for dimension {number} "
                f"bin edges in {path!r}"
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
    if "cycle" in header and int(cycle) != header["cycle"]:
        raise ValueError(
            f"PDF V2 cycle mismatch in {path!r}: {int(cycle)} != {header['cycle']}"
        )
    _require_finite(time_value, f"PDF V2 payload time in {path!r}")
    return int(rank), int(count), float(time_value), int(cycle), _V2_PREAMBLE.size


def _read_sparse_file(path, header, retained_bytes):
    payload_bytes = _require_file_size(path, _MAX_PAYLOAD_READ_BYTES, "PDF payload")
    _require_retained_bytes(
        retained_bytes + payload_bytes,
        "PDF sparse reconstruction peak",
    )
    with open(path, "rb") as handle:
        prefix = handle.read(_V2_PREAMBLE.size)
        if header.get("binary_magic") == "AKPDFV2" and not prefix.startswith(
            _V2_MAGIC
        ):
            raise ValueError(f"PDF V2 payload {path!r} is missing its AKPDFV2 preamble")
        if prefix.startswith(_V2_MAGIC):
            rank, nnz, time_value, cycle, offset = _read_v2_preamble(
                prefix, path, header, 1
            )
            index_itemsize = np.dtype(np.uint64).itemsize
            expected_size = offset + _checked_product(
                (nnz, index_itemsize + np.dtype(np.float64).itemsize),
                "sparse PDF V2 payload",
            )
        else:
            if len(prefix) < 12:
                raise ValueError(
                    f"truncated sparse PDF shard {path!r}: missing time or nnz"
                )
            time_value = float(np.frombuffer(prefix, dtype=np.float64, count=1)[0])
            _require_finite(time_value, f"sparse PDF payload time in {path!r}")
            nnz = int(np.frombuffer(prefix, dtype=np.uint32, count=1, offset=8)[0])
            rank = None
            cycle = None
            offset = 12
            index_itemsize = np.dtype(np.uint32).itemsize
            expected_size = offset + _checked_product(
                (nnz, index_itemsize + np.dtype(np.float64).itemsize),
                "sparse PDF payload",
            )
        if payload_bytes != expected_size:
            qualifier = (
                "truncated" if payload_bytes < expected_size else "has trailing bytes"
            )
            raise ValueError(
                f"sparse PDF shard {path!r} {qualifier}: nnz={nnz}, "
                f"expected {expected_size} bytes, found {payload_bytes}"
            )
        index_bytes = _checked_product(
            (nnz, index_itemsize), "sparse PDF index copy"
        )
        value_bytes = _checked_product(
            (nnz, np.dtype(np.float64).itemsize), "sparse PDF value copy"
        )
        _require_retained_bytes(
            retained_bytes + payload_bytes + 2 * index_bytes + value_bytes,
            "PDF sparse reconstruction peak",
        )
        unique_bytes = 2 * index_bytes + nnz * np.dtype(bool).itemsize
        _require_retained_bytes(
            retained_bytes
            + payload_bytes
            + index_bytes
            + value_bytes
            + unique_bytes,
            "PDF sparse duplicate-validation peak",
        )
        del prefix
        handle.seek(0)
        payload = handle.read()
    if payload.startswith(_V2_MAGIC):
        records = np.frombuffer(
            payload,
            dtype=np.dtype([("index", np.uint64), ("value", np.float64)]),
            count=nnz,
            offset=offset,
        )
        indices = records["index"].copy()
        values = records["value"].copy()
        total_bins = header["total_bins"]
        if nnz and np.max(indices) >= total_bins:
            bad = int(np.max(indices))
            raise ValueError(
                f"sparse PDF shard {path!r} has out-of-range index {bad}; "
                f"total_bins={total_bins}"
            )
        if nnz and np.unique(indices).size != nnz:
            raise ValueError(f"sparse PDF shard {path!r} contains duplicate COO indices")
        return time_value, cycle, indices, values, rank

    total_bins = header["total_bins"]
    indices = np.frombuffer(payload, dtype=np.uint32, count=nnz, offset=12).copy()
    values = np.frombuffer(
        payload, dtype=np.float64, count=nnz, offset=12 + 4 * nnz
    ).copy()
    if nnz and np.max(indices) >= total_bins:
        bad = int(np.max(indices))
        raise ValueError(
            f"sparse PDF shard {path!r} has out-of-range index {bad}; "
            f"total_bins={total_bins}"
        )
    if nnz and np.unique(indices).size != nnz:
        raise ValueError(f"sparse PDF shard {path!r} contains duplicate COO indices")
    return time_value, None, indices, values, None


def _read_dense_binary(path, header):
    payload_bytes = _require_file_size(path, _MAX_PAYLOAD_READ_BYTES, "PDF payload")
    _require_retained_bytes(
        _header_retained_bytes(header)
        + payload_bytes
        + header["total_bins"] * np.dtype(np.float64).itemsize,
        "PDF dense reconstruction peak",
    )
    with open(path, "rb") as handle:
        payload = handle.read()
    total_bins = header["total_bins"]
    if header.get("binary_magic") == "AKPDFV2" and not payload.startswith(_V2_MAGIC):
        raise ValueError(f"PDF V2 payload {path!r} is missing its AKPDFV2 preamble")
    if payload.startswith(_V2_MAGIC):
        _, count, time_value, cycle, offset = _read_v2_preamble(
            payload, path, header, 0
        )
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
    time_value = float(raw[0])
    _require_finite(time_value, f"dense PDF payload time in {path!r}")
    return time_value, None, raw[1:].copy()


def _validate_sparse_shard_metadata(path, header, payload_rank):
    kind, shard_id = _partition_info(path)
    distribution = header.get("distribution")
    if distribution != kind:
        raise ValueError(
            f"PDF shard {path!r} declares distribution={distribution!r}, "
            f"but resides in a {kind!r} layout"
        )
    if header.get("binary_magic") == "AKPDFV2":
        id_key = "rank" if kind == "rank" else "node"
        count_key = "number_of_ranks" if kind == "rank" else "number_of_nodes"
        if id_key not in header or count_key not in header:
            raise ValueError(
                f"PDF V2 {kind} shard {path!r} is missing required inventory metadata"
            )
    if kind == "rank":
        if "rank" in header and header["rank"] != shard_id:
            raise ValueError(
                f"PDF shard {path!r} declares rank={header['rank']}, "
                f"but its directory identifies {shard_id}"
            )
        if payload_rank is not None and payload_rank != shard_id:
            raise ValueError(
                f"PDF shard {path!r} payload declares rank={payload_rank}, "
                f"but its directory identifies {shard_id}"
            )
    elif kind == "node":
        if "node" in header and header["node"] != shard_id:
            raise ValueError(
                f"PDF shard {path!r} declares node={header['node']}, "
                f"but its directory identifies {shard_id}"
            )
        if header.get("binary_magic") == "AKPDFV2" and "node" not in header:
            raise ValueError(f"PDF node shard {path!r} is missing node metadata")
    for count_key, expected_distribution in (
        ("number_of_ranks", "rank"),
        ("number_of_nodes", "node"),
    ):
        if count_key in header:
            if header[count_key] <= 0:
                raise ValueError(
                    f"PDF shard {path!r} has invalid {count_key}={header[count_key]}"
                )
            if distribution != expected_distribution:
                raise ValueError(
                    f"PDF shard {path!r} defines {count_key} for "
                    f"distribution={distribution!r}"
                )


def _validate_sparse_sibling_inventory(files, headers):
    distribution = headers[0]["distribution"]
    count_key = "number_of_ranks" if distribution == "rank" else "number_of_nodes"
    count_values = [header.get(count_key) for header in headers]
    if all(value is None for value in count_values):
        return
    if any(value is None for value in count_values) or len(set(count_values)) != 1:
        raise ValueError(
            f"PDF {distribution} shard inventory has inconsistent {count_key} metadata"
        )
    expected_count = count_values[0]
    if expected_count != len(files):
        raise ValueError(
            f"PDF {distribution} shard inventory is incomplete: "
            f"expected {expected_count} shards, found {len(files)}"
        )
    expected_ids = set(range(len(files)))
    actual_ids = {_partition_info(path)[1] for path in files}
    if actual_ids != expected_ids:
        raise ValueError(
            f"PDF {distribution} shard inventory is incomplete: "
            f"expected IDs {sorted(expected_ids)!r}, found {sorted(actual_ids)!r}"
        )


def _read_legacy_data(data_path, header):
    time_value = None
    rows = []
    retained_bytes = _header_retained_bytes(header)
    _require_file_size(data_path, _MAX_PAYLOAD_READ_BYTES, "legacy PDF payload")
    with open(data_path, "r") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            time_match = _LEGACY_TIME_RE.match(line)
            if time_match is not None:
                if time_value is not None:
                    raise ValueError(
                        f"legacy PDF data {data_path!r} contains multiple times"
                    )
                time_value = float(time_match.group(1))
                _require_finite(time_value, f"legacy PDF payload time in {data_path!r}")
                continue
            if line.startswith("#"):
                continue
            row = _parse_ascii_floats(
                line,
                "legacy PDF payload row",
                externally_retained_bytes=retained_bytes,
            )
            retained_bytes += row.nbytes
            rows.append(row)
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
    _require_retained_bytes(
        retained_bytes + header["total_bins"] * np.dtype(np.float64).itemsize,
        "legacy PDF stacked payload peak",
    )
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
            raise ValueError(
                "legacy PDF data do not define rank/node shard reconstruction"
            )
        time_value, data = _read_legacy_data(data_path, header)
        cycle = None
    elif fmt == "dense":
        if _shard_kind(data_path) != "shared":
            raise ValueError(
                f"dense PDF data {data_path!r} unexpectedly resides in a shard directory"
            )
        time_value, cycle, data = _read_dense_binary(data_path, header)
    else:
        kind = _shard_kind(data_path)
        distribution = header.get("distribution", kind if kind != "shared" else None)
        if kind == "shared" or distribution != kind:
            raise ValueError(
                f"sparse PDF data {data_path!r} cannot be resolved as a "
                f"{distribution!r} shard family"
            )
        header["distribution"] = distribution
        retained_bytes = _header_retained_bytes(header)
        _require_retained_bytes(
            retained_bytes + header["total_bins"] * np.dtype(np.float64).itemsize,
            "PDF sparse reconstruction",
        )
        data = np.zeros(header["total_bins"], dtype=np.float64)
        time_value = None
        cycle = None
        shard_files = _glob_partition_files(data_path)
        shard_headers = []
        for shard in shard_files:
            shard_header_path = _header_for_shard(shard, explicit_header, header_path)
            shard_header = _read_pdf_header(
                shard_header_path, retained_bytes + data.nbytes
            )
            if (
                shard_header.get("binary_magic") == "AKPDFV2"
                and "distribution" not in shard_header
            ):
                raise ValueError(
                    f"PDF V2 shard header {shard_header_path!r} is missing "
                    "distribution metadata"
                )
            shard_header["distribution"] = shard_header.get("distribution", kind)
            _compare_headers(header, shard_header, shard_header_path)
            shard_time, shard_cycle, indices, values, payload_rank = _read_sparse_file(
                shard,
                shard_header,
                retained_bytes + _header_retained_bytes(shard_header) + data.nbytes,
            )
            _validate_sparse_shard_metadata(shard, shard_header, payload_rank)
            shard_headers.append(
                {
                    key: shard_header[key]
                    for key in (
                        "distribution",
                        "number_of_ranks",
                        "number_of_nodes",
                        "rank",
                        "node",
                    )
                    if key in shard_header
                }
            )
            if time_value is None:
                time_value = shard_time
                cycle = shard_cycle
            elif shard_time != time_value:
                raise ValueError(
                    f"PDF shard time mismatch in {shard!r}: "
                    f"{shard_time!r} != {time_value!r}"
                )
            if shard_cycle != cycle:
                raise ValueError(
                    f"PDF shard cycle mismatch in {shard!r}: "
                    f"{shard_cycle!r} != {cycle!r}"
                )
            np.add.at(data, indices, values)
            del indices, values, shard_header
        if time_value is None:
            raise ValueError(f"no PDF shards found for {data_path!r}")
        _validate_sparse_sibling_inventory(shard_files, shard_headers)
        header.pop("rank", None)
        header.pop("node", None)

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
