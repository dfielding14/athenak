"""
Functions to:
  (1) convert bin --> Python dictionary
  (2) convert Python dictionary --> athdf(xdmf) files

This module contains a collection of helper functions for readng and
writing athena file data formats. More information is provided in the
function docstrings.

----

In order to translate a binary file into athdf and corresponding xdmf
files, you could do the following:

  import bin_convert
  import os

  binary_fname = "path/to/file.bin"
  athdf_fname = os.path.splitext(binary_fname)[0] + ".athdf"
  xdmf_fname = athdf_fname + ".xdmf"
  filedata = bin_convert.read_binary(binary_fname)
  bin_convert.write_athdf(athdf_fname, filedata)
  bin_convert.write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)

Notice that write_xdmf_for(...) function expects the relative path to
the athdf file from the xdmf, so please be aware of this requirement.

----

The read_*(...) functions return a filedata dictionary-like object with

    filedata['header'] = array of strings
        ordered array of header, including all the header information
    filedata['time'] = float
        time from input file
    filedata['cycle'] = int
        cycle from input file
    filedata['var_names'] = array of strings
        ordered array of variable names, like ['dens', 'eint', ...]
    filedata['n_mbs'] = int
        total number of meshblocks in the file
    filedata['nx1_mb'] = int
        number of cells in x1 direction in MeshBlock
    filedata['nx2_mb'] = int
        number of cells in x2 direction in MeshBlock
    filedata['nx3_mb'] = int
        number of cells in x3 direction in MeshBlock
    filedata['nx1_out_mb'] = int
        number of output cells in x1 direction in MeshBlock (useful for slicing)
    filedata['nx2_out_mb'] = int
        number of output cells in x2 direction in MeshBlock (useful for slicing)
    filedata['nx3_out_mb'] = int
        number of output cells in x3 direction in MeshBlock (useful for slicing)
    filedata['Nx1'] = int
        total number of cell in x1 direction in root grid
    filedata['Nx2'] = int
        total number of cell in x2 direction in root grid
    filedata['Nx3'] = int
        total number of cell in x3 direction in root grid
    filedata['x1min'] = float
        coordinate minimum of root grid in x1 direction
    filedata['x1max'] = float
        coordinate maximum of root grid in x1 direction
    filedata['x2min'] = float
        coordinate minimum of root grid in x2 direction
    filedata['x2max'] = float
        coordinate maximum of root grid in x2 direction
    filedata['x3min'] = float
        coordinate minimum of root grid in x3 direction
    filedata['x3max'] = float
        coordinate maximum of root grid in x3 direction
    filedata['nvars'] = int
        number of output variables (including magnetic field if it exists)
    filedata['mb_index'] = array with shape [n_mbs, 6]
        is,ie,js,je,ks,ke range for output MeshBlock indexing (useful for slicing)
    filedata['mb_logical'] = array with shape [n_mbs, 4]
        i,j,k,level coordinates for each MeshBlock
    filedata['mb_geometry'] = array with shape [n_mbs, 6]
        x1i,x2i,x3i,dx1,dx2,dx3 including cell-centered location of left-most
        cell and offsets between cells
    filedata['mb_data'] = dict of arrays with shape [n_mbs, nx3, nx2, nx1]
        {'var1':var1_array, 'var2':var2_array, ...} dictionary of fluid data arrays
        for each variable in var_names
"""

import numpy as np
import os
import h5py
from numbers import Integral
import re
import sys

if __package__:
    from .io_reader_common import (
        ReaderLimits,
        add_reader_limit_arguments,
        bounded_partition_files,
        checked_product,
        checked_sum,
        conservative_text_metadata_bytes,
        count_ascii_tokens,
        metadata_record_bytes,
        normalize_reader_limits,
        reader_limits_from_args,
        require_live_bytes,
        require_read_bytes,
    )
else:
    from io_reader_common import (
        ReaderLimits,
        add_reader_limit_arguments,
        bounded_partition_files,
        checked_product,
        checked_sum,
        conservative_text_metadata_bytes,
        count_ascii_tokens,
        metadata_record_bytes,
        normalize_reader_limits,
        reader_limits_from_args,
        require_live_bytes,
        require_read_bytes,
    )


_SHARD_DIRECTORY_RE = re.compile(r"^(rank|node)_([0-9]{8})$")
_MAX_MESHBLOCK_ALLOCATION_BYTES = 512 * 1024 * 1024
_MAX_ATHDF_ALLOCATION_BYTES = 512 * 1024 * 1024
_MAX_BINARY_HEADER_BYTES = 16 * 1024 * 1024
_MAX_BINARY_PREHEADER_LINES = 4096
_MAX_LOGICAL_LEVEL = 30


def _partition_info(filename):
    """Return a strict ``(kind, integer ID)`` pair for a shard path."""
    shard_name = os.path.basename(os.path.dirname(os.path.abspath(filename)))
    match = _SHARD_DIRECTORY_RE.fullmatch(shard_name)
    if match is not None:
        return match.group(1), int(match.group(2))
    if shard_name.startswith(("rank_", "node_")):
        raise ValueError(f"invalid binary shard directory {shard_name!r}")
    return "shared", None


def _binary_reader_limits(limits):
    return normalize_reader_limits(
        limits,
        max_live_bytes=_MAX_MESHBLOCK_ALLOCATION_BYTES,
        max_header_read_bytes=_MAX_BINARY_HEADER_BYTES,
    )


def _athdf_reader_limits(limits):
    return normalize_reader_limits(
        limits,
        max_live_bytes=_MAX_ATHDF_ALLOCATION_BYTES,
        max_header_read_bytes=_MAX_BINARY_HEADER_BYTES,
    )


def _require_binary_bytes(nbytes, label, limits=None):
    """Reject one reader allocation before materializing unreasonable data."""
    return require_live_bytes(nbytes, label, _binary_reader_limits(limits))


def _checked_binary_product(values, label):
    """Multiply non-negative integer extents without permitting malformed values."""
    return checked_product(values, label)


def _require_athdf_bytes(nbytes, label, limits=None):
    """Reject athdf-like helper allocations before calling NumPy."""
    return require_live_bytes(nbytes, label, _athdf_reader_limits(limits))


def _preflight_athdf_coordinates(
    nx_vals, dtype, limits=None, retained_bytes=0
):
    """Bound face and center coordinate arrays for one athdf-like conversion."""
    itemsize = np.dtype(dtype).itemsize
    total = sum(
        _checked_binary_product((nx + 1 + nx, itemsize), "athdf-like coordinates")
        for nx in nx_vals
    )
    total = checked_sum(
        (retained_bytes, total), "athdf-like coordinate arrays"
    )
    _require_athdf_bytes(total, "athdf-like coordinate arrays", limits)
    linspace_temporary = max(nx + 1 for nx in nx_vals) * np.dtype(np.float64).itemsize
    _require_athdf_bytes(
        total + linspace_temporary,
        "athdf-like coordinate generation peak",
        limits,
    )
    return total


def _preflight_athdf_outputs(
    shape,
    quantities,
    dtype,
    return_levels,
    restricted_shape,
    retained_bytes=0,
    limits=None,
):
    """Bound dense fields, optional levels, and optional restriction bookkeeping."""
    cells = _checked_binary_product(shape, "athdf-like output shape")
    total = _checked_binary_product(
        (len(quantities), cells, np.dtype(dtype).itemsize),
        "athdf-like output fields",
    )
    if return_levels:
        total += _checked_binary_product(
            (cells, np.dtype(np.int32).itemsize), "athdf-like level map"
        )
    if restricted_shape is not None:
        total += _checked_binary_product(
            (*restricted_shape, np.dtype(bool).itemsize),
            "athdf-like restriction map",
        )
    total += retained_bytes
    _require_athdf_bytes(total, "athdf-like output arrays", limits)
    return total


def _preflight_athdf_temporary(
    retained_bytes, shape, dtype, label, limits=None
):
    """Bound one temporary array while the dense athdf-like result is retained."""
    temporary = _checked_binary_product(
        (*shape, np.dtype(dtype).itemsize), label
    )
    _require_athdf_bytes(retained_bytes + temporary, label + " peak", limits)


def _validate_binary_grid_metadata(
    root_grid_size, block_size, nghost, bounds, family, filename
):
    """Reject malformed grid dimensions before reconstructing dense products."""
    if any(size <= 0 for size in root_grid_size):
        raise ValueError(f"{family} file {filename!r} has invalid root-grid dimensions")
    if any(size <= 0 for size in block_size):
        raise ValueError(f"{family} file {filename!r} has invalid MeshBlock dimensions")
    if any(root % block != 0 for root, block in zip(root_grid_size, block_size)):
        raise ValueError(
            f"{family} file {filename!r} has root-grid dimensions that are not "
            "divisible by its MeshBlock dimensions"
        )
    if nghost < 0:
        raise ValueError(f"{family} file {filename!r} has a negative ghost-zone count")
    if not np.all(np.isfinite(bounds)):
        raise ValueError(f"{family} file {filename!r} has non-finite metadata")
    if any(lower >= upper for lower, upper in zip(bounds[::2], bounds[1::2])):
        raise ValueError(f"{family} file {filename!r} has invalid coordinate bounds")
    return tuple(root // block for root, block in zip(root_grid_size, block_size))


def _validate_meshblock_metadata(
    logical,
    geometry,
    family,
    filename,
    root_meshblock_counts,
    root_grid_size,
    root_bounds,
    nghost,
    limits=None,
    externally_retained_bytes=0,
):
    """Reject malformed logical locations before athdf-like exponent arithmetic."""
    _require_binary_bytes(
        checked_sum(
            (
                externally_retained_bytes,
                _logical_owner_dictionary_bytes(len(logical)),
            ),
            f"{family} logical MeshBlock validation peak in {filename!r}",
        ),
        f"{family} logical MeshBlock validation peak in {filename!r}",
        limits,
    )
    seen = set()
    for block_num, (location, bounds) in enumerate(zip(logical, geometry)):
        values = [int(value) for value in location]
        if any(value < 0 for value in values[:3]):
            raise ValueError(
                f"{family} meshblock {block_num} in {filename!r} has a negative "
                "logical location"
            )
        if values[3] < 0 or values[3] > _MAX_LOGICAL_LEVEL:
            raise ValueError(
                f"{family} meshblock {block_num} in {filename!r} has invalid "
                f"logical level {values[3]}"
            )
        logical_key = tuple(values)
        if logical_key in seen:
            raise ValueError(
                f"{family} file {filename!r} contains duplicate logical MeshBlocks"
            )
        seen.add(logical_key)
        for axis, value in enumerate(values[:3]):
            upper = root_meshblock_counts[axis] * 2 ** values[3]
            if value >= upper:
                raise ValueError(
                    f"{family} meshblock {block_num} in {filename!r} has "
                    f"out-of-range logical location {value} for axis {axis + 1}"
                )
        if not np.all(np.isfinite(bounds)):
            raise ValueError(
                f"{family} meshblock {block_num} in {filename!r} has non-finite "
                "geometry"
            )
        if any(lower >= upper for lower, upper in zip(bounds[::2], bounds[1::2])):
            raise ValueError(
                f"{family} meshblock {block_num} in {filename!r} has invalid geometry"
            )
        expected_bounds = []
        for axis, (root_blocks, logical_index) in enumerate(
            zip(root_meshblock_counts, values[:3])
        ):
            root_lower, root_upper = root_bounds[2 * axis:2 * axis + 2]
            refinement = 2 ** values[3]
            block_width = (root_upper - root_lower) / (root_blocks * refinement)
            expected_lower = root_lower + logical_index * block_width
            expected_upper = expected_lower + block_width
            expected_bounds.append((expected_lower, expected_upper))
        if any(
            not np.isclose(
                lower,
                expected_lower,
                rtol=0.0,
                atol=min(
                    abs(expected_upper - expected_lower) / 8.0,
                    2.0
                    * np.finfo(np.asarray(bounds).dtype).eps
                    * max(abs(expected_lower), abs(expected_upper), 1.0),
                ),
            )
            or not np.isclose(
                upper,
                expected_upper,
                rtol=0.0,
                atol=min(
                    abs(expected_upper - expected_lower) / 8.0,
                    2.0
                    * np.finfo(np.asarray(bounds).dtype).eps
                    * max(abs(expected_lower), abs(expected_upper), 1.0),
                ),
            )
            for (lower, upper), (expected_lower, expected_upper) in zip(
                zip(bounds[::2], bounds[1::2]), expected_bounds
            )
        ):
            raise ValueError(
                f"{family} meshblock {block_num} in {filename!r} has geometry "
                "outside its logical location"
            )


def _validate_coarsening_factor(factor, grid_sizes, filename):
    """Require the producer's documented coarsened-binary factor contract."""
    if (
        factor < 2
        or factor & (factor - 1)
        or any(size <= 0 or size % factor != 0 for size in grid_sizes)
    ):
        raise ValueError(
            f"coarsened binary file {filename!r} has invalid coarsening factor"
        )


def _validate_coarsened_moments(filename, pheader, nvars, var_list):
    """Validate the producer's scalar or grouped four-moment metadata contract."""
    number_of_moments = int(pheader["number of moments"])
    if number_of_moments not in (1, 4):
        raise ValueError(
            f"coarsened binary file {filename!r} has invalid number of moments"
        )
    if nvars % number_of_moments != 0:
        raise ValueError(
            f"coarsened binary file {filename!r} has incomplete moment groups"
        )
    if number_of_moments == 4:
        suffixes = ("_1st", "_2nd", "_3rd", "_4th")
        for group_start in range(0, nvars, number_of_moments):
            group = var_list[group_start:group_start + number_of_moments]
            if any(not label.endswith(suffix) for label, suffix in zip(group, suffixes)):
                raise ValueError(
                    f"coarsened binary file {filename!r} has malformed moment labels"
                )
            roots = [label[: -len(suffix)] for label, suffix in zip(group, suffixes)]
            if len(set(roots)) != 1:
                raise ValueError(
                    f"coarsened binary file {filename!r} has malformed moment labels"
                )
    return number_of_moments


def _read_limited_binary_line(fp, family, label, budget=None, limits=None):
    """Read one binary metadata line without accepting an unbounded record."""
    limits = _binary_reader_limits(limits)
    retained_raw_bytes = 0 if budget is None else budget[0]
    remaining_header_bytes = limits.max_header_read_bytes - retained_raw_bytes
    remaining_text_bytes = limits.max_live_bytes // 4 - retained_raw_bytes
    read_cap = max(0, min(remaining_header_bytes, remaining_text_bytes))
    line = fp.readline(read_cap + 1)
    if len(line) > read_cap:
        _require_binary_bytes(
            conservative_text_metadata_bytes(retained_raw_bytes + len(line)),
            f"{family} metadata text peak in {fp.name!r}",
            limits,
        )
        if budget is not None and retained_raw_bytes:
            raise ValueError(
                f"{family} metadata records in {fp.name!r} require at least "
                f"{retained_raw_bytes + len(line)} bytes, exceeding the practical "
                f"header limit of {limits.max_header_read_bytes} bytes"
            )
        raise ValueError(
            f"{family} {label} in {fp.name!r} exceeds the practical header limit "
            f"of {limits.max_header_read_bytes} bytes"
        )
    if budget is not None:
        budget[0] += len(line)
        if budget[0] > limits.max_header_read_bytes:
            raise ValueError(
                f"{family} metadata records in {fp.name!r} require {budget[0]} bytes, "
                f"exceeding the practical header limit of "
                f"{limits.max_header_read_bytes} bytes"
            )
        _require_binary_bytes(
            conservative_text_metadata_bytes(budget[0]),
            f"{family} metadata text peak in {fp.name!r}",
            limits,
        )
    return line


def _parse_binary_metadata_bytes(raw_line, family, label, filename):
    """Parse one fixed ``key=value`` metadata record without token expansion."""
    key, separator, value = raw_line.strip().partition(b"=")
    if not separator or not key or not value or b"=" in value:
        raise ValueError(f"invalid {family} {label} in {filename!r}")
    return key.strip(), value.strip()


def _parse_binary_metadata_int(raw_line, family, label, filename):
    """Parse one bounded integer metadata record without materializing tokens."""
    _, value = _parse_binary_metadata_bytes(raw_line, family, label, filename)
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(f"invalid {family} {label} in {filename!r}") from exc


def _parse_binary_preheader_text(raw_line, family, filename):
    """Decode one bounded preheader pair after the cumulative text preflight."""
    key, value = _parse_binary_metadata_bytes(
        raw_line, family, "preheader metadata", filename
    )
    try:
        return key.decode("utf-8"), value.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"invalid {family} preheader text in {filename!r}") from exc


def _validate_binary_format_header(raw_line, family, filename):
    """Require the exact versioned binary signature without splitting attacker text."""
    signature = raw_line.rstrip(b"\r\n")
    if not signature.startswith(b"Athena"):
        raise TypeError(f"bad file format in {filename!r} (should start with \"Athena\")")
    if signature != b"Athena binary output version=1.1":
        raise TypeError(f"unsupported {family} file format version in {filename!r}")


def _read_binary_parameter_dump(
    fp, header_size, family, limits=None, externally_retained_bytes=0
):
    """Read one embedded athinput dump after bounding its declared byte count."""
    limits = _binary_reader_limits(limits)
    if header_size < 0 or header_size > limits.max_header_read_bytes:
        raise ValueError(
            f"{family} parameter header in {fp.name!r} declares {header_size} bytes, "
            "exceeding the practical header limit of "
            f"{limits.max_header_read_bytes} bytes"
        )
    retained_bytes = checked_sum(
        (conservative_text_metadata_bytes(header_size), sys.getsizeof([])),
        f"{family} retained parameter header",
    )
    _require_binary_bytes(
        checked_sum(
            (externally_retained_bytes, retained_bytes),
            f"{family} parameter header parsing peak",
        ),
        f"{family} parameter header parsing peak",
        limits,
    )
    lines = []
    retained_string_bytes = 0
    bytes_read = 0
    while bytes_read < header_size:
        raw_line = fp.readline(header_size - bytes_read)
        if not raw_line:
            raise ValueError(
                f"truncated {family} parameter header in {fp.name!r}: "
                f"expected {header_size} bytes, found {bytes_read}"
            )
        bytes_read += len(raw_line)
        try:
            line = raw_line.decode("utf-8").split("#", 1)[0].strip()
        except UnicodeDecodeError as exc:
            raise ValueError(
                f"invalid {family} parameter header text in {fp.name!r}"
            ) from exc
        if not line:
            continue
        lines.append(line)
        retained_string_bytes = checked_sum(
            (retained_string_bytes, sys.getsizeof(line)),
            f"{family} retained parameter-header strings",
        )
        retained_bytes = checked_sum(
            (
                conservative_text_metadata_bytes(header_size),
                sys.getsizeof(lines),
                retained_string_bytes,
            ),
            f"{family} retained parameter header",
        )
        _require_binary_bytes(
            checked_sum(
                (externally_retained_bytes, retained_bytes),
                f"{family} parameter header parsing peak",
            ),
            f"{family} parameter header parsing peak",
            limits,
        )
    return lines, retained_bytes


def _preflight_binary_preheader_text(
    fp, family, budget, limits=None, retained_record_bytes=0
):
    """Charge decoded preheader text and retained records before materialization."""
    retained_bytes = checked_sum(
        (
            conservative_text_metadata_bytes(budget[0]),
            retained_record_bytes,
        ),
        f"{family} retained preheader metadata",
    )
    _require_binary_bytes(
        retained_bytes,
        f"{family} preheader metadata peak in {fp.name!r}",
        limits,
    )
    return retained_bytes


def _variable_token_peak_bytes(text, token_count):
    """Conservatively bound split variable strings and their pointer lists."""
    pointer_bytes = np.dtype(np.intp).itemsize
    return (
        sys.getsizeof([])
        + 2 * token_count * pointer_bytes
        + token_count * sys.getsizeof("")
        + 4 * len(text)
    )


def _variable_dictionary_bytes(variable_count):
    """Conservatively bound the retained per-variable dictionary and lists."""
    pointer_bytes = np.dtype(np.intp).itemsize
    return (
        sys.getsizeof({})
        + variable_count * (8 * pointer_bytes + sys.getsizeof([]))
    )


def _variable_name_set_bytes(variable_names):
    """Conservatively bound transient duplicate-name set construction."""
    return sys.getsizeof(set()) + sum(
        metadata_record_bytes(variable, charge_value=False)
        for variable in variable_names
    )


def _logical_owner_dictionary_bytes(meshblock_count):
    """Conservatively bound transient duplicate-detection dictionary entries."""
    pointer_bytes = np.dtype(np.intp).itemsize
    per_entry_bytes = (
        sys.getsizeof((0, 0, 0, 0))
        + 4 * sys.getsizeof(0)
        + 32 * pointer_bytes
    )
    return sys.getsizeof({}) + meshblock_count * per_entry_bytes


def _meshblock_arrays_bytes(
    mb_index, mb_logical, mb_geometry, mb_data, retained_variable_bytes
):
    """Return retained binary arrays before they are attached to filedata."""
    return checked_sum(
        (
            retained_variable_bytes,
            mb_index.nbytes,
            mb_logical.nbytes,
            mb_geometry.nbytes,
            *(values.nbytes for values in mb_data.values()),
        ),
        "retained binary MeshBlock arrays",
    )


def _is_effective_slice(
    levels, logical_locations, axis, retained_binary_bytes, limits
):
    """Classify one reduced axis after charging tuple-set materialization."""
    _require_athdf_bytes(
        checked_sum(
            (
                retained_binary_bytes,
                _logical_owner_dictionary_bytes(len(levels)),
            ),
            "athdf-like slice classification peak",
        ),
        "athdf-like slice classification peak",
        limits,
    )
    other_locations = {
        (
            level,
            logical_locations[index, (axis + 1) % 3],
            logical_locations[index, (axis + 2) % 3],
        )
        for index, level in enumerate(levels)
    }
    return len(other_locations) == len(levels)


def _read_binary_variable_list(
    fp, family, budget, limits=None, externally_retained_bytes=0
):
    """Read and preflight one variable-name record before splitting it."""
    limits = _binary_reader_limits(limits)
    raw_line = _read_limited_binary_line(
        fp, family, "variable list", budget, limits
    )
    try:
        text = raw_line.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(
            f"invalid {family} variable list in {fp.name!r}"
        ) from exc
    token_count = count_ascii_tokens(text)
    retained_bytes = _variable_token_peak_bytes(text, token_count)
    _require_binary_bytes(
        checked_sum(
            (externally_retained_bytes, retained_bytes),
            f"{family} variable tokenization peak in {fp.name!r}",
        ),
        f"{family} variable tokenization peak in {fp.name!r}",
        limits,
    )
    return text.split()[1:], retained_bytes


def _is_partitioned_path(filename):
    """Return whether *filename* is under a rank_*/node_* shard directory."""
    return _partition_info(filename)[0] != "shared"


def _glob_partition_files(shard_filename, limits=None):
    """Return all sibling files belonging to a rank- or node-sharded output."""
    shard_filename = os.path.abspath(shard_filename)
    shard_kind, _ = _partition_info(shard_filename)
    return bounded_partition_files(
        shard_filename,
        shard_kind,
        "binary",
        _binary_reader_limits(limits),
        _partition_info,
    )


def _optional_pheader_int(pheader, key, filename, family):
    """Parse one additive binary preheader integer without requiring it in old files."""
    if key not in pheader:
        return None
    try:
        return int(pheader[key])
    except ValueError as exc:
        raise ValueError(
            f"{family} file {filename!r} has invalid {key!r} metadata"
        ) from exc


def _binary_partition_metadata(filename, pheader, family):
    """Validate optional additive shard metadata while preserving old files."""
    kind, shard_id = _partition_info(filename)
    distribution = pheader.get("distribution")
    if distribution is not None:
        if distribution not in ("shared", "rank", "node"):
            raise ValueError(
                f"{family} file {filename!r} has invalid distribution {distribution!r}"
            )
        if distribution != kind:
            raise ValueError(
                f"{family} file {filename!r} declares distribution={distribution!r} "
                f"but resides in a {kind!r} layout"
            )
    rank = _optional_pheader_int(pheader, "rank", filename, family)
    node = _optional_pheader_int(pheader, "node", filename, family)
    number_of_ranks = _optional_pheader_int(
        pheader, "number of ranks", filename, family
    )
    number_of_nodes = _optional_pheader_int(
        pheader, "number of nodes", filename, family
    )
    number_of_meshblocks = _optional_pheader_int(
        pheader, "number of meshblocks", filename, family
    )
    if rank is not None:
        if kind != "rank" or rank != shard_id:
            raise ValueError(
                f"{family} file {filename!r} declares rank={rank}, "
                f"but resides in {os.path.basename(os.path.dirname(filename))!r}"
            )
    if node is not None:
        if kind != "node" or node != shard_id:
            raise ValueError(
                f"{family} file {filename!r} declares node={node}, "
                f"but resides in {os.path.basename(os.path.dirname(filename))!r}"
            )
    if number_of_ranks is not None:
        if kind != "rank" or number_of_ranks <= 0:
            raise ValueError(
                f"{family} file {filename!r} has invalid 'number of ranks' metadata"
            )
    if number_of_nodes is not None:
        if kind != "node" or number_of_nodes <= 0:
            raise ValueError(
                f"{family} file {filename!r} has invalid 'number of nodes' metadata"
            )
    if number_of_meshblocks is not None and number_of_meshblocks < 0:
        raise ValueError(
            f"{family} file {filename!r} has invalid 'number of meshblocks' metadata"
        )
    if kind == "node":
        required = {
            "distribution": distribution,
            "node": node,
            "number of nodes": number_of_nodes,
            "number of meshblocks": number_of_meshblocks,
        }
        missing = [key for key, value in required.items() if value is None]
        if missing:
            raise ValueError(
                f"{family} node shard {filename!r} is missing required "
                f"inventory metadata: {', '.join(missing)}"
            )
    return {
        "distribution": distribution if distribution is not None else kind,
        "rank": rank,
        "node": node,
        "number_of_ranks": number_of_ranks,
        "number_of_nodes": number_of_nodes,
        "number_of_meshblocks": number_of_meshblocks,
    }


def _validate_binary_sibling_inventory(shard_files, shard_data, family):
    """Require a complete unique ID inventory when additive count metadata exists."""
    distribution = shard_data[0]["distribution"]
    if distribution == "shared":
        return
    count_key = "number_of_ranks" if distribution == "rank" else "number_of_nodes"
    id_key = "rank" if distribution == "rank" else "node"
    count_label = count_key.replace("_", " ")
    expected_count = shard_data[0].get(count_key)
    if expected_count is None and all(
        item.get(count_key) is None for item in shard_data
    ):
        return
    if expected_count is None or any(
        item.get(count_key) != expected_count for item in shard_data
    ):
        raise ValueError(
            f"{family} {distribution} shard inventory has inconsistent "
            f"{count_label!r} metadata"
        )
    for expected_id, (path, item) in enumerate(zip(shard_files, shard_data)):
        kind, shard_id = _partition_info(path)
        if kind != distribution or item.get("distribution") != distribution:
            raise ValueError(
                f"{family} {distribution} shard inventory includes "
                f"non-{distribution} file {path!r}"
            )
        if item.get(id_key) is None:
            raise ValueError(
                f"{family} {distribution} shard {path!r} is missing {id_key} metadata"
            )
        if item[id_key] != shard_id:
            raise ValueError(
                f"{family} {distribution} shard {path!r} declares "
                f"{id_key}={item[id_key]}, "
                f"but its directory identifies {shard_id}"
            )
        if shard_id != expected_id:
            raise ValueError(
                f"{family} {distribution} shard inventory is incomplete: "
                f"expected ID {expected_id}, found {shard_id}"
            )
    if expected_count != len(shard_files):
        raise ValueError(
            f"{family} {distribution} shard inventory is incomplete: "
            f"expected {expected_count} shards, found {len(shard_files)}"
        )


def _validate_binary_shard(reference, candidate, path, family):
    """Validate metadata that must agree before binary shards can be combined."""
    keys = (
        "time",
        "cycle",
        "var_names",
        "nvars",
        "Nx1",
        "Nx2",
        "Nx3",
        "x1min",
        "x1max",
        "x2min",
        "x2max",
        "x3min",
        "x3max",
        "nx1_mb",
        "nx2_mb",
        "nx3_mb",
    )
    if family == "coarsened binary":
        keys += ("number_of_moments",)
    for key in keys:
        if candidate[key] != reference[key]:
            raise ValueError(
                f"{family} shard metadata mismatch for {key!r} in {path!r}: "
                f"{candidate[key]!r} != {reference[key]!r}"
            )
    if candidate["header"] != reference["header"]:
        raise ValueError(f"{family} shard metadata mismatch for 'header' in {path!r}")


def _binary_payload_bytes(filedata):
    """Return the retained variable-array bytes for one parsed binary shard."""
    return sum(
        values.nbytes
        for variable in filedata["var_names"]
        for values in filedata["mb_data"][variable]
    )


def _binary_metadata_bytes(filedata):
    """Return the retained fixed MeshBlock-array bytes for one parsed shard."""
    return sum(
        filedata[key].nbytes
        for key in ("mb_index", "mb_logical", "mb_geometry")
    )


def _binary_variable_metadata_bytes(filedata):
    """Return the retained variable-name and per-variable container bytes."""
    return filedata.get("_retained_variable_metadata_bytes", 0)


def _binary_retained_bytes(filedata):
    """Return retained binary payload, MeshBlock arrays, and variable metadata."""
    return checked_sum(
        (
            _binary_payload_bytes(filedata),
            _binary_metadata_bytes(filedata),
            _binary_variable_metadata_bytes(filedata),
        ),
        "retained binary input",
    )


def _binary_reader_materialization_peak_bytes(filedata):
    """Bound the parser's final-array and one-record materialization peak."""
    retained_bytes = _binary_retained_bytes(filedata)
    if filedata["n_mbs"] == 0:
        return retained_bytes
    geometry_itemsize = filedata["mb_geometry"].dtype.itemsize
    serialized_fixed_bytes = 24 + 16 + 6 * geometry_itemsize
    value_bytes = sum(
        filedata["mb_data"][variable][0].nbytes
        for variable in filedata["var_names"]
    )
    return checked_sum(
        (retained_bytes, 4096, serialized_fixed_bytes, 2 * value_bytes),
        "binary reader materialization peak",
    )


def _remaining_binary_limits(limits, retained_bytes, label):
    """Reserve retained sibling arrays before parsing the next binary shard."""
    if retained_bytes >= limits.max_live_bytes:
        _require_binary_bytes(retained_bytes, label, limits)
        raise ValueError(
            f"{label} leaves no live-memory budget for another shard"
        )
    return ReaderLimits(
        max_live_bytes=limits.max_live_bytes - retained_bytes,
        max_header_read_bytes=limits.max_header_read_bytes,
        max_payload_read_bytes=limits.max_payload_read_bytes,
    )


def _combine_partitioned_binary(shard_filename, reader, family, limits=None):
    """Combine one binary shard family, retaining valid empty sliced shards."""
    limits = _binary_reader_limits(limits)
    shard_files, inventory_bytes = _glob_partition_files(shard_filename, limits)
    shard_data = []
    reference = None
    aggregate_payload_bytes = 0
    aggregate_metadata_bytes = 0
    aggregate_variable_metadata_bytes = 0
    for path in shard_files:
        retained_bytes = checked_sum(
            (
                inventory_bytes,
                aggregate_payload_bytes,
                aggregate_metadata_bytes,
                aggregate_variable_metadata_bytes,
            ),
            f"{family} assembled retained arrays",
        )
        candidate_limits = _remaining_binary_limits(
            limits, retained_bytes, f"{family} assembled retained arrays"
        )
        try:
            candidate = reader(path, limits=candidate_limits)
        except ValueError as exc:
            if (
                retained_bytes
                and "exceeding the practical allocation limit" in str(exc)
            ):
                raise ValueError(
                    f"{family} assembled materialization peak requires more than "
                    f"{limits.max_live_bytes} bytes while parsing {path!r}"
                ) from exc
            raise
        if reference is None:
            reference = candidate
        else:
            _validate_binary_shard(reference, candidate, path, family)
        shard_data.append(candidate)
        aggregate_payload_bytes += _binary_payload_bytes(candidate)
        _require_binary_bytes(
            aggregate_payload_bytes,
            f"{family} assembled shard payload",
            limits,
        )
        aggregate_metadata_bytes += _binary_metadata_bytes(candidate)
        _require_binary_bytes(
            aggregate_metadata_bytes,
            f"{family} assembled shard metadata",
            limits,
        )
        aggregate_variable_metadata_bytes += _binary_variable_metadata_bytes(
            candidate
        )
        _require_binary_bytes(
            aggregate_variable_metadata_bytes,
            f"{family} assembled variable metadata",
            limits,
        )
        _require_binary_bytes(
            checked_sum(
                (
                    inventory_bytes,
                    aggregate_payload_bytes,
                    aggregate_metadata_bytes,
                    aggregate_variable_metadata_bytes,
                ),
                f"{family} assembled retained arrays",
            ),
            f"{family} assembled retained arrays",
            limits,
        )
    _validate_binary_sibling_inventory(shard_files, shard_data, family)
    _require_binary_bytes(
        2 * aggregate_payload_bytes,
        f"{family} assembled shard payload",
        limits,
    )
    _require_binary_bytes(
        2 * aggregate_metadata_bytes,
        f"{family} assembled shard metadata",
        limits,
    )
    _require_binary_bytes(
        2 * aggregate_variable_metadata_bytes,
        f"{family} assembled variable metadata",
        limits,
    )
    _require_binary_bytes(
        checked_sum(
            (
                inventory_bytes,
                2 * aggregate_payload_bytes,
                2 * aggregate_metadata_bytes,
                2 * aggregate_variable_metadata_bytes,
            ),
            f"{family} assembled materialization peak",
        ),
        f"{family} assembled materialization peak",
        limits,
    )

    nonempty = next((item for item in shard_data if item["n_mbs"] > 0), reference)
    nonempty_shape = tuple(nonempty[f"nx{axis}_out_mb"] for axis in (1, 2, 3))
    for path, candidate in zip(shard_files, shard_data):
        if candidate["n_mbs"] == 0:
            continue
        candidate_shape = tuple(candidate[f"nx{axis}_out_mb"] for axis in (1, 2, 3))
        if candidate_shape != nonempty_shape:
            raise ValueError(
                f"{family} shard output-shape mismatch in {path!r}: "
                f"{candidate_shape!r} != {nonempty_shape!r}"
            )

    total_meshblocks = sum(item["n_mbs"] for item in shard_data)
    logical_owner_bytes = _logical_owner_dictionary_bytes(total_meshblocks)
    _require_binary_bytes(
        checked_sum(
            (
                inventory_bytes,
                aggregate_payload_bytes,
                aggregate_metadata_bytes,
                aggregate_variable_metadata_bytes,
                logical_owner_bytes,
            ),
            f"{family} assembled logical ownership peak",
        ),
        f"{family} assembled logical ownership peak",
        limits,
    )
    logical_owners = {}
    for path, item in zip(shard_files, shard_data):
        for logical in item["mb_logical"]:
            logical_key = tuple(int(value) for value in logical)
            if logical_key in logical_owners:
                raise ValueError(
                    f"{family} shards contain duplicate logical MeshBlock "
                    f"{logical_key!r} in {path!r} and {logical_owners[logical_key]!r}"
                )
            logical_owners[logical_key] = path
    del logical_owners

    combined = nonempty.copy()
    combined["mb_index"] = np.empty((total_meshblocks, 6), dtype=np.int64)
    combined["mb_logical"] = np.empty((total_meshblocks, 4), dtype=np.int32)
    combined["mb_geometry"] = np.empty(
        (total_meshblocks, 6), dtype=nonempty["mb_geometry"].dtype
    )
    combined["mb_data"] = {}
    for var in reference["var_names"]:
        if total_meshblocks:
            sample = next(iter(nonempty["mb_data"][var]))
            combined["mb_data"][var] = np.empty(
                (total_meshblocks, *sample.shape), dtype=sample.dtype
            )
        else:
            combined["mb_data"][var] = np.asarray([])

    cursor = 0
    for path, item in zip(shard_files, shard_data):
        next_cursor = cursor + item["n_mbs"]
        if next_cursor == cursor:
            continue
        combined["mb_index"][cursor:next_cursor] = item["mb_index"]
        combined["mb_logical"][cursor:next_cursor] = item["mb_logical"]
        combined["mb_geometry"][cursor:next_cursor] = item["mb_geometry"]
        for var in reference["var_names"]:
            for offset, values in enumerate(item["mb_data"][var]):
                combined["mb_data"][var][cursor + offset] = values
        cursor = next_cursor

    combined["n_mbs"] = total_meshblocks
    meshblock_counts = [item["number_of_meshblocks"] for item in shard_data]
    combined["number_of_meshblocks"] = (
        sum(meshblock_counts)
        if all(count is not None for count in meshblock_counts)
        else None
    )
    combined["shard_files"] = shard_files
    combined["distribution"] = _partition_info(shard_filename)[0]
    combined["rank"] = None
    combined["node"] = None
    combined["_retained_variable_metadata_bytes"] = (
        aggregate_variable_metadata_bytes + inventory_bytes
    )
    return combined


def _require_meshblocks_for_athdf(filedata, filename):
    """Reject athdf-like conversion when a valid binary input has no MeshBlocks."""
    if filedata["n_mbs"] == 0:
        raise ValueError(
            f"cannot convert {filename!r} to athdf-like data: "
            "binary output contains no meshblocks"
        )


def _validate_requested_level(level):
    if (
        isinstance(level, (bool, np.bool_))
        or not isinstance(level, Integral)
        or level < 0
        or level > _MAX_LOGICAL_LEVEL
    ):
        raise ValueError(f"athdf-like conversion has invalid logical level {level!r}")
    return int(level)


def _crop_athdf_coordinates(data, selection):
    """Return coordinate arrays that describe the selected dense result only."""
    for axis, (lower, upper) in enumerate(selection, start=1):
        face = f"x{axis}f"
        center = f"x{axis}v"
        data[face] = data[face][lower:upper + 1]
        data[center] = data[center][lower:upper]


def _athdf_selection(data, nx_vals, lower_bounds, upper_bounds):
    """Return cell-index selections whose cells intersect requested bounds."""
    selection = []
    for axis, (nx, lower, upper) in enumerate(
        zip(nx_vals, lower_bounds, upper_bounds), start=1
    ):
        faces = data[f"x{axis}f"]
        low_index = 0
        high_index = nx
        if lower is not None:
            low_index = np.searchsorted(faces, lower, side="right") - 1
            low_index = min(max(low_index, 0), nx)
        if upper is not None:
            high_index = np.searchsorted(faces, upper, side="left")
            high_index = min(max(high_index, 0), nx)
        selection.extend((low_index, high_index))
    return tuple(selection)


def _populate_root_athdf_coordinates(
    data, filedata, nx_vals, root_grid_size, level, num_ghost, dtype, center_funcs
):
    """Populate root-grid coordinates, including requested exterior ghost cells."""
    for axis, (nx, root_nx, center_func) in enumerate(
        zip(nx_vals, root_grid_size, center_funcs), start=1
    ):
        face = f"x{axis}f"
        center = f"x{axis}v"
        lower = filedata[f"x{axis}min"]
        upper = filedata[f"x{axis}max"]
        if nx > 1 and num_ghost:
            spacing = (upper - lower) / (root_nx * 2**level)
            lower -= num_ghost * spacing
            upper += num_ghost * spacing
        data[face] = np.linspace(lower, upper, nx + 1, dtype=dtype)
        data[center] = np.empty(nx, dtype=dtype)
        for index in range(nx):
            data[center][index] = center_func(
                data[face][index], data[face][index + 1]
            )


def _validate_athdf_num_ghost(
    filedata, block_size, root_grid_size, num_ghost, max_level, level, bounds
):
    """Require a coherent, uniform ghost-zone reconstruction request."""
    if (
        isinstance(num_ghost, (bool, np.bool_))
        or not isinstance(num_ghost, Integral)
        or num_ghost < 0
    ):
        raise ValueError(f"athdf-like conversion has invalid num_ghost={num_ghost!r}")
    num_ghost = int(num_ghost)
    base_size = [filedata["nx1_mb"], filedata["nx2_mb"], filedata["nx3_mb"]]
    if not num_ghost:
        if any(output > base for output, base in zip(block_size, base_size)):
            raise ValueError(
                "athdf-like conversion detected ghost zones but num_ghost is zero"
            )
        return 0
    if level != max_level or not np.all(filedata["mb_logical"][:, 3] == max_level):
        raise ValueError(
            "athdf-like conversion cannot use ghost zones with different "
            "refinement levels"
        )
    if any(bound is not None for bound in bounds):
        raise ValueError("athdf-like conversion cannot select subsets with ghost zones")
    for output, base, root in zip(block_size, base_size, root_grid_size):
        expected = base + 2 * num_ghost if root > 1 else base
        if output != expected:
            raise ValueError(
                "athdf-like conversion num_ghost does not match emitted MeshBlock extents"
            )
    return num_ghost


def _copy_meshblock_to_athdf(
    data,
    filedata,
    quantities,
    block_num,
    level,
    block_size,
    nx_vals,
    selection,
    output_bytes,
    return_levels,
    subsample,
    fast_restrict,
    restricted_data,
    vol_func,
    num_ghost,
    limits=None,
):
    """Place one in-memory binary MeshBlock into an athdf-like dense result."""
    nx1, nx2, nx3 = nx_vals
    i_min, i_max, j_min, j_max, k_min, k_max = selection
    block_level = int(filedata["mb_logical"][block_num, 3])
    block_location = [
        int(value) for value in filedata["mb_logical"][block_num, :3]
    ]
    active = (nx1 > 1, nx2 > 1, nx3 > 1)

    if block_level <= level:
        s = 2 ** (level - block_level)
        interior = [
            size - 2 * num_ghost if is_active else size
            for size, is_active in zip(block_size, active)
        ]
        il_d = block_location[0] * interior[0] * s if active[0] else 0
        jl_d = block_location[1] * interior[1] * s if active[1] else 0
        kl_d = block_location[2] * interior[2] * s if active[2] else 0
        iu_d = il_d + block_size[0] * s if active[0] else 1
        ju_d = jl_d + block_size[1] * s if active[1] else 1
        ku_d = kl_d + block_size[2] * s if active[2] else 1

        il_s, iu_s = max(il_d, i_min) - il_d, min(iu_d, i_max) - il_d
        jl_s, ju_s = max(jl_d, j_min) - jl_d, min(ju_d, j_max) - jl_d
        kl_s, ku_s = max(kl_d, k_min) - kl_d, min(ku_d, k_max) - kl_d
        if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
            return

        il_d, iu_d = max(il_d, i_min) - i_min, min(iu_d, i_max) - i_min
        jl_d, ju_d = max(jl_d, j_min) - j_min, min(ju_d, j_max) - j_min
        kl_d, ku_d = max(kl_d, k_min) - k_min, min(ku_d, k_max) - k_min
        scale1, scale2, scale3 = (
            s if active[0] else 1,
            s if active[1] else 1,
            s if active[2] else 1,
        )
        index_lengths = (iu_s - il_s, ju_s - jl_s, ku_s - kl_s)
        index_bytes = sum(index_lengths) * np.dtype(np.int64).itemsize
        index_source_bytes = max(index_lengths) * np.dtype(np.int64).itemsize
        _require_athdf_bytes(
            output_bytes + index_bytes + index_source_bytes,
            "athdf-like prolongation index peak",
            limits,
        )
        i_indices = np.arange(il_s, iu_s) // scale1
        j_indices = np.arange(jl_s, ju_s) // scale2
        k_indices = np.arange(kl_s, ku_s) // scale3
        temporary_shape = (k_indices.size, j_indices.size, i_indices.size)
        for q in quantities:
            _preflight_athdf_temporary(
                output_bytes + index_bytes,
                temporary_shape,
                data[q].dtype,
                "athdf-like prolongation",
                limits,
            )
            block_data = filedata["mb_data"][q][block_num]
            data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                np.ix_(k_indices, j_indices, i_indices)
            ]
    else:
        s = 2 ** (block_level - level)
        for size, is_active in zip(block_size, active):
            if is_active and size % s != 0:
                raise ValueError(
                    "athdf-like restriction requires MeshBlock dimensions divisible "
                    "by the refinement ratio"
                )
        il_d = block_location[0] * block_size[0] // s if active[0] else 0
        jl_d = block_location[1] * block_size[1] // s if active[1] else 0
        kl_d = block_location[2] * block_size[2] // s if active[2] else 0
        iu_d = il_d + block_size[0] // s if active[0] else 1
        ju_d = jl_d + block_size[1] // s if active[1] else 1
        ku_d = kl_d + block_size[2] // s if active[2] else 1

        il_s, iu_s = max(il_d, i_min) - il_d, min(iu_d, i_max) - il_d
        jl_s, ju_s = max(jl_d, j_min) - jl_d, min(ju_d, j_max) - jl_d
        kl_s, ku_s = max(kl_d, k_min) - kl_d, min(ku_d, k_max) - kl_d
        if il_s >= iu_s or jl_s >= ju_s or kl_s >= ku_s:
            return

        il_d, iu_d = max(il_d, i_min) - i_min, min(iu_d, i_max) - i_min
        jl_d, ju_d = max(jl_d, j_min) - j_min, min(ju_d, j_max) - j_min
        kl_d, ku_d = max(kl_d, k_min) - k_min, min(ku_d, k_max) - k_min
        if active[0]:
            il_s, iu_s = il_s * s, iu_s * s
        if active[1]:
            jl_s, ju_s = jl_s * s, ju_s * s
        if active[2]:
            kl_s, ku_s = kl_s * s, ku_s * s

        if subsample:
            offset1 = s // 2 - 1 if active[0] else 0
            offset2 = s // 2 - 1 if active[1] else 0
            offset3 = s // 2 - 1 if active[2] else 0
            for q in quantities:
                block_data = filedata["mb_data"][q][block_num]
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_data[
                    kl_s + offset3:ku_s:s,
                    jl_s + offset2:ju_s:s,
                    il_s + offset1:iu_s:s,
                ]
        elif fast_restrict:
            offsets1 = range(s) if active[0] else (0,)
            offsets2 = range(s) if active[1] else (0,)
            offsets3 = range(s) if active[2] else (0,)
            divisor = (s if active[0] else 1) * (s if active[1] else 1)
            divisor *= s if active[2] else 1
            for q in quantities:
                block_data = filedata["mb_data"][q][block_num]
                for offset3 in offsets3:
                    for offset2 in offsets2:
                        for offset1 in offsets1:
                            data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] += block_data[
                                kl_s + offset3:ku_s:s,
                                jl_s + offset2:ju_s:s,
                                il_s + offset1:iu_s:s,
                            ]
                data[q][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] /= divisor
        else:
            bounds = filedata["mb_geometry"][block_num]
            exact_temporary = sum(
                (size + 1) * np.dtype(np.float64).itemsize
                for size in block_size
            )
            exact_temporary += sum(
                (upper - lower) * (s if is_active else 1) * np.dtype(np.int64).itemsize
                for lower, upper, is_active in (
                    (il_d, iu_d, active[0]),
                    (jl_d, ju_d, active[1]),
                    (kl_d, ku_d, active[2]),
                )
            )
            exact_temporary += max(
                (upper - lower) * np.dtype(np.int64).itemsize
                for lower, upper in (
                    (il_d, iu_d),
                    (jl_d, ju_d),
                    (kl_d, ku_d),
                )
            )
            _require_athdf_bytes(
                output_bytes + exact_temporary,
                "athdf-like exact restriction peak",
                limits,
            )
            faces = [
                np.linspace(bounds[2 * axis], bounds[2 * axis + 1], size + 1)
                if is_active
                else np.asarray([bounds[2 * axis], bounds[2 * axis + 1]])
                for axis, (size, is_active) in enumerate(zip(block_size, active))
            ]
            i_dest = np.repeat(range(il_d, iu_d), s if active[0] else 1)
            j_dest = np.repeat(range(jl_d, ju_d), s if active[1] else 1)
            k_dest = np.repeat(range(kl_d, ku_d), s if active[2] else 1)
            for k_source, k_dest_value in zip(range(kl_s, ku_s), k_dest):
                for j_source, j_dest_value in zip(range(jl_s, ju_s), j_dest):
                    for i_source, i_dest_value in zip(range(il_s, iu_s), i_dest):
                        volume = vol_func(
                            faces[0][i_source],
                            faces[0][i_source + 1],
                            faces[1][j_source],
                            faces[1][j_source + 1],
                            faces[2][k_source],
                            faces[2][k_source + 1],
                        )
                        for q in quantities:
                            data[q][k_dest_value, j_dest_value, i_dest_value] += (
                                volume * filedata["mb_data"][q][block_num][
                                    k_source, j_source, i_source
                                ]
                            )
            loc1 = block_location[0] // s if active[0] else 0
            loc2 = block_location[1] // s if active[1] else 0
            loc3 = block_location[2] // s if active[2] else 0
            restricted_data[loc3, loc2, loc1] = True

    if return_levels:
        data["Levels"][kl_d:ku_d, jl_d:ju_d, il_d:iu_d] = block_level


def _read_meshblocks(
    fp,
    filesize,
    nghost,
    locsizebytes,
    varsizebytes,
    var_list,
    family,
    limits=None,
    retained_variable_bytes=0,
):
    """Read uniform MeshBlock records directly into bounded retained arrays."""
    limits = _binary_reader_limits(limits)
    dtype_loc = np.float64 if locsizebytes == 8 else np.float32
    dtype_var = np.float64 if varsizebytes == 8 else np.float32
    nvars = len(var_list)
    retained_variable_bytes = checked_sum(
        (retained_variable_bytes, _variable_dictionary_bytes(nvars)),
        f"{family} retained variable metadata in {fp.name!r}",
    )
    _require_binary_bytes(
        retained_variable_bytes,
        f"{family} variable dictionary materialization peak in {fp.name!r}",
        limits,
    )
    serialized_fixed_bytes = 24 + 16 + 6 * locsizebytes
    retained_fixed_bytes = 6 * np.dtype(np.int64).itemsize
    retained_fixed_bytes += 4 * np.dtype(np.int32).itemsize
    retained_fixed_bytes += 6 * np.dtype(dtype_loc).itemsize

    def read_fixed_record():
        if filesize - fp.tell() < serialized_fixed_bytes:
            raise ValueError(
                f"truncated {family} meshblock metadata in {fp.name!r}"
            )
        index = np.frombuffer(fp.read(24), dtype=np.int32).astype(np.int64) - nghost
        logical = np.frombuffer(fp.read(16), dtype=np.int32).copy()
        geometry = np.frombuffer(fp.read(6 * locsizebytes), dtype=dtype_loc).copy()
        shape = (
            int(index[5] - index[4] + 1),
            int(index[3] - index[2] + 1),
            int(index[1] - index[0] + 1),
        )
        if any(length <= 0 for length in shape):
            raise ValueError(
                f"invalid {family} meshblock extent {shape!r} in {fp.name!r}"
            )
        return index, logical, geometry, shape

    if fp.tell() == filesize:
        return (
            np.empty((0, 6), dtype=np.int64),
            np.empty((0, 4), dtype=np.int32),
            np.empty((0, 6), dtype=dtype_loc),
            {var: np.empty((0,), dtype=dtype_var) for var in var_list},
            retained_variable_bytes,
        )

    first_record_offset = fp.tell()
    shape = None
    value_bytes = None
    meshblock_count = 0
    while fp.tell() < filesize:
        _, _, _, current_shape = read_fixed_record()
        if shape is None:
            shape = current_shape
            value_count = _checked_binary_product(
                (nvars, *shape), f"{family} meshblock values in {fp.name!r}"
            )
            value_bytes = _checked_binary_product(
                (value_count, varsizebytes),
                f"{family} meshblock payload in {fp.name!r}",
            )
            require_read_bytes(
                value_bytes,
                limits.max_payload_read_bytes,
                f"{family} meshblock payload in {fp.name!r}",
            )
        elif current_shape != shape:
            raise ValueError(
                f"{family} file {fp.name!r} has nonuniform MeshBlock extents: "
                f"{tuple(reversed(current_shape))!r} != {tuple(reversed(shape))!r}"
            )
        if filesize - fp.tell() < value_bytes:
            raise ValueError(f"truncated {family} meshblock values in {fp.name!r}")
        fp.seek(value_bytes, 1)
        meshblock_count += 1
    fp.seek(first_record_offset)
    metadata_bytes = _checked_binary_product(
        (meshblock_count, retained_fixed_bytes),
        f"{family} retained meshblock metadata in {fp.name!r}",
    )
    payload_bytes = _checked_binary_product(
        (meshblock_count, value_bytes),
        f"{family} retained meshblock payload in {fp.name!r}",
    )
    retained_bytes = checked_sum(
        (retained_variable_bytes, metadata_bytes, payload_bytes),
        f"{family} retained arrays in {fp.name!r}",
    )
    _require_binary_bytes(
        retained_variable_bytes + metadata_bytes,
        f"{family} meshblock metadata arrays in {fp.name!r}",
        limits,
    )
    _require_binary_bytes(
        retained_bytes, f"{family} retained arrays in {fp.name!r}", limits
    )
    _require_binary_bytes(
        checked_sum(
            (
                retained_bytes,
                4096,
                serialized_fixed_bytes,
                2 * value_bytes,
            ),
            f"{family} meshblock materialization peak in {fp.name!r}",
        ),
        f"{family} meshblock materialization peak in {fp.name!r}",
        limits,
    )
    mb_index = np.empty((meshblock_count, 6), dtype=np.int64)
    mb_logical = np.empty((meshblock_count, 4), dtype=np.int32)
    mb_geometry = np.empty((meshblock_count, 6), dtype=dtype_loc)
    mb_data = {
        var: np.empty((meshblock_count, *shape), dtype=dtype_var)
        for var in var_list
    }

    def store_record(record, index, logical, geometry, current_shape):
        if current_shape != shape:
            raise ValueError(
                f"{family} file {fp.name!r} has nonuniform MeshBlock extents: "
                f"{tuple(reversed(current_shape))!r} != {tuple(reversed(shape))!r}"
            )
        raw_values = fp.read(value_bytes)
        if len(raw_values) != value_bytes:
            raise ValueError(
                f"truncated {family} meshblock values in {fp.name!r}: "
                f"expected {value_bytes} bytes, found {len(raw_values)}"
            )
        values = np.frombuffer(raw_values, dtype=dtype_var).reshape((nvars,) + shape)
        mb_index[record] = index
        mb_logical[record] = logical
        mb_geometry[record] = geometry
        for vari, var in enumerate(var_list):
            mb_data[var][record] = values[vari]

    for record in range(meshblock_count):
        store_record(record, *read_fixed_record())
    return (
        mb_index,
        mb_logical,
        mb_geometry,
        mb_data,
        retained_variable_bytes,
    )


def _read_binary(filename, assemble_shards=False, *, limits=None):
    """
    Reads a bin file from filename to dictionary.

    Originally written by Lev Arzamasskiy (leva@ias.edu) on 11/15/2021
    Updated to support mesh refinement by George Wong (gnwong@ias.edu) on 01/27/2022
    Made faster by Drummond Fielding on 09/09/2024

    args:
      filename - string
          filename of bin file to read
      assemble_shards - bool, optional
          when True and filename is in a rank_* or node_* directory, discover
          sibling shards and return their combined meshblocks

    returns:
      filedata - dict
          dictionary of fluid file data
    """

    limits = _binary_reader_limits(limits)
    if assemble_shards and _is_partitioned_path(filename):
        return _combine_partitioned_binary(
            filename, _read_binary, "binary", limits
        )

    filedata = {}

    with open(filename, "rb") as fp:
        # load file and get size
        fp.seek(0, 2)
        filesize = fp.tell()
        fp.seek(0, 0)

        # load header information and validate file format
        header_budget = [0]
        code_header = _read_limited_binary_line(
            fp, "binary", "format header", header_budget, limits
        )
        _validate_binary_format_header(code_header, "binary", filename)

        pheader_count = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp, "binary", "preheader count", header_budget, limits
            ),
            "binary",
            "preheader count",
            filename,
        )
        if not 1 <= pheader_count <= _MAX_BINARY_PREHEADER_LINES:
            raise ValueError(f"binary file {filename!r} has invalid preheader count")
        pheader = {}
        retained_preheader_record_bytes = sys.getsizeof(pheader)
        retained_preheader_bytes = _preflight_binary_preheader_text(
            fp, "binary", header_budget, limits, retained_preheader_record_bytes
        )
        for _ in range(pheader_count - 1):
            line = _read_limited_binary_line(
                fp, "binary", "preheader line", header_budget, limits
            )
            key, val = _parse_binary_preheader_text(line, "binary", filename)
            if key in pheader:
                raise ValueError(
                    f"binary file {filename!r} has duplicate preheader metadata "
                    f"for {key!r}"
                )
            retained_preheader_record_bytes = checked_sum(
                (
                    retained_preheader_record_bytes,
                    metadata_record_bytes(key, val),
                ),
                f"binary retained preheader records in {filename!r}",
            )
            retained_preheader_bytes = _preflight_binary_preheader_text(
                fp,
                "binary",
                header_budget,
                limits,
                retained_preheader_record_bytes,
            )
            pheader[key] = val
        time = float(pheader["time"])
        cycle = int(pheader["cycle"])
        locsizebytes = int(pheader["size of location"])
        varsizebytes = int(pheader["size of variable"])
        partition_metadata = _binary_partition_metadata(filename, pheader, "binary")

        nvars = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp, "binary", "variable count", header_budget, limits
            ),
            "binary",
            "variable count",
            filename,
        )
        if nvars <= 0:
            raise ValueError(f"binary file {filename!r} has invalid variable count")
        var_list, retained_variable_bytes = _read_binary_variable_list(
            fp, "binary", header_budget, limits, retained_preheader_bytes
        )
        header_size = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp, "binary", "parameter-header size", header_budget, limits
            ),
            "binary",
            "parameter-header size",
            filename,
        )
        header, retained_parameter_bytes = _read_binary_parameter_dump(
            fp,
            header_size,
            "binary",
            limits,
            retained_preheader_bytes + retained_variable_bytes,
        )
        retained_variable_bytes = checked_sum(
            (
                retained_preheader_bytes,
                retained_variable_bytes,
                retained_parameter_bytes,
            ),
            f"binary retained parameter header in {filename!r}",
        )

        if locsizebytes not in [4, 8]:
            raise ValueError(f"unsupported location size (in bytes) {locsizebytes}")
        if varsizebytes not in [4, 8]:
            raise ValueError(f"unsupported variable size (in bytes) {varsizebytes}")

        # load grid information from header and validate
        def get_from_header(header, blockname, keyname):
            blockname = blockname.strip()
            keyname = keyname.strip()
            if not blockname.startswith("<"):
                blockname = "<" + blockname
            if blockname[-1] != ">":
                blockname += ">"
            block = "<none>"
            for line in [entry for entry in header]:
                if line.startswith("<"):
                    block = line
                    continue
                key, separator, value = line.partition("=")
                if not separator or "=" in value:
                    raise ValueError(
                        f"invalid binary parameter header record {line!r}"
                    )
                if block == blockname and key.strip() == keyname:
                    return value
            raise KeyError(f"no parameter called {blockname}/{keyname}")

        Nx1 = int(get_from_header(header, "<mesh>", "nx1"))
        Nx2 = int(get_from_header(header, "<mesh>", "nx2"))
        Nx3 = int(get_from_header(header, "<mesh>", "nx3"))
        nx1 = int(get_from_header(header, "<meshblock>", "nx1"))
        nx2 = int(get_from_header(header, "<meshblock>", "nx2"))
        nx3 = int(get_from_header(header, "<meshblock>", "nx3"))

        nghost = int(get_from_header(header, "<mesh>", "nghost"))

        x1min = float(get_from_header(header, "<mesh>", "x1min"))
        x1max = float(get_from_header(header, "<mesh>", "x1max"))
        x2min = float(get_from_header(header, "<mesh>", "x2min"))
        x2max = float(get_from_header(header, "<mesh>", "x2max"))
        x3min = float(get_from_header(header, "<mesh>", "x3min"))
        x3max = float(get_from_header(header, "<mesh>", "x3max"))

        bounds = (x1min, x1max, x2min, x2max, x3min, x3max)
        if not np.isfinite(time):
            raise ValueError(f"binary file {filename!r} has non-finite metadata")
        root_meshblock_counts = _validate_binary_grid_metadata(
            (Nx1, Nx2, Nx3), (nx1, nx2, nx3), nghost, bounds, "binary", filename
        )
        if len(var_list) != nvars:
            raise ValueError(
                f"binary variable count mismatch in {filename!r}: "
                f"declared {nvars}, listed {len(var_list)}"
            )
        _require_binary_bytes(
            checked_sum(
                (
                    retained_variable_bytes,
                    _variable_name_set_bytes(var_list),
                ),
                f"binary variable-name duplicate check in {filename!r}",
            ),
            f"binary variable-name duplicate check in {filename!r}",
            limits,
        )
        if len(set(var_list)) != len(var_list):
            raise ValueError(f"binary file {filename!r} has duplicate variable names")
        (
            mb_index,
            mb_logical,
            mb_geometry,
            mb_data,
            retained_variable_bytes,
        ) = _read_meshblocks(
            fp,
            filesize,
            nghost,
            locsizebytes,
            varsizebytes,
            var_list,
            "binary",
            limits,
            retained_variable_bytes,
        )
        _validate_meshblock_metadata(
            mb_logical,
            mb_geometry,
            "binary",
            filename,
            root_meshblock_counts,
            (Nx1, Nx2, Nx3),
            bounds,
            nghost,
            limits,
            _meshblock_arrays_bytes(
                mb_index,
                mb_logical,
                mb_geometry,
                mb_data,
                retained_variable_bytes,
            ),
        )
        mb_count = len(mb_index)
        if (
            partition_metadata["number_of_meshblocks"] is not None
            and partition_metadata["number_of_meshblocks"] != mb_count
        ):
            raise ValueError(
                f"binary file {filename!r} declares "
                f"{partition_metadata['number_of_meshblocks']} meshblocks, "
                f"found {mb_count}"
            )

    filedata["header"] = header
    filedata["time"] = time
    filedata["cycle"] = cycle
    filedata["var_names"] = var_list

    filedata["Nx1"] = Nx1
    filedata["Nx2"] = Nx2
    filedata["Nx3"] = Nx3
    filedata["nvars"] = nvars

    filedata["x1min"] = x1min
    filedata["x1max"] = x1max
    filedata["x2min"] = x2min
    filedata["x2max"] = x2max
    filedata["x3min"] = x3min
    filedata["x3max"] = x3max

    filedata["n_mbs"] = mb_count
    filedata["nx1_mb"] = nx1
    filedata["nx2_mb"] = nx2
    filedata["nx3_mb"] = nx3
    filedata["nx1_out_mb"] = (mb_index[0][1] - mb_index[0][0]) + 1 if mb_count else 0
    filedata["nx2_out_mb"] = (mb_index[0][3] - mb_index[0][2]) + 1 if mb_count else 0
    filedata["nx3_out_mb"] = (mb_index[0][5] - mb_index[0][4]) + 1 if mb_count else 0

    filedata["mb_index"] = mb_index
    filedata["mb_logical"] = mb_logical
    filedata["mb_geometry"] = mb_geometry
    filedata["mb_data"] = mb_data
    filedata["_retained_variable_metadata_bytes"] = retained_variable_bytes
    filedata.update(partition_metadata)

    return filedata


def _read_coarsened_binary(filename, assemble_shards=False, *, limits=None):
    """
    Reads a coarsened bin file from filename to dictionary.
    Originally written by Lev Arzamasskiy (leva@ias.edu) on 11/15/2021
    Updated to support mesh refinement by George Wong (gnwong@ias.edu) on 01/27/2022
    Updated to support coarsened outputs and for speed by Drummond Fielding on 09/09/2024

    args:
      filename - string
          filename of bin file to read
      assemble_shards - bool, optional
          when True and filename is in a rank_* or node_* directory, discover
          sibling shards and return their combined meshblocks

    returns:
      filedata - dict
          dictionary of fluid file data
    """

    limits = _binary_reader_limits(limits)
    if assemble_shards and _is_partitioned_path(filename):
        return _combine_partitioned_binary(
            filename, _read_coarsened_binary, "coarsened binary", limits
        )

    filedata = {}

    with open(filename, "rb") as fp:
        # load file and get size
        fp.seek(0, 2)
        filesize = fp.tell()
        fp.seek(0, 0)

        # load header information and validate file format
        header_budget = [0]
        code_header = _read_limited_binary_line(
            fp, "coarsened binary", "format header", header_budget, limits
        )
        _validate_binary_format_header(code_header, "coarsened binary", filename)

        pheader_count = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp,
                "coarsened binary",
                "preheader count",
                header_budget,
                limits,
            ),
            "coarsened binary",
            "preheader count",
            filename,
        )
        if not 1 <= pheader_count <= _MAX_BINARY_PREHEADER_LINES:
            raise ValueError(
                f"coarsened binary file {filename!r} has invalid preheader count"
            )
        pheader = {}
        retained_preheader_record_bytes = sys.getsizeof(pheader)
        retained_preheader_bytes = _preflight_binary_preheader_text(
            fp,
            "coarsened binary",
            header_budget,
            limits,
            retained_preheader_record_bytes,
        )
        for _ in range(pheader_count - 1):
            line = _read_limited_binary_line(
                fp, "coarsened binary", "preheader line", header_budget, limits
            )
            key, val = _parse_binary_preheader_text(
                line, "coarsened binary", filename
            )
            if key in pheader:
                raise ValueError(
                    f"coarsened binary file {filename!r} has duplicate preheader "
                    f"metadata for {key!r}"
                )
            retained_preheader_record_bytes = checked_sum(
                (
                    retained_preheader_record_bytes,
                    metadata_record_bytes(key, val),
                ),
                f"coarsened binary retained preheader records in {filename!r}",
            )
            retained_preheader_bytes = _preflight_binary_preheader_text(
                fp,
                "coarsened binary",
                header_budget,
                limits,
                retained_preheader_record_bytes,
            )
            pheader[key] = val
        time = float(pheader["time"])
        cycle = int(pheader["cycle"])
        locsizebytes = int(pheader["size of location"])
        varsizebytes = int(pheader["size of variable"])
        coarsen_factor = int(pheader["coarsening factor"])
        if coarsen_factor <= 0:
            raise ValueError(
                f"coarsened binary file {filename!r} has invalid coarsening factor"
            )
        partition_metadata = _binary_partition_metadata(
            filename, pheader, "coarsened binary"
        )

        nvars = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp, "coarsened binary", "variable count", header_budget, limits
            ),
            "coarsened binary",
            "variable count",
            filename,
        )
        if nvars <= 0:
            raise ValueError(
                f"coarsened binary file {filename!r} has invalid variable count"
            )
        var_list, retained_variable_bytes = _read_binary_variable_list(
            fp,
            "coarsened binary",
            header_budget,
            limits,
            retained_preheader_bytes,
        )
        if len(var_list) != nvars:
            raise ValueError(
                f"coarsened binary variable count mismatch in {filename!r}: "
                f"declared {nvars}, listed {len(var_list)}"
            )
        _require_binary_bytes(
            checked_sum(
                (
                    retained_preheader_bytes,
                    retained_variable_bytes,
                    _variable_name_set_bytes(var_list),
                ),
                f"coarsened binary variable-name duplicate check in {filename!r}",
            ),
            f"coarsened binary variable-name duplicate check in {filename!r}",
            limits,
        )
        if len(set(var_list)) != len(var_list):
            raise ValueError(
                f"coarsened binary file {filename!r} has duplicate variable names"
            )
        number_of_moments = _validate_coarsened_moments(
            filename, pheader, nvars, var_list
        )
        header_size = _parse_binary_metadata_int(
            _read_limited_binary_line(
                fp,
                "coarsened binary",
                "parameter-header size",
                header_budget,
                limits,
            ),
            "coarsened binary",
            "parameter-header size",
            filename,
        )
        header, retained_parameter_bytes = _read_binary_parameter_dump(
            fp,
            header_size,
            "coarsened binary",
            limits,
            retained_preheader_bytes + retained_variable_bytes,
        )
        retained_variable_bytes = checked_sum(
            (
                retained_preheader_bytes,
                retained_variable_bytes,
                retained_parameter_bytes,
            ),
            f"coarsened binary retained parameter header in {filename!r}",
        )

        if locsizebytes not in [4, 8]:
            raise ValueError(f"unsupported location size (in bytes) {locsizebytes}")
        if varsizebytes not in [4, 8]:
            raise ValueError(f"unsupported variable size (in bytes) {varsizebytes}")

        # load grid information from header and validate
        def get_from_header(header, blockname, keyname):
            blockname = blockname.strip()
            keyname = keyname.strip()
            if not blockname.startswith("<"):
                blockname = "<" + blockname
            if blockname[-1] != ">":
                blockname += ">"
            block = "<none>"
            for line in [entry for entry in header]:
                if line.startswith("<"):
                    block = line
                    continue
                key, separator, value = line.partition("=")
                if not separator or "=" in value:
                    raise ValueError(
                        f"invalid coarsened binary parameter header record {line!r}"
                    )
                if block == blockname and key.strip() == keyname:
                    return value
            raise KeyError(f"no parameter called {blockname}/{keyname}")

        Nx1 = int(get_from_header(header, "<mesh>", "nx1"))
        Nx2 = int(get_from_header(header, "<mesh>", "nx2"))
        Nx3 = int(get_from_header(header, "<mesh>", "nx3"))
        nx1 = int(get_from_header(header, "<meshblock>", "nx1"))
        nx2 = int(get_from_header(header, "<meshblock>", "nx2"))
        nx3 = int(get_from_header(header, "<meshblock>", "nx3"))

        nghost = int(get_from_header(header, "<mesh>", "nghost"))

        x1min = float(get_from_header(header, "<mesh>", "x1min"))
        x1max = float(get_from_header(header, "<mesh>", "x1max"))
        x2min = float(get_from_header(header, "<mesh>", "x2min"))
        x2max = float(get_from_header(header, "<mesh>", "x2max"))
        x3min = float(get_from_header(header, "<mesh>", "x3min"))
        x3max = float(get_from_header(header, "<mesh>", "x3max"))

        bounds = (x1min, x1max, x2min, x2max, x3min, x3max)
        if not np.isfinite(time):
            raise ValueError(
                f"coarsened binary file {filename!r} has non-finite metadata"
            )
        root_meshblock_counts = _validate_binary_grid_metadata(
            (Nx1, Nx2, Nx3),
            (nx1, nx2, nx3),
            nghost,
            bounds,
            "coarsened binary",
            filename,
        )
        _validate_coarsening_factor(
            coarsen_factor, (Nx1, Nx2, Nx3, nx1, nx2, nx3), filename
        )
        (
            mb_index,
            mb_logical,
            mb_geometry,
            mb_data,
            retained_variable_bytes,
        ) = _read_meshblocks(
            fp,
            filesize,
            nghost,
            locsizebytes,
            varsizebytes,
            var_list,
            "coarsened binary",
            limits,
            retained_variable_bytes,
        )
        _validate_meshblock_metadata(
            mb_logical,
            mb_geometry,
            "coarsened binary",
            filename,
            root_meshblock_counts,
            (Nx1, Nx2, Nx3),
            bounds,
            nghost,
            limits,
            _meshblock_arrays_bytes(
                mb_index,
                mb_logical,
                mb_geometry,
                mb_data,
                retained_variable_bytes,
            ),
        )
        mb_count = len(mb_index)
        if (
            partition_metadata["number_of_meshblocks"] is not None
            and partition_metadata["number_of_meshblocks"] != mb_count
        ):
            raise ValueError(
                f"coarsened binary file {filename!r} declares "
                f"{partition_metadata['number_of_meshblocks']} meshblocks, "
                f"found {mb_count}"
            )

    filedata["header"] = header
    filedata["time"] = time
    filedata["cycle"] = cycle
    filedata["var_names"] = var_list

    filedata["Nx1"] = Nx1 // coarsen_factor
    filedata["Nx2"] = Nx2 // coarsen_factor
    filedata["Nx3"] = Nx3 // coarsen_factor
    filedata["nvars"] = nvars
    filedata["number_of_moments"] = number_of_moments

    filedata["x1min"] = x1min
    filedata["x1max"] = x1max
    filedata["x2min"] = x2min
    filedata["x2max"] = x2max
    filedata["x3min"] = x3min
    filedata["x3max"] = x3max

    filedata["n_mbs"] = mb_count
    filedata["nx1_mb"] = nx1 // coarsen_factor
    filedata["nx2_mb"] = nx2 // coarsen_factor
    filedata["nx3_mb"] = nx3 // coarsen_factor
    filedata["nx1_out_mb"] = (mb_index[0][1] - mb_index[0][0]) + 1 if mb_count else 0
    filedata["nx2_out_mb"] = (mb_index[0][3] - mb_index[0][2]) + 1 if mb_count else 0
    filedata["nx3_out_mb"] = (mb_index[0][5] - mb_index[0][4]) + 1 if mb_count else 0

    filedata["mb_index"] = mb_index
    filedata["mb_logical"] = mb_logical
    filedata["mb_geometry"] = mb_geometry
    filedata["mb_data"] = mb_data
    filedata["_retained_variable_metadata_bytes"] = retained_variable_bytes
    filedata.update(partition_metadata)

    return filedata


def _public_binary_data(filedata):
    """Remove parser bookkeeping before returning one public reader result."""
    filedata.pop("_retained_variable_metadata_bytes", None)
    return filedata


def read_binary(filename, assemble_shards=False, *, limits=None):
    """Read one binary output without exposing parser accounting fields."""
    return _public_binary_data(
        _read_binary(filename, assemble_shards=assemble_shards, limits=limits)
    )


def read_coarsened_binary(filename, assemble_shards=False, *, limits=None):
    """Read one coarsened-binary output without parser accounting fields."""
    return _public_binary_data(
        _read_coarsened_binary(filename, assemble_shards=assemble_shards, limits=limits)
    )


def _read_all_ranks_binary(rank0_filename, *, limits=None):
    return _combine_partitioned_binary(
        rank0_filename, _read_binary, "binary", limits
    )


def _read_all_ranks_coarsened_binary(rank0_filename, *, limits=None):
    return _combine_partitioned_binary(
        rank0_filename, _read_coarsened_binary, "coarsened binary", limits
    )


def read_all_ranks_binary(rank0_filename, *, limits=None):
    """
    Reads binary files from all rank or node shards into a single dictionary.

    args:
      rank0_filename - string
          filename of any rank/node shard, or a shared binary file

    returns:
      combined_filedata - dict
          dictionary of combined fluid file data from all ranks
    """
    return _public_binary_data(_read_all_ranks_binary(rank0_filename, limits=limits))


def read_all_ranks_coarsened_binary(rank0_filename, *, limits=None):
    """
    Reads coarsened binary files from all rank or node shards.

    args:
      rank0_filename - string
          filename of any rank/node shard, or a shared coarsened binary file

    returns:
      combined_filedata - dict
          dictionary of combined fluid file data from all ranks
    """
    return _public_binary_data(
        _read_all_ranks_coarsened_binary(rank0_filename, limits=limits)
    )


def read_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
    *,
    limits=None,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    limits = _athdf_reader_limits(limits)
    filedata = _read_binary(filename, limits=limits)

    # Step 2: Organize data similar to athdf
    if raw:
        return _public_binary_data(filedata)

    _require_meshblocks_for_athdf(filedata, filename)
    retained_binary_bytes = _binary_retained_bytes(filedata)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = int(max(filedata["mb_logical"][:, 3]))
    if level is None:
        level = max_level
    level = _validate_requested_level(level)
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32
    num_ghost = _validate_athdf_num_ghost(
        filedata,
        block_size,
        root_grid_size,
        num_ghost,
        max_level,
        level,
        (x1_min, x1_max, x2_min, x2_max, x3_min, x3_max),
    )

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
            if _is_effective_slice(
                levels, logical_locations, d, retained_binary_bytes, limits
            ):
                nx_vals.append(1)
            else:  # nontrivial sum
                num_blocks_this_dim = 0
                for level_this_dim, loc_this_dim in zip(
                    levels, logical_locations[:, d]
                ):
                    if level_this_dim <= level:
                        possible_max = (loc_this_dim + 1) * 2 ** (
                            level - level_this_dim
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    else:
                        possible_max = (loc_this_dim + 1) // 2 ** (
                            level_this_dim - level
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                nx_vals.append(num_blocks_this_dim)
        elif block_size[d] == 1:  # singleton dimension
            nx_vals.append(1)
        else:  # normal case
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]
    coordinate_bytes = _preflight_athdf_coordinates(
        nx_vals, dtype, limits, retained_binary_bytes
    )

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]
    _populate_root_athdf_coordinates(
        data, filedata, nx_vals, root_grid_size, level, num_ghost, dtype, center_funcs
    )

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max, j_min, j_max, k_min, k_max = _athdf_selection(
        data,
        (nx1, nx2, nx3),
        (x1_min, x2_min, x3_min),
        (x1_max, x2_max, x3_max),
    )

    _crop_athdf_coordinates(data, ((i_min, i_max), (j_min, j_max), (k_min, k_max)))

    # Prepare arrays for data and bookkeeping
    restricted_shape = (
        (lx3, lx2, lx1)
        if not subsample and not fast_restrict and max_level > level
        else None
    )
    output_bytes = _preflight_athdf_outputs(
        (k_max - k_min, j_max - j_min, i_max - i_min),
        quantities,
        dtype,
        return_levels,
        restricted_shape,
        coordinate_bytes,
        limits,
    )
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.full(
                (k_max - k_min, j_max - j_min, i_max - i_min), -1, dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
        if return_levels:
            data["Levels"].fill(-1)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        _copy_meshblock_to_athdf(
            data,
            filedata,
            quantities,
            block_num,
            level,
            block_size,
            (nx1, nx2, nx3),
            (i_min, i_max, j_min, j_max, k_min, k_max),
            output_bytes,
            return_levels,
            subsample,
            fast_restrict,
            restricted_data if restricted_shape is not None else None,
            vol_func,
            num_ghost,
            limits,
        )

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
        for loc3 in range(lx3):
            for loc2 in range(lx2):
                for loc1 in range(lx1):
                    if restricted_data[loc3, loc2, loc1]:
                        il = loc1 * block_size[0]
                        jl = loc2 * block_size[1]
                        kl = loc3 * block_size[2]
                        iu = il + block_size[0]
                        ju = jl + block_size[1]
                        ku = kl + block_size[2]
                        il = max(il, i_min) - i_min
                        jl = max(jl, j_min) - j_min
                        kl = max(kl, k_min) - k_min
                        iu = min(iu, i_max) - i_min
                        ju = min(ju, j_max) - j_min
                        ku = min(ku, k_max) - k_min
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_rank_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
    *,
    limits=None,
):
    """Read one binary shard through the canonical logical-location mapper."""
    return read_binary_as_athdf(**locals())


def read_all_ranks_binary_as_athdf(
    rank0_filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
    *,
    limits=None,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    limits = _athdf_reader_limits(limits)
    filedata = _read_all_ranks_binary(rank0_filename, limits=limits)

    # Step 2: Organize data similar to athdf
    if raw:
        return _public_binary_data(filedata)

    _require_meshblocks_for_athdf(filedata, rank0_filename)
    retained_binary_bytes = _binary_retained_bytes(filedata)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = int(max(filedata["mb_logical"][:, 3]))
    if level is None:
        level = max_level
    level = _validate_requested_level(level)
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    levels = filedata["mb_logical"][:, 3]
    logical_locations = filedata["mb_logical"][:, :3]
    if dtype is None:
        dtype = np.float32
    num_ghost = _validate_athdf_num_ghost(
        filedata,
        block_size,
        root_grid_size,
        num_ghost,
        max_level,
        level,
        (x1_min, x1_max, x2_min, x2_max, x3_min, x3_max),
    )

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:  # sum or slice
            if _is_effective_slice(
                levels, logical_locations, d, retained_binary_bytes, limits
            ):
                nx_vals.append(1)
            else:  # nontrivial sum
                num_blocks_this_dim = 0
                for level_this_dim, loc_this_dim in zip(
                    levels, logical_locations[:, d]
                ):
                    if level_this_dim <= level:
                        possible_max = (loc_this_dim + 1) * 2 ** (
                            level - level_this_dim
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                    else:
                        possible_max = (loc_this_dim + 1) // 2 ** (
                            level_this_dim - level
                        )
                        num_blocks_this_dim = max(num_blocks_this_dim, possible_max)
                nx_vals.append(num_blocks_this_dim)
        elif block_size[d] == 1:  # singleton dimension
            nx_vals.append(1)
        else:  # normal case
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]
    coordinate_bytes = _preflight_athdf_coordinates(
        nx_vals, dtype, limits, retained_binary_bytes
    )
    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]
    _populate_root_athdf_coordinates(
        data, filedata, nx_vals, root_grid_size, level, num_ghost, dtype, center_funcs
    )

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max, j_min, j_max, k_min, k_max = _athdf_selection(
        data,
        (nx1, nx2, nx3),
        (x1_min, x2_min, x3_min),
        (x1_max, x2_max, x3_max),
    )

    _crop_athdf_coordinates(data, ((i_min, i_max), (j_min, j_max), (k_min, k_max)))

    # Prepare arrays for data and bookkeeping
    restricted_shape = (
        (lx3, lx2, lx1)
        if not subsample and not fast_restrict and max_level > level
        else None
    )
    output_bytes = _preflight_athdf_outputs(
        (k_max - k_min, j_max - j_min, i_max - i_min),
        quantities,
        dtype,
        return_levels,
        restricted_shape,
        coordinate_bytes,
        limits,
    )
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.full(
                (k_max - k_min, j_max - j_min, i_max - i_min), -1, dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
        if return_levels:
            data["Levels"].fill(-1)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        _copy_meshblock_to_athdf(
            data,
            filedata,
            quantities,
            block_num,
            level,
            block_size,
            (nx1, nx2, nx3),
            (i_min, i_max, j_min, j_max, k_min, k_max),
            output_bytes,
            return_levels,
            subsample,
            fast_restrict,
            restricted_data if restricted_shape is not None else None,
            vol_func,
            num_ghost,
            limits,
        )

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
        for loc3 in range(lx3):
            for loc2 in range(lx2):
                for loc1 in range(lx1):
                    if restricted_data[loc3, loc2, loc1]:
                        il = loc1 * block_size[0]
                        jl = loc2 * block_size[1]
                        kl = loc3 * block_size[2]
                        iu = il + block_size[0]
                        ju = jl + block_size[1]
                        ku = kl + block_size[2]
                        il = max(il, i_min) - i_min
                        jl = max(jl, j_min) - j_min
                        kl = max(kl, k_min) - k_min
                        iu = min(iu, i_max) - i_min
                        ju = min(ju, j_max) - j_min
                        ku = min(ku, k_max) - k_min
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_all_ranks_coarsened_binary_as_athdf(
    rank0_filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
    *,
    limits=None,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    limits = _athdf_reader_limits(limits)
    filedata = _read_all_ranks_coarsened_binary(rank0_filename, limits=limits)

    # Step 2: Organize data similar to athdf
    if raw:
        return _public_binary_data(filedata)

    _require_meshblocks_for_athdf(filedata, rank0_filename)
    retained_binary_bytes = _binary_retained_bytes(filedata)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = int(max(filedata["mb_logical"][:, 3]))
    if level is None:
        level = max_level
    level = _validate_requested_level(level)
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    if dtype is None:
        dtype = np.float32
    num_ghost = _validate_athdf_num_ghost(
        filedata,
        block_size,
        root_grid_size,
        num_ghost,
        max_level,
        level,
        (x1_min, x1_max, x2_min, x2_max, x3_min, x3_max),
    )

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:
            # Implement logic for sum or slice as in athdf
            nx_vals.append(root_grid_size[d] * 2**level)
        elif block_size[d] == 1:
            nx_vals.append(1)
        else:
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]
    coordinate_bytes = _preflight_athdf_coordinates(
        nx_vals, dtype, limits, retained_binary_bytes
    )

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]
    _populate_root_athdf_coordinates(
        data, filedata, nx_vals, root_grid_size, level, num_ghost, dtype, center_funcs
    )

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max, j_min, j_max, k_min, k_max = _athdf_selection(
        data,
        (nx1, nx2, nx3),
        (x1_min, x2_min, x3_min),
        (x1_max, x2_max, x3_max),
    )

    _crop_athdf_coordinates(data, ((i_min, i_max), (j_min, j_max), (k_min, k_max)))

    # Prepare arrays for data and bookkeeping
    restricted_shape = (
        (lx3, lx2, lx1)
        if not subsample and not fast_restrict and max_level > level
        else None
    )
    output_bytes = _preflight_athdf_outputs(
        (k_max - k_min, j_max - j_min, i_max - i_min),
        quantities,
        dtype,
        return_levels,
        restricted_shape,
        coordinate_bytes,
        limits,
    )
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.full(
                (k_max - k_min, j_max - j_min, i_max - i_min), -1, dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
        if return_levels:
            data["Levels"].fill(-1)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        _copy_meshblock_to_athdf(
            data,
            filedata,
            quantities,
            block_num,
            level,
            block_size,
            (nx1, nx2, nx3),
            (i_min, i_max, j_min, j_max, k_min, k_max),
            output_bytes,
            return_levels,
            subsample,
            fast_restrict,
            restricted_data if restricted_shape is not None else None,
            vol_func,
            num_ghost,
            limits,
        )

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
        for loc3 in range(lx3):
            for loc2 in range(lx2):
                for loc1 in range(lx1):
                    if restricted_data[loc3, loc2, loc1]:
                        il = loc1 * block_size[0]
                        jl = loc2 * block_size[1]
                        kl = loc3 * block_size[2]
                        iu = il + block_size[0]
                        ju = jl + block_size[1]
                        ku = kl + block_size[2]
                        il = max(il, i_min) - i_min
                        jl = max(jl, j_min) - j_min
                        kl = max(kl, k_min) - k_min
                        iu = min(iu, i_max) - i_min
                        ju = min(ju, j_max) - j_min
                        ku = min(ku, k_max) - k_min
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def read_single_rank_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    return_levels=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    *,
    meshblock_index_in_file=0,
    limits=None,
):
    """
    Reads a single rank binary file and organizes data similar to
    athdf format without writing to file.
    """
    # Step 1: Read binary data for a single rank
    limits = _athdf_reader_limits(limits)
    filedata = _read_binary(filename, limits=limits)

    if isinstance(meshblock_index_in_file, (bool, np.bool_)) or not isinstance(
        meshblock_index_in_file, Integral
    ):
        raise TypeError("meshblock_index_in_file must be an integer")
    meshblock_index_in_file = int(meshblock_index_in_file)

    if raw:
        if meshblock_index_in_file != 0:
            raise ValueError("meshblock_index_in_file must be 0 when raw=True")
        return _public_binary_data(filedata)

    _require_meshblocks_for_athdf(filedata, filename)
    if not 0 <= meshblock_index_in_file < filedata["n_mbs"]:
        raise IndexError(
            f"meshblock_index_in_file {meshblock_index_in_file} is out of range "
            f"for {filedata['n_mbs']} meshblocks"
        )

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    block_size = [filedata["nx1_mb"], filedata["nx2_mb"], filedata["nx3_mb"]]
    emitted_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    if emitted_size != block_size:
        raise ValueError(
            "single-meshblock athdf-like conversion does not support ghost-bearing "
            "or sliced binary extents"
        )
    retained_binary_bytes = _binary_retained_bytes(filedata)
    if dtype is None:
        dtype = np.float32
    coordinate_bytes = _preflight_athdf_coordinates(
        block_size, dtype, limits, retained_binary_bytes
    )

    # Set coordinate system and related functions
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    # Populate coordinate arrays
    center_funcs = [center_func_1, center_func_2, center_func_3]
    for d in range(1, 4):
        xf = f"x{d}f"
        xv = f"x{d}v"
        nx = block_size[d - 1]

        # Use the meshblock geometry for local min and max
        xmin = filedata["mb_geometry"][meshblock_index_in_file, (d - 1) * 2]
        xmax = filedata["mb_geometry"][meshblock_index_in_file, (d - 1) * 2 + 1]

        data[xf] = np.linspace(xmin, xmax, nx + 1, dtype=dtype)
        data[xv] = np.empty(nx, dtype=dtype)
        for i in range(nx):
            data[xv][i] = center_funcs[d - 1](data[xf][i], data[xf][i + 1])

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max, j_min, j_max, k_min, k_max = _athdf_selection(
        data,
        block_size,
        (x1_min, x2_min, x3_min),
        (x1_max, x2_max, x3_max),
    )

    _crop_athdf_coordinates(data, ((i_min, i_max), (j_min, j_max), (k_min, k_max)))

    # Prepare arrays for data
    _preflight_athdf_outputs(
        (k_max - k_min, j_max - j_min, i_max - i_min),
        quantities,
        dtype,
        return_levels,
        None,
        coordinate_bytes,
        limits,
    )
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.full(
                (k_max - k_min, j_max - j_min, i_max - i_min), -1, dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)

    # Process the single block
    for q in quantities:
        block_data = filedata["mb_data"][q][meshblock_index_in_file]
        data[q][...] = block_data[k_min:k_max, j_min:j_max, i_min:i_max]

    if return_levels:
        data["Levels"].fill(filedata["mb_logical"][meshblock_index_in_file, 3])

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = filedata["mb_logical"][meshblock_index_in_file, 3]

    return data


def read_coarsened_binary_as_athdf(
    filename,
    raw=False,
    data=None,
    quantities=None,
    dtype=None,
    level=None,
    return_levels=False,
    subsample=False,
    fast_restrict=False,
    x1_min=None,
    x1_max=None,
    x2_min=None,
    x2_max=None,
    x3_min=None,
    x3_max=None,
    vol_func=None,
    vol_params=None,
    face_func_1=None,
    face_func_2=None,
    face_func_3=None,
    center_func_1=None,
    center_func_2=None,
    center_func_3=None,
    num_ghost=0,
    *,
    limits=None,
):
    """
    Reads a bin file and organizes data similar to athdf format without writing to file.
    """
    # Step 1: Read binary data
    limits = _athdf_reader_limits(limits)
    filedata = _read_coarsened_binary(filename, limits=limits)

    # Step 2: Organize data similar to athdf
    if raw:
        return _public_binary_data(filedata)

    _require_meshblocks_for_athdf(filedata, filename)
    retained_binary_bytes = _binary_retained_bytes(filedata)

    # Prepare dictionary for results
    if data is None:
        data = {}
        new_data = True
    else:
        new_data = False

    # Extract size information
    max_level = int(max(filedata["mb_logical"][:, 3]))
    if level is None:
        level = max_level
    level = _validate_requested_level(level)
    block_size = [
        filedata["nx1_out_mb"],
        filedata["nx2_out_mb"],
        filedata["nx3_out_mb"],
    ]
    root_grid_size = [filedata["Nx1"], filedata["Nx2"], filedata["Nx3"]]
    if dtype is None:
        dtype = np.float32
    num_ghost = _validate_athdf_num_ghost(
        filedata,
        block_size,
        root_grid_size,
        num_ghost,
        max_level,
        level,
        (x1_min, x1_max, x2_min, x2_max, x3_min, x3_max),
    )

    # Calculate nx_vals
    nx_vals = []
    for d in range(3):
        if block_size[d] == 1 and root_grid_size[d] > 1:
            # Implement logic for sum or slice as in athdf
            nx_vals.append(root_grid_size[d] * 2**level)
        elif block_size[d] == 1:
            nx_vals.append(1)
        else:
            nx_vals.append(root_grid_size[d] * 2**level + 2 * num_ghost)
    nx1, nx2, nx3 = nx_vals
    lx1, lx2, lx3 = [nx // bs for nx, bs in zip(nx_vals, block_size)]
    coordinate_bytes = _preflight_athdf_coordinates(
        nx_vals, dtype, limits, retained_binary_bytes
    )

    # Set coordinate system and related functions
    # coord = "cartesian"  # Adjust based on your data
    if vol_func is None:

        def vol_func(xm, xp, ym, yp, zm, zp):
            return (xp - xm) * (yp - ym) * (zp - zm)

    # Define center functions if not provided
    if center_func_1 is None:

        def center_func_1(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_2 is None:

        def center_func_2(xm, xp):
            return 0.5 * (xm + xp)

    if center_func_3 is None:

        def center_func_3(xm, xp):
            return 0.5 * (xm + xp)

    center_funcs = [center_func_1, center_func_2, center_func_3]
    _populate_root_athdf_coordinates(
        data, filedata, nx_vals, root_grid_size, level, num_ghost, dtype, center_funcs
    )

    # Create list of quantities
    if quantities is None:
        quantities = filedata["var_names"]

    # Account for selection
    i_min, i_max, j_min, j_max, k_min, k_max = _athdf_selection(
        data,
        (nx1, nx2, nx3),
        (x1_min, x2_min, x3_min),
        (x1_max, x2_max, x3_max),
    )

    _crop_athdf_coordinates(data, ((i_min, i_max), (j_min, j_max), (k_min, k_max)))

    # Prepare arrays for data and bookkeeping
    restricted_shape = (
        (lx3, lx2, lx1)
        if not subsample and not fast_restrict and max_level > level
        else None
    )
    output_bytes = _preflight_athdf_outputs(
        (k_max - k_min, j_max - j_min, i_max - i_min),
        quantities,
        dtype,
        return_levels,
        restricted_shape,
        coordinate_bytes,
        limits,
    )
    if new_data:
        for q in quantities:
            data[q] = np.zeros(
                (k_max - k_min, j_max - j_min, i_max - i_min), dtype=dtype
            )
        if return_levels:
            data["Levels"] = np.full(
                (k_max - k_min, j_max - j_min, i_max - i_min), -1, dtype=np.int32
            )
    else:
        for q in quantities:
            data[q].fill(0.0)
        if return_levels:
            data["Levels"].fill(-1)
    if not subsample and not fast_restrict and max_level > level:
        restricted_data = np.zeros((lx3, lx2, lx1), dtype=bool)

    # Step 3: Process each block
    for block_num in range(filedata["n_mbs"]):
        _copy_meshblock_to_athdf(
            data,
            filedata,
            quantities,
            block_num,
            level,
            block_size,
            (nx1, nx2, nx3),
            (i_min, i_max, j_min, j_max, k_min, k_max),
            output_bytes,
            return_levels,
            subsample,
            fast_restrict,
            restricted_data if restricted_shape is not None else None,
            vol_func,
            num_ghost,
            limits,
        )

    # Step 4: Finalize data
    if level < max_level and not subsample and not fast_restrict:
        # Remove volume factors from restricted data
        for loc3 in range(lx3):
            for loc2 in range(lx2):
                for loc1 in range(lx1):
                    if restricted_data[loc3, loc2, loc1]:
                        il = loc1 * block_size[0]
                        jl = loc2 * block_size[1]
                        kl = loc3 * block_size[2]
                        iu = il + block_size[0]
                        ju = jl + block_size[1]
                        ku = kl + block_size[2]
                        il = max(il, i_min) - i_min
                        jl = max(jl, j_min) - j_min
                        kl = max(kl, k_min) - k_min
                        iu = min(iu, i_max) - i_min
                        ju = min(ju, j_max) - j_min
                        ku = min(ku, k_max) - k_min
                        for k in range(kl, ku):
                            for j in range(jl, ju):
                                for i in range(il, iu):
                                    x1m, x1p = data["x1f"][i], data["x1f"][i + 1]
                                    x2m, x2p = data["x2f"][j], data["x2f"][j + 1]
                                    x3m, x3p = data["x3f"][k], data["x3f"][k + 1]
                                    vol = vol_func(x1m, x1p, x2m, x2p, x3m, x3p)
                                    for q in quantities:
                                        data[q][k, j, i] /= vol

    # Add metadata
    data["Time"] = filedata["time"]
    data["NumCycles"] = filedata["cycle"]
    data["MaxLevel"] = max_level

    return data


def _preflight_write_athdf(
    fdata, vars_without_b, vars_only_b, number_of_moments, limits
):
    """Bound converter outputs and one coordinate-generation working set."""
    float_bytes = np.dtype(np.float64).itemsize
    nmb = fdata["n_mbs"]
    nx1 = fdata["nx1_out_mb"]
    nx2 = fdata["nx2_out_mb"]
    nx3 = fdata["nx3_out_mb"]
    output_counts = [
        nmb,
        _checked_binary_product((nmb, 3), "ATHDF logical locations"),
        _checked_binary_product(
            (len(vars_without_b), nmb, nx3, nx2, nx1), "ATHDF uov"
        ),
    ]
    if vars_only_b:
        output_counts.append(
            _checked_binary_product(
                (3 * number_of_moments, nmb, nx3, nx2, nx1),
                "ATHDF magnetic fields",
            )
        )
    output_counts.extend(
        _checked_binary_product(
            (nmb, 2 * size + 1), f"ATHDF x{axis} coordinates"
        )
        for axis, size in enumerate((nx1, nx2, nx3), start=1)
    )
    output_bytes = _checked_binary_product(
        (checked_sum(output_counts, "ATHDF output arrays"), float_bytes),
        "ATHDF output arrays",
    )
    retained_input_bytes = _binary_retained_bytes(fdata)
    _require_athdf_bytes(
        checked_sum(
            (retained_input_bytes, output_bytes),
            "ATHDF converter materialization peak",
        ),
        "ATHDF converter materialization peak",
        limits,
    )
    temporary_count = checked_sum(
        (
            3 * fdata["nx1_mb"] + 1,
            3 * fdata["nx2_mb"] + 1,
            3 * fdata["nx3_mb"] + 1,
        ),
        "ATHDF coordinate working set",
    )
    _require_athdf_bytes(
        checked_sum(
            (
                retained_input_bytes,
                output_bytes,
                temporary_count * float_bytes,
            ),
            "ATHDF coordinate generation peak",
        ),
        "ATHDF coordinate generation peak",
        limits,
    )


def write_athdf(
    filename, fdata, varsize_bytes=4, locsize_bytes=8, *, limits=None
):
    """
    Writes an athdf (hdf5) file from a loaded python filedata object.

    args:
      filename      - string
          filename for output athdf (hdf5) file
      fdata         - dict
          dictionary of fluid file data, e.g., as loaded from read_binary(...)
      varsize_bytes - int (default=4, options=4,8)
          number of bytes to use for output variable data
      locsize_bytes - int (default=8, options=4,8)
          number of bytes to use for output location data
    """

    limits = _athdf_reader_limits(limits)
    if varsize_bytes not in [4, 8]:
        raise ValueError(f"varsizebytes must be 4 or 8, not {varsize_bytes}")
    if locsize_bytes not in [4, 8]:
        raise ValueError(f"locsizebytes must be 4 or 8, not {locsize_bytes}")
    locfmt = "<f4" if locsize_bytes == 4 else "<f8"
    varfmt = "<f4" if varsize_bytes == 4 else "<f8"

    # extract Mesh/MeshBlock parameters
    nmb = fdata["n_mbs"]
    Nx1 = fdata["Nx1"]  # noqa: F841
    Nx2 = fdata["Nx2"]
    Nx3 = fdata["Nx3"]
    nx1 = fdata["nx1_mb"]
    nx2 = fdata["nx2_mb"]
    nx3 = fdata["nx3_mb"]
    nx1_out = fdata["nx1_out_mb"]
    nx2_out = fdata["nx2_out_mb"]
    nx3_out = fdata["nx3_out_mb"]

    number_of_moments = fdata.get("number_of_moments", 1)

    # check dimensionality/slicing
    two_d = Nx2 != 1 and Nx3 == 1
    three_d = Nx3 != 1
    x1slice = nx1_out == 1
    x2slice = nx2_out == 1 and (two_d or three_d)
    x3slice = nx3_out == 1 and three_d

    # keep variable order but separate out magnetic field
    vars_without_b = [v for v in fdata["var_names"] if "bcc" not in v]
    vars_only_b = [v for v in fdata["var_names"] if v not in vars_without_b]
    _preflight_write_athdf(
        fdata, vars_without_b, vars_only_b, number_of_moments, limits
    )

    if len(vars_only_b) > 0:
        B = np.zeros((3 * number_of_moments, nmb, nx3_out, nx2_out, nx1_out))
    Levels = np.zeros(nmb)
    LogicalLocations = np.zeros((nmb, 3))
    uov = np.zeros((len(vars_without_b), nmb, nx3_out, nx2_out, nx1_out))
    x1f = np.zeros((nmb, nx1_out + 1))
    x1v = np.zeros((nmb, nx1_out))
    x2f = np.zeros((nmb, nx2_out + 1))
    x2v = np.zeros((nmb, nx2_out))
    x3f = np.zeros((nmb, nx3_out + 1))
    x3v = np.zeros((nmb, nx3_out))

    for ivar, var in enumerate(vars_without_b):
        uov[ivar] = fdata["mb_data"][var]
    for ibvar, bvar in enumerate(vars_only_b):
        B[ibvar] = fdata["mb_data"][bvar]

    for mb in range(nmb):
        logical = fdata["mb_logical"][mb]
        LogicalLocations[mb] = logical[:3]
        Levels[mb] = logical[-1]
        geometry = fdata["mb_geometry"][mb]
        mb_x1f = np.linspace(geometry[0], geometry[1], nx1 + 1)
        mb_x1v = 0.5 * (mb_x1f[1:] + mb_x1f[:-1])
        mb_x2f = np.linspace(geometry[2], geometry[3], nx2 + 1)
        mb_x2v = 0.5 * (mb_x2f[1:] + mb_x2f[:-1])
        mb_x3f = np.linspace(geometry[4], geometry[5], nx3 + 1)
        mb_x3v = 0.5 * (mb_x3f[1:] + mb_x3f[:-1])
        if x1slice:
            x1f[mb] = np.array(
                mb_x1f[(fdata["mb_index"][mb][0]): (fdata["mb_index"][mb][0] + 2)]
            )
            x1v[mb] = np.array([np.average(mb_x1f)])
        else:
            x1f[mb] = mb_x1f
            x1v[mb] = mb_x1v
        if x2slice:
            x2f[mb] = np.array(
                mb_x2f[(fdata["mb_index"][mb][2]): (fdata["mb_index"][mb][2] + 2)]
            )
            x2v[mb] = np.array([np.average(x2f[mb])])
        else:
            x2f[mb] = mb_x2f
            x2v[mb] = mb_x2v
        if x3slice:
            x3f[mb] = np.array(
                mb_x3f[(fdata["mb_index"][mb][4]): (fdata["mb_index"][mb][4] + 2)]
            )
            x3v[mb] = np.array([np.average(x3f[mb])])
        else:
            x3f[mb] = mb_x3f
            x3v[mb] = mb_x3v

    # set dataset names and number of variables
    dataset_names = [np.array("uov", dtype="|S21")]
    dataset_nvars = [len(vars_without_b)]
    if len(vars_only_b) > 0:
        dataset_names.append(np.array("B", dtype="|S21"))
        dataset_nvars.append(len(vars_only_b))

    # Set Attributes
    hfp = h5py.File(filename, "w")
    hfp.attrs["Header"] = fdata["header"]
    hfp.attrs["Time"] = fdata["time"]
    hfp.attrs["NumCycles"] = fdata["cycle"]
    hfp.attrs["Coordinates"] = np.array("cartesian", dtype="|S11")
    hfp.attrs["NumMeshBlocks"] = fdata["n_mbs"]
    hfp.attrs["MaxLevel"] = int(max(Levels))
    hfp.attrs["MeshBlockSize"] = [
        fdata["nx1_out_mb"],
        fdata["nx2_out_mb"],
        fdata["nx3_out_mb"],
    ]
    hfp.attrs["RootGridSize"] = [fdata["Nx1"], fdata["Nx2"], fdata["Nx3"]]
    hfp.attrs["RootGridX1"] = [fdata["x1min"], fdata["x1max"], 1.0]
    hfp.attrs["RootGridX2"] = [fdata["x2min"], fdata["x2max"], 1.0]
    hfp.attrs["RootGridX3"] = [fdata["x3min"], fdata["x3max"], 1.0]
    hfp.attrs["DatasetNames"] = dataset_names
    hfp.attrs["NumVariables"] = dataset_nvars
    hfp.attrs["VariableNames"] = [
        np.array(i, dtype="|S21") for i in (vars_without_b + vars_only_b)
    ]

    # Create Datasets
    if len(vars_only_b) > 0:
        hfp.create_dataset("B", data=B, dtype=varfmt)
    hfp.create_dataset("Levels", data=Levels, dtype=">i4")
    hfp.create_dataset("LogicalLocations", data=LogicalLocations, dtype=">i8")
    hfp.create_dataset("uov", data=uov, dtype=varfmt)
    hfp.create_dataset("x1f", data=x1f, dtype=locfmt)
    hfp.create_dataset("x1v", data=x1v, dtype=locfmt)
    hfp.create_dataset("x2f", data=x2f, dtype=locfmt)
    hfp.create_dataset("x2v", data=x2v, dtype=locfmt)
    hfp.create_dataset("x3f", data=x3f, dtype=locfmt)
    hfp.create_dataset("x3v", data=x3v, dtype=locfmt)
    hfp.close()


def write_xdmf_for(xdmfname, dumpname, fdata, mode="auto"):
    """
    Writes an xdmf file for a fluid snapshot file.

    args:
      xdmfname - string
          name of xdmf file
      dumpname - string
          location of fluid data file relative to xdmfname directory
      fdata    - dict
          dictionary of fluid file data, e.g., as loaded from read_binary(...)
      mode     - string (unimplemented)
          force xdmf for format (auto sets by extension)
    """

    fp = open(xdmfname, "w")

    def write_meshblock(fp, mb, nx1, nx2, nx3, nmb, dumpname, vars_no_b, vars_w_b):
        fp.write(f"""  <Grid Name="MeshBlock{mb}" GridType="Uniform">\n""")
        fp.write("""   <Topology TopologyType="3DRectMesh" """)
        fp.write(f""" NumberOfElements="{nx3 + 1} {nx2 + 1} {nx1 + 1}"/>\n""")
        fp.write("""   <Geometry GeometryType="VXVYVZ">\n""")
        fp.write(
            f"""    <DataItem ItemType="HyperSlab" Dimensions="{nx1 + 1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx1 + 1} </DataItem>
     <DataItem Dimensions="{nmb} {nx1 + 1}" Format="HDF"> {dumpname}:/x1f </DataItem>
    </DataItem>
    <DataItem ItemType="HyperSlab" Dimensions="{nx2 + 1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx2 + 1} </DataItem>
     <DataItem Dimensions="{nmb} {nx2 + 1}" Format="HDF"> {dumpname}:/x2f </DataItem>
    </DataItem>
    <DataItem ItemType="HyperSlab" Dimensions="{nx3 + 1}">
     <DataItem Dimensions="3 2" NumberType="Int"> {mb} 0 1 1 1 {nx3 + 1} </DataItem>
     <DataItem Dimensions="{nmb} {nx3 + 1}" Format="HDF"> {dumpname}:/x3f </DataItem>
    </DataItem>
   </Geometry>\n"""
        )

        nvar_no_b = len(vars_no_b)
        for vi, var_name in enumerate(vars_no_b):
            fp.write(
                f"""   <Attribute Name="{var_name}" Center="Cell">
    <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">
     <DataItem Dimensions="3 5" NumberType="Int">
      {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}
     </DataItem>
     <DataItem Dimensions="{nvar_no_b} {nmb} {nx3} {nx2} {nx1}" Format="HDF">
      {dumpname}:/uov
     </DataItem>
    </DataItem>
   </Attribute>\n"""
            )

        nvar_w_b = len(vars_w_b)
        if nvar_w_b > 0:
            for vi, var_name in enumerate(vars_w_b):
                fp.write(
                    f"""   <Attribute Name="{var_name}" Center="Cell">
        <DataItem ItemType="HyperSlab" Dimensions="{nx3} {nx2} {nx1}">
         <DataItem Dimensions="3 5" NumberType="Int">
          {vi} {mb} 0 0 0 1 1 1 1 1 1 1 {nx3} {nx2} {nx1}
         </DataItem>
         <DataItem Dimensions="{nvar_w_b} {nmb} {nx3} {nx2} {nx1}" Format="HDF">
          {dumpname}:/B
         </DataItem>
        </DataItem>
       </Attribute>\n"""
                )

        fp.write("""  </Grid>\n""")

    fp.write(
        """<?xml version="1.0" ?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="2.0">
<Information Name="TimeVaryingMetaData" Value="True"/>\n"""
    )
    fp.write("""<Domain>\n""")
    fp.write("""<Grid Name="Mesh" GridType="Collection">\n""")
    fp.write(f""" <Time Value="{fdata['time']}"/>\n""")

    vars_without_b = [v for v in fdata["var_names"] if "bcc" not in v]
    vars_only_b = [v for v in fdata["var_names"] if v not in vars_without_b]

    nx1 = fdata["nx1_out_mb"]
    nx2 = fdata["nx2_out_mb"]
    nx3 = fdata["nx3_out_mb"]
    nmb = fdata["n_mbs"]

    for mb in range(nmb):
        write_meshblock(
            fp, mb, nx1, nx2, nx3, nmb, dumpname, vars_without_b, vars_only_b
        )

    fp.write("""</Grid>\n""")
    fp.write("""</Domain>\n""")
    fp.write("""</Xdmf>\n""")

    fp.close()


def convert_file(
    binary_fname, assemble_shards=False, coarsened=None, *, limits=None
):
    """
    Converts a binary file, optionally assembling its rank/node shard family.

    args:
      binary_filename - string
        filename of bin file to convert
      assemble_shards - bool, optional
        assemble sibling rank_* or node_* files before writing output
      coarsened - bool or None, optional
        use the coarsened-binary reader; when None, infer from a .cbin suffix

    This will create new files "binary_data.bin" -> "binary_data.athdf" and
    "binary_data.athdf.xdmf"
    """
    limits = _binary_reader_limits(limits)
    athdf_fname = os.path.splitext(binary_fname)[0] + ".athdf"
    xdmf_fname = athdf_fname + ".xdmf"
    if coarsened is None:
        coarsened = binary_fname.endswith(".cbin")
    if coarsened:
        filedata = _read_coarsened_binary(
            binary_fname, assemble_shards=assemble_shards, limits=limits
        )
    else:
        filedata = _read_binary(
            binary_fname, assemble_shards=assemble_shards, limits=limits
        )
    if filedata["n_mbs"] == 0:
        raise ValueError(
            f"cannot convert {binary_fname!r}: binary output contains no meshblocks"
        )
    write_athdf(athdf_fname, filedata, limits=limits)
    write_xdmf_for(xdmf_fname, os.path.basename(athdf_fname), filedata)


__all__ = [
    "read_binary",
    "read_coarsened_binary",
    "read_all_ranks_binary",
    "read_all_ranks_coarsened_binary",
    "read_binary_as_athdf",
    "read_all_ranks_binary_as_athdf",
    "read_all_ranks_coarsened_binary_as_athdf",
    "read_single_rank_binary_as_athdf",
    "read_rank_binary_as_athdf",
    "read_coarsened_binary_as_athdf",
    "write_athdf",
    "write_xdmf_for",
    "convert_file",
]


if __name__ == "__main__":
    import argparse

    try:
        from tqdm import tqdm
    except ModuleNotFoundError:

        def tqdm(L):
            for x in L:
                print(x)
                yield x

    parser = argparse.ArgumentParser(
        description="Convert AthenaK binary output to ATHDF/XDMF."
    )
    parser.add_argument(
        "--assemble-shards",
        action="store_true",
        help="assemble matching rank_* or node_* sibling shards before conversion",
    )
    parser.add_argument(
        "--coarsened",
        action="store_true",
        default=None,
        help="force the coarsened-binary reader (otherwise inferred from .cbin)",
    )
    add_reader_limit_arguments(parser)
    parser.add_argument("binary_files", nargs="+", help="binary files to convert")
    args = parser.parse_args(sys.argv[1:])
    limits = reader_limits_from_args(args)

    for binary_fname in tqdm(args.binary_files):
        convert_file(
            binary_fname,
            assemble_shards=args.assemble_shards,
            coarsened=args.coarsened,
            limits=limits,
        )
