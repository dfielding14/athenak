#!/usr/bin/env python3
"""Strict reusable arithmetic foundation for Q-011 Section 5.4 output analysis.

This module intentionally does not bind artifacts, choose campaign windows, or
set qualification thresholds.  Callers must supply those policy decisions.
"""

from __future__ import annotations

from dataclasses import dataclass
import io
import math
from pathlib import Path
import re
from typing import Mapping, Sequence

import numpy as np


_MAGIC = b"Athena binary output version=1.1\n"
_LINE_PATTERNS = {
    "preheader": re.compile(rb"  size of preheader=([0-9]+)\n"),
    "nvars": re.compile(rb"  number of variables=([0-9]+)\n"),
    "variables": re.compile(rb"  variables:(?:  ([^\n]*))?\n"),
    "header_offset": re.compile(rb"  header offset=([0-9]+)\n"),
}
_REQUIRED_PREHEADER_KEYS = {
    "time",
    "cycle",
    "size of location",
    "size of variable",
}
_REQUIRED_INPUT_PARAMETERS = {
    ("mesh", "nx1"),
    ("mesh", "nx2"),
    ("mesh", "nx3"),
    ("mesh", "nghost"),
    ("mesh", "x1min"),
    ("mesh", "x1max"),
    ("mesh", "x2min"),
    ("mesh", "x2max"),
    ("mesh", "x3min"),
    ("mesh", "x3max"),
    ("meshblock", "nx1"),
    ("meshblock", "nx2"),
    ("meshblock", "nx3"),
}
_MAX_ASCII_HEADER_BYTES = 64 * 1024 * 1024
_MAX_VARIABLES = 4096
_MAX_BLOCKS = 1 << 20
_MAX_CELLS_PER_BLOCK = 1 << 31
_MAX_PHYSICAL_REFINEMENT_LEVEL = 30


class AnalysisError(ValueError):
    """Raised when retained output or requested arithmetic fails closed."""


@dataclass(frozen=True)
class AthenaBinaryBlock:
    """One cell-centered AthenaK mesh-block payload."""

    index_bounds: tuple[int, int, int, int, int, int]
    logical_location: tuple[int, int, int]
    level: int
    geometry: tuple[float, float, float, float, float, float]
    fields: Mapping[str, np.ndarray]

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return the payload shape in ``(x3, x2, x1)`` order."""
        i0, i1, j0, j1, k0, k1 = self.index_bounds
        return (k1 - k0 + 1, j1 - j0 + 1, i1 - i0 + 1)

    @property
    def leaf_key(self) -> tuple[int, int, int, int]:
        """Return the physical-level leaf identity."""
        return (*self.logical_location, self.level)


@dataclass(frozen=True)
class AthenaBinaryDataset:
    """Strictly parsed metadata and mesh-block payloads from AthenaK ``*.bin``."""

    source: str
    time: float
    cycle: int
    location_size: int
    variable_size: int
    variable_names: tuple[str, ...]
    input_parameters: Mapping[str, Mapping[str, str]]
    root_grid_shape: tuple[int, int, int]
    meshblock_shape: tuple[int, int, int]
    nghost: int
    domain_bounds: tuple[float, float, float, float, float, float]
    blocks: tuple[AthenaBinaryBlock, ...]

    @property
    def root_block_shape(self) -> tuple[int, int, int]:
        """Return root-grid MeshBlock counts in ``(x1, x2, x3)`` order."""
        return tuple(
            root // block
            for root, block in zip(self.root_grid_shape, self.meshblock_shape)
        )


@dataclass(frozen=True)
class CompositeGrid:
    """A validated finest-level Cartesian composite grid."""

    x1_faces: np.ndarray
    x2_faces: np.ndarray
    x3_faces: np.ndarray
    values: np.ndarray
    source_levels: np.ndarray
    target_level: int


@dataclass(frozen=True)
class AreaRestriction:
    """Conservative area restriction result."""

    values: np.ndarray
    areas: np.ndarray


@dataclass(frozen=True)
class FixedHistogram:
    """Histogram arithmetic on caller-supplied fixed bins."""

    bin_edges: np.ndarray
    counts: np.ndarray
    weighted_counts: np.ndarray
    underflow_count: int
    overflow_count: int
    underflow_weight: float
    overflow_weight: float

    @property
    def bin_centers(self) -> np.ndarray:
        """Return arithmetic bin centers."""
        return 0.5 * (self.bin_edges[:-1] + self.bin_edges[1:])

    @property
    def geometric_bin_centers(self) -> np.ndarray:
        """Return geometric centers for positive bins."""
        _require(np.all(self.bin_edges > 0.0), "logarithmic bin edges must be positive")
        return np.sqrt(self.bin_edges[:-1] * self.bin_edges[1:])


@dataclass(frozen=True)
class LogLogSlopeFit:
    """Least-squares fit of fixed-bin values against geometric bin centers."""

    slope: float
    intercept: float
    selected_bin_indices: tuple[int, ...]


@dataclass(frozen=True)
class ShockFront:
    """Strongest caller-constrained one-dimensional density jump."""

    index: int
    x1: float
    density_gradient: float


@dataclass(frozen=True)
class UpstreamAmplification:
    """Area-weighted upstream magnetic-field magnitude relative to a reference."""

    selected_cell_count: int
    selected_area: float
    mean_magnetic_magnitude: float
    amplification: float


@dataclass(frozen=True)
class MatchedGridResiduals:
    """Area-weighted residual norms for arrays already placed on one grid."""

    mean_absolute: float
    root_mean_square: float
    maximum_absolute: float
    relative_mean_absolute: float | None
    relative_root_mean_square: float | None
    relative_maximum_absolute: float | None


class _Cursor:
    """Small exact-read wrapper that turns truncation into contextual failures."""

    def __init__(self, payload: bytes, source: str):
        self._stream = io.BytesIO(payload)
        self._size = len(payload)
        self._source = source

    @property
    def offset(self) -> int:
        return self._stream.tell()

    @property
    def remaining(self) -> int:
        return self._size - self.offset

    def read_exact(self, size: int, label: str) -> bytes:
        _require(size >= 0, f"{self._source}: negative read size for {label}")
        payload = self._stream.read(size)
        if len(payload) != size:
            raise AnalysisError(
                f"{self._source}: truncated {label} at byte {self.offset - len(payload)}"
            )
        return payload

    def readline(self, label: str) -> bytes:
        line = self._stream.readline()
        if not line:
            raise AnalysisError(f"{self._source}: truncated {label} at byte {self.offset}")
        if not line.endswith(b"\n"):
            raise AnalysisError(f"{self._source}: unterminated {label}")
        return line


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AnalysisError(message)


def _parse_uint_line(line: bytes, kind: str, source: str) -> int:
    match = _LINE_PATTERNS[kind].fullmatch(line)
    _require(match is not None, f"{source}: malformed {kind} line")
    return int(match.group(1))


def _parse_key_value_line(line: bytes, label: str) -> tuple[str, str]:
    try:
        decoded = line.decode("utf-8").rstrip("\n")
    except UnicodeDecodeError as exc:
        raise AnalysisError(f"{label}: preheader is not UTF-8") from exc
    _require(decoded.count("=") == 1, f"{label}: malformed preheader key/value line")
    key, value = (item.strip() for item in decoded.split("=", 1))
    _require(bool(key) and bool(value), f"{label}: empty preheader key or value")
    return key, value


def _parse_input_parameters(payload: bytes, source: str) -> dict[str, dict[str, str]]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise AnalysisError(f"{source}: parameter header is not UTF-8") from exc
    blocks: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    for lineno, original in enumerate(text.splitlines(), start=1):
        line = original.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            _require(
                line.endswith(">") and line.count("<") == 1 and line.count(">") == 1,
                f"{source}: parameter header line {lineno}: malformed block",
            )
            name = line[1:-1].strip()
            _require(bool(name), f"{source}: parameter header line {lineno}: empty block")
            _require(
                name not in blocks,
                f"{source}: parameter header line {lineno}: duplicate block {name}",
            )
            current = {}
            blocks[name] = current
            continue
        _require(
            current is not None,
            f"{source}: parameter header line {lineno}: parameter precedes block",
        )
        _require(
            line.count("=") == 1,
            f"{source}: parameter header line {lineno}: malformed parameter",
        )
        key, value = (item.strip() for item in line.split("=", 1))
        # AthenaK can emit runtime-added optional parameters with an empty value,
        # for example particles/pic_deltaf_f0 when delta-f mode is disabled.
        _require(
            bool(key),
            f"{source}: parameter header line {lineno}: empty parameter key",
        )
        _require(
            key not in current,
            f"{source}: parameter header line {lineno}: duplicate parameter {key}",
        )
        current[key] = value
    missing = sorted(
        f"{block}/{key}"
        for block, key in _REQUIRED_INPUT_PARAMETERS
        if block not in blocks or key not in blocks[block]
    )
    _require(not missing, f"{source}: missing parameter header entries: {missing}")
    return blocks


def _parameter_int(
    parameters: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> int:
    try:
        value = int(parameters[block][key])
    except ValueError as exc:
        raise AnalysisError(f"{source}: {block}/{key} is not an integer") from exc
    return value


def _parameter_float(
    parameters: Mapping[str, Mapping[str, str]], block: str, key: str, source: str
) -> float:
    try:
        value = float(parameters[block][key])
    except ValueError as exc:
        raise AnalysisError(f"{source}: {block}/{key} is not a float") from exc
    _require(math.isfinite(value), f"{source}: {block}/{key} is not finite")
    return value


def _shape_from_bounds(
    bounds: Sequence[int], source: str, block_index: int
) -> tuple[int, int, int]:
    i0, i1, j0, j1, k0, k1 = bounds
    _require(
        i1 >= i0 and j1 >= j0 and k1 >= k0,
        f"{source}: block {block_index}: reversed index bounds",
    )
    shape = (k1 - k0 + 1, j1 - j0 + 1, i1 - i0 + 1)
    cells = math.prod(shape)
    _require(
        0 < cells <= _MAX_CELLS_PER_BLOCK,
        f"{source}: block {block_index}: unreasonable cell count",
    )
    return shape


def _decode_array(
    cursor: _Cursor, count: int, dtype: np.dtype, label: str
) -> np.ndarray:
    payload = cursor.read_exact(count * dtype.itemsize, label)
    return np.frombuffer(payload, dtype=dtype, count=count).copy()


def parse_athenak_binary_bytes(
    payload: bytes, *, source: str = "<bytes>"
) -> AthenaBinaryDataset:
    """Parse and structurally validate one AthenaK binary ``version=1.1`` payload."""
    _require(isinstance(payload, bytes), f"{source}: binary payload must be bytes")
    cursor = _Cursor(payload, source)
    _require(cursor.readline("magic header") == _MAGIC, f"{source}: bad magic header")

    preheader_count = _parse_uint_line(cursor.readline("preheader count"), "preheader", source)
    _require(preheader_count >= 1, f"{source}: invalid preheader count")
    preheader: dict[str, str] = {}
    for _ in range(preheader_count - 1):
        key, value = _parse_key_value_line(cursor.readline("preheader field"), source)
        _require(key not in preheader, f"{source}: duplicate preheader field {key}")
        preheader[key] = value
    _require(
        set(preheader) == _REQUIRED_PREHEADER_KEYS,
        f"{source}: unexpected preheader fields {sorted(preheader)}",
    )
    try:
        time = float(preheader["time"])
        cycle = int(preheader["cycle"])
        location_size = int(preheader["size of location"])
        variable_size = int(preheader["size of variable"])
    except ValueError as exc:
        raise AnalysisError(f"{source}: malformed preheader value") from exc
    _require(math.isfinite(time), f"{source}: time is not finite")
    _require(cycle >= 0, f"{source}: cycle must be non-negative")
    _require(location_size in (4, 8), f"{source}: unsupported location size")
    _require(variable_size in (4, 8), f"{source}: unsupported variable size")

    nvars = _parse_uint_line(cursor.readline("variable count"), "nvars", source)
    _require(0 < nvars <= _MAX_VARIABLES, f"{source}: invalid variable count")
    variable_line = cursor.readline("variable names")
    match = _LINE_PATTERNS["variables"].fullmatch(variable_line)
    _require(match is not None, f"{source}: malformed variables line")
    encoded_names = (
        [] if match.group(1) is None else match.group(1).strip().split(b"  ")
    )
    try:
        variable_names = tuple(name.decode("utf-8") for name in encoded_names)
    except UnicodeDecodeError as exc:
        raise AnalysisError(f"{source}: variable name is not UTF-8") from exc
    _require(
        len(variable_names) == nvars and all(variable_names),
        f"{source}: variable count does not match names",
    )
    _require(
        len(set(variable_names)) == len(variable_names),
        f"{source}: duplicate variable name",
    )

    header_size = _parse_uint_line(cursor.readline("header offset"), "header_offset", source)
    _require(
        header_size <= _MAX_ASCII_HEADER_BYTES,
        f"{source}: unreasonable parameter header size",
    )
    parameters = _parse_input_parameters(
        cursor.read_exact(header_size, "parameter header"), source
    )
    root_grid_shape = tuple(
        _parameter_int(parameters, "mesh", f"nx{axis}", source) for axis in range(1, 4)
    )
    meshblock_shape = tuple(
        _parameter_int(parameters, "meshblock", f"nx{axis}", source)
        for axis in range(1, 4)
    )
    nghost = _parameter_int(parameters, "mesh", "nghost", source)
    domain_bounds = tuple(
        value
        for axis in range(1, 4)
        for value in (
            _parameter_float(parameters, "mesh", f"x{axis}min", source),
            _parameter_float(parameters, "mesh", f"x{axis}max", source),
        )
    )
    _require(all(size > 0 for size in root_grid_shape), f"{source}: invalid root-grid shape")
    _require(all(size > 0 for size in meshblock_shape), f"{source}: invalid MeshBlock shape")
    _require(
        math.prod(meshblock_shape) <= _MAX_CELLS_PER_BLOCK,
        f"{source}: unreasonable MeshBlock cell count",
    )
    _require(nghost >= 0, f"{source}: nghost must be non-negative")
    _require(
        all(root % block == 0 for root, block in zip(root_grid_shape, meshblock_shape)),
        f"{source}: root grid is not divisible by MeshBlock shape",
    )
    _require(
        math.prod(root // block for root, block in zip(root_grid_shape, meshblock_shape))
        <= _MAX_BLOCKS,
        f"{source}: too many root-grid MeshBlocks",
    )
    _require(
        all(upper > lower for lower, upper in zip(domain_bounds[::2], domain_bounds[1::2])),
        f"{source}: invalid domain bounds",
    )

    location_dtype = np.dtype("=f4" if location_size == 4 else "=f8")
    variable_dtype = np.dtype("=f4" if variable_size == 4 else "=f8")
    blocks: list[AthenaBinaryBlock] = []
    while cursor.remaining:
        block_index = len(blocks)
        _require(block_index < _MAX_BLOCKS, f"{source}: too many MeshBlocks")
        bounds_array = _decode_array(
            cursor, 6, np.dtype("=i4"), f"block {block_index} bounds"
        )
        bounds = tuple(int(value) for value in bounds_array)
        shape = _shape_from_bounds(bounds, source, block_index)
        logical = _decode_array(
            cursor, 4, np.dtype("=i4"), f"block {block_index} logical location"
        )
        logical_location = tuple(int(value) for value in logical[:3])
        level = int(logical[3])
        geometry_array = _decode_array(
            cursor, 6, location_dtype, f"block {block_index} geometry"
        )
        _require(
            np.all(np.isfinite(geometry_array)),
            f"{source}: block {block_index}: geometry is not finite",
        )
        geometry = tuple(float(value) for value in geometry_array)
        _require(
            all(upper > lower for lower, upper in zip(geometry[::2], geometry[1::2])),
            f"{source}: block {block_index}: invalid geometry bounds",
        )
        cell_count = math.prod(shape)
        field_array = _decode_array(
            cursor,
            nvars * cell_count,
            variable_dtype,
            f"block {block_index} field data",
        ).reshape((nvars, *shape))
        _require(
            np.all(np.isfinite(field_array)),
            f"{source}: block {block_index}: field data is not finite",
        )
        fields = {
            name: field_array[index].copy() for index, name in enumerate(variable_names)
        }
        blocks.append(
            AthenaBinaryBlock(
                index_bounds=bounds,
                logical_location=logical_location,
                level=level,
                geometry=geometry,
                fields=fields,
            )
        )
    _require(bool(blocks), f"{source}: no MeshBlocks")
    dataset = AthenaBinaryDataset(
        source=source,
        time=time,
        cycle=cycle,
        location_size=location_size,
        variable_size=variable_size,
        variable_names=variable_names,
        input_parameters=parameters,
        root_grid_shape=root_grid_shape,
        meshblock_shape=meshblock_shape,
        nghost=nghost,
        domain_bounds=domain_bounds,
        blocks=tuple(blocks),
    )
    validate_leaf_blocks(dataset)
    return dataset


def read_athenak_binary(path: str | Path) -> AthenaBinaryDataset:
    """Read and strictly parse one AthenaK binary file."""
    binary_path = Path(path)
    try:
        payload = binary_path.read_bytes()
    except OSError as exc:
        raise AnalysisError(f"{binary_path}: unable to read binary output") from exc
    return parse_athenak_binary_bytes(payload, source=str(binary_path))


def _logical_axis_extent(dataset: AthenaBinaryDataset, axis: int, level: int) -> int:
    root_blocks = dataset.root_block_shape[axis]
    return root_blocks * (2**level if dataset.root_grid_shape[axis] > 1 else 1)


def _validate_block_geometry(dataset: AthenaBinaryDataset, block: AthenaBinaryBlock) -> None:
    epsilon = np.finfo(np.float32 if dataset.location_size == 4 else np.float64).eps
    for axis, logical_index in enumerate(block.logical_location):
        domain_min = dataset.domain_bounds[2 * axis]
        domain_max = dataset.domain_bounds[2 * axis + 1]
        extent = _logical_axis_extent(dataset, axis, block.level)
        width = (domain_max - domain_min) / extent
        expected_min = domain_min + logical_index * width
        expected_max = expected_min + width
        actual_min = block.geometry[2 * axis]
        actual_max = block.geometry[2 * axis + 1]
        tolerance = 32.0 * epsilon * max(1.0, abs(domain_min), abs(domain_max))
        _require(
            math.isclose(actual_min, expected_min, rel_tol=32.0 * epsilon, abs_tol=tolerance)
            and math.isclose(actual_max, expected_max, rel_tol=32.0 * epsilon, abs_tol=tolerance),
            f"{dataset.source}: leaf {block.leaf_key}: geometry disagrees with logical location",
        )


def validate_leaf_blocks(dataset: AthenaBinaryDataset) -> None:
    """Validate levels, logical extents, geometry, unique leaves, and field shapes."""
    seen: set[tuple[int, int, int, int]] = set()
    for block in dataset.blocks:
        _require(
            0 <= block.level <= _MAX_PHYSICAL_REFINEMENT_LEVEL,
            f"{dataset.source}: invalid physical refinement level",
        )
        _require(
            all(index >= 0 for index in block.logical_location),
            f"{dataset.source}: negative logical location",
        )
        for axis, logical_index in enumerate(block.logical_location):
            extent = _logical_axis_extent(dataset, axis, block.level)
            _require(
                logical_index < extent,
                f"{dataset.source}: logical location outside level extent",
            )
        _require(
            block.leaf_key not in seen,
            f"{dataset.source}: duplicate leaf {block.leaf_key}",
        )
        seen.add(block.leaf_key)
        if block.shape == tuple(reversed(dataset.meshblock_shape)):
            _validate_block_geometry(dataset, block)
        _require(
            set(block.fields) == set(dataset.variable_names),
            f"{dataset.source}: leaf {block.leaf_key}: field names drifted",
        )
        for name, values in block.fields.items():
            _require(
                values.shape == block.shape,
                f"{dataset.source}: leaf {block.leaf_key}: {name} shape drifted",
            )
            _require(
                np.all(np.isfinite(values)),
                f"{dataset.source}: leaf {block.leaf_key}: {name} is not finite",
            )


def _compatible_metadata(dataset: AthenaBinaryDataset) -> tuple[object, ...]:
    return (
        dataset.time,
        dataset.cycle,
        dataset.location_size,
        dataset.variable_size,
        dataset.variable_names,
        dataset.root_grid_shape,
        dataset.meshblock_shape,
        dataset.nghost,
        dataset.domain_bounds,
        dataset.input_parameters,
    )


def merge_athenak_binary_datasets(
    datasets: Sequence[AthenaBinaryDataset],
) -> AthenaBinaryDataset:
    """Merge strict same-snapshot rank datasets before composite-grid validation."""
    _require(bool(datasets), "at least one binary dataset is required")
    first = datasets[0]
    metadata = _compatible_metadata(first)
    blocks: list[AthenaBinaryBlock] = []
    sources: list[str] = []
    for dataset in datasets:
        _require(
            _compatible_metadata(dataset) == metadata,
            f"{dataset.source}: rank binary metadata disagrees with {first.source}",
        )
        sources.append(dataset.source)
        blocks.extend(dataset.blocks)
    merged = AthenaBinaryDataset(
        source=", ".join(sources),
        time=first.time,
        cycle=first.cycle,
        location_size=first.location_size,
        variable_size=first.variable_size,
        variable_names=first.variable_names,
        input_parameters=first.input_parameters,
        root_grid_shape=first.root_grid_shape,
        meshblock_shape=first.meshblock_shape,
        nghost=first.nghost,
        domain_bounds=first.domain_bounds,
        blocks=tuple(blocks),
    )
    validate_leaf_blocks(merged)
    return merged


def _faces(dataset: AthenaBinaryDataset, axis: int, cell_count: int) -> np.ndarray:
    lower = dataset.domain_bounds[2 * axis]
    upper = dataset.domain_bounds[2 * axis + 1]
    return np.linspace(lower, upper, cell_count + 1, dtype=np.float64)


def compose_leaf_field(
    dataset: AthenaBinaryDataset, field: str, *, target_level: int | None = None
) -> CompositeGrid:
    """Compose a complete Cartesian leaf field and reject overlap or holes."""
    validate_leaf_blocks(dataset)
    _require(field in dataset.variable_names, f"{dataset.source}: unknown field {field}")
    maximum_level = max(block.level for block in dataset.blocks)
    if target_level is None:
        target_level = maximum_level
    _require(
        isinstance(target_level, int) and not isinstance(target_level, bool),
        "target level must be an integer",
    )
    _require(
        maximum_level <= target_level <= _MAX_PHYSICAL_REFINEMENT_LEVEL,
        "composite target level must span input leaves and remain supported",
    )
    output_shape_xyz = tuple(
        root * (2**target_level if root > 1 else 1)
        for root in dataset.root_grid_shape
    )
    output_shape = tuple(reversed(output_shape_xyz))
    values = np.empty(output_shape, dtype=np.float64)
    source_levels = np.empty(output_shape, dtype=np.int32)
    occupancy = np.zeros(output_shape, dtype=np.uint8)
    nominal_shape = tuple(reversed(dataset.meshblock_shape))
    for block in dataset.blocks:
        _require(
            block.shape == nominal_shape,
            f"{dataset.source}: leaf {block.leaf_key}: sliced payload cannot form full composite",
        )
        _validate_block_geometry(dataset, block)
        scale = 2 ** (target_level - block.level)
        starts_xyz = []
        repeats_xyz = []
        for axis, logical_index in enumerate(block.logical_location):
            active = dataset.root_grid_shape[axis] > 1
            starts_xyz.append(
                logical_index * dataset.meshblock_shape[axis] * scale if active else 0
            )
            repeats_xyz.append(scale if active else 1)
        repeated = np.asarray(block.fields[field], dtype=np.float64)
        for array_axis, repeat in enumerate(reversed(repeats_xyz)):
            if repeat > 1:
                repeated = np.repeat(repeated, repeat, axis=array_axis)
        starts = tuple(reversed(starts_xyz))
        slices = tuple(
            slice(start, start + size) for start, size in zip(starts, repeated.shape)
        )
        _require(
            all(region.stop <= limit for region, limit in zip(slices, output_shape)),
            f"{dataset.source}: leaf {block.leaf_key}: composite placement exceeds domain",
        )
        _require(
            not np.any(occupancy[slices]),
            f"{dataset.source}: composite overlap at leaf {block.leaf_key}",
        )
        values[slices] = repeated
        source_levels[slices] = block.level
        occupancy[slices] = 1
    _require(np.all(occupancy), f"{dataset.source}: composite grid contains holes")
    return CompositeGrid(
        x1_faces=_faces(dataset, 0, output_shape_xyz[0]),
        x2_faces=_faces(dataset, 1, output_shape_xyz[1]),
        x3_faces=_faces(dataset, 2, output_shape_xyz[2]),
        values=values,
        source_levels=source_levels,
        target_level=target_level,
    )


def _restriction_factors(factor: int | tuple[int, int]) -> tuple[int, int]:
    if isinstance(factor, int) and not isinstance(factor, bool):
        factors = (factor, factor)
    else:
        _require(
            isinstance(factor, tuple) and len(factor) == 2,
            "restriction factor must be an integer or a (y, x) tuple",
        )
        factors = factor
    _require(
        all(
            isinstance(value, int) and not isinstance(value, bool) and value > 0
            for value in factors
        ),
        "restriction factors must be positive integers",
    )
    return factors


def _finite_array(values: object, label: str, *, allow_empty: bool = False) -> np.ndarray:
    array = np.asarray(values, dtype=np.float64)
    _require(allow_empty or array.size > 0, f"{label} must not be empty")
    _require(np.all(np.isfinite(array)), f"{label} must be finite")
    return array


def restrict_area_average(
    values: object,
    *,
    factor: int | tuple[int, int] = 2,
    cell_areas: object | None = None,
) -> AreaRestriction:
    """Conservatively restrict the last two axes by area-weighted averaging."""
    array = _finite_array(values, "restriction values")
    _require(array.ndim >= 2, "restriction values must have at least two dimensions")
    fy, fx = _restriction_factors(factor)
    ny, nx = array.shape[-2:]
    _require(
        ny % fy == 0 and nx % fx == 0,
        "restriction shape must be divisible by factor",
    )
    if cell_areas is None:
        areas = np.ones_like(array)
    else:
        areas = _finite_array(cell_areas, "cell areas")
        _require(areas.shape == array.shape, "cell areas must match restriction values")
        _require(np.all(areas > 0.0), "cell areas must be positive")
    leading = array.shape[:-2]
    reshape = (*leading, ny // fy, fy, nx // fx, fx)
    restricted_areas = areas.reshape(reshape).sum(axis=(-3, -1))
    restricted_values = (array * areas).reshape(reshape).sum(axis=(-3, -1))
    restricted_values /= restricted_areas
    return AreaRestriction(values=restricted_values, areas=restricted_areas)


def fixed_histogram(
    values: object, bin_edges: object, *, weights: object | None = None
) -> FixedHistogram:
    """Accumulate deterministic counts and non-negative weights on fixed bins."""
    samples = _finite_array(values, "histogram values", allow_empty=True).reshape(-1)
    edges = _finite_array(bin_edges, "histogram bin edges").reshape(-1)
    _require(edges.size >= 2, "histogram requires at least one bin")
    _require(np.all(np.diff(edges) > 0.0), "histogram bin edges must increase")
    if weights is None:
        sample_weights = np.ones_like(samples)
    else:
        sample_weights = _finite_array(
            weights, "histogram weights", allow_empty=True
        ).reshape(-1)
        _require(sample_weights.shape == samples.shape, "histogram weights must match values")
        _require(np.all(sample_weights >= 0.0), "histogram weights must be non-negative")
    underflow = samples < edges[0]
    overflow = samples > edges[-1]
    admitted = ~(underflow | overflow)
    bin_indices = np.searchsorted(edges, samples[admitted], side="right") - 1
    bin_indices[bin_indices == edges.size - 1] -= 1
    counts = np.bincount(bin_indices, minlength=edges.size - 1).astype(np.int64)
    weighted_counts = np.bincount(
        bin_indices, weights=sample_weights[admitted], minlength=edges.size - 1
    ).astype(np.float64)
    return FixedHistogram(
        bin_edges=edges.copy(),
        counts=counts,
        weighted_counts=weighted_counts,
        underflow_count=int(np.count_nonzero(underflow)),
        overflow_count=int(np.count_nonzero(overflow)),
        underflow_weight=float(np.sum(sample_weights[underflow])),
        overflow_weight=float(np.sum(sample_weights[overflow])),
    )


def fit_fixed_bin_loglog_slope(
    bin_edges: object,
    bin_values: object,
    *,
    fit_window: tuple[float, float],
    minimum_bins: int = 3,
) -> LogLogSlopeFit:
    """Fit ``log(bin_values)`` against positive fixed-bin geometric centers."""
    edges = _finite_array(bin_edges, "fit bin edges").reshape(-1)
    values = _finite_array(bin_values, "fit bin values").reshape(-1)
    _require(edges.size == values.size + 1, "fit bin edge/value sizes disagree")
    _require(np.all(edges > 0.0), "fit bin edges must be positive")
    _require(np.all(np.diff(edges) > 0.0), "fit bin edges must increase")
    _require(np.all(values >= 0.0), "fit bin values must be non-negative")
    _require(
        isinstance(minimum_bins, int)
        and not isinstance(minimum_bins, bool)
        and minimum_bins >= 2,
        "minimum fit-bin count must be an integer of at least two",
    )
    lower, upper = fit_window
    _require(
        math.isfinite(lower) and math.isfinite(upper) and upper > lower,
        "fit window must be finite and increasing",
    )
    centers = np.sqrt(edges[:-1] * edges[1:])
    selected = (centers >= lower) & (centers <= upper) & (values > 0.0)
    indices = np.flatnonzero(selected)
    _require(indices.size >= minimum_bins, "insufficient positive fit bins")
    design = np.column_stack((np.log(centers[selected]), np.ones(indices.size)))
    slope, intercept = np.linalg.lstsq(design, np.log(values[selected]), rcond=None)[0]
    return LogLogSlopeFit(
        slope=float(slope),
        intercept=float(intercept),
        selected_bin_indices=tuple(int(index) for index in indices),
    )


def _validated_window(window: tuple[float, float], label: str) -> tuple[float, float]:
    _require(isinstance(window, tuple) and len(window) == 2, f"{label} must be a tuple")
    lower, upper = window
    _require(
        math.isfinite(lower) and math.isfinite(upper) and upper > lower,
        f"{label} must be finite and increasing",
    )
    return lower, upper


def detect_shock_front(
    x1_centers: object,
    density_profile: object,
    *,
    search_window: tuple[float, float],
    gradient_sign: str = "either",
) -> ShockFront:
    """Select the strongest density jump within an explicit physical window."""
    x1 = _finite_array(x1_centers, "shock x1 centers").reshape(-1)
    density = _finite_array(density_profile, "shock density profile").reshape(-1)
    _require(x1.size == density.size and x1.size >= 3, "shock profile shape is invalid")
    _require(np.all(np.diff(x1) > 0.0), "shock x1 centers must increase")
    lower, upper = _validated_window(search_window, "shock search window")
    _require(
        gradient_sign in {"either", "positive", "negative"},
        "shock gradient sign must be either, positive, or negative",
    )
    gradient = np.gradient(density, x1, edge_order=2)
    if gradient_sign == "positive":
        score = gradient
    elif gradient_sign == "negative":
        score = -gradient
    else:
        score = np.abs(gradient)
    selected = (x1 >= lower) & (x1 <= upper)
    indices = np.flatnonzero(selected)
    _require(indices.size > 0, "shock search window selects no cells")
    selected_scores = score[indices]
    strongest = float(np.max(selected_scores))
    _require(strongest > 0.0, "shock search window contains no requested density jump")
    maxima = indices[selected_scores == strongest]
    _require(maxima.size == 1, "shock search window has an ambiguous strongest jump")
    index = int(maxima[0])
    return ShockFront(index=index, x1=float(x1[index]), density_gradient=float(gradient[index]))


def _broadcast_positive_areas(shape: tuple[int, ...], cell_areas: object | None) -> np.ndarray:
    if cell_areas is None:
        return np.ones(shape, dtype=np.float64)
    areas = _finite_array(cell_areas, "cell areas")
    try:
        broadcast = np.broadcast_to(areas, shape)
    except ValueError as exc:
        raise AnalysisError("cell areas do not broadcast to values") from exc
    _require(np.all(broadcast > 0.0), "cell areas must be positive")
    return np.asarray(broadcast, dtype=np.float64)


def estimate_upstream_amplification(
    x1_centers: object,
    magnetic_magnitude: object,
    *,
    upstream_window: tuple[float, float],
    reference_magnetic_field: float,
    cell_areas: object | None = None,
) -> UpstreamAmplification:
    """Return an area-weighted ``|B| / B0`` estimate inside a supplied window."""
    x1 = _finite_array(x1_centers, "upstream x1 centers").reshape(-1)
    field = _finite_array(magnetic_magnitude, "magnetic magnitude")
    _require(field.ndim >= 1 and field.shape[-1] == x1.size, "magnetic x1 shape disagrees")
    _require(np.all(np.diff(x1) > 0.0), "upstream x1 centers must increase")
    _require(np.all(field >= 0.0), "magnetic magnitude must be non-negative")
    _require(
        math.isfinite(reference_magnetic_field) and reference_magnetic_field > 0.0,
        "reference magnetic field must be finite and positive",
    )
    lower, upper = _validated_window(upstream_window, "upstream window")
    selected_x1 = (x1 >= lower) & (x1 <= upper)
    _require(np.any(selected_x1), "upstream window selects no cells")
    selection = np.broadcast_to(selected_x1, field.shape)
    areas = _broadcast_positive_areas(field.shape, cell_areas)
    selected_area = float(np.sum(areas[selection]))
    mean_field = float(np.sum(field[selection] * areas[selection]) / selected_area)
    return UpstreamAmplification(
        selected_cell_count=int(np.count_nonzero(selection)),
        selected_area=selected_area,
        mean_magnetic_magnitude=mean_field,
        amplification=mean_field / reference_magnetic_field,
    )


def matched_grid_residual_metrics(
    reference: object, candidate: object, *, cell_areas: object | None = None
) -> MatchedGridResiduals:
    """Compute absolute and unfloored relative norms on an already matched grid."""
    expected = _finite_array(reference, "residual reference")
    actual = _finite_array(candidate, "residual candidate")
    _require(expected.shape == actual.shape, "matched residual shapes disagree")
    areas = _broadcast_positive_areas(expected.shape, cell_areas)
    total_area = float(np.sum(areas))
    residual = np.abs(actual - expected)
    mean_absolute = float(np.sum(residual * areas) / total_area)
    root_mean_square = float(np.sqrt(np.sum(residual**2 * areas) / total_area))
    maximum_absolute = float(np.max(residual))
    reference_absolute = float(np.sum(np.abs(expected) * areas) / total_area)
    reference_rms = float(np.sqrt(np.sum(expected**2 * areas) / total_area))
    reference_maximum = float(np.max(np.abs(expected)))
    return MatchedGridResiduals(
        mean_absolute=mean_absolute,
        root_mean_square=root_mean_square,
        maximum_absolute=maximum_absolute,
        relative_mean_absolute=(
            mean_absolute / reference_absolute if reference_absolute > 0.0 else None
        ),
        relative_root_mean_square=(
            root_mean_square / reference_rms if reference_rms > 0.0 else None
        ),
        relative_maximum_absolute=(
            maximum_absolute / reference_maximum if reference_maximum > 0.0 else None
        ),
    )


__all__ = [
    "AnalysisError",
    "AreaRestriction",
    "AthenaBinaryBlock",
    "AthenaBinaryDataset",
    "CompositeGrid",
    "FixedHistogram",
    "LogLogSlopeFit",
    "MatchedGridResiduals",
    "ShockFront",
    "UpstreamAmplification",
    "compose_leaf_field",
    "detect_shock_front",
    "estimate_upstream_amplification",
    "fit_fixed_bin_loglog_slope",
    "fixed_histogram",
    "matched_grid_residual_metrics",
    "merge_athenak_binary_datasets",
    "parse_athenak_binary_bytes",
    "read_athenak_binary",
    "restrict_area_average",
    "validate_leaf_blocks",
]
