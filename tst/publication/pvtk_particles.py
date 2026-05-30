"""Minimal reader for AthenaK legacy binary particle VTK files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np

_VTK_SCALAR_DTYPES = {
    "float": np.dtype(">f4"),
    "int": np.dtype(">i4"),
}


@dataclass
class ParticleVTKData:
    points: np.ndarray
    scalars: Dict[str, np.ndarray]
    vectors: Dict[str, np.ndarray]


def _read_line(blob: bytes, idx: int) -> tuple[str, int]:
    end = blob.find(b"\n", idx)
    if end < 0:
        raise ValueError("Malformed VTK file: missing line terminator")
    line = blob[idx:end].decode("ascii").strip()
    return line, end + 1


def read_particle_vtk(path: str | Path) -> ParticleVTKData:
    """Read one AthenaK `pvtk/*.part.vtk` file.

    The parser assumes the legacy binary layout produced by `vtk_prtcl.cpp`.
    """
    blob = Path(path).read_bytes()
    idx = 0

    npoints = None
    while idx < len(blob):
        line, idx = _read_line(blob, idx)
        if line.startswith("POINTS "):
            tokens = line.split()
            if len(tokens) != 3 or tokens[2] != "float":
                raise ValueError("Malformed POINTS header")
            npoints = int(tokens[1])
            if npoints < 0:
                raise ValueError("POINTS count must be non-negative")
            break
    if npoints is None:
        raise ValueError("POINTS section not found")

    npt_floats = 3 * npoints
    npt_bytes = 4 * npt_floats
    if idx + npt_bytes > len(blob):
        raise ValueError("Unexpected EOF in POINTS binary block")
    points = np.frombuffer(blob, dtype=">f4", count=npt_floats, offset=idx)
    points = points.astype(np.float64).reshape(npoints, 3)
    idx += npt_bytes

    scalars: Dict[str, np.ndarray] = {}
    vectors: Dict[str, np.ndarray] = {}

    point_data_count = None
    while idx < len(blob):
        line, idx = _read_line(blob, idx)
        if not line:
            continue
        if line.startswith("POINT_DATA "):
            tokens = line.split()
            if point_data_count is not None:
                raise ValueError("Duplicate POINT_DATA section")
            if len(tokens) != 2:
                raise ValueError("Malformed POINT_DATA header")
            point_data_count = int(tokens[1])
            if point_data_count != npoints:
                raise ValueError("POINT_DATA count does not match POINTS count")
            continue
        if line.startswith("SCALARS "):
            if point_data_count is None:
                raise ValueError("SCALARS section appears before POINT_DATA")
            tokens = line.split()
            if len(tokens) != 3:
                raise ValueError("Malformed SCALARS header")
            name = tokens[1]
            scalar_type = tokens[2]
            if name in scalars or name in vectors:
                raise ValueError(f"Duplicate SCALARS section: {name}")
            if scalar_type not in _VTK_SCALAR_DTYPES:
                raise ValueError(f"Unsupported SCALARS type: {scalar_type}")
            dtype = _VTK_SCALAR_DTYPES[scalar_type]
            look, idx = _read_line(blob, idx)
            if look != "LOOKUP_TABLE default":
                raise ValueError("Expected LOOKUP_TABLE after SCALARS header")
            nbytes = dtype.itemsize * point_data_count
            if idx + nbytes > len(blob):
                raise ValueError("Unexpected EOF in SCALARS data block")
            values = np.frombuffer(
                blob,
                dtype=dtype,
                count=point_data_count,
                offset=idx,
            )
            if scalar_type == "float":
                values = values.astype(np.float64)
            else:
                values = values.astype(np.int64)
            scalars[name] = values
            idx += nbytes
            continue
        if line.startswith("VECTORS "):
            if point_data_count is None:
                raise ValueError("VECTORS section appears before POINT_DATA")
            tokens = line.split()
            if len(tokens) != 3:
                raise ValueError("Malformed VECTORS header")
            name = tokens[1]
            vector_type = tokens[2]
            if name in vectors or name in scalars:
                raise ValueError(f"Duplicate VECTORS section: {name}")
            if vector_type != "float":
                raise ValueError(f"Unsupported VECTORS type: {vector_type}")
            nvec_floats = 3 * point_data_count
            nbytes = 4 * nvec_floats
            if idx + nbytes > len(blob):
                raise ValueError("Unexpected EOF in VECTORS data block")
            values = np.frombuffer(
                blob,
                dtype=">f4",
                count=nvec_floats,
                offset=idx,
            ).astype(np.float64).reshape(point_data_count, 3)
            vectors[name] = values
            idx += nbytes
            continue
        raise ValueError(f"Unexpected VTK content: {line}")

    if point_data_count is None:
        raise ValueError("POINT_DATA section not found")
    return ParticleVTKData(points=points, scalars=scalars, vectors=vectors)
