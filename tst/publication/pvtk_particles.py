"""Minimal reader for AthenaK legacy binary particle VTK files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict

import numpy as np


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
            if len(tokens) < 3:
                raise ValueError("Malformed POINTS header")
            npoints = int(tokens[1])
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

    point_data_count = npoints
    while idx < len(blob):
        line, idx = _read_line(blob, idx)
        if not line:
            continue
        if line.startswith("POINT_DATA "):
            tokens = line.split()
            if len(tokens) >= 2:
                point_data_count = int(tokens[1])
            continue
        if line.startswith("SCALARS "):
            tokens = line.split()
            if len(tokens) < 2:
                raise ValueError("Malformed SCALARS header")
            name = tokens[1]
            look, idx = _read_line(blob, idx)
            if not look.startswith("LOOKUP_TABLE"):
                raise ValueError("Expected LOOKUP_TABLE after SCALARS header")
            nbytes = 4 * point_data_count
            if idx + nbytes > len(blob):
                raise ValueError("Unexpected EOF in SCALARS data block")
            values = np.frombuffer(
                blob,
                dtype=">f4",
                count=point_data_count,
                offset=idx,
            ).astype(np.float64)
            scalars[name] = values
            idx += nbytes
            continue
        if line.startswith("VECTORS "):
            tokens = line.split()
            if len(tokens) < 2:
                raise ValueError("Malformed VECTORS header")
            name = tokens[1]
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

    return ParticleVTKData(points=points, scalars=scalars, vectors=vectors)
