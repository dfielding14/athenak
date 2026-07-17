#!/usr/bin/env python3
"""Write AthenaK MHD blocks and particle tracks as XDMF/HDF5."""

from __future__ import annotations

import os
from pathlib import Path
import re
from typing import Sequence
import xml.etree.ElementTree as ET

import h5py
import numpy as np


MHD_VECTORS = {
    "fluid_velocity": ("velx", "vely", "velz"),
    "magnetic_field": ("bcc1", "bcc2", "bcc3"),
}
TRACK_VECTORS = {
    "particle_velocity": ("vx", "vy", "vz"),
    "magnetic_field": ("bx", "by", "bz"),
    "curvature": ("k1", "k2", "k3"),
    "magnetic_field_gradient": ("db1", "db2", "db3"),
}
COORDINATE_FIELDS = ("x", "y", "z")


def xmf_path(path: str | Path) -> Path:
    """Return an XMF path, adding the conventional suffix when needed."""

    path = Path(path)
    return path if path.suffix.lower() == ".xmf" else path.with_suffix(".xmf")


def piece_path(xmf_filename: str | Path, rank: int | None = None) -> Path:
    """Return the HDF5 payload path paired with an XMF collection."""

    path = xmf_path(xmf_filename)
    suffix = ".xdmf" if rank is None else f".rank{rank:05d}"
    return path.parent / f"{path.stem}{suffix}.h5"


def visualization_track_fields(fields: Sequence[str] | None) -> list[str] | None:
    """Ensure a requested track-field subset still contains its geometry."""

    if fields is None:
        return None
    return list(dict.fromkeys((*COORDINATE_FIELDS, *fields)))


def _hdf_name(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", name)


def _vector_fields(
    fields: Sequence[str], definitions: dict[str, tuple[str, str, str]]
) -> tuple[dict[str, tuple[str, str, str]], list[str]]:
    available = set(fields)
    vectors = {
        name: components
        for name, components in definitions.items()
        if set(components) <= available
    }
    consumed = {component for components in vectors.values() for component in components}
    scalars = [
        name for name in fields if name not in consumed and name not in vectors
    ]
    return vectors, scalars


def _array_meta(
    name: str,
    dataset: h5py.Dataset,
    attribute_type: str = "Scalar",
) -> dict:
    return {
        "name": name,
        "path": dataset.name,
        "shape": tuple(dataset.shape),
        "dtype": dataset.dtype.str,
        "attribute_type": attribute_type,
    }


def _write_meshblock(parent: h5py.Group, block: dict, index: int) -> dict:
    group = parent.create_group(f"block_{index:06d}")
    coordinates = group.create_group("coordinates")
    cell_data = group.create_group("cell_data")

    coordinate_meta = []
    for name in ("x1f", "x2f", "x3f"):
        dataset = coordinates.create_dataset(name, data=block[name])
        coordinate_meta.append(_array_meta(name, dataset))

    fields = tuple(block["VariableNames"])
    vectors, scalars = _vector_fields(fields, MHD_VECTORS)
    arrays = []
    for name in scalars:
        dataset = cell_data.create_dataset(_hdf_name(name), data=block[name])
        arrays.append(_array_meta(name, dataset))
    for name, components in vectors.items():
        shape = block[components[0]].shape + (3,)
        dataset = cell_data.create_dataset(
            _hdf_name(name), shape=shape, dtype=block[components[0]].dtype
        )
        for component, source in enumerate(components):
            dataset[..., component] = block[source]
        dataset.attrs["components"] = ",".join(components)
        arrays.append(_array_meta(name, dataset, "Vector"))

    source = str(block.get("SourceFile", ""))
    logical = np.asarray(block.get("LogicalLocation", ()), dtype=np.int64)
    domain_bounds = np.asarray(block.get("DomainBounds", block["Bounds"]))
    periodic_axes = tuple(int(axis) for axis in block.get("PeriodicAxes", ()))
    group.attrs["source_file"] = source
    group.attrs["time"] = block.get("Time", 0.0)
    group.attrs["cycle"] = block.get("NumCycles", 0)
    group.attrs["meshblock_index"] = block.get("MeshBlockIndex", index)
    group.attrs["logical_location"] = logical
    group.attrs["domain_bounds"] = domain_bounds
    group.attrs["periodic_axes"] = periodic_axes

    shape = tuple(block[fields[0]].shape)
    return {
        "name": f"block_{index:06d}",
        "shape": shape,
        "coordinates": coordinate_meta,
        "arrays": arrays,
        "source_file": source,
        "time": float(block.get("Time", 0.0)),
        "cycle": int(block.get("NumCycles", 0)),
        "meshblock_index": int(block.get("MeshBlockIndex", index)),
        "logical_location": logical.tolist(),
        "domain_bounds": domain_bounds.tolist(),
        "periodic_axes": list(periodic_axes),
    }


def _track_arrays(tracks: dict) -> tuple[dict[str, int], dict, list[str]]:
    fields = tuple(tracks["fields"])
    field = {name: index for index, name in enumerate(fields)}
    missing = sorted(set(COORDINATE_FIELDS) - set(field))
    if missing:
        raise ValueError(f"track geometry fields are missing: {missing}")
    vectors, scalars = _vector_fields(fields, TRACK_VECTORS)
    scalars = [
        name
        for name in scalars
        if name not in COORDINATE_FIELDS and name not in ("time", "cycle")
    ]
    return field, vectors, scalars


def _derived_track_names(field: dict[str, int]) -> list[str]:
    names = []
    if {"bx", "by", "bz"} <= set(field) and "magnetic_field_magnitude" not in field:
        names.append("magnetic_field_magnitude")
    if {"k1", "k2", "k3"} <= set(field) and "curvature_magnitude" not in field:
        names.append("curvature_magnitude")
    if (
        {"vx", "vy", "vz", "bx", "by", "bz"} <= set(field)
        and "mu_M" not in field
    ):
        names.append("mu_M")
    return names


def _create_track_point_data(
    group: h5py.Group,
    npoints: int,
    values_dtype: np.dtype,
    vectors: dict[str, tuple[str, str, str]],
    scalars: Sequence[str],
    derived: Sequence[str],
    include_inside: bool,
) -> tuple[dict[str, h5py.Dataset], list[dict]]:
    datasets = {
        "time": group.create_dataset("time", shape=(npoints,), dtype=np.float64),
        "cycle": group.create_dataset("cycle", shape=(npoints,), dtype=np.int64),
    }
    metadata = [
        _array_meta("time", datasets["time"]),
        _array_meta("cycle", datasets["cycle"]),
    ]
    for name in scalars:
        dataset = group.create_dataset(
            _hdf_name(name), shape=(npoints,), dtype=values_dtype
        )
        datasets[name] = dataset
        metadata.append(_array_meta(name, dataset))
    for name, components in vectors.items():
        dataset = group.create_dataset(
            _hdf_name(name), shape=(npoints, 3), dtype=values_dtype
        )
        dataset.attrs["components"] = ",".join(components)
        datasets[name] = dataset
        metadata.append(_array_meta(name, dataset, "Vector"))
    for name in derived:
        dataset = group.create_dataset(
            _hdf_name(name), shape=(npoints,), dtype=values_dtype
        )
        datasets[name] = dataset
        metadata.append(_array_meta(name, dataset))
    if include_inside:
        dataset = group.create_dataset(
            "inside_meshblock", shape=(npoints,), dtype=np.uint8
        )
        datasets["inside_meshblock"] = dataset
        metadata.append(_array_meta("inside_meshblock", dataset))
    return datasets, metadata


def _create_track_cell_data(
    group: h5py.Group, tracks: dict, ncells: int
) -> tuple[dict[str, h5py.Dataset], list[dict]]:
    datasets = {
        "source_row": group.create_dataset(
            "source_row", shape=(ncells,), dtype=np.int64
        )
    }
    metadata = [_array_meta("source_row", datasets["source_row"])]
    particles = tracks["particles"]
    for name in particles.dtype.names or ():
        dataset = group.create_dataset(
            _hdf_name(name), shape=(ncells,), dtype=particles.dtype[name]
        )
        datasets[name] = dataset
        metadata.append(_array_meta(name, dataset))
    return datasets, metadata


def _write_track_point_values(
    datasets: dict[str, h5py.Dataset],
    destination: slice,
    values: np.ndarray,
    sample_times: np.ndarray,
    sample_cycles: np.ndarray,
    field: dict[str, int],
    vectors: dict[str, tuple[str, str, str]],
    scalars: Sequence[str],
    derived: Sequence[str],
    inside: np.ndarray | None = None,
) -> None:
    flat = values.reshape(-1, values.shape[-1])
    datasets["time"][destination] = sample_times.reshape(-1)
    datasets["cycle"][destination] = sample_cycles.reshape(-1)
    for name in scalars:
        datasets[name][destination] = flat[:, field[name]]
    for name, components in vectors.items():
        indices = [field[component] for component in components]
        datasets[name][destination] = flat[:, indices]

    if "magnetic_field_magnitude" in derived or "mu_M" in derived:
        magnetic = flat[:, [field[name] for name in ("bx", "by", "bz")]]
        bmag = np.linalg.norm(magnetic, axis=1)
        if "magnetic_field_magnitude" in derived:
            datasets["magnetic_field_magnitude"][destination] = bmag
    if "curvature_magnitude" in derived:
        curvature = flat[:, [field[name] for name in ("k1", "k2", "k3")]]
        datasets["curvature_magnitude"][destination] = np.linalg.norm(
            curvature, axis=1
        )
    if "mu_M" in derived:
        velocity = flat[:, [field[name] for name in ("vx", "vy", "vz")]]
        floor = np.array(1.0e-30, dtype=flat.dtype)
        safe_bmag = np.maximum(bmag, floor)
        v_parallel = np.sum(velocity * magnetic, axis=1) / safe_bmag
        v_perp2 = np.maximum(
            np.sum(velocity * velocity, axis=1) - v_parallel * v_parallel,
            floor,
        )
        datasets["mu_M"][destination] = v_perp2 / (2.0 * safe_bmag)
    if inside is not None:
        datasets["inside_meshblock"][destination] = inside.reshape(-1)


def _write_complete_tracks(
    group: h5py.Group,
    tracks: dict,
    particle_batch: int,
    domain_bounds: np.ndarray | None,
    periodic_axes: Sequence[int],
) -> dict | None:
    values = np.asarray(tracks["values"])
    nparticles, ntimes, _ = values.shape
    if nparticles == 0 or ntimes < 2:
        return None

    field, vectors, scalars = _track_arrays(tracks)
    if periodic_axes and domain_bounds is not None:
        coordinates = [field[name] for name in COORDINATE_FIELDS]
        runs = [
            _periodic_runs(
                values[particle][..., coordinates], domain_bounds, periodic_axes
            )
            for particle in range(nparticles)
        ]
        return _write_segmented_tracks(
            group, tracks, runs, "complete_periodic"
        )

    npoints = nparticles * ntimes
    derived = _derived_track_names(field)
    points = group.create_dataset(
        "points", shape=(npoints, 3), dtype=values.dtype
    )
    index_dtype = np.int32 if npoints <= np.iinfo(np.int32).max else np.int64
    connectivity = group.create_dataset(
        "connectivity", shape=(nparticles, ntimes), dtype=index_dtype
    )
    point_data, point_meta = _create_track_point_data(
        group.create_group("point_data"),
        npoints,
        values.dtype,
        vectors,
        scalars,
        derived,
        include_inside="inside_meshblock" in tracks,
    )
    cell_data, cell_meta = _create_track_cell_data(
        group.create_group("cell_data"), tracks, nparticles
    )

    coordinate_indices = [field[name] for name in COORDINATE_FIELDS]
    times = np.asarray(tracks["times"])
    cycles = np.asarray(tracks["cycles"])
    inside = tracks.get("inside_meshblock")
    for begin in range(0, nparticles, particle_batch):
        end = min(begin + particle_batch, nparticles)
        destination = slice(begin * ntimes, end * ntimes)
        batch = values[begin:end]
        points[destination] = batch[..., coordinate_indices].reshape(-1, 3)
        sample_times = np.broadcast_to(times, (end - begin, ntimes))
        sample_cycles = np.broadcast_to(cycles, (end - begin, ntimes))
        batch_inside = None if inside is None else np.asarray(inside[begin:end])
        _write_track_point_values(
            point_data,
            destination,
            batch,
            sample_times,
            sample_cycles,
            field,
            vectors,
            scalars,
            derived,
            batch_inside,
        )
        first = begin * ntimes
        last = end * ntimes
        connectivity[begin:end] = np.arange(
            first, last, dtype=index_dtype
        ).reshape(end - begin, ntimes)

    cell_data["source_row"][:] = tracks["source_rows"]
    for name in tracks["particles"].dtype.names or ():
        cell_data[name][:] = tracks["particles"][name]

    return {
        "geometry": "complete",
        "topology": "Polyline",
        "npoints": npoints,
        "ncells": nparticles,
        "nodes_per_cell": ntimes,
        "points": _array_meta("points", points),
        "connectivity": _array_meta("connectivity", connectivity),
        "point_arrays": point_meta,
        "cell_arrays": cell_meta,
    }


def _runs_from_breaks(breaks: np.ndarray, stop: int) -> np.ndarray:
    starts = np.concatenate((np.array([0]), breaks + 1))
    stops = np.concatenate((breaks + 1, np.array([stop])))
    return np.column_stack((starts, stops)).astype(np.int64, copy=False)


def _periodic_runs(
    points: np.ndarray,
    domain_bounds: np.ndarray,
    periodic_axes: Sequence[int],
) -> np.ndarray:
    bounds = np.asarray(domain_bounds)
    axes = np.asarray(periodic_axes, dtype=int)
    widths = bounds[axes, 1] - bounds[axes, 0]
    jumps = np.any(np.abs(np.diff(points[:, axes], axis=0)) > 0.5 * widths, axis=1)
    return _runs_from_breaks(np.flatnonzero(jumps), points.shape[0])


def _inside_runs(mask: np.ndarray) -> np.ndarray:
    changes = np.diff(np.pad(np.asarray(mask, dtype=np.int8), (1, 1)))
    starts = np.flatnonzero(changes == 1)
    stops = np.flatnonzero(changes == -1)
    return np.column_stack((starts, stops)).astype(np.int64, copy=False)


def _samples_for_runs(values: np.ndarray, runs: np.ndarray) -> np.ndarray:
    if runs.shape[0] == 1:
        return values[runs[0, 0]:runs[0, 1]]
    return np.concatenate([values[start:stop] for start, stop in runs])


def _write_segmented_tracks(
    group: h5py.Group,
    tracks: dict,
    particle_runs: Sequence[np.ndarray],
    geometry_name: str,
) -> dict | None:
    values = np.asarray(tracks["values"])
    run_lengths = [runs[:, 1] - runs[:, 0] for runs in particle_runs]
    npoints = int(sum(lengths.sum() for lengths in run_lengths))
    ncells = sum(len(lengths) for lengths in run_lengths)
    connectivity_size = npoints + 2 * ncells
    if npoints == 0 or ncells == 0:
        return None

    field, vectors, scalars = _track_arrays(tracks)
    derived = _derived_track_names(field)
    points = group.create_dataset(
        "points", shape=(npoints, 3), dtype=values.dtype
    )
    index_dtype = np.int32 if npoints <= np.iinfo(np.int32).max else np.int64
    connectivity = group.create_dataset(
        "connectivity", shape=(connectivity_size,), dtype=index_dtype
    )
    point_data, point_meta = _create_track_point_data(
        group.create_group("point_data"),
        npoints,
        values.dtype,
        vectors,
        scalars,
        derived,
        include_inside=False,
    )
    cell_data, cell_meta = _create_track_cell_data(
        group.create_group("cell_data"), tracks, ncells
    )

    coordinate_indices = [field[name] for name in COORDINATE_FIELDS]
    times = np.asarray(tracks["times"])
    cycles = np.asarray(tracks["cycles"])
    point_cursor = 0
    cell_cursor = 0
    connectivity_cursor = 0
    for particle, runs in enumerate(particle_runs):
        if not len(runs):
            continue
        sample = _samples_for_runs(values[particle], runs)
        sample_times = _samples_for_runs(times, runs)
        sample_cycles = _samples_for_runs(cycles, runs)
        point_slice = slice(point_cursor, point_cursor + sample.shape[0])
        points[point_slice] = sample[:, coordinate_indices]
        _write_track_point_values(
            point_data,
            point_slice,
            sample,
            sample_times,
            sample_cycles,
            field,
            vectors,
            scalars,
            derived,
        )

        local_point = point_cursor
        lengths = run_lengths[particle]
        mixed = np.empty(int(lengths.sum() + 2 * len(lengths)), dtype=index_dtype)
        local_connectivity = 0
        for length in lengths:
            length = int(length)
            cell_type = 1 if length == 1 else 2
            mixed[local_connectivity:local_connectivity + 2] = (cell_type, length)
            start = local_connectivity + 2
            mixed[start:start + length] = np.arange(
                local_point, local_point + length, dtype=index_dtype
            )
            local_point += int(length)
            local_connectivity = start + length
        mixed_slice = slice(connectivity_cursor, connectivity_cursor + mixed.size)
        connectivity[mixed_slice] = mixed

        count = len(runs)
        cell_slice = slice(cell_cursor, cell_cursor + count)
        cell_data["source_row"][cell_slice] = tracks["source_rows"][particle]
        for name in tracks["particles"].dtype.names or ():
            cell_data[name][cell_slice] = tracks["particles"][particle][name]
        point_cursor += sample.shape[0]
        cell_cursor += count
        connectivity_cursor += mixed.size

    return {
        "geometry": geometry_name,
        "topology": "Mixed",
        "npoints": npoints,
        "ncells": ncells,
        "nodes_per_cell": None,
        "points": _array_meta("points", points),
        "connectivity": _array_meta("connectivity", connectivity),
        "point_arrays": point_meta,
        "cell_arrays": cell_meta,
    }


def _write_inside_tracks(group: h5py.Group, tracks: dict) -> dict | None:
    values = np.asarray(tracks["values"])
    mask = np.asarray(tracks.get("inside_meshblock"))
    if mask.shape != values.shape[:2]:
        raise ValueError("inside track geometry requires inside_meshblock samples")
    return _write_segmented_tracks(
        group, tracks, [_inside_runs(row) for row in mask], "inside"
    )


def _write_tracks(
    parent: h5py.Group,
    tracks: dict,
    geometry: str,
    particle_batch: int,
    domain_bounds: np.ndarray | None,
    periodic_axes: Sequence[int],
) -> dict | None:
    if particle_batch < 1:
        raise ValueError("particle_batch must be at least one")
    group = parent.create_group("tracks")
    group.attrs["source_file"] = str(tracks.get("source_file", ""))
    group.attrs["fields"] = ",".join(tracks["fields"])
    group.attrs["periodic_axes"] = periodic_axes
    if domain_bounds is not None:
        group.attrs["domain_bounds"] = domain_bounds
    if geometry == "complete":
        return _write_complete_tracks(
            group, tracks, particle_batch, domain_bounds, periodic_axes
        )
    if geometry == "inside":
        return _write_inside_tracks(group, tracks)
    raise ValueError("track geometry must be 'complete' or 'inside'")


def write_visualization_piece(
    filename: str | Path,
    meshblocks: Sequence[dict],
    tracks: dict | None = None,
    track_geometry: str = "complete",
    particle_batch: int = 4,
    domain_bounds: np.ndarray | None = None,
    periodic_axes: Sequence[int] | None = None,
) -> dict:
    """Write one independent HDF5 piece and return its small XDMF metadata."""

    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "x") as handle:
        handle.attrs["format"] = "athenak_xdmf_payload_v1"
        mhd = handle.create_group("mhd")
        blocks = [
            _write_meshblock(mhd, block, index)
            for index, block in enumerate(meshblocks)
        ]
        track_meta = None
        if tracks is not None:
            first_block = meshblocks[0] if meshblocks else {}
            if domain_bounds is None:
                domain_bounds = first_block.get("DomainBounds")
            if periodic_axes is None:
                periodic_axes = first_block.get("PeriodicAxes", ())
            resolved_axes = tuple(periodic_axes)
            track_meta = _write_tracks(
                handle,
                tracks,
                track_geometry,
                particle_batch,
                domain_bounds,
                resolved_axes,
            )
    return {"filename": str(path.resolve()), "blocks": blocks, "tracks": track_meta}


def _dimensions(shape: Sequence[int]) -> str:
    return " ".join(str(value) for value in shape)


def _number_type(dtype_string: str) -> tuple[str, str]:
    dtype = np.dtype(dtype_string)
    if dtype.kind == "f":
        return "Float", str(dtype.itemsize)
    if dtype.kind == "i":
        return "Int", str(dtype.itemsize)
    if dtype.kind in ("u", "b"):
        return "UInt", str(dtype.itemsize)
    raise TypeError(f"unsupported XDMF dtype {dtype}")


def _hdf_reference(xmf_filename: Path, piece: dict, array: dict) -> str:
    relative = os.path.relpath(piece["filename"], xmf_filename.parent)
    return f"{Path(relative).as_posix()}:{array['path']}"


def _data_item(
    parent: ET.Element,
    xmf_filename: Path,
    piece: dict,
    array: dict,
) -> ET.Element:
    number_type, precision = _number_type(array["dtype"])
    item = ET.SubElement(
        parent,
        "DataItem",
        {
            "Dimensions": _dimensions(array["shape"]),
            "NumberType": number_type,
            "Precision": precision,
            "Format": "HDF",
        },
    )
    item.text = _hdf_reference(xmf_filename, piece, array)
    return item


def _add_mhd_grid(
    collection: ET.Element,
    xmf_filename: Path,
    piece: dict,
    block: dict,
    piece_index: int,
    block_index: int,
) -> None:
    name = f"piece_{piece_index:05d}_{block['name']}_{block_index:06d}"
    grid = ET.SubElement(
        collection, "Grid", {"Name": name, "GridType": "Uniform"}
    )
    ET.SubElement(
        grid,
        "Information",
        {"Name": "source_file", "Value": block["source_file"]},
    )
    ET.SubElement(
        grid,
        "Information",
        {"Name": "meshblock_index", "Value": str(block["meshblock_index"])},
    )
    ET.SubElement(
        grid,
        "Information",
        {
            "Name": "logical_location",
            "Value": _dimensions(block["logical_location"]),
        },
    )
    nz, ny, nx = block["shape"]
    ET.SubElement(
        grid,
        "Topology",
        {
            "TopologyType": "3DRectMesh",
            "Dimensions": _dimensions((nz + 1, ny + 1, nx + 1)),
        },
    )
    geometry = ET.SubElement(grid, "Geometry", {"GeometryType": "VXVYVZ"})
    for coordinate in block["coordinates"]:
        _data_item(geometry, xmf_filename, piece, coordinate)
    for array in block["arrays"]:
        attribute = ET.SubElement(
            grid,
            "Attribute",
            {
                "Name": array["name"],
                "AttributeType": array["attribute_type"],
                "Center": "Cell",
            },
        )
        _data_item(attribute, xmf_filename, piece, array)


def _add_track_grid(
    collection: ET.Element,
    xmf_filename: Path,
    piece: dict,
    tracks: dict,
    piece_index: int,
) -> None:
    grid = ET.SubElement(
        collection,
        "Grid",
        {"Name": f"tracks_piece_{piece_index:05d}", "GridType": "Uniform"},
    )
    ET.SubElement(
        grid,
        "Information",
        {"Name": "track_geometry", "Value": tracks["geometry"]},
    )
    topology_attributes = {
        "TopologyType": tracks["topology"],
        "NumberOfElements": str(tracks["ncells"]),
    }
    if tracks["nodes_per_cell"] is not None:
        topology_attributes["NodesPerElement"] = str(tracks["nodes_per_cell"])
    topology = ET.SubElement(grid, "Topology", topology_attributes)
    _data_item(topology, xmf_filename, piece, tracks["connectivity"])
    geometry = ET.SubElement(grid, "Geometry", {"GeometryType": "XYZ"})
    _data_item(geometry, xmf_filename, piece, tracks["points"])
    for center, arrays in (
        ("Node", tracks["point_arrays"]),
        ("Cell", tracks["cell_arrays"]),
    ):
        for array in arrays:
            attribute = ET.SubElement(
                grid,
                "Attribute",
                {
                    "Name": array["name"],
                    "AttributeType": array["attribute_type"],
                    "Center": center,
                },
            )
            _data_item(attribute, xmf_filename, piece, array)


def write_xdmf_collection(filename: str | Path, pieces: Sequence[dict]) -> Path:
    """Write one XMF file referencing independent HDF5 piece payloads."""

    path = xmf_path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    root = ET.Element("Xdmf", {"Version": "3.0"})
    domain = ET.SubElement(root, "Domain")
    mhd_collection = ET.SubElement(
        domain,
        "Grid",
        {"Name": "MHD", "GridType": "Collection", "CollectionType": "Spatial"},
    )
    track_collection = ET.SubElement(
        domain,
        "Grid",
        {
            "Name": "ParticleTracks",
            "GridType": "Collection",
            "CollectionType": "Spatial",
        },
    )
    for piece_index, piece in enumerate(pieces):
        for block_index, block in enumerate(piece["blocks"]):
            _add_mhd_grid(
                mhd_collection,
                path,
                piece,
                block,
                piece_index,
                block_index,
            )
        if piece["tracks"] is not None:
            _add_track_grid(
                track_collection,
                path,
                piece,
                piece["tracks"],
                piece_index,
            )

    ET.indent(root, space="  ")
    document = '<?xml version="1.0" ?>\n' + ET.tostring(
        root, encoding="unicode"
    )
    with path.open("x") as handle:
        handle.write(document)
        handle.write("\n")
    return path
