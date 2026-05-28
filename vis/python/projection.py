"""Read and compose AthenaK uniform or native-AMR projection output."""

import glob
from pathlib import Path
import struct

import numpy as np


def _partition_files(path):
    shard = path.parent.name
    if shard.startswith("rank_"):
        pattern = str(path.parent.parent / "rank_*" / path.name)
    elif shard.startswith("node_"):
        pattern = str(path.parent.parent / "node_*" / path.name)
    else:
        return [path]
    files = [Path(item) for item in sorted(glob.glob(pattern))]
    if not files:
        raise ValueError(f"no projection shard files match {pattern}")
    return files


def _parse_uniform(path):
    metadata = {}
    columns = None
    with path.open("r", encoding="ascii") as handle:
        for line in handle:
            if not line.startswith("#"):
                break
            content = line[1:].strip()
            if content.startswith("columns:"):
                columns = content.split(":", 1)[1].split()
                continue
            for item in content.split():
                if "=" in item:
                    key, value = item.split("=", 1)
                    metadata[key] = value
    if columns is None:
        raise ValueError(f"{path} does not contain a projection column header")
    raw = np.loadtxt(path, comments="#", ndmin=2)
    nx = int(metadata["nx"])
    ny = int(metadata["ny"])
    if raw.shape[0] != nx * ny:
        raise ValueError(f"{path} has {raw.shape[0]} samples, expected {nx}x{ny}")
    values = raw.reshape(ny, nx, raw.shape[1])
    return {
        "metadata": metadata,
        "columns": columns,
        "x": values[0, :, 2],
        "y": values[:, 0, 3],
        "fields": {
            name: values[:, :, index] for index, name in enumerate(columns[4:], start=4)
        },
        "patches": [],
    }


def _read_native_header(handle, path):
    first = handle.readline().decode("ascii").strip()
    if first != "AthenaK projection output version=3.0":
        raise ValueError(f"{path} is not a native-AMR projection file")
    metadata = {"version": "3.0"}
    variables = None
    while True:
        line = handle.readline()
        if not line:
            raise ValueError(f"unexpected EOF in header for {path}")
        text = line.decode("ascii").strip()
        if text.startswith("variables:"):
            variables = text.split(":", 1)[1].split()
            continue
        if "=" not in text:
            continue
        key, value = text.split("=", 1)
        metadata[key.strip()] = value.strip()
        if key.strip() == "header_offset":
            break
    if metadata.get("layout") != "native_amr" or variables is None:
        raise ValueError(f"invalid native-AMR projection header in {path}")
    handle.seek(int(metadata["header_offset"]), 1)
    return metadata, variables


def _read_native_file(path):
    patches = []
    with path.open("rb") as handle:
        metadata, variables = _read_native_header(handle, path)
        nvars = int(metadata["number_of_variables"])
        nmoments = int(metadata["number_of_moments"])
        if nvars != len(variables) or nmoments not in (2, 3):
            raise ValueError(f"invalid moment dimensions in {path}")
        for _ in range(int(metadata["npatches"])):
            record = handle.read(struct.calcsize("=4i4d"))
            if len(record) != struct.calcsize("=4i4d"):
                raise ValueError(f"truncated patch metadata in {path}")
            gid, level, nx, ny, xmin, xmax, ymin, ymax = struct.unpack("=4i4d", record)
            npixels = nx * ny
            weight = np.fromfile(handle, dtype=np.float64, count=npixels)
            first = np.fromfile(handle, dtype=np.float64, count=nvars * npixels)
            second = (
                np.fromfile(handle, dtype=np.float64, count=nvars * npixels)
                if nmoments == 3
                else None
            )
            if weight.size != npixels or first.size != nvars * npixels or (
                second is not None and second.size != nvars * npixels
            ):
                raise ValueError(f"truncated patch payload in {path}")
            patches.append(
                {
                    "gid": gid,
                    "level": level,
                    "nx": nx,
                    "ny": ny,
                    "xmin": xmin,
                    "xmax": xmax,
                    "ymin": ymin,
                    "ymax": ymax,
                    "weight": weight.reshape(ny, nx),
                    "first": first.reshape(nvars, ny, nx),
                    "second": None if second is None else second.reshape(nvars, ny, nx),
                }
            )
    return metadata, variables, patches


def _validate_native_headers(reference, candidate, path):
    keys = (
        "time",
        "cycle",
        "projection_axes",
        "weighting",
        "statistic",
        "image_dimension",
        "image_axis0",
        "image_axis1",
        "root_nx",
        "root_ny",
        "number_of_variables",
        "number_of_moments",
    )
    for key in keys:
        if reference.get(key) != candidate.get(key):
            raise ValueError(f"projection shard header mismatch for {key} in {path}")


def _compose_native(metadata, variables, patches, level):
    finest = max((patch["level"] for patch in patches), default=0)
    if level is None:
        level = finest
    if level < 0:
        raise ValueError("projection composition level must be non-negative")
    nx = int(metadata["root_nx"]) * (2**level)
    ny = (
        int(metadata["root_ny"]) * (2**level)
        if int(metadata["image_dimension"]) == 2
        else 1
    )
    xmin = float(metadata["domain_xmin"])
    xmax = float(metadata["domain_xmax"])
    ymin = float(metadata["domain_ymin"])
    ymax = float(metadata["domain_ymax"])
    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    nvars = len(variables)
    weight = np.zeros((ny, nx), dtype=np.float64)
    first = np.zeros((nvars, ny, nx), dtype=np.float64)
    second = (
        np.zeros((nvars, ny, nx), dtype=np.float64)
        if metadata["number_of_moments"] == "3"
        else None
    )
    for patch in patches:
        pdx = (patch["xmax"] - patch["xmin"]) / patch["nx"]
        pdy = (patch["ymax"] - patch["ymin"]) / patch["ny"]
        for py in range(patch["ny"]):
            plo_y = patch["ymin"] + py * pdy
            phi_y = plo_y + pdy
            iy0 = max(0, int(np.floor((plo_y - ymin) / dy + 1.0e-10)))
            iy1 = min(ny - 1, int(np.ceil((phi_y - ymin) / dy - 1.0e-10)) - 1)
            for px in range(patch["nx"]):
                plo_x = patch["xmin"] + px * pdx
                phi_x = plo_x + pdx
                ix0 = max(0, int(np.floor((plo_x - xmin) / dx + 1.0e-10)))
                ix1 = min(nx - 1, int(np.ceil((phi_x - xmin) / dx - 1.0e-10)) - 1)
                for iy in range(iy0, iy1 + 1):
                    fy = min(phi_y, ymin + (iy + 1) * dy) - max(plo_y, ymin + iy * dy)
                    fy /= dy
                    for ix in range(ix0, ix1 + 1):
                        fx = min(phi_x, xmin + (ix + 1) * dx) - max(
                            plo_x, xmin + ix * dx
                        )
                        fraction = fy * fx / dx
                        weight[iy, ix] += patch["weight"][py, px] * fraction
                        first[:, iy, ix] += patch["first"][:, py, px] * fraction
                        if second is not None:
                            second[:, iy, ix] += patch["second"][:, py, px] * fraction
    weighting = metadata["weighting"]
    statistic = metadata["statistic"]
    fields = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        for n, name in enumerate(variables):
            if weighting == "integral":
                fields[name] = first[n]
            elif statistic == "stddev":
                mean = first[n] / weight
                fields[name] = np.sqrt(np.maximum(second[n] / weight - mean * mean, 0.0))
            else:
                fields[name] = first[n] / weight
    composed = dict(metadata)
    composed.update(
        {
            "nx": str(nx),
            "ny": str(ny),
            "xmin": str(xmin),
            "xmax": str(xmax),
            "ymin": str(ymin),
            "ymax": str(ymax),
            "projection_level": str(level),
        }
    )
    x = xmin + (np.arange(nx, dtype=np.float64) + 0.5) * dx
    y = ymin + (np.arange(ny, dtype=np.float64) + 0.5) * dy
    return {
        "metadata": composed,
        "columns": ["i", "j", "x", "y"] + variables,
        "x": x,
        "y": y,
        "fields": fields,
        "moments": {"weight": weight, "first": first, "second": second},
        "patches": patches,
    }


def read_projection(filename, level=None):
    """Read a projection and compose native-AMR shards at ``level`` when requested."""
    path = Path(filename)
    with path.open("rb") as handle:
        first = handle.readline()
    if first.startswith(b"#"):
        if level is not None:
            raise ValueError(
                "level selection is only supported for native-AMR projections"
            )
        return _parse_uniform(path)

    files = _partition_files(path)
    metadata, variables, patches = _read_native_file(files[0])
    for other in files[1:]:
        other_metadata, other_variables, other_patches = _read_native_file(other)
        _validate_native_headers(metadata, other_metadata, other)
        if variables != other_variables:
            raise ValueError(f"projection variable mismatch in {other}")
        patches.extend(other_patches)
    return _compose_native(metadata, variables, patches, level)


def plot_projection(filename, variable=None, output=None, cmap="magma", level=None):
    """Plot one reduced field and optionally write the figure to ``output``."""
    import matplotlib.pyplot as plt

    projection = read_projection(filename, level=level)
    fields = projection["fields"]
    if variable is None:
        variable = next(iter(fields))
    if variable not in fields:
        raise KeyError(f"{variable!r} is not present; choices are {list(fields)}")
    image = fields[variable]
    dimension = int(projection["metadata"].get("image_dimension", "2"))
    fig, ax = plt.subplots(figsize=(6.0, 5.0), constrained_layout=True)
    if dimension == 1:
        ax.plot(projection["x"], image[0], color="tab:blue")
        ax.set_xlabel(projection["metadata"]["image_axis0"])
        ax.set_ylabel(variable)
    else:
        extent = [
            float(projection["metadata"]["xmin"]),
            float(projection["metadata"]["xmax"]),
            float(projection["metadata"]["ymin"]),
            float(projection["metadata"]["ymax"]),
        ]
        rendered = ax.imshow(
            image, origin="lower", extent=extent, cmap=cmap, aspect="auto"
        )
        ax.set_xlabel(projection["metadata"]["image_axis0"])
        ax.set_ylabel(projection["metadata"]["image_axis1"])
        fig.colorbar(rendered, ax=ax, label=variable)
    axes = projection["metadata"].get("projection_axes", "")
    ax.set_title(
        f"{variable}: {projection['metadata'].get('statistic', 'value')} "
        f"({projection['metadata']['weighting']}, axes={axes})"
    )
    if output is not None:
        fig.savefig(output, dpi=180)
    return fig, ax
