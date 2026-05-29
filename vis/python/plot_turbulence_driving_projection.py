#!/usr/bin/env python3
"""Plot projected turbulence forcing from AthenaK binary `turb_force` output.

Example:
  python vis/python/plot_turbulence_driving_projection.py \
      --panel "Constant Edot=run/bin/turb_edot_fixed.force.00004.bin" \
      --panel "Tiled include=run/bin/turb_tiled_include.force.00004.bin" \
      --panel "AMR exclude=run/bin/turb_amr_exclude.force.00004.bin" \
      --output turbulence_driving_projection_gallery.png
"""

import argparse
import struct

import matplotlib
import numpy as np

matplotlib.use("agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_parameter_dump(file_obj, header_offset):
    """Read the parameter text embedded in an AthenaK binary header."""
    parameter_text = file_obj.read(header_offset).decode("ascii")
    sections = {}
    section = None
    for raw_line in parameter_text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            section = line[1:-1]
            sections[section] = {}
        elif section is not None and "=" in line:
            key, value = line.split("=", 1)
            sections[section][key.strip()] = value.strip()
    return sections


def read_force_blocks(filename):
    """Read force arrays, geometry, and metadata from a binary output file."""
    with open(filename, "rb") as file_obj:
        if file_obj.readline().decode("ascii") != "Athena binary output version=1.1\n":
            raise RuntimeError(f"{filename}: not an AthenaK binary output")
        file_obj.readline()
        time = float(file_obj.readline().decode("ascii").split("=", 1)[1])
        cycle = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        location_size = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        variable_size = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        file_obj.readline()
        variables = file_obj.readline().decode("ascii").split(":", 1)[1].split()
        header_offset = int(file_obj.readline().decode("ascii").split("=", 1)[1])
        parameters = read_parameter_dump(file_obj, header_offset)

        if location_size not in (4, 8) or variable_size not in (4, 8):
            raise RuntimeError(f"{filename}: unsupported floating-point size")
        location_dtype = np.dtype("=f4" if location_size == 4 else "=f8")
        variable_dtype = np.dtype("=f4" if variable_size == 4 else "=f8")
        required = ("force1", "force2", "force3")
        try:
            force_indices = [variables.index(name) for name in required]
        except ValueError as error:
            raise RuntimeError(
                f"{filename}: output must contain turb_force"
            ) from error

        blocks = []
        while True:
            raw_indices = file_obj.read(6 * 4)
            if not raw_indices:
                break
            if len(raw_indices) != 6 * 4:
                raise RuntimeError(f"{filename}: truncated MeshBlock record")
            indices = struct.unpack("=6i", raw_indices)
            logical = struct.unpack("=4i", file_obj.read(4 * 4))
            limits = np.frombuffer(file_obj.read(6 * location_size), dtype=location_dtype)
            nx = indices[1] - indices[0] + 1
            ny = indices[3] - indices[2] + 1
            nz = indices[5] - indices[4] + 1
            count = nx * ny * nz
            data = np.frombuffer(
                file_obj.read(len(variables) * count * variable_size),
                dtype=variable_dtype,
            ).reshape((len(variables), nz, ny, nx))
            force = data[force_indices, :, :, :]
            blocks.append({"level": logical[3], "limits": limits, "force": force})

    return {"time": time, "cycle": cycle, "parameters": parameters, "blocks": blocks}


def project_force_magnitude(filename):
    """Integrate |a| along x3 and return a finest-level x1-x2 image."""
    output = read_force_blocks(filename)
    mesh = output["parameters"]["mesh"]
    nx_root = int(mesh["nx1"])
    ny_root = int(mesh["nx2"])
    x1min = float(mesh["x1min"])
    x1max = float(mesh["x1max"])
    x2min = float(mesh["x2min"])
    x2max = float(mesh["x2max"])
    x3min = float(mesh["x3min"])
    x3max = float(mesh["x3max"])
    max_level = max(block["level"] for block in output["blocks"])
    nx = nx_root * 2**max_level
    ny = ny_root * 2**max_level
    dx = (x1max - x1min) / nx
    dy = (x2max - x2min) / ny
    projected = np.zeros((ny, nx))
    levels = np.full((ny, nx), -1, dtype=int)

    for block in output["blocks"]:
        force = block["force"]
        level = block["level"]
        limits = block["limits"]
        dz = (limits[5] - limits[4]) / force.shape[1]
        magnitude = np.sqrt(np.sum(force * force, axis=0))
        plane = np.sum(magnitude, axis=0) * dz
        repeat = 2 ** (max_level - level)
        plane = np.repeat(np.repeat(plane, repeat, axis=0), repeat, axis=1)
        i0 = int(round((limits[0] - x1min) / dx))
        j0 = int(round((limits[2] - x2min) / dy))
        i1 = i0 + plane.shape[1]
        j1 = j0 + plane.shape[0]
        projected[j0:j1, i0:i1] += plane
        levels[j0:j1, i0:i1] = np.maximum(levels[j0:j1, i0:i1], level)

    projected /= x3max - x3min
    return projected, levels, (x1min, x1max, x2min, x2max), output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        action="append",
        required=True,
        metavar="LABEL=FILE",
        help="panel title and `turb_force` binary file; repeat for multiple panels",
    )
    parser.add_argument("--output", required=True, help="output PNG path")
    parser.add_argument("--dpi", type=int, default=180)
    args = parser.parse_args()

    panels = []
    for item in args.panel:
        if "=" not in item:
            parser.error("--panel must be written LABEL=FILE")
        label, filename = item.split("=", 1)
        panels.append((label, filename))

    figure, axes = plt.subplots(1, len(panels), figsize=(5.0 * len(panels), 4.2),
                                constrained_layout=True, squeeze=False)
    for axis, (label, filename) in zip(axes[0], panels):
        projection, levels, extent, output = project_force_magnitude(filename)
        image = axis.imshow(projection, origin="lower", extent=extent, cmap="magma")
        if levels.max() > 0:
            axis.contour(levels, levels=[0.5], origin="lower", extent=extent,
                         colors="cyan", linewidths=0.7)
        axis.set_title(f"{label}\ncycle {output['cycle']}")
        axis.set_xlabel("$x_1$")
        axis.set_ylabel("$x_2$")
        axis.set_aspect("equal")
        figure.colorbar(image, ax=axis, label=r"$L_3^{-1}\int |\mathbf{a}|\,dx_3$")
    figure.savefig(args.output, dpi=args.dpi)


if __name__ == "__main__":
    main()
