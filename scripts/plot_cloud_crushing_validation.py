#!/usr/bin/env python3
"""Plot low-resolution cloud-crushing frame-tracking diagnostics.

This script reads an AthenaK run directory produced by
`inputs/hydro/cloud_crushing_snr.athinput`, parses the `FrameTracker` diagnostic
blocks from an Athena stdout log, and reads native binary `hydro_w` snapshots.
It writes one time-series validation plot and one x-y midplane density plot.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputs" / "hydro" / "cloud_crushing_snr.athinput"
ATOMIC_MASS_UNIT_CGS = 1.66053906660e-24
MSUN_CGS = 1.98847e33
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new  # noqa: E402


def strip_comment(line: str) -> str:
  return line.split("#", 1)[0].strip()


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
  blocks: dict[str, dict[str, str]] = {}
  block: str | None = None
  for raw_line in path.read_text().splitlines():
    line = strip_comment(raw_line)
    if not line:
      continue
    if line.startswith("<") and line.endswith(">"):
      block = line[1:-1].strip()
      blocks.setdefault(block, {})
      continue
    if block is None or "=" not in line:
      continue
    key, value = line.split("=", 1)
    blocks.setdefault(block, {})[key.strip()] = value.strip()
  return blocks


def parse_float(value: str) -> float:
  return float(value.split()[0])


def block_float(
    blocks: dict[str, dict[str, str]],
    block: str,
    key: str,
    default: float,
) -> float:
  if key in blocks.get(block, {}):
    return parse_float(blocks[block][key])
  return default


def parse_frame_log(path: Path) -> list[dict[str, Any]]:
  cycle_re = re.compile(
      r"cycle=\s*([0-9]+)\s+time=([+\-0-9.eE]+)\s+dt=([+\-0-9.eE]+)"
  )
  field_re = re.compile(r"\s*([A-Za-z0-9_]+)=([+\-0-9.eE]+)\s*$")
  rows: list[dict[str, Any]] = []
  current: dict[str, Any] | None = None
  cycle: int | None = None
  time = math.nan
  dt = math.nan

  for line in path.read_text().splitlines():
    cycle_match = cycle_re.search(line)
    if cycle_match:
      cycle = int(cycle_match.group(1))
      time = float(cycle_match.group(2))
      dt = float(cycle_match.group(3))
      continue

    if line.startswith("FrameTracker "):
      if current is not None:
        rows.append(current)
      current = {
          "status": line[len("FrameTracker ") :].strip(),
          "cycle": cycle,
          "time": time,
          "dt": dt,
      }
      continue

    field_match = field_re.match(line)
    if current is not None and field_match:
      current[field_match.group(1)] = float(field_match.group(2))

  if current is not None:
    rows.append(current)
  return rows


def read_binary_snapshot(path: Path) -> dict[str, Any]:
  if "rank_00000000" in path.parts:
    data = bin_convert_new.read_all_ranks_binary_as_athdf(str(path), quantities=["dens"])
  else:
    data = bin_convert_new.read_binary_as_athdf(str(path), quantities=["dens"])

  rho = data.get("dens")
  if rho is None:
    raise ValueError(f"snapshot {path} does not contain dens")

  x = np.asarray(data["x1v"], dtype=float)
  y = np.asarray(data["x2v"], dtype=float)
  z = np.asarray(data["x3v"], dtype=float)
  zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
  spacing = (
      float(np.mean(np.diff(np.asarray(data["x1f"], dtype=float)))),
      float(np.mean(np.diff(np.asarray(data["x2f"], dtype=float)))),
      float(np.mean(np.diff(np.asarray(data["x3f"], dtype=float)))),
  )
  return {
      "path": path,
      "time": float(data["Time"]),
      "shape": tuple(int(n) for n in rho.shape),
      "spacing": spacing,
      "x1v": x,
      "x2v": y,
      "x3v": z,
      "x1f": np.asarray(data["x1f"], dtype=float),
      "x2f": np.asarray(data["x2f"], dtype=float),
      "x3f": np.asarray(data["x3f"], dtype=float),
      "x": xx.ravel(),
      "y": yy.ravel(),
      "z": zz.ravel(),
      "fields": {"dens": np.asarray(rho, dtype=float).ravel()},
  }


def read_binary_series(run_dir: Path) -> list[dict[str, Any]]:
  files = sorted((run_dir / "bin").glob("*.bin"))
  if not files:
    files = sorted((run_dir / "bin" / "rank_00000000").glob("*.bin"))
  if not files:
    files = sorted(run_dir.glob("*.bin"))
  return [read_binary_snapshot(path) for path in files]


def cold_centroids(
    snapshots: list[dict[str, Any]],
    rho_min: float,
    rho_max: float,
) -> dict[str, np.ndarray]:
  times: list[float] = []
  xs: list[float] = []
  masses: list[float] = []
  for snap in snapshots:
    rho = snap["fields"].get("dens")
    if rho is None:
      continue
    mask = (rho >= rho_min) & (rho <= rho_max)
    dvol = float(np.prod(np.asarray(snap["spacing"])))
    weight = np.where(mask, rho*dvol, 0.0)
    total = float(np.sum(weight))
    if total <= 0.0:
      continue
    times.append(float(snap["time"]))
    xs.append(float(np.sum(weight*snap["x"])/total))
    masses.append(total)
  return {
      "time": np.asarray(times),
      "x1_centroid": np.asarray(xs),
      "weight": np.asarray(masses),
  }


def dense_mass_history(
    snapshots: list[dict[str, Any]],
    rho_min: float,
    mass_unit_cgs: float,
) -> dict[str, np.ndarray]:
  times: list[float] = []
  mass_code: list[float] = []
  for snap in snapshots:
    rho = snap["fields"].get("dens")
    if rho is None:
      continue
    mask = rho > rho_min
    dvol = float(np.prod(np.asarray(snap["spacing"])))
    mass = float(np.sum(rho[mask])*dvol)
    times.append(float(snap["time"]))
    mass_code.append(mass)
  mass_code_array = np.asarray(mass_code)
  mass_cgs = mass_code_array*mass_unit_cgs
  return {
      "time": np.asarray(times),
      "mass_code": mass_code_array,
      "mass_g": mass_cgs,
      "mass_msun": mass_cgs/MSUN_CGS,
  }


def write_dense_mass_csv(path: Path, history: dict[str, np.ndarray]) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w") as stream:
    stream.write("time_Myr,mass_code,mass_g,mass_Msun\n")
    for time, mass_code, mass_g, mass_msun in zip(
        history["time"],
        history["mass_code"],
        history["mass_g"],
        history["mass_msun"],
    ):
      stream.write(
          f"{time:.16e},{mass_code:.16e},{mass_g:.16e},{mass_msun:.16e}\n"
      )


def plot_dense_mass(path: Path, history: dict[str, np.ndarray], rho_min: float) -> None:
  path.parent.mkdir(parents=True, exist_ok=True)
  fig, ax = plt.subplots(figsize=(7.2, 4.4), constrained_layout=True)
  ax.plot(history["time"], history["mass_msun"], marker="o", ms=3.0, lw=1.5)
  ax.set_xlabel("time [Myr]")
  ax.set_ylabel(r"$M(\rho > \rho_{\rm cloud,0}/3)$ [$M_\odot$]")
  ax.set_title(f"Dense gas mass, threshold rho > {rho_min:.6g} [code]")
  ax.grid(True, alpha=0.3)
  fig.savefig(path, dpi=180)
  plt.close(fig)


def select_evenly_spaced(items: list[dict[str, Any]], count: int) -> list[dict[str, Any]]:
  if count >= len(items):
    return items
  indices = np.linspace(0, len(items) - 1, count)
  return [items[int(round(index))] for index in indices]


def density_plane(snapshot: dict[str, Any]) -> np.ndarray:
  rho = snapshot["fields"].get("dens")
  if rho is None:
    raise ValueError(f"snapshot {snapshot['path']} does not contain dens")
  rho_3d = np.asarray(rho, dtype=float).reshape(snapshot["shape"])
  k = int(np.argmin(np.abs(snapshot["x3v"])))
  return rho_3d[k, :, :]


def plot_density_slices_vertical(
    path: Path,
    snapshots: list[dict[str, Any]],
    target_x1: float,
    number_density_unit: float,
    slice_count: int,
) -> None:
  if not snapshots:
    raise ValueError("no binary snapshots found")
  selected = select_evenly_spaced(snapshots, slice_count)
  panels = [
      np.log10(np.maximum(density_plane(snapshot)*number_density_unit, 1.0e-99))
      for snapshot in selected
  ]
  vmin = min(float(np.nanmin(panel)) for panel in panels)
  vmax = max(float(np.nanmax(panel)) for panel in panels)
  extent = [
      float(selected[0]["x1f"][0]),
      float(selected[0]["x1f"][-1]),
      float(selected[0]["x2f"][0]),
      float(selected[0]["x2f"][-1]),
  ]

  path.parent.mkdir(parents=True, exist_ok=True)
  fig_height = max(7.0, 1.65*len(selected) + 1.2)
  fig, axes = plt.subplots(
      len(selected),
      1,
      figsize=(14.0, fig_height),
      constrained_layout=True,
      squeeze=False,
  )
  image = None
  for ax, panel, snapshot in zip(axes[:, 0], panels, selected):
    image = ax.imshow(
        panel,
        origin="lower",
        extent=extent,
        interpolation="nearest",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    ax.axvline(target_x1, color="white", lw=1.0, ls="--", alpha=0.75)
    ax.set_aspect("equal", adjustable="box")
    ax.set_ylabel("x2 [pc]")
    ax.set_title(f"t = {float(snapshot['time']):.4f} Myr   {snapshot['path'].name}")
  axes[-1, 0].set_xlabel("x1 [pc]")
  if image is not None:
    fig.colorbar(
        image,
        ax=axes[:, 0].tolist(),
        label=r"$\log_{10} n$ [cm$^{-3}$]",
        fraction=0.025,
        pad=0.015,
    )
  fig.savefig(path, dpi=180)
  plt.close(fig)


def finite_values(rows: list[dict[str, Any]], key: str) -> tuple[np.ndarray, np.ndarray]:
  x: list[float] = []
  y: list[float] = []
  for row in rows:
    value = row.get(key)
    time = row.get("time")
    if value is None or time is None:
      continue
    if math.isfinite(float(value)) and math.isfinite(float(time)):
      x.append(float(time))
      y.append(float(value))
  return np.asarray(x), np.asarray(y)


def plot_tracker_validation(
    path: Path,
    rows: list[dict[str, Any]],
    snapshot_centroids: dict[str, np.ndarray],
    target_x1: float,
) -> None:
  if not rows:
    raise ValueError("no FrameTracker diagnostic rows found")
  path.parent.mkdir(parents=True, exist_ok=True)
  fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)

  ax = axes[0, 0]
  for key, label, style in (
      ("x1_centroid", "tracker centroid", "-"),
      ("x1_control_position", "control position", "--"),
      ("x1_filtered_position", "filtered position", ":"),
  ):
    time, values = finite_values(rows, key)
    if values.size:
      ax.plot(time, values, style, label=label)
  if snapshot_centroids["time"].size:
    ax.scatter(
        snapshot_centroids["time"],
        snapshot_centroids["x1_centroid"],
        s=30,
        color="black",
        label="binary cold centroid",
        zorder=4,
    )
  ax.axhline(target_x1, color="0.35", lw=1.2, label="target")
  ax.set_xlabel("time [code]")
  ax.set_ylabel("x1 position [code]")
  ax.legend(loc="best")

  ax = axes[0, 1]
  time, values = finite_values(rows, "x1_frame_velocity")
  if values.size:
    ax.plot(time, values, color="tab:blue", label="frame velocity")
  ax.axhline(0.0, color="0.45", lw=0.8)
  ax.set_xlabel("time [code]")
  ax.set_ylabel("frame velocity [code]")
  ax2 = ax.twinx()
  for key, label, color, style in (
      ("x1_frame_displacement", "frame displacement", "tab:orange", "-"),
      ("x1_boost", "applied boost", "tab:green", "--"),
  ):
    time, values = finite_values(rows, key)
    if values.size:
      ax2.plot(time, values, color=color, ls=style, label=label)
  ax2.set_ylabel("displacement / boost [code]")
  lines, labels = ax.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax.legend(lines + lines2, labels + labels2, loc="best")

  ax = axes[1, 0]
  time, values = finite_values(rows, "x1_position_error")
  if values.size:
    ax.plot(time, values, color="tab:blue", label="position error")
  ax.axhline(0.0, color="0.45", lw=0.8)
  ax.set_xlabel("time [code]")
  ax.set_ylabel("position error [code]")
  ax2 = ax.twinx()
  for key, label, color, style in (
      ("x1_filtered_velocity", "filtered velocity", "tab:orange", "-"),
      ("x1_mean_velocity", "tracked mean velocity", "tab:red", "--"),
  ):
    time, values = finite_values(rows, key)
    if values.size:
      ax2.plot(time, values, color=color, ls=style, label=label)
  ax2.set_ylabel("velocity [code]")
  lines, labels = ax.get_legend_handles_labels()
  lines2, labels2 = ax2.get_legend_handles_labels()
  ax.legend(lines + lines2, labels + labels2, loc="best")

  ax = axes[1, 1]
  time, weight = finite_values(rows, "x1_global_weight")
  if weight.size:
    ax.plot(time, weight, label="tracker weight")
  if snapshot_centroids["time"].size:
    ax.scatter(
        snapshot_centroids["time"],
        snapshot_centroids["weight"],
        s=30,
        color="black",
        label="binary cold weight",
        zorder=4,
    )
  ax.set_xlabel("time [code]")
  ax.set_ylabel("selected density weight [code]")
  ax.set_yscale("log")
  ax.legend(loc="best")

  fig.suptitle("Cloud-crushing frame-tracking validation")
  fig.savefig(path, dpi=180)
  plt.close(fig)


def midplane_points(snapshot: dict[str, Any], density_unit: float) -> dict[str, np.ndarray]:
  rho = snapshot["fields"].get("dens")
  if rho is None:
    raise ValueError(f"snapshot {snapshot['path']} does not contain dens")
  z = snapshot["z"]
  z0 = z[np.argmin(np.abs(z))]
  mask = np.isclose(z, z0)
  rho_cgs = np.maximum(rho[mask]*density_unit, 1.0e-99)
  return {
      "x": snapshot["x"][mask],
      "y": snapshot["y"][mask],
      "logrho": np.log10(rho_cgs),
      "time": np.asarray([snapshot["time"]]),
  }


def plot_midplane_density(
    path: Path,
    snapshots: list[dict[str, Any]],
    snapshot_centroids: dict[str, np.ndarray],
    target_x1: float,
    density_unit: float,
) -> None:
  if not snapshots:
    raise ValueError("no binary snapshots found")
  selected = [snapshots[0]]
  if len(snapshots) > 1:
    selected.append(snapshots[-1])

  panels = [midplane_points(snapshot, density_unit) for snapshot in selected]
  vmin = min(float(np.min(panel["logrho"])) for panel in panels)
  vmax = max(float(np.max(panel["logrho"])) for panel in panels)
  path.parent.mkdir(parents=True, exist_ok=True)
  fig, axes = plt.subplots(1, len(panels), figsize=(5.6*len(panels), 4.8), squeeze=False)
  image = None
  for ax, panel, snapshot in zip(axes[0], panels, selected):
    image = ax.scatter(
        panel["x"],
        panel["y"],
        c=panel["logrho"],
        s=18,
        marker="s",
        cmap="viridis",
        vmin=vmin,
        vmax=vmax,
    )
    ax.axvline(target_x1, color="white", lw=1.2, ls="--", label="target")
    if snapshot_centroids["time"].size:
      centroid = np.interp(
          float(snapshot["time"]),
          snapshot_centroids["time"],
          snapshot_centroids["x1_centroid"],
      )
      ax.axvline(centroid, color="black", lw=1.2, ls="-", label="cold centroid")
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x1 [code]")
    ax.set_ylabel("x2 [code]")
    ax.set_title(f"t={float(snapshot['time']):.4e}")
    ax.legend(loc="upper right")
  if image is not None:
    fig.colorbar(image, ax=axes.ravel().tolist(), label=r"$\log_{10}\rho$ [g cm$^{-3}$]")
  fig.suptitle("Cloud-crushing midplane density")
  fig.savefig(path, dpi=180)
  plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("run_dir", type=Path, help="AthenaK run directory")
  parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
  parser.add_argument("--log", type=Path, default=None, help="Athena stdout log")
  parser.add_argument(
      "--output-prefix",
      type=Path,
      default=Path("cloud_crushing_lowres"),
      help="Prefix for output PNG files",
  )
  parser.add_argument("--target-min", type=float, default=None)
  parser.add_argument("--target-max", type=float, default=None)
  parser.add_argument("--slice-count", type=int, default=6)
  return parser


def main() -> None:
  args = build_parser().parse_args()
  blocks = parse_athinput(args.input)
  log_path = args.log if args.log is not None else args.run_dir / "athena.log"
  rows = parse_frame_log(log_path)
  snapshots = read_binary_series(args.run_dir)

  target_min = args.target_min
  if target_min is None:
    target_min = block_float(blocks, "frame_tracking", "target_min", 5.0)
  target_max = args.target_max
  if target_max is None:
    target_max = block_float(blocks, "frame_tracking", "target_max", 200.0)
  target_x1 = block_float(blocks, "frame_tracking", "x1_target", 3.0)
  length_cgs = block_float(blocks, "units", "length_cgs", 1.0)
  mass_cgs = block_float(blocks, "units", "mass_cgs", 1.0)
  mu = block_float(blocks, "units", "mu", 1.0)
  density_unit = mass_cgs/(length_cgs**3)
  number_density_unit = density_unit/(mu*ATOMIC_MASS_UNIT_CGS)

  snapshot_centroids = cold_centroids(snapshots, target_min, target_max)
  dense_mass = dense_mass_history(snapshots, target_min, mass_cgs)
  validation_path = args.output_prefix.with_name(args.output_prefix.name + "_validation.png")
  midplane_path = args.output_prefix.with_name(args.output_prefix.name + "_midplane.png")
  dense_mass_path = args.output_prefix.with_name(args.output_prefix.name + "_dense_mass.png")
  dense_mass_csv_path = args.output_prefix.with_name(args.output_prefix.name + "_dense_mass.csv")
  slices_path = args.output_prefix.with_name(
      args.output_prefix.name
      + f"_density_slices_{args.slice_count}_vertical_equal_aspect.png"
  )
  plot_tracker_validation(validation_path, rows, snapshot_centroids, target_x1)
  plot_midplane_density(midplane_path, snapshots, snapshot_centroids, target_x1, density_unit)
  plot_dense_mass(dense_mass_path, dense_mass, target_min)
  write_dense_mass_csv(dense_mass_csv_path, dense_mass)
  plot_density_slices_vertical(
      slices_path,
      snapshots,
      target_x1,
      number_density_unit,
      args.slice_count,
  )

  print(f"read {len(rows)} frame-tracking diagnostic block(s)")
  print(f"read {len(snapshots)} binary snapshot(s)")
  print(f"wrote {validation_path}")
  print(f"wrote {midplane_path}")
  print(f"wrote {dense_mass_path}")
  print(f"wrote {dense_mass_csv_path}")
  print(f"wrote {slices_path}")
  if rows:
    last = rows[-1]
    print(
        "final frame state: "
        f"time={last.get('time', math.nan):.6g}, "
        f"v_frame={last.get('x1_frame_velocity', math.nan):.6g}, "
        f"x_frame={last.get('x1_frame_displacement', math.nan):.6g}"
    )
  if snapshot_centroids["time"].size:
    print(
        "final binary cold centroid: "
        f"time={snapshot_centroids['time'][-1]:.6g}, "
        f"x1={snapshot_centroids['x1_centroid'][-1]:.6g}"
    )
  if dense_mass["time"].size:
    print(
        "final dense mass: "
        f"time={dense_mass['time'][-1]:.6g}, "
        f"M={dense_mass['mass_msun'][-1]:.6g} Msun"
    )


if __name__ == "__main__":
  main()
