#!/usr/bin/env python3
"""Plot the Sedov-Taylor inner_x1 boundary state used by cloud_crushing.cpp.

By default this reads `inputs/hydro/cloud_crushing_snr.athinput`, mirrors the
Sedov-boundary formulas in `src/pgen/cloud_crushing.cpp`, and plots the primitive
state and fluxes imposed at a selected inner_x1 ghost-cell location.
The boundary implementation injects the full xi-dependent
Taylor-von Neumann-Sedov interior profile inside the shock radius.

Examples:

  python3 scripts/plot_sedov_boundary.py

  python3 scripts/plot_sedov_boundary.py \
      --sedov-energy-cgs 3e50 --sedov-origin-distance 50 \
      --sedov-radius-at-start 51 --x2 1.0 --output sedov_boundary.png
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = REPO_ROOT / "inputs" / "hydro" / "cloud_crushing_snr.athinput"
ISM_COOLING = REPO_ROOT / "src" / "srcterms" / "ismcooling.hpp"

AMU_CGS = 1.67262192369e-24
KBOLTZ_CGS = 1.3806488e-16


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


def value_from(
    args: argparse.Namespace,
    attr: str,
    blocks: dict[str, dict[str, str]],
    block: str,
    key: str,
    default: float,
) -> float:
  cli_value = getattr(args, attr)
  if cli_value is not None:
    return float(cli_value)
  if key in blocks.get(block, {}):
    return parse_float(blocks[block][key])
  return default


def load_ism_lhd() -> np.ndarray:
  text = ISM_COOLING.read_text()
  match = re.search(r"lhd\s*\[[^\]]*\]\s*=\s*\{(.*?)\};", text, re.S)
  if not match:
    raise ValueError(f"could not locate lhd array in {ISM_COOLING}")
  data = np.fromstring(match.group(1), sep=",")
  if data.size != 102:
    raise ValueError(f"unexpected lhd array size {data.size}")
  return data


def ism_cooling(temp: np.ndarray | float, lhd: np.ndarray) -> np.ndarray | float:
  scalar = np.isscalar(temp)
  t = np.asarray(temp, dtype=float)
  logt = np.log10(t)
  rate = np.empty_like(t)

  low = logt <= 4.2
  if np.any(low):
    tlow = t[low]
    rate[low] = (
        2.0e-19 * np.exp(-1.184e5 / (tlow + 1.0e3))
        + 2.8e-28 * np.sqrt(tlow) * np.exp(-92.0 / tlow)
    )

  high = logt > 8.15
  if np.any(high):
    rate[high] = 10.0 ** (0.45 * logt[high] - 26.065)

  mid = ~(low | high)
  if np.any(mid):
    lt = logt[mid]
    ipps = np.clip((25.0 * lt - 103).astype(int), 0, 100)
    x0 = 4.12 + 0.04 * ipps.astype(float)
    dx = lt - x0
    logcool = (lhd[ipps + 1] * dx - lhd[ipps] * (dx - 0.04)) * 25.0
    rate[mid] = 10.0 ** logcool

  return float(rate) if scalar else rate


def equilibrium_function(temp: float, pressure_over_k: float, hrate: float, lhd: np.ndarray) -> float:
  return pressure_over_k * float(ism_cooling(temp, lhd)) / temp - hrate


def bisect_root(lo: float, hi: float, pressure_over_k: float, hrate: float, lhd: np.ndarray) -> float:
  flo = equilibrium_function(lo, pressure_over_k, hrate, lhd)
  for _ in range(80):
    mid = math.sqrt(lo * hi)
    fmid = equilibrium_function(mid, pressure_over_k, hrate, lhd)
    if flo * fmid <= 0.0:
      hi = mid
    else:
      lo = mid
      flo = fmid
  return math.sqrt(lo * hi)


def stable_equilibria(pressure_over_k: float, hrate: float, lhd: np.ndarray) -> tuple[float, float]:
  temp_min = 1.0
  temp_max = 1.0e7
  roots: list[float] = []
  prev_t = temp_min
  prev_f = equilibrium_function(prev_t, pressure_over_k, hrate, lhd)
  for n in range(1, 4097):
    frac = n / 4096.0
    temp = temp_min * (temp_max / temp_min) ** frac
    f = equilibrium_function(temp, pressure_over_k, hrate, lhd)
    if prev_f == 0.0 or f == 0.0 or prev_f * f < 0.0:
      root = bisect_root(prev_t, temp, pressure_over_k, hrate, lhd)
      if prev_f < 0.0 and f > 0.0:
        roots.append(root)
    prev_t = temp
    prev_f = f
  if len(roots) < 2:
    raise ValueError("could not find both cold and warm stable equilibria")
  return roots[0], roots[-1]


def cell_centers(n: int, xmin: float, xmax: float) -> np.ndarray:
  return xmin + (np.arange(n, dtype=float) + 0.5) * (xmax - xmin) / n


def tvns_log_xi5_from_log_eta(log_eta: np.ndarray, gamma: float) -> np.ndarray:
  nu1 = -(13.0 * gamma**2 - 7.0 * gamma + 12.0) / ((3.0 * gamma - 1.0) * (2.0 * gamma + 1.0))
  nu2 = 5.0 * (gamma - 1.0) / (2.0 * gamma + 1.0)
  eta = np.exp(log_eta)
  similarity_v = (1.0 + eta) / gamma
  a = 0.5 * (gamma + 1.0) * similarity_v
  b = ((gamma + 1.0) / (7.0 - gamma)) * (5.0 - (3.0 * gamma - 1.0) * similarity_v)
  c = ((gamma + 1.0) / (gamma - 1.0)) * eta
  return -2.0 * np.log(a) + nu1 * np.log(b) + nu2 * np.log(c)


def tvns_profile(xi: np.ndarray, gamma: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
  if gamma <= 1.0 or gamma >= 2.0:
    raise ValueError("the analytic TVNS Sedov profile requires 1 < gamma < 2")

  xi_floor = 1.0e-12
  xi_limited = np.clip(np.asarray(xi, dtype=float), xi_floor, 1.0)
  target = 5.0 * np.log(xi_limited)
  eta_shock = (gamma - 1.0) / (gamma + 1.0)
  log_eta_shock = math.log(eta_shock)
  lo = np.full_like(xi_limited, -700.0)
  hi = np.full_like(xi_limited, log_eta_shock)

  for _ in range(80):
    mid = 0.5 * (lo + hi)
    move_lo = tvns_log_xi5_from_log_eta(mid, gamma) < target
    lo = np.where(move_lo, mid, lo)
    hi = np.where(move_lo, hi, mid)

  log_eta = 0.5 * (lo + hi)
  eta = np.exp(log_eta)
  similarity_v = (1.0 + eta) / gamma
  nu1 = -(13.0 * gamma**2 - 7.0 * gamma + 12.0) / ((3.0 * gamma - 1.0) * (2.0 * gamma + 1.0))
  nu3 = 3.0 / (2.0 * gamma + 1.0)
  nu4 = -nu1 / (2.0 - gamma)
  nu5 = -2.0 / (2.0 - gamma)
  compression = (gamma + 1.0) / (gamma - 1.0)
  a = ((gamma + 1.0) / (gamma - 1.0)) * eta
  b = ((gamma + 1.0) / (7.0 - gamma)) * (5.0 - (3.0 * gamma - 1.0) * similarity_v)
  c = ((gamma + 1.0) / (gamma - 1.0)) * (1.0 - similarity_v)
  density_ratio = compression * a**nu3 * b**nu4 * c**nu5
  sound_speed_shape = (
      gamma
      * (gamma - 1.0)
      * (1.0 - similarity_v)
      * similarity_v**2
      / (2.0 * eta)
  )
  return density_ratio, similarity_v, sound_speed_shape


def build_series(args: argparse.Namespace) -> dict[str, Any]:
  blocks = parse_athinput(args.input) if args.input else {}
  lhd = load_ism_lhd()

  x1min = value_from(args, "x1min", blocks, "mesh", "x1min", -2.0)
  x1max = value_from(args, "x1max", blocks, "mesh", "x1max", 10.0)
  x2min = value_from(args, "x2min", blocks, "mesh", "x2min", -3.0)
  x2max = value_from(args, "x2max", blocks, "mesh", "x2max", 3.0)
  x3min = value_from(args, "x3min", blocks, "mesh", "x3min", -3.0)
  x3max = value_from(args, "x3max", blocks, "mesh", "x3max", 3.0)
  nx1 = int(value_from(args, "nx1", blocks, "mesh", "nx1", 192))
  nx2 = int(value_from(args, "nx2", blocks, "mesh", "nx2", 64))
  nx3 = int(value_from(args, "nx3", blocks, "mesh", "nx3", 64))
  nghost = int(value_from(args, "nghost", blocks, "mesh", "nghost", 2))

  length_cgs = value_from(args, "length_cgs", blocks, "units", "length_cgs", 1.0)
  time_cgs = value_from(args, "time_cgs", blocks, "units", "time_cgs", 1.0)
  mass_cgs = value_from(args, "mass_cgs", blocks, "units", "mass_cgs", 1.0)
  mu = value_from(args, "mu", blocks, "units", "mu", 1.0)
  gamma = value_from(args, "gamma", blocks, "hydro", "gamma", 5.0 / 3.0)
  velocity_unit = length_cgs / time_cgs

  pressure_over_k = value_from(
      args, "pressure_over_k", blocks, "problem", "pressure_over_k", 3162.277660168379
  )
  hrate = value_from(args, "hrate", blocks, "hydro", "hrate", 2.0e-26)

  cold_temp, warm_temp = stable_equilibria(pressure_over_k, hrate, lhd)
  cold_number_density = pressure_over_k / cold_temp
  warm_number_density = pressure_over_k / warm_temp

  cold_mass_density_cgs = cold_number_density * mu * AMU_CGS
  ambient_mass_density_cgs = warm_number_density * mu * AMU_CGS
  ambient_pressure_cgs = pressure_over_k * KBOLTZ_CGS

  energy_cgs = value_from(args, "sedov_energy_cgs", blocks, "problem", "sedov_energy_cgs", 1.0e51)
  beta = value_from(args, "sedov_beta", blocks, "problem", "sedov_beta", 1.15167)
  start_time = value_from(args, "sedov_start_time", blocks, "problem", "sedov_start_time", 0.0)
  origin_distance = value_from(
      args, "sedov_origin_distance", blocks, "problem", "sedov_origin_distance", 30.0
  )

  if args.sedov_origin_x1 is not None:
    origin_x1 = args.sedov_origin_x1
  elif "sedov_origin_x1" in blocks.get("problem", {}):
    origin_x1 = parse_float(blocks["problem"]["sedov_origin_x1"])
  else:
    origin_x1 = x1min - origin_distance
  origin_x2 = value_from(args, "sedov_origin_x2", blocks, "problem", "sedov_origin_x2", 0.0)
  origin_x3 = value_from(args, "sedov_origin_x3", blocks, "problem", "sedov_origin_x3", 0.0)
  radius_at_start = value_from(
      args, "sedov_radius_at_start", blocks, "problem", "sedov_radius_at_start", origin_distance
  )

  radius_start_cgs = radius_at_start * length_cgs
  age_start_cgs = math.sqrt((radius_start_cgs / beta) ** 5 * ambient_mass_density_cgs / energy_cgs)

  tmax = args.t_max
  if tmax is None:
    tmax = 10.0
  if tmax <= args.t_min:
    raise ValueError("--t-max must be greater than --t-min")
  if args.samples < 2:
    raise ValueError("--samples must be at least 2")
  times = np.linspace(args.t_min, tmax, args.samples)
  times_cgs = times * time_cgs
  ages = age_start_cgs + np.maximum(times - start_time, 0.0) * time_cgs
  radii_cgs = beta * (energy_cgs * ages**2 / ambient_mass_density_cgs) ** 0.2
  radii_code = radii_cgs / length_cgs
  shock_speed_cgs = 0.4 * radii_cgs / ages

  dx1 = (x1max - x1min) / nx1
  ghost_index = args.ghost_index
  if ghost_index < 1 or ghost_index > nghost:
    raise ValueError(f"--ghost-index must be between 1 and nghost={nghost}")
  x1 = x1min - (ghost_index - 0.5) * dx1
  x2 = args.x2
  x3 = args.x3
  frame_v1 = args.frame_velocity_x1
  frame_v2 = args.frame_velocity_x2
  frame_v3 = args.frame_velocity_x3
  frame_x1 = args.frame_displacement_x1 + frame_v1 * times
  frame_x2 = args.frame_displacement_x2 + frame_v2 * times
  frame_x3 = args.frame_displacement_x3 + frame_v3 * times
  frame_v1_cgs = frame_v1 * velocity_unit
  frame_v2_cgs = frame_v2 * velocity_unit
  frame_v3_cgs = frame_v3 * velocity_unit
  dx1_sample = x1 + frame_x1 - origin_x1
  dx2_sample = x2 + frame_x2 - origin_x2
  dx3_sample = x3 + frame_x3 - origin_x3
  radius_sample = np.sqrt(dx1_sample**2 + dx2_sample**2 + dx3_sample**2)
  radius_sample_cgs = radius_sample * length_cgs
  shocked = radii_code >= radius_sample
  xi = radius_sample / radii_code
  xi_profile = np.maximum(xi, 1.0e-12)
  density_ratio, velocity_shape, sound_speed_shape = tvns_profile(xi, gamma)
  radial_x = np.divide(dx1_sample, radius_sample, out=np.zeros_like(radius_sample), where=radius_sample > 0.0)
  radial_y = np.divide(dx2_sample, radius_sample, out=np.zeros_like(radius_sample), where=radius_sample > 0.0)
  radial_z = np.divide(dx3_sample, radius_sample, out=np.zeros_like(radius_sample), where=radius_sample > 0.0)

  shocked_density = ambient_mass_density_cgs * density_ratio
  shocked_pressure = shocked_density * (shock_speed_cgs * xi_profile) ** 2 * sound_speed_shape / gamma
  shocked_radial_speed = shock_speed_cgs * xi * velocity_shape
  rho = np.where(shocked, shocked_density, ambient_mass_density_cgs)
  pressure = np.where(shocked, shocked_pressure, ambient_pressure_cgs)
  vx = np.where(shocked, shocked_radial_speed * radial_x - frame_v1_cgs, -frame_v1_cgs)
  vy = np.where(shocked, shocked_radial_speed * radial_y - frame_v2_cgs, -frame_v2_cgs)
  vz = np.where(shocked, shocked_radial_speed * radial_z - frame_v3_cgs, -frame_v3_cgs)
  mass_flux = rho * vx
  momentum_flux = rho * vx**2 + pressure

  y_centers = cell_centers(nx2, x2min, x2max)
  z_centers = cell_centers(nx3, x3min, x3max)
  yy, zz = np.meshgrid(y_centers, z_centers, indexing="xy")
  shocked_fraction = np.array([
      np.mean(
          np.sqrt(
              (x1 + frame_x1[i] - origin_x1) ** 2
              + (yy + frame_x2[i] - origin_x2) ** 2
              + (zz + frame_x3[i] - origin_x3) ** 2
          )
          <= radii_code[i]
      )
      for i in range(args.samples)
  ])

  return {
      "times_cgs": times_cgs,
      "sedov_xi": xi,
      "sample_radius_cgs": radius_sample_cgs,
      "frame_x1_cgs": frame_x1 * length_cgs,
      "frame_x2_cgs": frame_x2 * length_cgs,
      "frame_x3_cgs": frame_x3 * length_cgs,
      "frame_v1_cgs": np.full_like(times_cgs, frame_v1_cgs),
      "frame_v2_cgs": np.full_like(times_cgs, frame_v2_cgs),
      "frame_v3_cgs": np.full_like(times_cgs, frame_v3_cgs),
      "radius_cgs": radii_cgs,
      "shock_speed_cgs": shock_speed_cgs,
      "rho": rho,
      "pressure": pressure,
      "vx": vx,
      "vy": vy,
      "vz": vz,
      "mass_flux": mass_flux,
      "momentum_flux": momentum_flux,
      "shocked": shocked.astype(float),
      "shocked_fraction": shocked_fraction,
      "metadata": {
          "cold_temp": cold_temp,
          "warm_temp": warm_temp,
          "cold_number_density": cold_number_density,
          "warm_number_density": warm_number_density,
          "cold_mass_density_cgs": cold_mass_density_cgs,
          "ambient_mass_density_cgs": ambient_mass_density_cgs,
          "ambient_pressure_cgs": ambient_pressure_cgs,
          "x1": x1,
          "x2": x2,
          "x3": x3,
          "x1_cgs": x1 * length_cgs,
          "x2_cgs": x2 * length_cgs,
          "x3_cgs": x3 * length_cgs,
          "radius_sample_initial": radius_sample[0],
          "radius_sample_final": radius_sample[-1],
          "radius_sample_initial_cgs": radius_sample_cgs[0],
          "radius_sample_final_cgs": radius_sample_cgs[-1],
          "frame_velocity_cgs": (frame_v1_cgs, frame_v2_cgs, frame_v3_cgs),
          "age_start_code": age_start_cgs / time_cgs,
          "age_start_cgs": age_start_cgs,
          "energy_cgs": energy_cgs,
          "gamma": gamma,
          "origin": (origin_x1, origin_x2, origin_x3),
          "origin_cgs": (origin_x1 * length_cgs, origin_x2 * length_cgs, origin_x3 * length_cgs),
          "radius_at_start": radius_at_start,
          "radius_at_start_cgs": radius_at_start * length_cgs,
      },
  }


def write_csv(path: Path, series: dict[str, Any]) -> None:
  fields = [
      "time_s",
      "sedov_xi",
      "sample_radius_cm",
      "frame_x1_cm",
      "frame_x2_cm",
      "frame_x3_cm",
      "frame_v1_cm_s",
      "frame_v2_cm_s",
      "frame_v3_cm_s",
      "shock_radius_cm",
      "shock_speed_cm_s",
      "rho_g_cm3",
      "pressure_dyn_cm2",
      "vx_cm_s",
      "vy_cm_s",
      "vz_cm_s",
      "mass_flux_g_cm2_s",
      "momentum_flux_dyn_cm2",
      "sample_is_shocked",
      "inner_x1_shocked_fraction",
  ]
  path.parent.mkdir(parents=True, exist_ok=True)
  with path.open("w", newline="") as stream:
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(fields)
    for i, time in enumerate(series["times_cgs"]):
      writer.writerow([
          time,
          series["sedov_xi"][i],
          series["sample_radius_cgs"][i],
          series["frame_x1_cgs"][i],
          series["frame_x2_cgs"][i],
          series["frame_x3_cgs"][i],
          series["frame_v1_cgs"][i],
          series["frame_v2_cgs"][i],
          series["frame_v3_cgs"][i],
          series["radius_cgs"][i],
          series["shock_speed_cgs"][i],
          series["rho"][i],
          series["pressure"][i],
          series["vx"][i],
          series["vy"][i],
          series["vz"][i],
          series["mass_flux"][i],
          series["momentum_flux"][i],
          series["shocked"][i],
          series["shocked_fraction"][i],
      ])


def plot_positive(ax: Any, x: np.ndarray, y: np.ndarray, label: str, **kwargs: Any) -> None:
  values = np.asarray(y, dtype=float)
  mask = values > 0.0
  if np.any(mask):
    ax.plot(x[mask], values[mask], label=label, **kwargs)


def plot_reference_line(ax: Any, y: float, label: str, **kwargs: Any) -> None:
  if y > 0.0:
    ax.axhline(y, label=label, **kwargs)


def plot_series(path: Path, series: dict[str, Any]) -> None:
  time = series["times_cgs"]
  meta = series["metadata"]

  fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
  ax = axes[0, 0]
  plot_positive(ax, time, series["radius_cgs"], label="shock radius")
  plot_positive(ax, time, series["sample_radius_cgs"], color="0.35", ls="--", label="sample radius")
  ax.set_xlabel("time [s]")
  ax.set_ylabel("radius [cm]")
  ax.set_yscale("log")
  ax.legend(loc="best")

  ax = axes[0, 1]
  plot_positive(ax, time, series["shock_speed_cgs"], label="shock speed")
  plot_positive(ax, time, series["vx"], label="boundary vx")
  ax.set_xlabel("time [s]")
  ax.set_ylabel("speed [cm s^-1]")
  ax.set_yscale("log")
  ax.legend(loc="best")

  ax = axes[1, 0]
  plot_positive(ax, time, series["rho"], label="TVNS boundary density [g cm^-3]")
  plot_positive(ax, time, series["pressure"], label="TVNS boundary pressure [dyn cm^-2]")
  plot_reference_line(
      ax,
      meta["ambient_mass_density_cgs"],
      color="tab:blue",
      ls=":",
      label="ambient ISM density [g cm^-3]",
  )
  plot_reference_line(
      ax,
      meta["cold_mass_density_cgs"],
      color="tab:cyan",
      ls=":",
      label="cold cloud density [g cm^-3]",
  )
  plot_reference_line(
      ax,
      meta["ambient_pressure_cgs"],
      color="tab:orange",
      ls=":",
      label="ambient ISM pressure [dyn cm^-2]",
  )
  ax.set_xlabel("time [s]")
  ax.set_ylabel("boundary primitive/reference values [CGS]")
  ax.set_yscale("log")
  ax.legend(loc="best")

  ax = axes[1, 1]
  plot_positive(ax, time, series["mass_flux"], label="rho vx [g cm^-2 s^-1]")
  plot_positive(ax, time, series["momentum_flux"], label="rho vx^2 + P [dyn cm^-2]")
  plot_positive(ax, time, series["shocked_fraction"], label="shocked area fraction")
  ax.set_xlabel("time [s]")
  ax.set_ylabel("flux / fraction [CGS or dimensionless]")
  ax.set_yscale("log")
  ax.legend(loc="best")

  fig.suptitle(
      "Sedov inner_x1 boundary at "
      f"(x1,x2,x3)=({meta['x1_cgs']:.3e},{meta['x2_cgs']:.3e},{meta['x3_cgs']:.3e}) cm; "
      f"E={meta['energy_cgs']:.2e} erg, R_start={meta['radius_at_start_cgs']:.3e} cm"
  )
  path.parent.mkdir(parents=True, exist_ok=True)
  fig.savefig(path, dpi=180)
  plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="AthenaK input file")
  parser.add_argument("--output", type=Path, default=Path("sedov_inner_x1_boundary.png"))
  parser.add_argument("--csv-output", type=Path, default=None)

  parser.add_argument("--t-min", type=float, default=0.0)
  parser.add_argument("--t-max", type=float, default=10.0)
  parser.add_argument("--samples", type=int, default=512)
  parser.add_argument("--ghost-index", type=int, default=1, help="1 is the innermost ghost cell")
  parser.add_argument("--x2", type=float, default=0.0)
  parser.add_argument("--x3", type=float, default=0.0)
  parser.add_argument("--frame-displacement-x1", type=float, default=0.0)
  parser.add_argument("--frame-displacement-x2", type=float, default=0.0)
  parser.add_argument("--frame-displacement-x3", type=float, default=0.0)
  parser.add_argument("--frame-velocity-x1", type=float, default=0.0)
  parser.add_argument("--frame-velocity-x2", type=float, default=0.0)
  parser.add_argument("--frame-velocity-x3", type=float, default=0.0)

  for name in (
      "x1min",
      "x1max",
      "x2min",
      "x2max",
      "x3min",
      "x3max",
      "nx1",
      "nx2",
      "nx3",
      "nghost",
      "length_cgs",
      "time_cgs",
      "mass_cgs",
      "mu",
      "gamma",
      "pressure_over_k",
      "hrate",
      "sedov_energy_cgs",
      "sedov_beta",
      "sedov_start_time",
      "sedov_origin_distance",
      "sedov_origin_x1",
      "sedov_origin_x2",
      "sedov_origin_x3",
      "sedov_radius_at_start",
  ):
    parser.add_argument(f"--{name.replace('_', '-')}", type=float, default=None)

  return parser


def main() -> None:
  args = build_parser().parse_args()
  series = build_series(args)
  plot_series(args.output, series)
  if args.csv_output is not None:
    write_csv(args.csv_output, series)

  meta = series["metadata"]
  print(f"wrote {args.output}")
  if args.csv_output is not None:
    print(f"wrote {args.csv_output}")
  print(
      "equilibria: "
      f"T_cold={meta['cold_temp']:.6g} K, "
      f"T_warm={meta['warm_temp']:.6g} K, "
      f"n_cold={meta['cold_number_density']:.6g} cm^-3, "
      f"n_warm={meta['warm_number_density']:.6g} cm^-3"
  )
  print(
      "ambient/cold CGS: "
      f"P_ISM={meta['ambient_pressure_cgs']:.6g} dyn cm^-2, "
      f"rho_ISM={meta['ambient_mass_density_cgs']:.6g} g cm^-3, "
      f"rho_cold={meta['cold_mass_density_cgs']:.6g} g cm^-3"
  )
  print(
      "sample: "
      f"x=({meta['x1_cgs']:.6g}, {meta['x2_cgs']:.6g}, {meta['x3_cgs']:.6g}) cm, "
      f"radius_from_origin={meta['radius_sample_initial_cgs']:.6g} -> "
      f"{meta['radius_sample_final_cgs']:.6g} cm, "
      f"age_start={meta['age_start_cgs']:.6g} s"
  )
  print(
      "frame: "
      f"v=({meta['frame_velocity_cgs'][0]:.6g}, "
      f"{meta['frame_velocity_cgs'][1]:.6g}, "
      f"{meta['frame_velocity_cgs'][2]:.6g}) cm s^-1"
  )


if __name__ == "__main__":
  main()
