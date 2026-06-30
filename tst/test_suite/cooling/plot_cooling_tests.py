#!/usr/bin/env python3
"""
Generate analytic-vs-numerical plots for the cooling regression tests.

The script runs the built-in cooling_test pgen with the already-built AthenaK
binary and writes documentation figures under docs/source/modules/figures/.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


ROOT = Path(__file__).resolve().parents[3]
ATHENA = ROOT / "build" / "src" / "athena"
RUN_DIR = ROOT / "build" / "cooling_plot_runs"
FIG_DIR = ROOT / "docs" / "source" / "modules" / "figures"
HYDRO_INPUT = "inputs/cooling_test_hydro.athinput"
MHD_INPUT = "inputs/cooling_test_mhd.athinput"
COMPOSITION_INPUT = "inputs/cooling_test_hydro_composition.athinput"
CONVERGENCE_TABLE_DIR = RUN_DIR / "convergence_tables"

sys.path.insert(0, str(ROOT / "vis" / "python"))
import athena_read  # noqa: E402


@dataclass(frozen=True)
class RateCase:
  label: str
  basename: str
  input_file: str
  physics: str
  flags: tuple[str, ...]
  gross: float
  net: float


def prepare_run_dir() -> None:
  if not ATHENA.exists():
    raise FileNotFoundError(f"Build AthenaK first; missing {ATHENA}")
  RUN_DIR.mkdir(parents=True, exist_ok=True)
  FIG_DIR.mkdir(parents=True, exist_ok=True)
  inputs_link = RUN_DIR / "inputs"
  target = ROOT / "tst" / "inputs"
  if inputs_link.is_symlink() or not inputs_link.exists():
    if inputs_link.exists() or inputs_link.is_symlink():
      inputs_link.unlink()
    inputs_link.symlink_to(target, target_is_directory=True)
  elif inputs_link.resolve() != target:
    raise RuntimeError(f"{inputs_link} exists and does not point to {target}")


def run_case(basename: str, input_file: str, physics: str,
             flags: tuple[str, ...]) -> dict[str, np.ndarray]:
  for path in RUN_DIR.glob(f"{basename}.*.hst"):
    path.unlink()
  command = [str(ATHENA), "-i", input_file, f"job/basename={basename}", *flags]
  result = subprocess.run(command, cwd=RUN_DIR, text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, check=False)
  if result.returncode != 0:
    print(result.stdout)
    print(result.stderr, file=sys.stderr)
    raise RuntimeError(f"AthenaK failed for {basename}")
  return athena_read.hst(str(RUN_DIR / f"{basename}.{physics}.hst"), raw=True)


def strictly_evolved_rows(data: dict[str, np.ndarray]) -> np.ndarray:
  rows = [0]
  last_time = data["time"][0]
  for n in range(1, data["time"].size):
    if data["time"][n] > last_time:
      rows.append(n)
      last_time = data["time"][n]
  return np.array(rows, dtype=int)


def plot_constant_solution() -> Path:
  rate = 1.0e-2
  fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.0), sharex=True,
                           constrained_layout=True)
  colors = {"rk1": "#2a6f97", "rk2": "#b35c00", "rk3": "#5a7f2b"}

  for integrator in ("rk1", "rk2", "rk3"):
    data = run_case(
        f"cooling_plot_constant_{integrator}",
        HYDRO_INPUT,
        "hydro",
        (
            f"time/integrator={integrator}",
            "time/tlim=2.0e-2",
            "time/nlim=64",
            "cooling/timestep=true",
            "cooling/timestep_factor=1.6666666666666667e-5",
        ),
    )
    rows = strictly_evolved_rows(data)
    time = data["time"][rows]
    energy = data["tot-E"][rows]
    loss = energy[0] - energy
    exact_loss = rate*time
    axes[0].plot(time, exact_loss, color="0.18", lw=1.4, alpha=0.55)
    axes[0].plot(time, loss, "o", ms=3.2, color=colors[integrator],
                 label=integrator)

    interval_rows = rows[1:]
    axes[1].plot(data["time"][interval_rows], data["cool_net"][interval_rows] - rate,
                 "o-", ms=3.0, lw=1.1, color=colors[integrator],
                 label=integrator)
    axes[2].plot(time, loss - exact_loss, "o-", ms=3.0, lw=1.1,
                 color=colors[integrator], label=integrator)

  axes[0].set_ylabel("energy lost")
  axes[0].set_title("Constant cooling: numerical solution lies on analytic loss")
  axes[0].legend(title="integrator", ncols=3, frameon=False)
  axes[1].axhline(0.0, color="0.18", lw=1.2, label="analytic")
  axes[1].set_ylabel("net rate error")
  axes[1].legend(frameon=False)
  axes[2].axhline(0.0, color="0.18", lw=1.0)
  axes[2].set_ylabel("loss error")
  axes[2].set_xlabel("time")
  for ax in axes:
    ax.grid(True, color="0.88", lw=0.8)

  out = FIG_DIR / "cooling_constant_solution.png"
  fig.savefig(out, dpi=180)
  plt.close(fig)
  return out


def rate_cases() -> tuple[RateCase, ...]:
  tiny_step = ("time/tlim=1.0e-8", "time/nlim=1")
  return (
      RateCase("power law", "cooling_plot_powerlaw", HYDRO_INPUT, "hydro",
               tiny_step, 1.0e-2, 1.0e-2),
      RateCase("table 1D", "cooling_plot_table_1d", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=table",
                            "cooling/cooling_table=inputs/tables/cooling_lambda_1d.tbl"),
               1.0e-2, 1.0e-2),
      RateCase("table 2D", "cooling_plot_table_2d", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=table",
                            "cooling/cooling_table=inputs/tables/cooling_lambda_2d.tbl"),
               1.0e-2, 1.0e-2),
      RateCase("table 3D", "cooling_plot_table_3d", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=table",
                            "cooling/cooling_table=inputs/tables/cooling_lambda_3d.tbl"),
               1.0e-2, 1.0e-2),
      RateCase("default log", "cooling_plot_default_log", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=table",
                            "cooling/cooling_table=inputs/tables/"
                            "cooling_lambda_default_log.tbl",
                            "problem/pressure=10.0"),
               3.0e-2, 3.0e-2),
      RateCase("mixed axes", "cooling_plot_mixed_axes", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=table",
                            "cooling/cooling_table=inputs/tables/"
                            "cooling_lambda_mixed_axes.tbl",
                            "problem/pressure=10.0"),
               3.0e-2, 3.0e-2),
      RateCase("piecewise", "cooling_plot_piecewise", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=piecewise_powerlaw",
                            "cooling/cooling_breaks=0.5,1.5",
                            "cooling/cooling_slopes=0.0,0.0,0.0"),
               1.0e-2, 1.0e-2),
      RateCase("heat const", "cooling_plot_heat_const", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/heating_model=constant",
                            "cooling/heating_rate=2.0e-3"),
               1.0e-2, 8.0e-3),
      RateCase("heat table", "cooling_plot_heat_table", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/heating_model=table",
                            "cooling/heating_table=inputs/tables/cooling_gamma_1d.tbl"),
               1.0e-2, 8.0e-3),
      RateCase("modifier", "cooling_plot_modifier", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_modifier_model=constant",
                            "cooling/cooling_modifier=2.0"),
               2.0e-2, 2.0e-2),
      RateCase("user hook", "cooling_plot_user", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=user",
                            "cooling/heating_model=user",
                            "problem/register_user_cooling=true",
                            "problem/user_lambda=1.0e-2",
                            "problem/user_gamma=2.0e-3"),
               1.0e-2, 8.0e-3),
      RateCase("CGM", "cooling_plot_cgm", HYDRO_INPUT, "hydro",
               tiny_step + ("cooling/cooling_model=cgm",
                            "cooling/cgm_pie_table=inputs/tables/cooling_cgm_pie.tbl",
                            "cooling/cgm_cie_table=inputs/tables/cooling_cgm_cie.tbl",
                            "cooling/cooling_density=mass_density",
                            "cooling/shielding_model=none"),
               1.0e-2, 1.0e-2),
      RateCase("CGM shield", "cooling_plot_cgm_shield", COMPOSITION_INPUT, "hydro",
               tiny_step, 3.0e-2, 3.0e-2),
      RateCase("MHD", "cooling_plot_mhd", MHD_INPUT, "mhd",
               tiny_step, 1.0e-2, 1.0e-2),
  )


def plot_rate_sweep() -> Path:
  cases = rate_cases()
  labels = [case.label for case in cases]
  gross_exact = np.array([case.gross for case in cases])
  net_exact = np.array([case.net for case in cases])
  gross_num = []
  net_num = []

  for case in cases:
    data = run_case(case.basename, case.input_file, case.physics, case.flags)
    row = strictly_evolved_rows(data)[1]
    gross_num.append(data["cool_gross"][row])
    net_num.append(data["cool_net"][row])
  gross_num = np.array(gross_num)
  net_num = np.array(net_num)

  x = np.arange(len(cases))
  width = 0.36
  fig, axes = plt.subplots(2, 1, figsize=(10.5, 7.2), sharex=True,
                           constrained_layout=True,
                           gridspec_kw={"height_ratios": [2.2, 1.0]})
  axes[0].bar(x - width/2, gross_exact, width, color="0.82",
              edgecolor="0.25", linewidth=0.8, label="analytic gross")
  axes[0].bar(x + width/2, net_exact, width, color="0.95",
              edgecolor="0.25", linewidth=0.8, label="analytic net")
  axes[0].plot(x - width/2, gross_num, "o", color="#2a6f97",
               label="numerical gross")
  axes[0].plot(x + width/2, net_num, "s", color="#b35c00",
               label="numerical net")
  axes[0].set_ylabel("rate")
  axes[0].set_title("Cooling model sweep: numerical rates match analytic expectations")
  axes[0].legend(ncols=2, frameon=False)
  axes[0].grid(True, axis="y", color="0.88", lw=0.8)

  axes[1].plot(x, gross_num - gross_exact, "o-", color="#2a6f97",
               label="gross error")
  axes[1].plot(x, net_num - net_exact, "s-", color="#b35c00",
               label="net error")
  axes[1].axhline(0.0, color="0.18", lw=1.0)
  axes[1].set_ylabel("rate error")
  axes[1].set_xticks(x)
  axes[1].set_xticklabels(labels, rotation=35, ha="right")
  axes[1].grid(True, axis="y", color="0.88", lw=0.8)
  axes[1].legend(frameon=False)

  out = FIG_DIR / "cooling_model_rate_errors.png"
  fig.savefig(out, dpi=180)
  plt.close(fig)
  return out


def analytic_table_lambda(temperature: np.ndarray | float) -> np.ndarray | float:
  return 1.0e-2*(1.0 + 0.2*np.asarray(temperature)**2)


def write_convergence_table(naxis: int) -> Path:
  CONVERGENCE_TABLE_DIR.mkdir(parents=True, exist_ok=True)
  path = CONVERGENCE_TABLE_DIR / f"lambda_linear_n{naxis}.tbl"
  coords = np.linspace(0.0, 2.0, naxis)
  values = analytic_table_lambda(coords)
  lines = [
      "ATHENAK_COOLING_TABLE 1",
      "ndim 1",
      "value_kind Lambda",
      "value_units code",
      "value_scale linear",
      "bounds fatal",
      "data_order c_row_major_last_axis_fastest",
      f"axis0 temperature code linear 0.0 2.0 {naxis}",
      "data",
      *(f"{value:.17e}" for value in values),
  ]
  path.write_text("\n".join(lines) + "\n")
  return path


def plot_table_interpolation_convergence() -> Path:
  target_temperature = 0.73
  exact = float(analytic_table_lambda(target_temperature))
  resolutions = np.array([5, 9, 17, 33, 65])
  dx = 2.0/(resolutions - 1)
  errors = []

  for naxis in resolutions:
    table = write_convergence_table(int(naxis))
    data = run_case(
        f"cooling_plot_interp_n{naxis}",
        HYDRO_INPUT,
        "hydro",
        (
            "time/tlim=1.0e-8",
            "time/nlim=1",
            "cooling/cooling_model=table",
            f"cooling/cooling_table={table.relative_to(RUN_DIR)}",
            f"problem/pressure={target_temperature:.17e}",
        ),
    )
    row = strictly_evolved_rows(data)[1]
    errors.append(abs(data["cool_gross"][row] - exact))

  errors = np.array(errors)
  reference = errors[0]*(dx/dx[0])**2
  fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
  ax.loglog(dx, errors, "o-", color="#2a6f97", label="numerical error")
  ax.loglog(dx, reference, "--", color="0.35", label="second-order guide")
  ax.invert_xaxis()
  ax.set_xlabel("table spacing")
  ax.set_ylabel("absolute rate error")
  ax.set_title("Table interpolation convergence for a known analytic curve")
  ax.grid(True, which="both", color="0.88", lw=0.8)
  ax.legend(frameon=False)
  out = FIG_DIR / "cooling_table_interpolation_convergence.png"
  fig.savefig(out, dpi=180)
  plt.close(fig)
  return out


def main() -> None:
  prepare_run_dir()
  outputs = [
      plot_constant_solution(),
      plot_rate_sweep(),
      plot_table_interpolation_convergence(),
  ]
  for path in outputs:
    print(path.relative_to(ROOT))


if __name__ == "__main__":
  main()
