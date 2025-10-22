#!/usr/bin/env python3
"""One-off utility to visualise AthenaK cooling curves.

The script reads the lookup tables defined in
`src/srcterms/ismcooling.hpp` and `src/srcterms/cooling_tables.hpp`
and produces:

* `ism_cooling_curve.png` – ISM cooling coefficient vs. temperature
* `cgm_cooling_curves.png` – CGM PIE cooling for H+He and metals at
  representative densities (n_H = 10^{-6}, 10^{-4}, 10^{-2}, 10^{0} cm^{-3})

Outputs default to `docs/source/_static`, but the directory can be
overridden via `--output`.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = REPO_ROOT / "src" / "srcterms"


def extract_array(text: str, name: str) -> np.ndarray:
  pattern = rf"{name}\s*\[[^\]]*\]\s*=\s*\{{(.*?)\}};"
  match = re.search(pattern, text, re.S)
  if not match:
    raise ValueError(f"could not locate array {name}")
  data = np.fromstring(match.group(1), sep=",")
  if data.size == 0:
    raise ValueError(f"array {name} is empty")
  return data


def load_ism_cooling() -> np.ndarray:
  text = (SRC_DIR / "ismcooling.hpp").read_text()
  coeffs = extract_array(text, r"lhd")
  if coeffs.size != 102:
    raise ValueError("unexpected lhd array size")
  return coeffs


def ism_cooling(temp: np.ndarray, lhd: np.ndarray) -> np.ndarray:
  logt = np.log10(temp)
  rate = np.empty_like(temp)

  mask_low = logt <= 4.2
  if np.any(mask_low):
    t = temp[mask_low]
    rate[mask_low] = (
        2.0e-19 * np.exp(-1.184e5 / (t + 1.0e3)) + 2.8e-28 * np.sqrt(t) * np.exp(-92.0 / t)
    )

  mask_high = logt > 8.15
  if np.any(mask_high):
    lt = logt[mask_high]
    rate[mask_high] = 10.0 ** (0.45 * lt - 26.065)

  mask_mid = ~(mask_low | mask_high)
  if np.any(mask_mid):
    lt = logt[mask_mid]
    ipps = np.clip((25.0 * lt - 103).astype(int), 0, 100)
    x0 = 4.12 + 0.04 * ipps.astype(float)
    dx = lt - x0
    logcool = (lhd[ipps + 1] * dx - lhd[ipps] * (dx - 0.04)) * 25.0
    rate[mask_mid] = 10.0 ** logcool

  return rate


def load_cgm_tables():
  text = (SRC_DIR / "cooling_tables.hpp").read_text()

  h_dim = int(re.search(r"H_He_Cooling_DIM_0\s*=\s*(\d+);", text).group(1))
  n_dim = int(re.search(r"H_He_Cooling_DIM_1\s*=\s*(\d+);", text).group(1))
  metal_dim0 = int(re.search(r"Metal_Cooling_DIM_0\s*=\s*(\d+);", text).group(1))
  metal_dim1 = int(re.search(r"Metal_Cooling_DIM_1\s*=\s*(\d+);", text).group(1))
  if (h_dim, n_dim) != (metal_dim0, metal_dim1):
    raise ValueError("H+He and metal cooling tables shape mismatch")

  hhe = extract_array(text, r"H_He_Cooling_ARR").reshape(h_dim, n_dim)
  metals = extract_array(text, r"Metal_Cooling_ARR").reshape(h_dim, n_dim)
  tbins = extract_array(text, r"Tbins_ARR")
  nbins = extract_array(text, r"nHbins_ARR")

  return tbins, nbins, hhe, metals


def build_parser() -> argparse.ArgumentParser:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
      "--output",
      type=Path,
      default=REPO_ROOT / "docs" / "source" / "_static",
      help="directory for generated plots (default: %(default)s)",
  )
  return parser


def ensure_temp_axis(tbins: np.ndarray) -> np.ndarray:
  # Tbins array is stored as log10(T / K).
  if np.all((tbins > 0) & (tbins < 12)):
    return 10.0 ** tbins
  raise ValueError("unexpected T bin values")


def ensure_density_axis(nbins: np.ndarray) -> np.ndarray:
  # nH bins are log10 number density (cm^-3).
  return 10.0 ** nbins


def plot_ism_cooling(output_dir: Path, lhd: np.ndarray):
  temps = np.logspace(2, 9, 400)
  rates = ism_cooling(temps, lhd)

  fig, ax = plt.subplots(figsize=(6, 4))
  ax.loglog(temps, rates, color="#006699")
  ax.set_xlabel(r"Temperature $T$ [K]")
  ax.set_ylabel(r"$\Lambda_{\mathrm{ISM}}(T)$ [erg cm$^{3}$ s$^{-1}$]")
  ax.set_title("ISM Cooling Curve (Schure et al. 2009 + Koyama & Inutsuka 2002)")
  ax.grid(True, which="both", alpha=0.3)

  fig.tight_layout()
  path = output_dir / "ism_cooling_curve.png"
  fig.savefig(path, dpi=200)
  plt.close(fig)
  return path


def plot_cgm_cooling(output_dir: Path, temps: np.ndarray, densities: np.ndarray,
                     hhe: np.ndarray, metals: np.ndarray):
  log_t = np.log10(temps)
  fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
  selected = {
      1e-6: "n = 10$^{-6}$",
      1e-4: "n = 10$^{-4}$",
      1e-2: "n = 10$^{-2}$",
      1e0: "n = 10$^{0}$",
  }
  colors = ["#004488", "#4374B3", "#D95F02", "#A50F15"]

  for idx, (target_n, label) in enumerate(selected.items()):
    # Find nearest density bin
    density_idx = int(np.argmin(np.abs(densities - target_n)))
    axes[0].loglog(temps, hhe[:, density_idx], label=label, color=colors[idx])
    axes[1].loglog(temps, metals[:, density_idx], label=label, color=colors[idx])

  axes[0].set_xlabel(r"Temperature $T$ [K]")
  axes[0].set_ylabel(r"$\Lambda_{H+He}(T, n)$ [erg cm$^{3}$ s$^{-1}$]")
  axes[0].set_title("Primordial Cooling (PIE)")
  axes[1].set_xlabel(r"Temperature $T$ [K]")
  axes[1].set_ylabel(r"$\Lambda_{\mathrm{metals}}(T, n)$ [erg cm$^{3}$ s$^{-1}$]")
  axes[1].set_title("Metal-line Cooling (PIE)")

  for ax in axes:
    ax.grid(True, which="both", alpha=0.3)

  axes[1].legend(loc="lower left", frameon=False)
  fig.tight_layout()
  path = output_dir / "cgm_cooling_curves.png"
  fig.savefig(path, dpi=200)
  plt.close(fig)
  return path


def main() -> None:
  parser = build_parser()
  args = parser.parse_args()
  output_dir: Path = args.output
  output_dir.mkdir(parents=True, exist_ok=True)

  # Load data
  lhd = load_ism_cooling()
  tbins, nbins, hhe, metals = load_cgm_tables()
  temps = ensure_temp_axis(tbins)
  densities = ensure_density_axis(nbins)

  ism_path = plot_ism_cooling(output_dir, lhd)
  cgm_path = plot_cgm_cooling(output_dir, temps, densities, hhe, metals)

  print(f"Wrote {ism_path}")
  print(f"Wrote {cgm_path}")


if __name__ == "__main__":
  main()
