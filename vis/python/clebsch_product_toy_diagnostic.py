#!/usr/bin/env python3
"""Toy 1D Fourier-product diagnostic for the Clebsch regularity discussion."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _generate_real_series(rng: np.random.Generator, nmodes: int,
                          coeff_amp_exp: float) -> np.ndarray:
    ngrid = 2 * nmodes
    coeffs = np.zeros(ngrid, dtype=np.complex128)
    for k in range(1, nmodes):
        sigma = k ** (-coeff_amp_exp)
        zk = (rng.normal() + 1j * rng.normal()) * sigma / np.sqrt(2.0)
        coeffs[k] = zk
        coeffs[-k] = np.conj(zk)
    return np.fft.ifft(coeffs).real


def _mean_mode_energy(samples: list[np.ndarray]) -> np.ndarray:
    nmodes = samples[0].size // 2
    energy = np.zeros(nmodes + 1, dtype=np.float64)
    for field in samples:
        coeffs = np.fft.fft(field)
        for k in range(1, nmodes):
            energy[k] += 0.5 * (abs(coeffs[k]) ** 2 + abs(coeffs[-k]) ** 2)
    return energy / len(samples)


def _fit_slope(energy: np.ndarray, kmin: int, kmax: int) -> float:
    ks = np.arange(kmin, min(kmax, energy.size - 1) + 1, dtype=np.float64)
    mask = energy[ks.astype(int)] > 0.0
    xs = np.log(ks[mask])
    ys = np.log(energy[ks[mask].astype(int)])
    slope, _ = np.polyfit(xs, ys, 1)
    return float(slope)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nmodes", type=int, default=256,
                        help="Number of positive Fourier modes in each factor.")
    parser.add_argument("--samples", type=int, default=64,
                        help="Number of ensemble samples.")
    parser.add_argument("--seed", type=int, default=12345,
                        help="Random seed.")
    parser.add_argument(
        "--factor-exp",
        type=float,
        nargs="+",
        default=[0.75, 1.0],
        help="Coefficient-amplitude exponents to test. Brownian-like corresponds to 1.0.",
    )
    parser.add_argument("--output", type=Path, default=None,
                        help="Optional JSON output path.")
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    results = []
    for coeff_amp_exp in args.factor_exp:
        factor_a = []
        factor_b = []
        products = []
        for _ in range(args.samples):
            fa = _generate_real_series(rng, args.nmodes, coeff_amp_exp)
            fb = _generate_real_series(rng, args.nmodes, coeff_amp_exp)
            factor_a.append(fa)
            factor_b.append(fb)
            products.append(fa * fb)

        factor_energy = 0.5 * (_mean_mode_energy(factor_a) + _mean_mode_energy(factor_b))
        product_energy = _mean_mode_energy(products)
        kmin = max(8, args.nmodes // 32)
        kmax = max(kmin + 4, args.nmodes // 4)
        factor_slope = _fit_slope(factor_energy, kmin, kmax)
        product_slope = _fit_slope(product_energy, kmin, kmax)
        results.append(
            {
                "coeff_amp_exp": float(coeff_amp_exp),
                "fit_band": [int(kmin), int(kmax)],
                "factor_mode_energy_slope": factor_slope,
                "product_mode_energy_slope": product_slope,
                "slope_difference": product_slope - factor_slope,
            }
        )

    report = {
        "seed": args.seed,
        "nmodes": args.nmodes,
        "samples": args.samples,
        "results": results,
    }
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
