"""Shared analysis helpers for PIC regression scripts.

These helpers keep fitting and spectrum extraction logic consistent across the
Entity-mirroring test ladder.
"""

from __future__ import annotations

import numpy as np


def detrend_mean(series: np.ndarray) -> np.ndarray:
    """Return a copy with the mean removed."""
    values = np.asarray(series, dtype=float)
    return values - np.mean(values)


def dominant_frequency(time: np.ndarray, signal: np.ndarray) -> tuple[float, float]:
    """Return dominant positive frequency and amplitude from a real time series."""
    t = np.asarray(time, dtype=float)
    s = detrend_mean(np.asarray(signal, dtype=float))
    if t.size != s.size or t.size < 2:
        raise ValueError('time and signal must have matching length >= 2')
    dt = t[1] - t[0]
    if not np.allclose(np.diff(t), dt, rtol=0.0, atol=1.0e-12):
        raise ValueError('time grid must be uniform for FFT frequency extraction')

    fft = np.fft.rfft(s)
    freq = np.fft.rfftfreq(t.size, d=dt)
    if freq.size <= 1:
        return 0.0, 0.0

    amp = np.abs(fft)
    idx = int(np.argmax(amp[1:]) + 1)
    norm_amp = 2.0 * amp[idx] / t.size
    return float(freq[idx]), float(norm_amp)


def mode_amplitude(series: np.ndarray) -> float:
    """Return RMS amplitude for a real-valued mode signal."""
    s = detrend_mean(np.asarray(series, dtype=float))
    return float(np.sqrt(np.mean(s * s)))


def fit_exponential_growth(
    time: np.ndarray,
    amplitude: np.ndarray,
    tmin: float,
    tmax: float,
    floor: float = 1.0e-30,
) -> tuple[float, float, float]:
    """Fit amplitude ~ exp(gamma*t + c) over [tmin, tmax].

    Returns (gamma, intercept, r2).
    """
    t = np.asarray(time, dtype=float)
    a = np.asarray(amplitude, dtype=float)
    if t.size != a.size or t.size < 2:
        raise ValueError('time and amplitude must have matching length >= 2')

    mask = (t >= tmin) & (t <= tmax) & (a > floor)
    if np.count_nonzero(mask) < 2:
        raise ValueError('insufficient fit points in requested window')

    tx = t[mask]
    ly = np.log(a[mask])

    coeff = np.polyfit(tx, ly, 1)
    gamma = float(coeff[0])
    intercept = float(coeff[1])

    fit = np.polyval(coeff, tx)
    resid = ly - fit
    ss_res = float(np.sum(resid * resid))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = 1.0 if ss_tot == 0.0 else 1.0 - ss_res / ss_tot
    return gamma, intercept, float(r2)


def fit_exponential_growth_windowed(
    time: np.ndarray,
    amplitude: np.ndarray,
    min_points: int = 6,
    floor: float = 1.0e-30,
    min_growth_factor: float = 1.0,
) -> dict[str, float]:
    """Find the best positive-growth exponential fit over all time windows.

    The selected window maximizes ``gamma * max(r2, 0)`` subject to:
    1. at least ``min_points`` samples in the fit window,
    2. positive growth rate,
    3. end/start amplitude ratio >= ``min_growth_factor``.
    """
    t = np.asarray(time, dtype=float)
    a = np.asarray(amplitude, dtype=float)
    if t.size != a.size or t.size < min_points:
        raise ValueError('time and amplitude must have matching length >= min_points')

    best = None
    for i0 in range(0, t.size - min_points + 1):
        for i1 in range(i0 + min_points, t.size + 1):
            tx = t[i0:i1]
            ax = np.maximum(a[i0:i1], floor)
            growth_factor = float(ax[-1] / max(ax[0], floor))
            if growth_factor < min_growth_factor:
                continue

            gamma, intercept, r2 = fit_exponential_growth(
                tx, ax, float(tx[0]), float(tx[-1]), floor=floor
            )
            if gamma <= 0.0:
                continue

            score = float(gamma * max(r2, 0.0))
            if best is None or score > best['score']:
                best = {
                    'score': score,
                    'gamma': float(gamma),
                    'intercept': float(intercept),
                    'r2': float(r2),
                    'tmin': float(tx[0]),
                    'tmax': float(tx[-1]),
                    'npts': float(tx.size),
                    'growth_factor': growth_factor,
                }

    if best is None:
        gamma, intercept, r2 = fit_exponential_growth(
            t, np.maximum(a, floor), float(t[0]), float(t[-1]), floor=floor
        )
        growth_factor = float(
            np.maximum(a[-1], floor) / max(np.maximum(a[0], floor), floor)
        )
        best = {
            'score': float(gamma * max(r2, 0.0)),
            'gamma': float(gamma),
            'intercept': float(intercept),
            'r2': float(r2),
            'tmin': float(t[0]),
            'tmax': float(t[-1]),
            'npts': float(t.size),
            'growth_factor': growth_factor,
        }

    return best


def polarization_split(
    by_mode: np.ndarray,
    bz_mode: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return right/left circular polarization amplitudes from transverse modes."""
    by = np.asarray(by_mode, dtype=complex)
    bz = np.asarray(bz_mode, dtype=complex)
    if by.size != bz.size:
        raise ValueError('mode arrays must have matching shape')

    right = 0.5 * (by - 1.0j * bz)
    left = 0.5 * (by + 1.0j * bz)
    return right, left


def run(**kwargs):
    """Utility module: no standalone runtime stage."""
    return


def analyze():
    """Utility module: discovery compatibility for package-level test runs."""
    return True
