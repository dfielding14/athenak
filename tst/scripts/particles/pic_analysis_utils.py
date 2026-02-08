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
