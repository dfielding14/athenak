#!/usr/bin/env python3
"""Analyze one retained CGL-LF turbulence run; never launch a campaign.

Example: analyze_cgl_lf_physics_benchmark.py RUN --time-start 6 --time-end 10
The four figure groups are descriptive finite-window diagnostics, not a paper
image test, coefficient validation, convergence test, or active/passive control.
"""

from __future__ import annotations

import argparse
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import subprocess
import sys
import types

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import analyze_cgl_lf_paper as paper  # noqa: E402

FIELDS = paper.REQUIRED_FIELDS
FORCE_FIELDS = ("force1", "force2", "force3")
DEFINITIONS = {
    "X": "2*(p_perp-p_parallel)/B^2; AthenaK magnetic pressure is B^2/2",
    "pressure_input": "primitive binary 'eint' is p_parallel, not thermal energy",
    "PDF": "sum(cell volume in bin)/(domain volume * bin width); shared edges include every retained cell",
    "time_average": "trapezoidal integration of snapshot diagnostics over actual retained times, divided by covered duration",
    "blocks": "contiguous physical-time blocks; linearly interpolated diagnostic endpoints; sample SD/range of block means, not an iid confidence interval",
    "derivatives": "second-order centered periodic differences on the uniform full 3D mesh",
    "S_parallel": "b_i b_j partial_j u_i, b=B/|B| (local instantaneous direction)",
    "projected_gradients": "G_ij=partial_j u_i, P_ij=delta_ij-b_i*b_j; parallel/perpendicular components are b.G.b, P.G.b, b.G.P, P.G.P, projected after differentiating u; no derivatives of b enter these four components",
    "parallel_velocity": "u_parallel=u dot b; this includes the retained bulk velocity",
    "parallel_velocity_derivative": "b dot grad(u dot b) = S_parallel + u dot [(b dot grad)b] in the continuum; discretization need not obey the product rule exactly",
    "induction": "S_parallel-div(u) reconstructs ideal material D ln|B|/Dt; not a measured time derivative or resistive/numerical induction budget",
    "pressure_balance": "corr(delta p_perp, delta(B^2/2)); R=mean[(delta p_perp+delta(B^2/2))^2]/(var(p_perp)+var(B^2/2)); exact compensation gives corr=-1,R=0",
    "FFT": "F=fftn(field)/N; shell power=sum_shell |F|^2 / dk; integral over all shells (including k_perp=0) equals the stated real-space mean square",
    "shells": "k_perp=sqrt(kx^2+ky^2) relative to the initial z guide field; physical radians/length; bins [n*dk,(n+1)*dk), dk=min(2*pi/Lx,2*pi/Ly)",
    "kinetic_spectrum": "FFT of sqrt(rho/2)*(u-u_bulk), u_bulk=<rho*u>/<rho>; no further mean removal; integral=mean[rho*|u-u_bulk|^2/2]",
    "magnetic_spectrum": "FFT of (B-<B>)/sqrt(2); integral=mean[|B-<B>|^2/2]",
    "pressure_spectra": "FFT of each scalar minus its volume mean; unnormalized physical pressure units; integral=variance",
    "gradient_spectra": "FFT of each reconstructed scalar/vector minus its component means; signed fields are squared only by the power spectrum; spectra do not determine signs",
    "Mach": "u_rms about volume-mean velocity / sqrt(gamma*<p_iso>/<rho>), p_iso=(p_parallel+2*p_perp)/3; isotropic-pressure proxy, not a CGL characteristic-wave Mach number",
    "beta": "volume mean of 2*p_iso/B^2; also report 2*<p_iso>/<B^2>, which differs in general",
    "deltaB": "sqrt(<|B-<B>|^2>)/B0; mean field removed separately at each snapshot",
    "forcing": "actual accumulated applied forcing work is user-history force_work; force_pwr is instantaneous rho*u dot f and is not its exact quadrature",
    "force_decomposition": "nonzero Fourier acceleration modes: P_compressive=sum |k dot fhat|^2/k^2; P_solenoidal=P_total-P_compressive; this is acceleration power, not injected-energy partition",
    "nominal_force_mixture": "expected_solenoidal_power_fraction and expected_solenoidal_fraction retain the nominal ratio of expected isotropic innovation powers, 2*s^2/[2*s^2+(1-s)^2]; this is not the expectation of an instantaneous fraction or a target for one finite OU realization",
    "precision": "primitive/force snapshots may be float32; strict threshold crossings are descriptive and not tight full-precision admissibility tests",
    "X_rounding_envelope": "half the larger adjacent float32 gap for each stored pressure/B component; bound |delta X| <= [2*(delta p_perp+delta p_parallel)+|X|*delta(B^2)]/[B^2-delta(B^2)], delta(B^2)=sum(2*|Bi|*delta Bi+delta Bi^2); no dynamical-error claim",
    "limiter_history": "mirror_vol/fire_vol are inclusive threshold predicates; with nu_coll=0, both soft limiters enabled and backups disabled, nu_eff/(limiter_nu_coll*volume) measures the strict soft-rate fraction at history sampling; hard_vol counts strict physical firehose violation even with backups off",
    "LF_counters": "lf_nstage counts cumulative active-cell LF-stage checks; lf_hardbd counts cumulative hard-bound crossings in those unprojected checks, with repeated cells counted repeatedly. nonfin/nonpos also inspect active cells; dfloor/pfloor accumulate EOS refresh events including refreshed halo cells. All persist restart and do not reset at history output. In audited post-WO2 simulation revisions 7a37710f6 and 71ad25ebc, lf_hwproj is a reserved, uninstrumented field: its retained value does not measure wall projections or their absence. Physical stage crossings are not automatically numerical failure",
}


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_input(lines):
    result, block = {}, None
    for original in lines:
        line = original.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1]
            result.setdefault(block, {})
        elif block and "=" in line:
            key, value = line.split("=", 1)
            result[block][key.strip()] = value.strip()
    return result


def retained_file(path):
    return {"path": str(path.resolve()), "sha256": sha(path),
            "bytes": path.stat().st_size}


def read_uniform(path, quantities):
    """Reuse the existing exact-rank reader, with generic field selection."""
    rank_files = paper.snapshot_sibling_paths(path)
    local = path.parent.name == "rank_00000000"
    raw = (paper.read_exact_rank_set_binary(rank_files) if local
           else paper.bin_convert.read_binary(str(path)))
    paper.validate_exact_rank_set_meshblocks(raw)
    if np.any(np.asarray(raw["mb_logical"])[:, 3] != 0):
        raise ValueError("FFT diagnostics require an unrefined uniform full mesh")
    if any(raw[f"nx{i}_out_mb"] != raw[f"nx{i}_mb"] for i in (1, 2, 3)):
        raise ValueError("sliced/ghost-inclusive snapshots are not supported")
    missing = set(quantities) - set(raw["var_names"])
    if missing:
        raise ValueError(f"{path}: missing fields {sorted(missing)}")
    name = "read_all_ranks_binary_as_athdf" if local else "read_binary_as_athdf"
    original = getattr(paper.bin_convert, name)
    namespace = original.__globals__.copy()
    namespace["read_all_ranks_binary" if local else "read_binary"] = lambda _: raw
    converter = types.FunctionType(original.__code__, namespace,
                                   argdefs=original.__defaults__)
    converter.__kwdefaults__ = original.__kwdefaults__
    values = converter(str(path), quantities=list(quantities), dtype=np.float64)
    lengths = tuple(float(raw[f"x{i}max"] - raw[f"x{i}min"])
                    for i in (1, 2, 3))
    for i in (1, 2, 3):
        widths = np.diff(values[f"x{i}f"])
        if len(widths) < 4 or not np.allclose(widths, widths[0], rtol=1e-10, atol=0):
            raise ValueError("full 3D uniform Cartesian cells are required")
    fields = {name: np.asarray(values[name], dtype=float) for name in quantities}
    if any(not np.isfinite(value).all() for value in fields.values()):
        raise ValueError(f"nonfinite snapshot field: {path}")
    header = parse_input(raw["header"])
    if any(header.get("mesh", {}).get(f"{side}x{i}_bc") != "periodic"
           for i in (1, 2, 3) for side in ("i", "o")):
        raise ValueError("periodic boundaries are required for these FFT/derivative diagnostics")
    info = {"time": float(raw["time"]), "cycle": int(raw["cycle"]),
            "shape_zyx": list(next(iter(fields.values())).shape),
            "lengths_xyz": lengths, "cell_volume": float(np.prod(lengths)
                / next(iter(fields.values())).size),
            "payload_dtype": str(np.asarray(raw["mb_data"][quantities[0]]).dtype),
            "files": [retained_file(p) for p in rank_files]}
    return fields, lengths, header, info


def time_weights(times):
    times = np.asarray(times, dtype=float)
    if len(times) == 1:
        return np.ones(1)
    dt = np.diff(times)
    if np.any(dt <= 0):
        raise ValueError("times must be unique and increasing")
    weights = np.r_[dt[0], dt[:-1] + dt[1:], dt[-1]] * 0.5
    return weights / (times[-1] - times[0])


def interval_mean(times, values, start, end):
    times, values = np.asarray(times), np.asarray(values)
    selected = (times > start) & (times < end)
    grid = np.r_[start, times[selected], end]
    flat = values.reshape(len(times), -1)
    samples = np.stack([np.interp(grid, times, flat[:, j])
                        for j in range(flat.shape[1])], axis=1)
    return np.tensordot(time_weights(grid), samples, axes=1).reshape(values.shape[1:])


def temporal_summary(times, values, block_duration, start=None, end=None):
    times, values = np.asarray(times), np.asarray(values)
    lo = max(times[0], times[0] if start is None else start)
    hi = min(times[-1], times[-1] if end is None else end)
    if hi < lo:
        raise ValueError("no temporal overlap")
    if len(times) > 1:
        grid = np.unique(np.r_[lo, times[(times > lo) & (times < hi)], hi])
        flat = values.reshape(len(times), -1)
        values = np.stack([np.interp(grid, times, flat[:, j])
            for j in range(flat.shape[1])], axis=1).reshape((len(grid),)+values.shape[1:])
        times = grid
    mean = np.tensordot(time_weights(times), values, axes=1)
    blocks, starts = [], []
    for i in range(int(np.floor((times[-1] - times[0]) / block_duration + 1e-10))):
        start = times[0] + i * block_duration
        blocks.append(interval_mean(times, values, start, start + block_duration))
        starts.append(float(start))
    result = {"mean": np.asarray(mean).tolist(), "block_starts": starts,
              "block_duration": block_duration, "block_means": np.asarray(blocks).tolist(),
              "block_count": len(blocks), "block_sd": None, "block_min": None,
              "block_max": None, "effective_window": [float(lo), float(hi)]}
    if len(blocks) >= 2:
        result.update(block_sd=np.std(blocks, axis=0, ddof=1).tolist(),
                      block_min=np.min(blocks, axis=0).tolist(),
                      block_max=np.max(blocks, axis=0).tolist())
    if values.ndim == 1 and len(times) >= 2:
        slope = float(np.polyfit(times - times[0], values, 1)[0])
        result["linear_slope"] = slope
        result["fitted_change_over_window"] = slope * float(times[-1] - times[0])
    return result


def unavailable_summary(times, values, start, end, reason):
    """Explain undefined diagnostics without changing their integration window."""
    missing = [float(t) for t, value in zip(times, values) if value is None]
    return {"available": False, "reason": reason + "; required sampling or endpoint-bracketing values are undefined. "
            "No interpolation across undefined values or narrowing of the averaging window is performed.",
            "requested_window": [start, end], "undefined_sample_times": missing,
            "undefined_bracketing_times": [t for t in missing if t < start or t > end],
            "valid_samples": [{"time": float(t), "value": float(value)}
                              for t, value in zip(times, values) if value is not None]}


def autocorrelation(times, values):
    times, values = np.asarray(times), np.asarray(values)
    if len(times) < 8:
        return {"available": False, "reason": "need at least eight samples"}
    dt = float(np.median(np.diff(times)))
    if max(np.diff(times)) > 1.5*dt:
        return {"available": False, "reason": "sampling gap exceeds 1.5 median cadences"}
    resampled = not np.allclose(np.diff(times), dt, rtol=1e-6)
    if resampled:
        grid = np.arange(times[0], times[-1]+dt*1e-6, dt)
        values = np.interp(grid, times, values)
        times = grid
    centered = values - np.mean(values)
    variance = np.dot(centered, centered)
    if variance <= 0:
        return {"available": False, "reason": "zero sample variance"}
    acf = np.correlate(centered, centered, mode="full")[len(values) - 1:] / variance
    stop = next((i for i in range(1, len(acf)) if acf[i] <= 0), len(acf) // 2)
    tau = float(np.diff(times)[0] * (0.5 + np.sum(acf[1:stop])))
    return {"available": True, "resampled_to_median_cadence": resampled,
            "uniform_grid_start": float(times[0]), "uniform_grid_end": float(times[-1]),
            "uniform_grid_spacing": float(np.diff(times)[0]), "samples": len(times),
            "lag": (np.arange(stop + 1) * np.diff(times)[0]).tolist(),
            "acf": acf[:stop + 1].tolist(), "positive_sequence_integral_time": tau,
            "nominal_effective_samples": float((times[-1] - times[0]) / (2 * tau)),
            "caution": "mean-subtracted, not detrended; finite-window/drift-biased diagnostic, not proof of independent blocks"}


def volume_pdf(values, edges, cell_volumes):
    volumes = np.broadcast_to(cell_volumes, values.shape)
    count = np.histogram(values.ravel(), bins=edges, weights=volumes.ravel())[0]
    return count / (np.sum(volumes) * np.diff(edges))


def grids(shape, lengths):
    return np.meshgrid(*[2 * np.pi * np.fft.fftfreq(n, d=l / n)
        for n, l in zip(shape, lengths[::-1])], indexing="ij")[::-1]


def spectrum(fields, lengths, remove_mean=True):
    shape = fields[0].shape
    kx, ky, kz = grids(shape, lengths)
    dk = min(2 * np.pi / lengths[0], 2 * np.pi / lengths[1])
    shell = np.floor(np.hypot(kx, ky) / dk + 1e-12).astype(int)
    power, real_power = np.zeros(shape), 0.0
    for field in fields:
        centered = field - np.mean(field) if remove_mean else field
        real_power += float(np.mean(centered ** 2))
        power += abs(np.fft.fftn(centered) / centered.size) ** 2
    binned = np.bincount(shell.ravel(), weights=power.ravel())
    k_nyquist = min(np.pi*n/l for n, l in zip(shape, lengths[::-1]))
    k2 = kx*kx+ky*ky+kz*kz
    resolved = (k2 <= (.25*k_nyquist)**2) & (k2 > 0)
    resolved_power = np.bincount(shell.ravel(), weights=(power*resolved).ravel(),
                                minlength=len(binned))
    return {"dk": dk, "k": ((np.arange(len(binned)) + 0.5) * dk).tolist(),
            "power": (binned / dk).tolist(), "real_space_power": real_power,
            "full_k_resolved_power": (resolved_power/dk).tolist(),
            "spectral_integral": float(np.sum(binned)),
            "parseval_relative_error": float(abs(np.sum(binned) - real_power)
                / max(real_power, np.finfo(float).tiny))}


def pressure_balance(perp, magnetic):
    a, b = perp - np.mean(perp), magnetic - np.mean(magnetic)
    va, vb, cov = float(np.mean(a*a)), float(np.mean(b*b)), float(np.mean(a*b))
    return {"correlation": cov / np.sqrt(va*vb) if va*vb > 0 else None,
            "normalized_residual_variance": float(np.mean((a+b)**2)) / (va+vb)
                if va+vb > 0 else None,
            "covariance": cov, "variance_perp": va, "variance_magnetic": vb}


def float32_x_envelope(ppar, perp, magnetic):
    """Conservative rounding interval for X from independently stored floats.

    Half the larger adjacent float32 gap encloses rounding of each component.
    This is a representation sensitivity bound, not a dynamical error bound.
    """
    def half_ulp(value):
        stored = np.asarray(value, dtype=np.float32)
        up = np.nextafter(stored, np.float32(np.inf)).astype(float)
        down = np.nextafter(stored, np.float32(-np.inf)).astype(float)
        return .5*np.maximum(up-stored, stored-down)
    b2 = sum(v*v for v in magnetic)
    numerator = 2*(perp-ppar)
    dn = 2*(half_ulp(perp)+half_ulp(ppar))
    db2 = sum(2*abs(v)*half_ulp(v)+half_ulp(v)**2 for v in magnetic)
    return np.divide(dn+abs(numerator/b2)*db2, b2-db2,
                     out=np.full_like(b2, np.inf), where=b2 > db2)


def snapshot_products(fields, lengths, model, near_width):
    rho, ppar, perp = (fields[key] for key in ("dens", "eint", "p_perp"))
    if min(np.min(rho), np.min(ppar), np.min(perp)) <= 0:
        raise ValueError("positive density and both pressures required")
    u = [fields[key] for key in ("velx", "vely", "velz")]
    b = [fields[key] for key in ("bcc1", "bcc2", "bcc3")]
    b2 = sum(x*x for x in b)
    if np.min(b2) <= 0:
        raise ValueError("local-field diagnostics undefined where B=0")
    bhat = [x/np.sqrt(b2) for x in b]
    up = sum(x*y for x, y in zip(u, bhat))
    uperp = [x-y*up for x, y in zip(u, bhat)]
    grad = [paper.periodic_gradient(x, lengths) for x in u]
    strain = sum(bhat[i]*bhat[j]*grad[i][j] for i in range(3) for j in range(3))
    div = sum(grad[i][i] for i in range(3))
    gradup = paper.periodic_gradient(up, lengths)
    along_up = sum(bhat[i]*gradup[i] for i in range(3))
    # Project the velocity-gradient tensor itself, as in the four-component
    # decomposition. Differentiating u_perp instead would introduce grad(b).
    gb = [sum(grad[i][j]*bhat[j] for j in range(3)) for i in range(3)]
    bg = [sum(bhat[i]*grad[i][j] for i in range(3)) for j in range(3)]
    parallel_perp = [gb[i]-bhat[i]*strain for i in range(3)]
    perp_parallel = [bg[j]-bhat[j]*strain for j in range(3)]
    transverse = [grad[i][j]-bhat[i]*bg[j]-gb[i]*bhat[j]+bhat[i]*strain*bhat[j]
                  for i in range(3) for j in range(3)]
    piso, mag = (ppar + 2*perp)/3, b2/2
    x = 2*(perp-ppar)/b2
    bulk = [np.mean(rho*v)/np.mean(rho) for v in u]
    kinetic_fields = [np.sqrt(rho/2)*(v-v0) for v, v0 in zip(u, bulk)]
    magnetic_fields = [(v-np.mean(v))/np.sqrt(2) for v in b]
    spectral_fields = {
        "p_parallel": [ppar], "p_perp": [perp], "magnetic_pressure": [mag],
        "S_parallel": [strain], "div_u": [div], "induction": [strain-div],
        "b_grad_u_parallel": [along_up], "perp_gradient_u_perp": transverse,
        "parallel_gradient_u_perp": parallel_perp,
        "perp_gradient_u_parallel": perp_parallel,
        "u_parallel": [up], "u_perp": uperp,
    }
    spectra = {name: spectrum(values, lengths) for name, values in spectral_fields.items()}
    spectra["kinetic"] = spectrum(kinetic_fields, lengths, remove_mean=False)
    spectra["magnetic"] = spectrum(magnetic_fields, lengths, remove_mean=False)
    scalars = {"B_over_B0_mean": float(np.mean(np.sqrt(b2)/model["B0"])),
        "deltaB_rms_over_B0": float(np.sqrt(2*sum(np.mean(v*v) for v in magnetic_fields))/model["B0"]),
        "beta_volume_mean": float(np.mean(2*piso/b2)),
        "beta_ratio_of_means": float(2*np.mean(piso)/np.mean(b2)),
        "Mach_isotropic_proxy": float(np.sqrt(sum(np.mean((v-np.mean(v))**2) for v in u)
            / (model["gamma_sound"]*np.mean(piso)/np.mean(rho)))),
        "u_parallel_rms": float(np.sqrt(np.mean(up*up))),
        "u_parallel_fraction": float(np.mean(up*up)/max(sum(np.mean(v*v) for v in u), np.finfo(float).tiny)),
        "u_perp_rms": float(np.sqrt(sum(np.mean(v*v) for v in uperp))),
        "S_parallel_rms": float(np.sqrt(np.mean(strain*strain))),
        "div_u_rms": float(np.sqrt(np.mean(div*div))),
        "induction_rms": float(np.sqrt(np.mean((strain-div)**2))),
        "curvature_plus_discretization_rms": float(np.sqrt(np.mean((along_up-strain)**2))),
        "perp_gradient_u_perp_rms": float(np.sqrt(sum(np.mean(v*v) for v in transverse))),
        "parallel_gradient_u_perp_rms": float(np.sqrt(sum(np.mean(v*v) for v in parallel_perp))),
        "perp_gradient_u_parallel_rms": float(np.sqrt(sum(np.mean(v*v) for v in perp_parallel))),
        "kinetic_density": float(np.mean(0.5*rho*sum(v*v for v in u))),
        "magnetic_density": float(np.mean(mag)),
        "thermal_density": float(np.mean(perp+0.5*ppar)),
    }
    for name, value in pressure_balance(perp, mag).items():
        if value is not None:
            scalars[f"pressure_{name}"] = float(value)
    rounding = float32_x_envelope(ppar, perp, b)
    scalars["X_float32_rounding_envelope_max"] = float(np.max(rounding))
    for label, threshold in (("mirror", model["mirror_X"]), ("firehose", model["firehose_X"])):
        strict = x > threshold if label == "mirror" else x < threshold
        inclusive = x >= threshold if label == "mirror" else x <= threshold
        scalars[label+"_strict"] = float(np.mean(strict))
        scalars[label+"_inclusive"] = float(np.mean(inclusive))
        scalars[label+"_near"] = float(np.mean(abs(x-threshold) <= near_width))
        interior = (x <= threshold) if label == "mirror" else (x >= threshold)
        scalars[label+"_near_interior"] = float(np.mean((abs(x-threshold) <= near_width) & interior))
        beyond = x-rounding > threshold if label == "mirror" else x+rounding < threshold
        scalars[label+"_strict_beyond_rounding_envelope"] = float(np.mean(beyond))
        scalars[label+"_rounding_ambiguous"] = float(np.mean(abs(x-threshold) <= rounding))
    return scalars, spectra, {"B_over_B0": np.sqrt(b2)/model["B0"], "X": x}


def force_products(fields, lengths):
    modes = [np.fft.fftn(fields[key])/fields[key].size for key in FORCE_FIELDS]
    k = grids(modes[0].shape, lengths)
    k2 = sum(v*v for v in k)
    longitudinal = sum(v*f for v, f in zip(k, modes))
    nonzero = k2 > 0
    compressive = float(np.sum(abs(longitudinal[nonzero])**2/k2[nonzero]))
    total = float(sum(np.sum(abs(f[nonzero])**2) for f in modes))
    return {"total_acceleration_power": total, "compressive_acceleration_power": compressive,
            "solenoidal_acceleration_power": max(0., total-compressive),
            "solenoidal_fraction": (total-compressive)/total if total > 0 else None,
            "zero_mode_power": float(sum(abs(f[0, 0, 0])**2 for f in modes))}


def paths_from_metadata(run, outputs, name, default_glob=None):
    paths = outputs.get(name)
    if paths is None:
        pattern = outputs.get(name.removesuffix("s")+"_glob", default_glob)
        if name == "force_snapshots":
            pattern = outputs.get("forcing_glob", pattern)
        paths = sorted(run.glob(pattern)) if pattern else []
    elif isinstance(paths, str):
        paths = [paths]
    result = []
    for value in paths:
        path = Path(value)
        path = path if path.is_absolute() else run/path
        if path.parent.name.startswith("rank_") and path.parent.name != "rank_00000000":
            continue
        if not path.is_file():
            raise ValueError(f"retained file not found: {path}")
        result.append(path.resolve())
    return list(dict.fromkeys(result))


def merge_histories(paths, restart_boundaries=None):
    rows, provenance, duplicates, branches = {}, [], [], []
    boundaries = {}
    for event in restart_boundaries or []:
        boundaries.setdefault(event["before_file_index"], []).append(event)
    # A restart boundary belongs to the lineage, even when that continuation
    # produced no history file. Process the final boundary after the last file.
    for index in range(len(paths) + 1):
        for event in boundaries.get(index, []):
            boundary = event["restart_time"]
            boundary_key = float(f"{boundary:.12g}")
            discarded = sorted(t for t in rows if t > boundary_key)
            branches.append({"restart_time": boundary, "discarded_times": discarded,
                             "segment": event["segment"],
                             "basis": "retained restart header"})
            for t in discarded:
                del rows[t]
        if index == len(paths):
            break
        path = paths[index]
        data = paper.parse_history(path)
        provenance.append(retained_file(path))
        for i, time in enumerate(data["time"]):
            row = {name: float(value[i]) for name, value in data.items()}
            # Canonical printed-time key merges harmless representation noise.
            key = float(f"{time:.12g}")
            if rows and key < max(rows):
                discarded = [t for t in rows if t >= key]
                branches.append({"restart_time": key, "discarded_times": sorted(discarded),
                                 "path": str(path)})
                for t in discarded:
                    del rows[t]
            if key in rows:
                changes = [name for name in row if name in rows[key] and row[name] != rows[key][name]]
                duplicates.append({"time": float(time), "kept_path": str(path),
                                   "changed_columns": changes})
            rows[key] = row
    if not rows:
        return {}, {"available": False, "reason": "no retained history rows",
                    "files": provenance, "duplicate_rows": duplicates,
                    "discarded_restart_branches": branches}
    ordered = [rows[t] for t in sorted(rows)]
    keys = set.intersection(*(set(row) for row in ordered))
    data = {key: np.asarray([row[key] for row in ordered]) for key in keys}
    if any(not np.isfinite(v).all() for v in data.values()):
        raise ValueError("nonfinite history value")
    return data, {"available": True, "files": provenance, "duplicate_rows": duplicates,
                  "discarded_restart_branches": branches,
                  "dedup_rule": "chronological input order; a backward time jump removes the stale future branch, last duplicate row kept (12 significant-digit time key); changes audited"}


def load_segments(run, metadata, metadata_path):
    """A union manifest lists immutable segment directories in lineage order."""
    directories = metadata.get("segments")
    if directories is None:
        entries = [(run, metadata, metadata_path)]
    else:
        if not isinstance(directories, list) or not directories:
            raise ValueError("segments must be a nonempty chronological directory list")
        entries = []
        for value in directories:
            directory = (run/str(value)).resolve()
            path = directory/"benchmark_metadata.json"
            item = json.loads(path.read_text())
            if "segments" in item:
                raise ValueError("nested union manifests are not supported")
            entries.append((directory, item, path))
    result = []
    for directory, item, path in entries:
        simulation = item.get("simulation", {})
        record = {"directory": directory, "metadata": item, "metadata_file": retained_file(path),
                  "input": None, "input_record": None, "restart_time": None}
        input_name = simulation.get("input_path")
        if input_name:
            input_path = (directory/input_name).resolve()
            record["input_record"] = retained_file(input_path)
            expected = simulation.get("input_sha256")
            if expected and expected != record["input_record"]["sha256"]:
                raise ValueError(f"retained effective-input hash changed: {input_path}")
            record["input"] = parse_input(input_path.read_text().splitlines())
        restart_name = item.get("launch", {}).get("restart")
        if restart_name:
            restart = (directory/restart_name).resolve()
            text, _ = paper.restart_parameter_dump(restart)
            record["restart_time"] = paper.restart_marker_time(text, restart)
            record["restart_file"] = retained_file(restart)
            expected_restart = item.get("launch", {}).get("restart_sha256")
            if expected_restart and expected_restart != record["restart_file"]["sha256"]:
                raise ValueError(f"retained restart hash changed: {restart}")
            if result and not any(restart.is_relative_to(row["directory"]) for row in result):
                raise ValueError("segment restart is not retained beneath an earlier segment")
        elif result:
            raise ValueError("every continuation segment must identify its retained restart")
        result.append(record)
    return result


def segment_paths(segments, kind, pattern, snapshots=False):
    paths, boundaries, discarded = [], [], []
    for segment in segments:
        current = paths_from_metadata(segment["directory"], segment["metadata"].get("outputs", {}), kind, pattern)
        boundary = segment["restart_time"]
        if boundary is not None:
            if snapshots:
                stale = [path for path in paths if paper.snapshot_time(path) > boundary+1e-12]
                discarded += [{"path": str(path), "restart_time": boundary} for path in stale]
                paths = [path for path in paths if path not in stale]
            else:
                boundaries.append({"before_file_index": len(paths),
                                   "restart_time": boundary,
                                   "segment": str(segment["directory"])})
        paths.extend(current)
    return paths, boundaries, discarded


def bracketing_paths(paths, start, end):
    if not paths:
        return []
    dated = [(paper.snapshot_time(path), path) for path in paths]
    times = sorted(set(time for time, _ in dated))
    if times[-1] < start or times[0] > end:
        return []
    lo = times[max(0, np.searchsorted(times, start, side="right")-1)]
    hi = times[min(len(times)-1, np.searchsorted(times, end, side="left"))]
    return [path for time, path in dated if lo <= time <= hi]


def crosscheck_input(header, retained):
    for section in ("mhd", "problem", "turb_driving", "mesh"):
        for key, value in retained.get(section, {}).items():
            actual = header.get(section, {}).get(key)
            if actual is None:
                raise ValueError(f"retained input key absent from snapshot: {section}/{key}")
            try:
                equal = np.isclose(float(value), float(actual), rtol=1e-12, atol=0)
            except ValueError:
                equal = value == actual
            if not equal:
                raise ValueError(f"retained input differs from snapshot: {section}/{key}: {value} vs {actual}")


def infer_model(header, metadata):
    problem, mhd, driving = (header.get(k, {}) for k in ("problem", "mhd", "turb_driving"))
    required = {"problem": ("b0",),
        "mhd": ("mirror_threshold", "firehose_threshold", "gamma", "nu_coll", "limiter_nu_coll",
                "backup_limiters", "mirror_limiter", "firehose_limiter", "passive"),
        "turb_driving": ("tcorr", "sol_fraction", "projection_policy", "dedt", "nlow", "nhigh",
                         "k_shell_unit", "record_injected_work", "normalization",
                         "driving_type", "physical_k_shell")}
    missing = [f"{section}/{key}" for section, keys in required.items()
               for key in keys if key not in header.get(section, {})]
    if missing:
        raise ValueError("retained effective input lacks required model fields: "+", ".join(missing))

    def boolean(section, key):
        value = header[section][key].lower()
        if value not in ("true", "false"):
            raise ValueError(f"invalid retained boolean {section}/{key}: {value}")
        return value == "true"
    if (driving["normalization"] != "edot" or int(driving["driving_type"]) != 0
            or not boolean("turb_driving", "physical_k_shell")):
        raise ValueError("benchmark analyzer requires turb_driving normalization=edot, "
                         "driving_type=0, and physical_k_shell=true")
    model = {"B0": abs(float(problem["b0"])), "mirror_X": float(mhd["mirror_threshold"]),
             "firehose_X": -float(mhd["firehose_threshold"]),
             "gamma_sound": float(mhd["gamma"]),
             "forcing_correlation_time": float(driving["tcorr"]),
             "sol_fraction": float(driving["sol_fraction"]),
             "projection_policy": driving["projection_policy"],
             "dedt_per_volume": float(driving["dedt"]),
             "forcing_kmin": float(driving["nlow"])*float(driving["k_shell_unit"]),
             "forcing_kmax": float(driving["nhigh"])*float(driving["k_shell_unit"]),
             "record_injected_work": boolean("turb_driving", "record_injected_work"),
             "nu_coll": float(mhd["nu_coll"]),
             "limiter_nu_coll": float(mhd["limiter_nu_coll"]),
             "backup_limiters": boolean("mhd", "backup_limiters"),
             "mirror_limiter": boolean("mhd", "mirror_limiter"),
             "firehose_limiter": boolean("mhd", "firehose_limiter"),
             "passive": boolean("mhd", "passive")}
    if model["B0"] <= 0 or model["gamma_sound"] <= 0:
        raise ValueError("positive retained B0 and gamma required")
    for key, supplied in metadata.get("model", {}).items():
        if key in model and not (np.isclose(float(supplied), model[key], rtol=1e-12)
                                if isinstance(model[key], float) else supplied == model[key]):
            raise ValueError(f"retained input disagrees with metadata model.{key}")
    fraction = model["sol_fraction"]
    model["expected_solenoidal_power_fraction"] = (2*fraction*fraction
        / (2*fraction*fraction+(1-fraction)**2)) if model["projection_policy"] == "solenoidal_compressive" else None
    return model


def history_products(user, mhd, model, start, end, block_duration):
    result = {"available": bool(user), "series": {}, "window": {}, "energy_budget": {
        "available": False, "requested_interval": [start, end],
        "requested_window_covered": False}}
    if not user:
        result["reason"] = "user history missing; snapshots are not exact energy/source ledgers"
        return result
    t = user["time"]
    series = {name: value for name, value in user.items()}
    if "volume" in user:
        volume = user["volume"]
        for label in ("mirror", "fire", "hard"):
            if label+"_vol" in user:
                kind = "strict" if label == "hard" else "inclusive"
                series[label+"_"+kind+"_history_fraction"] = user[label+"_vol"]/volume
        if ("nu_eff" in user and model["nu_coll"] == 0 and not model["backup_limiters"]
                and model["limiter_nu_coll"] > 0 and model["mirror_limiter"] and model["firehose_limiter"]):
            series["strict_soft_rate_history_fraction"] = user["nu_eff"]/(model["limiter_nu_coll"]*volume)
        if "beta" in user:
            series["beta_volume_mean"] = user["beta"]/volume
    selected = (t >= start-1e-12) & (t <= end+1e-12)
    if np.count_nonzero(selected) >= 2:
        for name in ("kinetic", "magnetic", "therm_cgl", "beta_volume_mean", "force_work"):
            if name in series:
                result["window"][name] = temporal_summary(t, series[name], block_duration, start, end)
                result["window"][name]["autocorrelation"] = autocorrelation(t[selected], series[name][selected])
    if "tot-E" in mhd:
        result["conserved_energy_history"] = {"time": mhd["time"].tolist(), "tot-E": mhd["tot-E"].tolist()}
    if ("force_work" in user and "tot-E" in mhd and model["record_injected_work"]
            and model.get("passive") is False):
        lo, hi = max(start, t[0], mhd["time"][0]), min(end, t[-1], mhd["time"][-1])
        if hi > lo:
            work = np.interp([lo, hi], t, user["force_work"])
            energy = np.interp([lo, hi], mhd["time"], mhd["tot-E"])
            dw, de = float(work[1]-work[0]), float(energy[1]-energy[0])
            result["energy_budget"] = {"available": True, "interval": [lo, hi],
                "requested_interval": [start, end],
                "requested_window_covered": bool(lo <= start and hi >= end),
                "actual_applied_work": dw, "actual_mean_total_power": dw/(hi-lo),
                "conserved_total_energy_change": de, "residual_E_minus_work": de-dw,
                "relative_residual": abs(de-dw)/max(abs(de), abs(dw), np.finfo(float).tiny),
                "sampling": "linear interpolation only if explicit endpoints are between retained history rows"}
    if not result["energy_budget"]["available"]:
        result["energy_budget"]["reason"] = "need active CGL, enabled measured force_work, conserved tot-E, and overlapping interval"
    result["series"] = {name: value.tolist() for name, value in series.items()}
    if mhd:
        selection = (mhd["time"] >= start-1e-12) & (mhd["time"] <= end+1e-12)
        indices = np.flatnonzero(selection)
        if len(indices) >= 2:
            lo, hi = max(start, mhd["time"][0]), min(end, mhd["time"][-1])
            result["LF_counter_increments"] = {key: float(np.diff(np.interp([lo, hi], mhd["time"], values))[0])
                for key, values in mhd.items() if key.startswith("lf_")}
    return result


def common_energy_history(history):
    """Compare cumulative ledgers on their overlap with one physical baseline."""
    user = history.get("series", {})
    energy = history.get("conserved_energy_history", {})
    if "force_work" not in user or "tot-E" not in energy:
        return None
    lo = max(user["time"][0], energy["time"][0])
    hi = min(user["time"][-1], energy["time"][-1])
    if hi <= lo:
        return None
    result = {"interval": [lo, hi]}
    for name, record, column in (("work", user, "force_work"),
                                  ("energy", energy, "tot-E")):
        times = np.asarray(record["time"])
        grid = np.r_[lo, times[(times > lo) & (times < hi)], hi]
        values = np.interp(grid, times, record[column])
        result[name] = {"time": grid, "increment": values-values[0]}
    return result


def solver_health(segments, user, mhd, start, end):
    result = {"classification": "inconclusive", "segments": [], "reasons": [],
              "numerical_counters": {}, "physical_stage_activity": {},
              "post_operator_hard_volume": {"available": False}}
    for segment in segments:
        launch = segment["metadata"].get("launch", {})
        code = launch.get("returncode")
        result["segments"].append({"directory": str(segment["directory"]),
            "returncode": code, "completion_status": "unavailable" if code is None else "completed" if code == 0 else "failed"})
        if code is None:
            result["reasons"].append("segment completion returncode absent from retained metadata")
        elif code != 0:
            result["reasons"].append("nonzero simulation returncode")
    for group, names in (("numerical_counters", ("lf_dfloor", "lf_pfloor", "lf_nonfin", "lf_nonpos")),
                         ("physical_stage_activity", ("lf_hardbd", "lf_hwproj"))):
        for name in names:
            if name not in mhd:
                result[group][name] = {"available": False}
                continue
            values, times = mhd[name], mhd["time"]
            lo, hi = max(start, times[0]), min(end, times[-1])
            result[group][name] = {"available": True, "full_retained_min": float(min(values)),
                "full_retained_max": float(max(values)), "last": float(values[-1]),
                "window_increment": float(np.diff(np.interp([lo, hi], times, values))[0]) if hi > lo else None}
            if group == "numerical_counters" and max(values) > 0:
                result["reasons"].append(name+" records numerical floor/nonfinite/nonpositive activity")
    result["physical_stage_activity"]["lf_hardbd"]["meaning"] = (
        "cumulative active-cell hard-bound crossings at unprojected LF stages; repeated cell-stage events, not unique cells or a time-integrated volume fraction")
    revisions = [segment["metadata"].get("simulation", {}).get("revision") for segment in segments]
    audited_revisions = ("7a37710f6c224e24e7c7f364e7e0b812b3a9494c",
                         "71ad25ebce73d33db048defd8f585a7dce0528c2")
    audited = bool(revisions) and all(revision in audited_revisions for revision in revisions)
    result["physical_stage_activity"]["lf_hwproj"].update({
        "instrumentation": "reserved_uninstrumented" if audited else "not_verified_for_retained_revision",
        "audited_simulation_revision": "7a37710f6c224e24e7c7f364e7e0b812b3a9494c",
        "audited_simulation_revisions": list(audited_revisions),
        "meaning": "post-WO2 audited source initializes and serializes this field but never increments it; retained values are preserved for provenance and do not establish wall-projection activity or its absence"})
    result["LF_counter_lifecycle"] = (
        "audited LF stage/admissibility counters are cumulative and restart-persistent, not reset per history or process; shared-file restart restores the prior global sum on rank0 only, so later MPI history sums preserve it once")
    if mhd:
        result["achieved_history_range"] = [float(mhd["time"][0]), float(mhd["time"][-1])]
        result["requested_window_covered"] = bool(mhd["time"][0] <= start and mhd["time"][-1] >= end)
    else:
        result["achieved_history_range"] = None
        result["requested_window_covered"] = False
    if "hard_vol" in user and "volume" in user:
        fraction = user["hard_vol"]/user["volume"]
        result["post_operator_hard_volume"] = {"available": True,
            "min_fraction": float(min(fraction)), "max_fraction": float(max(fraction)),
            "meaning": "strict hard-bound predicate at retained history times; distinct from unprojected LF stage crossings"}
    complete = all(row["available"] for row in result["numerical_counters"].values())
    if any("nonzero" in reason or "records numerical" in reason for reason in result["reasons"]):
        result["classification"] = "concerning"
    elif complete and not result["reasons"] and result["requested_window_covered"]:
        result["classification"] = "consistent"
    result["scope"] = "completion, coverage and numerical counter evidence; no threshold on physical hard-bound/projection activity or energy residual"
    return result


def json_value(value):
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_value(v) for v in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--time-start", required=True, type=float)
    parser.add_argument("--time-end", required=True, type=float)
    parser.add_argument("--block-duration", default=2., type=float)
    parser.add_argument("--near-width", default=.05, type=float)
    parser.add_argument("--pdf-bins", default=128, type=int)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args(argv)
    if args.time_end <= args.time_start or args.block_duration <= 0 or args.near_width <= 0 or args.pdf_bins < 8:
        parser.error("require ordered interval, positive block/band widths, and >=8 PDF bins")
    run = args.run_dir.resolve()
    output = (args.output_dir or run/"analysis").resolve()
    metadata_path = (args.metadata or run/"benchmark_metadata.json").resolve()
    metadata = json.loads(metadata_path.read_text())
    segments = load_segments(run, metadata, metadata_path)
    paths, _, discarded_snapshots = segment_paths(segments, "snapshots", "**/*.mhd_w_bcc.*.bin", True)
    paths = bracketing_paths(paths, args.time_start, args.time_end)
    if not paths:
        return unavailable_output(output, metadata_path, metadata, args,
            "no retained primitive snapshots overlap the explicit averaging interval", segments)
    records, duplicates, limits, reference_header = {}, [], {}, None
    for path in paths:
        fields, lengths, header, info = read_uniform(path, FIELDS)
        if reference_header is None:
            reference_header = header
            model = infer_model(header, metadata)
            model["nominal_total_forcing_power"] = model["dedt_per_volume"]*float(np.prod(lengths))
            for segment in segments:
                if segment["input"] is not None:
                    crosscheck_input(header, segment["input"])
        elif any(header.get(key) != reference_header.get(key) for key in ("mhd", "problem", "turb_driving")):
            raise ValueError("model parameters changed across retained snapshots")
        digest = hashlib.sha256()
        for name in FIELDS:
            digest.update(np.ascontiguousarray(fields[name]).tobytes())
        info["field_sha256"] = digest.hexdigest()
        key = float(f"{info['time']:.12g}")
        if key in records:
            if records[key]["info"]["field_sha256"] != info["field_sha256"]:
                raise ValueError(f"conflicting physical states at duplicate time {key}")
            duplicates.append({"time": key, "omitted": str(path), "kept": str(records[key]["path"])})
            continue
        scalars, spectra, distributions = snapshot_products(fields, lengths, model, args.near_width)
        for name, value in distributions.items():
            low, high = float(np.min(value)), float(np.max(value))
            old = limits.get(name, (low, high))
            limits[name] = min(low, old[0]), max(high, old[1])
        records[key] = {"path": path, "info": info, "scalars": scalars, "spectra": spectra}
        print(f"analyzed t={key:g}: {path.name}", flush=True)
    ordered = [records[key] for key in sorted(records)]
    times = np.asarray(sorted(records))
    def average(values):
        return temporal_summary(times, values, args.block_duration, args.time_start, args.time_end)
    if any(row["info"]["shape_zyx"] != ordered[0]["info"]["shape_zyx"]
           or row["info"]["lengths_xyz"] != ordered[0]["info"]["lengths_xyz"] for row in ordered):
        raise ValueError("grid geometry changed across snapshots")
    edges = {}
    for name, (low, high) in limits.items():
        pad = max((high-low)*1e-8, 1e-12*max(1, abs(low), abs(high)))
        edges[name] = np.linspace(low-pad, high+pad, args.pdf_bins+1)
    pdf_rows = {name: [] for name in edges}
    # A second streaming pass gives every PDF identical edges without retaining
    # all cell arrays from every snapshot in memory or clipping physical tails.
    for row in ordered:
        fields, _, _, info = read_uniform(row["path"], FIELDS)
        if info["files"] != row["info"]["files"]:
            raise ValueError("snapshot files changed during analysis")
        b2 = sum(fields[name]**2 for name in ("bcc1", "bcc2", "bcc3"))
        distributions = {"B_over_B0": np.sqrt(b2)/model["B0"],
                         "X": 2*(fields["p_perp"]-fields["eint"])/b2}
        for name, values in distributions.items():
            pdf_rows[name].append(volume_pdf(values, edges[name], info["cell_volume"]))
    scalars = {name: average([row["scalars"][name] for row in ordered])
               for name in set.intersection(*(set(row["scalars"]) for row in ordered))}
    unavailable = {}
    scalar_names = set.union(*(set(row["scalars"]) for row in ordered)) | {
        "pressure_correlation", "pressure_normalized_residual_variance"}
    for name in sorted(scalar_names-set(scalars)):
        cause = ("pressure correlation requires nonzero variance in both pressure fluctuations"
                 if name == "pressure_correlation" else
                 "normalized pressure residual requires nonzero summed pressure variance"
                 if name == "pressure_normalized_residual_variance" else
                 "scalar diagnostic is undefined in one or more retained snapshots")
        unavailable["scalars."+name] = unavailable_summary(times,
            [row["scalars"].get(name) for row in ordered], args.time_start, args.time_end, cause)
    spectra = {}
    for name in ordered[0]["spectra"]:
        rows = [row["spectra"][name] for row in ordered]
        spectra[name] = {**average([row["power"] for row in rows]),
            "k": rows[0]["k"], "dk": rows[0]["dk"],
            "max_parseval_relative_error": max(row["parseval_relative_error"] for row in rows),
            "real_space_power": average([row["real_space_power"] for row in rows]),
            "full_k_resolved": average([row["full_k_resolved_power"] for row in rows])}
    pdfs = {name: {**average(rows), "edges": edges[name].tolist()}
            for name, rows in pdf_rows.items()}
    user_paths, user_boundaries, _ = segment_paths(segments, "user_history", "**/*.user.hst")
    mhd_paths, mhd_boundaries, _ = segment_paths(segments, "mhd_history", "**/*.mhd.hst")
    user, user_info = merge_histories(user_paths, user_boundaries)
    mhd, mhd_info = merge_histories(mhd_paths, mhd_boundaries)
    history = history_products(user, mhd, model, args.time_start, args.time_end, args.block_duration)
    health = solver_health(segments, user, mhd, args.time_start, args.time_end)
    force_paths, _, discarded_forcing = segment_paths(segments, "force_snapshots", "**/*.turb_force.*.bin", True)
    force_paths = bracketing_paths(force_paths, args.time_start, args.time_end)
    forces = {}
    for path in force_paths:
        time = paper.snapshot_time(path)
        fields, force_lengths, force_header, info = read_uniform(path, FORCE_FIELDS)
        products = force_products(fields, force_lengths)
        key = float(f"{time:.12g}")
        if key not in records:
            raise ValueError(f"forcing snapshot time {key} lacks a matching primitive snapshot")
        reference = records[key]["info"]
        if info["shape_zyx"] != reference["shape_zyx"] or force_lengths != reference["lengths_xyz"]:
            raise ValueError("forcing and primitive snapshot geometry differ")
        if any(force_header.get(section) != reference_header.get(section)
               for section in ("mesh", "mhd", "problem", "turb_driving")):
            raise ValueError("forcing and primitive snapshot physical input headers differ")
        digest = hashlib.sha256()
        for name in FORCE_FIELDS:
            digest.update(np.ascontiguousarray(fields[name]).tobytes())
        info["field_sha256"] = digest.hexdigest()
        if key in forces and forces[key]["info"]["field_sha256"] != info["field_sha256"]:
            raise ValueError("conflicting forcing fields at duplicate time")
        forces[key] = {"info": info, "values": products}
    force_summary = {"available": bool(forces), "snapshots": list(forces.values())}
    if forces:
        force_times = sorted(forces)
        force_summary["statistics"] = {name: temporal_summary(force_times,
            [forces[t]["values"][name] for t in force_times], args.block_duration,
            args.time_start, args.time_end)
            for name in next(iter(forces.values()))["values"]
            if all(forces[t]["values"][name] is not None for t in force_times)}
        force_summary["expected_solenoidal_fraction"] = model["expected_solenoidal_power_fraction"]
        statistics = force_summary["statistics"]
        for name in next(iter(forces.values()))["values"]:
            if name not in statistics:
                unavailable["forcing_decomposition."+name] = unavailable_summary(force_times,
                    [forces[t]["values"][name] for t in force_times], args.time_start, args.time_end,
                    "solenoidal acceleration-power fraction requires nonzero total nonzero-mode forcing power")
        total = statistics["total_acceleration_power"]["mean"]
        force_summary["ratio_of_time_mean_powers"] = (statistics["solenoidal_acceleration_power"]["mean"]/total if total > 0 else None)
        if total <= 0:
            unavailable["forcing_decomposition.ratio_of_time_mean_powers"] = {
                "available": False, "requested_window": [args.time_start, args.time_end],
                "reason": "time-mean nonzero-mode forcing power is zero; its solenoidal fraction is undefined"}
    else:
        force_summary["reason"] = "no retained turb_force snapshots; fixed-z history components cannot replace a Helmholtz decomposition"
    block_count = next(iter(scalars.values()))["block_count"]
    reasons = []
    if block_count < 4:
        reasons.append(f"only {block_count} complete time blocks; at least four requested for descriptive variability")
    if times[0] > args.time_start+1e-12 or times[-1] < args.time_end-1e-12:
        reasons.append("snapshot coverage does not bracket both requested endpoints; no extrapolation performed")
    if not history["energy_budget"]["available"]:
        reasons.append("actual forcing/conserved-energy budget unavailable")
    elif not history["energy_budget"]["requested_window_covered"]:
        reasons.append("actual forcing/conserved-energy budget covers only "
                       +str(history["energy_budget"]["interval"])
                       +", not the full requested window")
    if not forces:
        reasons.append("realized forcing-mode mixture unmeasured")
    elif set(forces) != set(records):
        reasons.append("forcing snapshots cover only part of the retained primitive sampling times")
    for label, info in (("user", user_info), ("mhd", mhd_info)):
        if any(row["changed_columns"] for row in info.get("duplicate_rows", [])):
            reasons.append(f"{label} duplicate history times contain changed columns; inspect retained audit")
    autocorr = {name: autocorrelation(times, [row["scalars"][name] for row in ordered])
        for name in ("deltaB_rms_over_B0", "beta_volume_mean", "S_parallel_rms", "u_parallel_rms")}
    nxyz = np.asarray(ordered[0]["info"]["shape_zyx"])[::-1]
    nyquist = min(np.pi*nxyz[0]/lengths[0], np.pi*nxyz[1]/lengths[1])
    full_nyquist = min(np.pi*n/l for n, l in zip(nxyz, lengths))
    scales = {"k_perp_nyquist_inscribed": float(nyquist), "resolved_kmax": float(.25*nyquist),
              "full_k_resolved_max": float(.25*full_nyquist),
              "cutoff_band_kmin": float(.75*nyquist),
              "definition": "resolved guide <=Nyquist/4 (>=8 cells per axis wavelength); scalar resolved ratios require 0<full |k|<=min(Nyquist_xyz)/4, including pure parallel modes at k_perp=0; plotted perpendicular shells sum all k_parallel, so low-k_perp alone does not guarantee resolved full wavevector; not a convergence claim"}
    ratio = {}
    kval = np.asarray(spectra["S_parallel"]["k"])
    for label, mask in (("resolved", np.ones(kval.shape, dtype=bool)),
                        ("cutoff", kval >= .75*nyquist)):
        numerator_row = spectra["S_parallel"]["full_k_resolved"] if label == "resolved" else spectra["S_parallel"]
        denominator_row = spectra["perp_gradient_u_perp"]["full_k_resolved"] if label == "resolved" else spectra["perp_gradient_u_perp"]
        numerator = np.sum(np.asarray(numerator_row["mean"])[mask])
        denominator = np.sum(np.asarray(denominator_row["mean"])[mask])
        ratio[label+"_strain_to_perpendicular_gradient_power"] = float(numerator/denominator) if denominator > 0 else None
        if denominator <= 0:
            unavailable["gradient_comparisons."+label+"_strain_to_perpendicular_gradient_power"] = {
                "available": False, "requested_window": [args.time_start, args.time_end],
                "reason": "perpendicular-gradient comparator has zero power in the declared "+label+" band"}
    revisions = [segment["metadata"].get("simulation", {}).get("revision") for segment in segments]
    if any(not value for value in revisions):
        reasons.append("simulation revision absent from retained launch provenance")
    try:
        revision = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        revision = None
    classification = "inconclusive"
    parseval_error = max(rec["max_parseval_relative_error"] for rec in spectra.values())
    pdf_error = max(abs(np.dot(rec["mean"], np.diff(rec["edges"]))-1) for rec in pdfs.values())
    pressure_corr = scalars.get("pressure_correlation", {}).get("mean")
    pressure_residual = scalars.get("pressure_normalized_residual_variance", {}).get("mean")
    pressure_direction = ("consistent" if pressure_corr is not None and pressure_corr < 0
                          and pressure_residual is not None and pressure_residual < 1 else "inconclusive")
    if pressure_corr is not None and pressure_corr > 0:
        pressure_direction = "concerning"
    findings = {
        "marginality": {"classification": "inconclusive", "reason": "report near-band residence and strict/rounding-robust crossings separately; snapshot precision and projection cadence do not establish full-time admissibility"},
        "pressure_balance": {"classification": pressure_direction, "reason": "negative covariance and R<1 indicate compensation relative to uncorrelated fields; extent and variability require review, with no paper amplitude tolerance", "correlation": pressure_corr, "residual": pressure_residual},
        "gradients": {"classification": "inconclusive", "reason": "resolved/cutoff tensor ratios and parallel velocity are measured; one active box cannot establish suppression caused by anisotropy feedback", **ratio},
        "spectra_energy": {"classification": "inconclusive", "reason": "Parseval and measured forcing ledger available as stated; thermal drift and cascade extent require finite-window review, without a target slope"}}
    result = {"schema_version": 1, "definitions": DEFINITIONS,
        "requested_window": [args.time_start, args.time_end],
        "retained_window": [max(float(times[0]), args.time_start), min(float(times[-1]), args.time_end)],
        "sampling": {"snapshots": len(times), "times": times.tolist(), "near_threshold_halfwidth_X": args.near_width,
            "bracketing_snapshot_range": [float(times[0]), float(times[-1])],
            "max_snapshot_gap": float(max(np.diff(times))) if len(times) > 1 else None,
            "endpoint_rule": "linearly interpolate diagnostic values/PDF bins/spectral bins at requested endpoints using bracketing snapshots; no field interpolation or extrapolation",
            "block_duration": args.block_duration, "complete_blocks": block_count,
            "block_anchor": max(float(times[0]), args.time_start),
            "instantaneous_output_caveat": "post-operator snapshots may miss transient threshold excursions; inclusive history switches, strict exceedance, and near bands have different meanings"},
        "model": model, "scalars": scalars, "spectra": spectra, "PDFs": pdfs,
        "scale_bands": scales, "gradient_comparisons": ratio, "history": history,
        "forcing_decomposition": force_summary, "empirical_autocorrelation": autocorr,
        "unavailable_metric_reasons": unavailable,
        "temporal_dependence_cautions": [f"{name}: block duration is less than twice the empirical positive-sequence integral correlation time"
            for name, row in autocorr.items() if row.get("available") and args.block_duration < 2*row["positive_sequence_integral_time"]],
        "numerical_integrity": {"classification": "consistent" if max(parseval_error, pdf_error) < 1e-10 else "concerning",
            "max_parseval_relative_error": parseval_error, "max_PDF_integral_error": pdf_error,
            "scope": "analysis normalization/finite-field checks, not simulation physics acceptance"},
        "simulation_integrity": health,
        "directional_findings": findings,
        "adequacy": {"classification": classification, "reasons": reasons,
            "sampling_prerequisites_satisfied": not reasons,
            "scope": "descriptive single-box finite-window comparison; no causal active/passive, asymptotic convergence, or LF-coefficient claim",
            "heating": "no cooling: secular thermal-energy/beta drift is expected and must be reported, not mistaken for a stationary thermal ensemble",
            "classification_rule": "overall scientific classification remains inconclusive pending joint scientific review; sampling prerequisites, numerical analysis checks and group directional findings are separate"},
        "snapshots": [{"info": row["info"], "scalars": row["scalars"]} for row in ordered],
        "provenance": {"metadata": retained_file(metadata_path), "retained_metadata": metadata,
            "simulation_revision": revisions[0] if len(set(revisions)) == 1 else revisions,
            "analysis_revision": revision,
            "analysis_script": retained_file(Path(__file__)), "reader_script": retained_file(Path(paper.__file__)),
            "binary_reader_script": retained_file(Path(paper.bin_convert.__file__)),
            "software": {"python": sys.version, "numpy": np.__version__, "matplotlib": version("matplotlib")},
            "embedded_input": reference_header, "duplicate_snapshots": duplicates,
            "discarded_snapshot_branches": discarded_snapshots, "discarded_forcing_branches": discarded_forcing,
            "segments": [{key: str(value) if isinstance(value, Path) else value for key, value in segment.items()} for segment in segments],
            "user_history": user_info, "mhd_history": mhd_info,
            "argv": sys.argv if argv is None else argv}}
    output.mkdir(parents=True, exist_ok=True)
    result = json_value(result)
    (output/"metrics.json").write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
    make_figures(result, output)
    write_report(result, output)
    print(f"{classification}: {output/'report.md'}", flush=True)
    return 0


def unavailable_output(output, metadata_path, metadata, args, reason, segments):
    """Retain a machine-readable missing-data result and shareable placeholders."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    output.mkdir(parents=True, exist_ok=True)
    revisions = [segment["metadata"].get("simulation", {}).get("revision") for segment in segments]
    result = {"schema_version": 1, "adequacy": {"classification": "inconclusive", "reasons": [reason]},
        "requested_window": [args.time_start, args.time_end], "definitions": DEFINITIONS,
        "provenance": {"metadata": retained_file(metadata_path), "retained_metadata": metadata,
                       "simulation_revision": revisions[0] if len(set(revisions)) == 1 else revisions,
                       "segments": [{key: str(value) if isinstance(value, Path) else value
                                     for key, value in segment.items()} for segment in segments],
                       "analysis_script": retained_file(Path(__file__))}}
    (output/"metrics.json").write_text(json.dumps(result, indent=2)+"\n")
    lines = ["# CGL-LF single-run physics benchmark", "", "**Inconclusive: missing data.** "+reason+".", ""]
    for name in ("marginality", "pressure_balance", "gradients", "spectra_energy"):
        fig, ax = plt.subplots(figsize=(9, 3))
        ax.axis("off")
        ax.text(.05, .7, name.replace("_", " "), fontsize=16)
        ax.text(.05, .4, reason, wrap=True)
        for suffix in ("png", "pdf"):
            fig.savefig(output/f"{name}.{suffix}", dpi=150)
        plt.close(fig)
        lines.append(f"- [{name}]({name}.png), [PDF]({name}.pdf).")
    lines += ["", "[metrics.json](metrics.json) retains provenance and missing-data reason; no physical comparison was made."]
    (output/"report.md").write_text("\n".join(lines)+"\n")
    return 0


def make_figures(data, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": .2})
    times = data["sampling"]["times"]
    snapshots = data["snapshots"]
    scalars = data["scalars"]

    def series(ax, key, label=None):
        ax.plot(times, [row["scalars"].get(key, np.nan) for row in snapshots], label=label or key)

    def spec(ax, key, label=None):
        record = data["spectra"][key]
        k, mean = np.asarray(record["k"]), np.asarray(record["mean"])
        valid = mean > 0
        ax.loglog(k[valid], mean[valid], label=label or key)
        if record["block_min"] is not None:
            low, high = np.asarray(record["block_min"]), np.asarray(record["block_max"])
            ax.fill_between(k[valid], np.maximum(low[valid], np.finfo(float).tiny), high[valid], alpha=.12)
        ax.set_xlabel(r"$k_\perp$ [radians / length]")
        if not getattr(ax, "_benchmark_scale_marked", False):
            ax.axvspan(k[0], data["model"]["forcing_kmax"], color="tab:green", alpha=.07,
                       label="projected forcing support")
            ax.axvline(data["scale_bands"]["resolved_kmax"], color="black", ls="--", lw=.8,
                       label="8-cell resolved guide")
            ax.axvspan(data["scale_bands"]["cutoff_band_kmin"], k[-1], color="grey", alpha=.12,
                       label="perpendicular cutoff band")
            ax._benchmark_scale_marked = True

    def finish(fig, name):
        for ax in fig.axes:
            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(fontsize=8)
        fig.suptitle(f"CGL-LF finite-window benchmark: t={data['retained_window']} | {data['sampling']['complete_blocks']} blocks")
        fig.tight_layout(rect=(0, 0, 1, .96))
        for extension in ("png", "pdf"):
            fig.savefig(output/f"{name}.{extension}", dpi=170)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for ax, key in zip(axes[0], ("B_over_B0", "X")):
        rec = data["PDFs"][key]
        centers = .5*(np.asarray(rec["edges"][:-1])+rec["edges"][1:])
        ax.plot(centers, rec["mean"], label="time-mean volume PDF")
        if rec["block_min"] is not None:
            ax.fill_between(centers, rec["block_min"], rec["block_max"], alpha=.2, label="block range")
        ax.set(xlabel=key, ylabel="probability density")
    for label in ("mirror", "firehose"):
        axes[0, 1].axvline(data["model"][label+"_X"], ls="--", color="black", lw=.8)
        series(axes[1, 0], label+"_strict", label+" strict exceedance")
        series(axes[1, 1], label+"_near", label+" symmetric near band")
        series(axes[1, 1], label+"_near_interior", label+" near interior")
    hist = data["history"]["series"]
    for label, name in (("mirror", "mirror"), ("firehose", "fire")):
        key = name+"_inclusive_history_fraction"
        if key in hist:
            axes[1, 0].plot(hist["time"], hist[key], alpha=.5, label=label+" inclusive history")
    for ax in axes[1]:
        ax.set(xlabel="time", ylabel="volume fraction")
    axes[1, 0].axvspan(*data["requested_window"], color="black", alpha=.08, label="averaging window")
    axes[1, 1].set_xlim(data["requested_window"])
    axes[1, 1].set_title(f"near halfwidth |X-Xthreshold| ≤ {data['sampling']['near_threshold_halfwidth_X']:g}")
    finish(fig, "marginality")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    for key in ("p_parallel", "p_perp", "magnetic_pressure"):
        spec(axes[0], key)
    axes[0].set_ylabel("pressure variance / dk")
    series(axes[1], "pressure_correlation", "signed pressure correlation")
    axes[1].axhline(0, color="black", lw=.7)
    axes[1].set(xlabel="time", ylabel="correlation", ylim=(-1.05, 1.05))
    series(axes[2], "pressure_normalized_residual_variance", "normalized residual variance")
    axes[2].set(xlabel="time", ylabel="R (0 = exact compensation)")
    finish(fig, "pressure_balance")

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for key in ("S_parallel", "parallel_gradient_u_perp", "perp_gradient_u_parallel", "perp_gradient_u_perp"):
        spec(axes[0, 0], key)
    for key in ("S_parallel", "div_u", "induction", "b_grad_u_parallel"):
        spec(axes[0, 1], key)
    for ax in axes[0]:
        ax.set_ylabel("gradient fluctuation power / dk")
    for key in ("u_parallel_rms", "u_perp_rms"):
        series(axes[1, 0], key)
    axes[1, 0].set(xlabel="time", ylabel="local-field velocity RMS")
    for key in ("S_parallel_rms", "induction_rms", "curvature_plus_discretization_rms"):
        series(axes[1, 1], key)
    axes[1, 1].set(xlabel="time", ylabel="gradient RMS")
    finish(fig, "gradients")

    fig, axes = plt.subplots(2, 3, figsize=(16, 8))
    for key in ("kinetic", "magnetic"):
        spec(axes[0, 0], key)
    axes[0, 0].set_ylabel("fluctuation energy density / dk")
    if hist:
        for key in ("kinetic", "magnetic", "therm_cgl", "force_work"):
            if key in hist:
                axes[0, 1].plot(hist["time"], hist[key], label=key)
        for ax in (axes[0, 1],):
            ax.axvspan(*data["requested_window"], alpha=.08, color="black")
        axes[0, 1].set(xlabel="time", ylabel="domain-integrated energy / accumulated work")
    else:
        axes[0, 1].text(.05, .5, "Missing energy histories", transform=axes[0, 1].transAxes)
    for key in ("Mach_isotropic_proxy", "deltaB_rms_over_B0", "u_parallel_fraction"):
        series(axes[1, 0], key)
    axes[1, 0].set(xlabel="time", ylabel="dimensionless measured amplitude")
    series(axes[1, 1], "beta_volume_mean", "volume mean beta")
    series(axes[1, 1], "beta_ratio_of_means", "ratio-of-means beta")
    axes[1, 1].set(xlabel="time", ylabel="beta (heating is retained)")
    comparison = common_energy_history(data["history"])
    if comparison:
        axes[0, 2].plot(comparison["energy"]["time"], comparison["energy"]["increment"], label="change in conserved E")
        axes[0, 2].plot(comparison["work"]["time"], comparison["work"]["increment"], ls="--", label="actual accumulated forcing work")
        axes[0, 2].set(xlabel="time", ylabel=f"energy increment from common t={comparison['interval'][0]:g}")
        axes[0, 2].axvspan(*data["requested_window"], color="black", alpha=.08)
    else:
        axes[0, 2].text(.05, .5, "Overlapping forcing/energy ledgers unavailable", transform=axes[0, 2].transAxes)
    forcing = data["forcing_decomposition"]
    if forcing["available"]:
        rows = sorted(forcing["snapshots"], key=lambda row: row["info"]["time"])
        axes[1, 2].plot([row["info"]["time"] for row in rows],
            [row["values"]["solenoidal_fraction"] for row in rows], label="realized acceleration solenoidal fraction")
        expected = forcing["expected_solenoidal_fraction"]
        if expected is not None:
            axes[1, 2].axhline(expected, color="black", ls="--", label="nominal ratio of expected innovation powers")
        axes[1, 2].set(xlabel="time", ylabel="Helmholtz acceleration-power fraction", ylim=(0, 1))
    else:
        axes[1, 2].text(.05, .5, "Forcing snapshots unavailable", transform=axes[1, 2].transAxes)
    finish(fig, "spectra_energy")


def write_report(data, output):
    sampling = data["sampling"]
    lines = ["# CGL-LF single-run physics benchmark", "",
        f"**{data['adequacy']['classification'].capitalize()}** for a finite-window descriptive comparison.", "",
        f"Requested interval {data['requested_window']}; retained interval {data['retained_window']}; "
        f"{sampling['snapshots']} snapshots and {sampling['complete_blocks']} complete blocks of duration {sampling['block_duration']}.", "",
        "Block bands show temporal variability, not independent-snapshot confidence intervals. No cooling is present: evolving thermal energy/beta are reported rather than assumed stationary.", "",
        f"Bracketing snapshot times span {sampling['bracketing_snapshot_range']}; maximum gap {sampling['max_snapshot_gap']}. "
        "Diagnostic values (not cell fields) are linearly interpolated to covered requested endpoints. "
        f"Snapshot blocks are anchored at t={sampling['block_anchor']:g}; this equals the requested start only when its coverage is available. No extrapolation is performed.", ""]
    for reason in data["adequacy"]["reasons"]:
        lines.append(f"- {reason}.")
    lines += ["", "| Figure group | PDF |", "| --- | --- |"]
    for name in ("marginality", "pressure_balance", "gradients", "spectra_energy"):
        lines.append(f"| [{name}]({name}.png) | [vector PDF]({name}.pdf) |")
    lines += ["", "| Directional finding | Classification | Reason |", "| --- | --- | --- |"]
    for name, finding in data["directional_findings"].items():
        lines.append(f"| {name} | {finding['classification']} | {finding['reason']} |")
    lines += ["", "Analysis normalization checks: `"+json.dumps(data["numerical_integrity"])+"`. "
        "These do not classify simulation health or replace scientific review.", ""]
    lines += ["**Simulation integrity: "+data["simulation_integrity"]["classification"]+".** "
        +"; ".join(data["simulation_integrity"]["reasons"]), "",
        "Numerical floor/nonfinite/nonpositive counters: `"+json.dumps(data["simulation_integrity"]["numerical_counters"])+"`.", "",
        "Unprojected LF stage crossings and reserved restart field (not automatically failures): `"+json.dumps(data["simulation_integrity"]["physical_stage_activity"])+"`.", "",
        "**Projection-count limitation:** in the audited post-WO2 simulation revision, `lf_hwproj` is reserved/uninstrumented. "
        "Its raw value and increment do not measure wall-projection activity, and zero does not demonstrate absence of projections. "
        "`lf_hardbd` and nonfinite/nonpositive counters count active-cell stage checks; density/pressure-floor counters count EOS refresh events including refreshed halo cells. "
        "All are instrumented cumulative counts, preserved across restart and history output.", "",
        "Post-operator hard-bound history: `"+json.dumps(data["simulation_integrity"]["post_operator_hard_volume"])+"`.", ""]
    lines += ["", "| Measurement | Time mean | Block SD |", "| --- | ---: | ---: |"]
    for name in ("Mach_isotropic_proxy", "deltaB_rms_over_B0", "beta_volume_mean", "u_parallel_fraction",
                 "mirror_strict", "firehose_strict", "mirror_near", "firehose_near",
                 "pressure_correlation", "pressure_normalized_residual_variance", "S_parallel_rms", "induction_rms"):
        rec = data["scalars"].get(name)
        if rec:
            sd = "unavailable" if rec["block_sd"] is None else f"{rec['block_sd']:.6g}"
            lines.append(f"| {name} | {rec['mean']:.6g} | {sd} |")
        elif "scalars."+name in data["unavailable_metric_reasons"]:
            lines.append(f"| {name} | unavailable (see reason below) | unavailable |")
    if data["unavailable_metric_reasons"]:
        lines += ["", "Unavailable metric summaries (valid instantaneous samples remain retained; the requested window is unchanged):", ""]
        for name, record in data["unavailable_metric_reasons"].items():
            times_note = (" Undefined sample times: "+str(record["undefined_sample_times"])+
                          "; required endpoint-bracketing times: "+str(record["undefined_bracketing_times"])+"."
                          if "undefined_sample_times" in record else "")
            lines.append(f"- **{name}**: {record['reason']}{times_note}")
    lines += ["", "Pressure spectra alone cannot establish compensation; use the signed correlation and normalized residual together. "
        "Local S_parallel is not b·grad(u_parallel). The latter includes field-direction curvature; S_parallel-div(u) is the ideal compressible induction proxy. "
        "Parallel velocity is retained explicitly. Compare the resolved and cutoff bands separately. "
        "The conservative resolved guide is Nyquist/4 (eight cells per wavelength); scalar resolved ratios also restrict full |k|. "
        "Plotted k_perp shells sum all k_parallel, so small k_perp alone does not imply fully resolved gradients. "
        "Green shading marks the projection of the physical forcing shell onto k_perp, including zero; no universal spectral slope is prescribed.", "",
        "Gradient band measurements: `"+json.dumps(data["gradient_comparisons"])+"`.", "",
        "Actual applied forcing energy budget: `"+json.dumps(data["history"]["energy_budget"])+"`.", ""]
    lines.append(f"Retained dedt={data['model']['dedt_per_volume']} is nominal power per domain volume; "
        f"nominal total power={data['model']['nominal_total_forcing_power']}. The measured accumulated work remains authoritative.")
    force = data["forcing_decomposition"]
    if force["available"]:
        fraction = force.get("statistics", {}).get("solenoidal_fraction")
        if fraction is None:
            reason = data["unavailable_metric_reasons"]["forcing_decomposition.solenoidal_fraction"]["reason"]
            lines.append("Time-mean instantaneous Helmholtz solenoidal acceleration-power fraction unavailable: "+reason)
        else:
            lines.append("Realized Helmholtz solenoidal acceleration-power fraction: `"+json.dumps(fraction)+"`.")
        lines.append(f"Nominal ratio of expected isotropic innovation powers: {force['expected_solenoidal_fraction']}. "
            "This is not the expectation of an instantaneous fraction or a target for this finite OU realization. "
            f"Ratio of time-mean solenoidal to total acceleration powers (a different statistic): {force['ratio_of_time_mean_powers']}. "
            "Neither statistic is an energy-injection partition.")
    else:
        lines.append("Forcing decomposition unavailable: "+force["reason"]+".")
    lines += ["", "| Temporal diagnostic | Integral correlation time | Caution |", "| --- | ---: | --- |"]
    for name, row in data["empirical_autocorrelation"].items():
        tau = row.get("positive_sequence_integral_time", "unavailable")
        lines.append(f"| {name} | {tau} | {row.get('reason', 'Finite-window estimate; inspect drift')} |")
    for warning in data["temporal_dependence_cautions"]:
        lines.append("\n- "+warning+".")
    lines += ["", "| History drift | Fitted change in window | Block means |", "| --- | ---: | --- |"]
    for name in ("kinetic", "magnetic", "therm_cgl", "beta_volume_mean"):
        row = data["history"]["window"].get(name)
        if row:
            lines.append(f"| {name} | {row.get('fitted_change_over_window')} | {row['block_means']} |")
    lines += ["", "Empirical autocorrelation and linear trends are retained in metrics.json. Positive-sequence correlation times are finite-window estimates, not evidence that OU-sized blocks are independent. "
        "Fewer than four complete blocks warrants extending this same realization before quoting variability; a longer heated run is still a finite-time ensemble.", "",
        "Thresholds: mirror X="+str(data["model"]["mirror_X"])+", firehose X="+str(data["model"]["firehose_X"])+
        f"; symmetric near halfwidth {sampling['near_threshold_halfwidth_X']}. Strict exceedance, inclusive solver history switches, and near-threshold bands are separate. "
        "Snapshot payload precision and post-operator/source/projection cadence limit admissibility claims; the dynamic half-ULP float32 sensitivity envelope is a representation bound, not altered physical thresholds or a bound on dynamical error.", "",
        "One box does not establish active/passive causality, spatial convergence, a universal spectral slope, or LF coefficient correctness. "
        "The directional classification above uses explicitly stated coverage prerequisites and signed compensation evidence, not a fitted paper-image tolerance.", "",
        "## Definitions", ""]
    lines += [f"- **{name}**: {definition}." for name, definition in data["definitions"].items()]
    lines += ["", "## Provenance", "",
        f"Simulation revision from retained launch metadata: `{data['provenance']['simulation_revision']}`. "
        f"Analysis checkout revision: `{data['provenance']['analysis_revision']}` (separate from simulation provenance).", "",
        "[metrics.json](metrics.json) includes embedded effective input, file hashes, analysis/reader hashes, runtime metadata, sampling, block means, Parseval errors, and restart-time deduplication audits.", ""]
    (output/"report.md").write_text("\n".join(lines))


if __name__ == "__main__":
    raise SystemExit(main())
