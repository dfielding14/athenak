#!/usr/bin/env python3
"""Extract shock-rich engineering diagnostics from AthenaK output artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from typing import Dict

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import bin_convert_new as bin_convert  # noqa: E402

from pvtk_particles import read_particle_vtk  # noqa: E402
from artifact_lineage import write_companion_manifest  # noqa: E402


_CYCLE_RE = re.compile(r"\.(\d+)\.bin$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as fobj:
        for chunk in iter(lambda: fobj.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _files_by_cycle(bin_dir: Path, basename: str, file_id: str) -> Dict[int, Path]:
    mapping: Dict[int, Path] = {}
    pattern = f"{basename}.{file_id}.*.bin"
    for path in sorted(bin_dir.glob(pattern)):
        match = _CYCLE_RE.search(path.name)
        if match is None:
            continue
        mapping[int(match.group(1))] = path
    return mapping


def _profile_and_shock(data: dict) -> tuple[np.ndarray, np.ndarray, float, int]:
    x = np.asarray(data["x1v"], dtype=float)
    rho = np.asarray(data["dens"], dtype=float)

    if rho.ndim == 3:
        rho_x = np.mean(rho, axis=(0, 1))
    elif rho.ndim == 2:
        rho_x = np.mean(rho, axis=0)
    else:
        rho_x = rho

    grad = np.gradient(rho_x, x)
    idx = int(np.argmax(np.abs(grad)))
    return x, rho_x, float(x[idx]), idx


def _periodic_masks(
    x: np.ndarray, x0: float, frac: float
) -> tuple[np.ndarray, np.ndarray]:
    xmin = float(np.min(x))
    xmax = float(np.max(x))
    length = xmax - xmin
    width = frac * length
    rel = np.mod(x - x0 + 0.5 * length, length) - 0.5 * length
    upstream = (rel > 0.0) & (rel <= width)
    downstream = (rel < 0.0) & (np.abs(rel) <= width)
    return upstream, downstream


def _reduce_to_xprofile(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 3:
        return np.mean(arr, axis=(0, 1))
    if arr.ndim == 2:
        return np.mean(arr, axis=0)
    return arr


def _safe_mean(arr: np.ndarray, mask: np.ndarray) -> float:
    if np.count_nonzero(mask) == 0:
        return float("nan")
    return float(np.mean(arr[mask]))


def _read_pvtk_velocity(pvtk_path: Path) -> tuple[np.ndarray, np.ndarray]:
    pdata = read_particle_vtk(pvtk_path)
    if "vel" in pdata.vectors:
        vel = np.asarray(pdata.vectors["vel"], dtype=float)
    elif pdata.vectors:
        first_key = sorted(pdata.vectors.keys())[0]
        vel = np.asarray(pdata.vectors[first_key], dtype=float)
    else:
        raise RuntimeError("No vector velocity field found in particle VTK file")

    pos = np.asarray(pdata.points, dtype=float)
    p_mag = np.linalg.norm(vel, axis=1)
    return pos[:, 0], p_mag


def _log_bins(values: np.ndarray, nbins: int) -> np.ndarray:
    vals = values[np.isfinite(values)]
    vals = vals[vals > 0.0]
    if vals.size == 0:
        return np.logspace(-6.0, 0.0, nbins + 1)

    vmin = float(np.min(vals))
    vmax = float(np.max(vals))
    if not math.isfinite(vmin) or not math.isfinite(vmax) or vmax <= 0.0:
        return np.logspace(-6.0, 0.0, nbins + 1)
    if vmin == vmax:
        vmin *= 0.8
        vmax *= 1.2
    return np.logspace(np.log10(vmin), np.log10(vmax), nbins + 1)


def _write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract shock profiles, x-p proxy, and spectra diagnostics."
    )
    parser.add_argument("--basename", required=True)
    parser.add_argument(
        "--bin-dir",
        default=str(REPO_ROOT / "tst" / "build" / "src" / "bin"),
    )
    parser.add_argument(
        "--pvtk-dir",
        default=str(REPO_ROOT / "tst" / "build" / "src" / "pvtk"),
    )
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--window-frac", type=float, default=0.15)
    parser.add_argument("--x-bins", type=int, default=128)
    parser.add_argument("--p-bins", type=int, default=96)
    parser.add_argument("--save-csv", action="store_true")
    args = parser.parse_args()

    bin_dir = Path(args.bin_dir).resolve()
    pvtk_dir = Path(args.pvtk_dir).resolve()
    out_dir = Path(args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    rho_map = _files_by_cycle(bin_dir, args.basename, "rho")
    bmag_map = _files_by_cycle(bin_dir, args.basename, "bmag")
    j2_map = _files_by_cycle(bin_dir, args.basename, "j2")

    shared = sorted(set(rho_map) & set(bmag_map) & set(j2_map))
    if not shared:
        raise RuntimeError("No shared rho/bmag/j2 cycles found for basename")

    series = []
    latest_payload = {}
    for cyc in shared:
        rho_data = bin_convert.read_binary_as_athdf(str(rho_map[cyc]))
        bmag_data = bin_convert.read_binary_as_athdf(str(bmag_map[cyc]))
        j2_data = bin_convert.read_binary_as_athdf(str(j2_map[cyc]))

        x, rho_x, xshock, ishock = _profile_and_shock(rho_data)
        bmag_x = _reduce_to_xprofile(np.asarray(bmag_data["bmag"], dtype=float))
        j2_x = _reduce_to_xprofile(np.asarray(j2_data["j2"], dtype=float))

        umask, dmask = _periodic_masks(x, xshock, args.window_frac)

        row = {
            "cycle": int(cyc),
            "time": float(rho_data["Time"]),
            "shock_index": int(ishock),
            "shock_x": float(xshock),
            "rho_upstream": _safe_mean(rho_x, umask),
            "rho_downstream": _safe_mean(rho_x, dmask),
            "bmag_upstream": _safe_mean(bmag_x, umask),
            "bmag_downstream": _safe_mean(bmag_x, dmask),
            "j2_upstream": _safe_mean(j2_x, umask),
            "j2_downstream": _safe_mean(j2_x, dmask),
        }
        series.append(row)

        if cyc == shared[-1]:
            latest_payload = {
                "x": x,
                "rho_x": rho_x,
                "bmag_x": bmag_x,
                "j2_x": j2_x,
                "umask": umask.astype(np.int8),
                "dmask": dmask.astype(np.int8),
                "shock_x": float(xshock),
                "shock_cycle": int(cyc),
                "shock_time": float(rho_data["Time"]),
            }

    if args.save_csv:
        csv_path = out_dir / "shock_series.csv"
        with csv_path.open("w", encoding="utf-8", newline="") as fobj:
            writer = csv.DictWriter(fobj, fieldnames=list(series[0].keys()))
            writer.writeheader()
            for row in series:
                writer.writerow(row)

    shock_series_json = out_dir / "shock_series.json"
    _write_json(
        shock_series_json,
        {
            "basename": args.basename,
            "n_samples": len(series),
            "window_frac": args.window_frac,
            "rows": series,
        },
    )

    profiles_npz = out_dir / "shock_profiles_latest.npz"
    np.savez(
        profiles_npz,
        x=latest_payload["x"],
        rho_x=latest_payload["rho_x"],
        bmag_x=latest_payload["bmag_x"],
        j2_x=latest_payload["j2_x"],
        upstream_mask=latest_payload["umask"],
        downstream_mask=latest_payload["dmask"],
        shock_x=np.array([latest_payload["shock_x"]], dtype=float),
        shock_cycle=np.array([latest_payload["shock_cycle"]], dtype=int),
        shock_time=np.array([latest_payload["shock_time"]], dtype=float),
    )

    pvtk_files = sorted(pvtk_dir.glob(f"{args.basename}.prtcl_all.*.part.vtk"))
    if not pvtk_files:
        raise RuntimeError("No particle VTK files found for basename")

    pvtk_path = pvtk_files[-1]
    xprt, p_mag = _read_pvtk_velocity(pvtk_path)

    shock_x = float(latest_payload["shock_x"])
    xbins = np.linspace(
        float(np.min(latest_payload["x"])),
        float(np.max(latest_payload["x"])),
        args.x_bins + 1,
    )
    pbins = _log_bins(p_mag, args.p_bins)

    hxp, xedges, pedges = np.histogram2d(xprt, p_mag, bins=(xbins, pbins))
    hglobal, _ = np.histogram(p_mag, bins=pbins)

    umask_p, dmask_p = _periodic_masks(xprt, shock_x, args.window_frac)
    hup, _ = np.histogram(p_mag[umask_p], bins=pbins)
    hdown, _ = np.histogram(p_mag[dmask_p], bins=pbins)

    phase_space_npz = out_dir / "phase_space_xp.npz"
    np.savez(
        phase_space_npz,
        hist=hxp,
        x_edges=xedges,
        p_edges=pedges,
        pvtk_file=np.array([str(pvtk_path)], dtype=object),
    )

    spectra_npz = out_dir / "spectra.npz"
    np.savez(
        spectra_npz,
        p_edges=pbins,
        dndp_global=hglobal,
        dndp_upstream=hup,
        dndp_downstream=hdown,
    )

    summary = {
        "basename": args.basename,
        "evidence_class": "engineering_proxy",
        "not_sun_bai_reproduction": True,
        "qualification_status": "unqualified_engineering_diagnostics",
        "physical_model": "orszag_tang_shock_rich_engineering_stress",
        "bin_dir": str(bin_dir),
        "pvtk_dir": str(pvtk_dir),
        "output_dir": str(out_dir),
        "n_profile_samples": len(shared),
        "n_particles_latest": int(p_mag.size),
        "shock_x_latest": float(shock_x),
        "pvtk_latest": str(pvtk_path),
        "pvtk_latest_sha256": _sha256(pvtk_path),
        "rho_latest_sha256": _sha256(rho_map[shared[-1]]),
        "bmag_latest_sha256": _sha256(bmag_map[shared[-1]]),
        "j2_latest_sha256": _sha256(j2_map[shared[-1]]),
        "profile_cycles": [int(c) for c in shared],
    }
    diagnostics_summary = out_dir / "diagnostics_summary.json"
    _write_json(diagnostics_summary, summary)
    consumed = [pvtk_path]
    for cyc in shared:
        consumed.extend((rho_map[cyc], bmag_map[cyc], j2_map[cyc]))
    outputs = [
        shock_series_json,
        profiles_npz,
        phase_space_npz,
        spectra_npz,
        diagnostics_summary,
    ]
    if args.save_csv:
        outputs.append(out_dir / "shock_series.csv")
    write_companion_manifest(
        out_dir / "diagnostics_artifact_manifest.json",
        generator="analyze_pic_shock_storyline.py",
        physical_model="orszag_tang_shock_rich_engineering_stress",
        inputs=consumed,
        outputs=outputs,
    )

    print("diagnostics_summary:", out_dir / "diagnostics_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
