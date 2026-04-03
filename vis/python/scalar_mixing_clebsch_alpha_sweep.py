#!/usr/bin/env python3
"""Run a 3D Clebsch alpha sweep and plot realized velocity power spectra."""

from __future__ import annotations

import argparse
import glob
import json
import math
import re
import shlex
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from bin_convert_new import read_binary_as_athdf


ALPHAS = (1.0 / 6.0, 1.0 / 3.0, 1.0 / 2.0)
COLORS = ("#275dad", "#8b3fb3", "#bd632f")


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_athena_path() -> Path:
    for relpath in (
        Path("build-codex/src/athena"),
        Path("build-mpi/src/athena"),
        Path("build/src/athena"),
    ):
        candidate = _repo_root() / relpath
        if candidate.exists():
            return candidate
    return _repo_root() / "build-codex" / "src" / "athena"


def _run_root(athena: Path) -> Path:
    return athena.parent.parent if athena.parent.name == "src" else athena.parent


def _parse_init_seconds(log_text: str) -> float | None:
    match = re.search(r"velocity_init_wall_seconds=([0-9.eE+-]+)", log_text)
    return None if match is None else float(match.group(1))


def _shell_spectrum(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ntot = float(vx.size)
    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    energy = (np.abs(uhat_x) ** 2 + np.abs(uhat_y) ** 2 + np.abs(uhat_z) ** 2) / (ntot * ntot)

    nz, ny, nx = vx.shape
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kx = np.fft.fftfreq(nx) * nx
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    shell = np.ceil(np.sqrt(kx3 * kx3 + ky3 * ky3 + kz3 * kz3) - 1.0e-12).astype(np.int64)
    shell_energy = np.bincount(shell.ravel(), weights=energy.ravel())
    kvals = np.arange(shell_energy.size, dtype=np.float64)
    return kvals, shell_energy


def _fit_slope(kvals: np.ndarray, spectrum: np.ndarray, kmin: int, kmax: int) -> float:
    mask = (kvals >= kmin) & (kvals <= kmax) & (spectrum > 0.0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    coeff = np.polyfit(np.log(kvals[mask]), np.log(spectrum[mask]), 1)
    return float(coeff[0])


def _shape_normalize(spectrum: np.ndarray, kref: int) -> np.ndarray:
    normalized = np.full_like(spectrum, np.nan, dtype=np.float64)
    if kref >= spectrum.size or spectrum[kref] <= 0.0:
        return normalized
    normalized[:] = spectrum / spectrum[kref]
    return normalized


def _format_alpha(alpha: float) -> str:
    if np.isclose(alpha, 1.0 / 6.0):
        return "1/6"
    if np.isclose(alpha, 1.0 / 3.0):
        return "1/3"
    if np.isclose(alpha, 1.0 / 2.0):
        return "1/2"
    return f"{alpha:.6g}"


def _velocity_slope(alpha: float) -> float:
    return 2.0 * alpha + 1.0


def _run_case(athena: Path, alpha: float, nres: int, nlow: int, nhigh: int,
              seed: int, k_crit: float, force: bool, launcher: str) -> dict[str, object]:
    run_root = _run_root(athena)
    output_dir = run_root / "bin"
    output_dir.mkdir(parents=True, exist_ok=True)
    slug = _format_alpha(alpha).replace("/", "o").replace(".", "p")
    kcrit_slug = f"{k_crit:.3f}".rstrip("0").rstrip(".").replace(".", "p")
    basename = f"scalar_mix_clebsch_alpha_{slug}_n{nres}_kmax{nhigh}_kcrit{kcrit_slug}"
    diag_path = run_root / f"{basename}.turb_init_diag.json"
    timing_path = run_root / f"{basename}.velocity_init.json"
    hydro_pattern = output_dir / f"{basename}.hydro_w.*.bin"
    hydro_matches = sorted(glob.glob(str(hydro_pattern)))

    if hydro_matches and diag_path.exists() and not force:
        init_seconds = None
        if timing_path.exists():
            timing = json.loads(timing_path.read_text(encoding="utf-8"))
            if timing.get("init_seconds") is not None:
                init_seconds = float(timing["init_seconds"])
        hydro_path = Path(hydro_matches[-1])
    else:
        for path in hydro_matches:
            Path(path).unlink()
        if diag_path.exists():
            diag_path.unlink()
        if timing_path.exists():
            timing_path.unlink()

        cmd = shlex.split(launcher) + [
            str(athena),
            "-i",
            str(_repo_root() / "inputs" / "tests" / "scalar_mixing_clebsch_3d_alpha_sweep.athinput"),
            f"job/basename={basename}",
            "problem/turb_velocity_method=clebsch",
            f"problem/turb_alpha={alpha:.16f}",
            f"problem/turb_rseed={seed}",
            f"problem/turb_nlow={nlow}",
            f"problem/turb_nhigh={nhigh}",
            f"problem/turb_k_crit={k_crit}",
            "problem/turb_dump_generator_diagnostics=true",
            f"mesh/nx1={nres}",
            f"mesh/nx2={nres}",
            f"mesh/nx3={nres}",
            f"meshblock/nx1={nres // 2}",
            f"meshblock/nx2={nres // 2}",
            f"meshblock/nx3={nres // 2}",
            "time/tlim=1.0e-6",
            "output1/dt=1.0e-6",
        ]
        completed = subprocess.run(
            cmd,
            cwd=run_root,
            capture_output=True,
            text=True,
            check=True,
        )
        init_seconds = _parse_init_seconds(completed.stdout + completed.stderr)
        timing_path.write_text(
            json.dumps(
                {
                    "alpha": alpha,
                    "k_crit": k_crit,
                    "basename": basename,
                    "command": cmd,
                    "init_seconds": init_seconds,
                },
                indent=2,
            ) + "\n",
            encoding="utf-8",
        )
        hydro_matches = sorted(glob.glob(str(hydro_pattern)))
        if not hydro_matches:
            raise RuntimeError(f"Missing hydro output for alpha={alpha}")
        hydro_path = Path(hydro_matches[-1])

    data = read_binary_as_athdf(str(hydro_path), quantities=["velx", "vely", "velz"], dtype=np.float64)
    vx = np.asarray(data["velx"], dtype=np.float64)
    vy = np.asarray(data["vely"], dtype=np.float64)
    vz = np.asarray(data["velz"], dtype=np.float64)
    kvals, shell_energy = _shell_spectrum(vx, vy, vz)
    fitted_slope = _fit_slope(kvals, shell_energy, nlow, nhigh)
    shape = _shape_normalize(shell_energy, kref=max(3, nlow))
    return {
        "alpha": alpha,
        "label": _format_alpha(alpha),
        "basename": basename,
        "hydro_path": str(hydro_path),
        "diag_path": str(diag_path),
        "init_seconds": init_seconds,
        "k_crit": k_crit,
        "target_velocity_slope": _velocity_slope(alpha),
        "fitted_fft_shell_slope": fitted_slope,
        "kvals": kvals,
        "shell_energy": shell_energy,
        "shell_shape": shape,
    }


def _plot(results: list[dict[str, object]], output: Path, nlow: int, nhigh: int,
          nres: int, show_max: int | None) -> None:
    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    kref = max(3, nlow)
    kplot_lo = max(1, nlow)
    max_available = max(int(np.asarray(result["kvals"]).size - 1) for result in results)
    kplot_hi = min(max_available, nhigh if show_max is None else show_max)

    ax.axhline(1.0, color="k", lw=1.2, ls="--", label="target-compensated flat law")

    for result, color in zip(results, COLORS):
        kvals = np.asarray(result["kvals"], dtype=np.float64)
        shape = np.asarray(result["shell_shape"], dtype=np.float64)
        slope = float(result["target_velocity_slope"])
        mask = kvals > 0.0
        compensated = np.full_like(shape, np.nan, dtype=np.float64)
        compensated[mask] = shape[mask] * (kvals[mask] / float(kref)) ** slope
        stop = min(len(kvals) - 1, kplot_hi)
        ax.loglog(
            kvals[kplot_lo:stop + 1],
            compensated[kplot_lo:stop + 1],
            color=color,
            lw=2.4,
            marker="o",
            ms=3.0,
            label=(
                rf"$\alpha={result['label']}$, "
                rf"$s_v={slope:.3f}$, "
                rf"fit$={-float(result['fitted_fft_shell_slope']):.3f}$"
            ),
        )

    ax.axvline(nlow, color="0.6", lw=1.0, ls=":")
    ax.axvline(nhigh, color="0.6", lw=1.0, ls=":")
    ax.set_xlabel("k")
    ax.set_ylabel(r"$[E_v(k) / E_v(k_{\rm ref})]\,(k/k_{\rm ref})^{s_v}$")
    ax.set_title(
        rf"3D Clebsch Compensated Velocity Spectrum Sweep ($N={nres}^3$, $k_{{\min}}={nlow}$, $k_{{\max}}={nhigh}$)"
    )
    ax.set_xlim(float(nlow), float(kplot_hi))
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, loc="lower left")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, default=_default_athena_path())
    parser.add_argument(
        "--output",
        type=Path,
        default=_repo_root() / "docs" / "figures" / "scalar_mixing_velocity_methods"
        / "clebsch3d_alpha_sweep.png",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=_repo_root() / "docs" / "figures" / "scalar_mixing_velocity_methods"
        / "clebsch3d_alpha_sweep.json",
    )
    parser.add_argument("--nres", type=int, default=128)
    parser.add_argument("--nlow", type=int, default=2)
    parser.add_argument("--nhigh", type=int, default=32)
    parser.add_argument("--seed", type=int, default=24680)
    parser.add_argument("--k-crit", type=float, default=16.0)
    parser.add_argument("--show-max", type=int, default=None)
    parser.add_argument("--launcher", default="")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    athena = args.athena.resolve()
    if not athena.exists():
        raise FileNotFoundError(f"Athena executable not found: {athena}")

    results = []
    for alpha in ALPHAS:
        print(f"[run] alpha={_format_alpha(alpha)}")
        results.append(
            _run_case(
                athena=athena,
                alpha=alpha,
                nres=args.nres,
                nlow=args.nlow,
                nhigh=args.nhigh,
                seed=args.seed,
                k_crit=args.k_crit,
                force=args.force,
                launcher=args.launcher,
            )
        )

    _plot(results, args.output.resolve(), args.nlow, args.nhigh, args.nres, args.show_max)
    args.summary.resolve().write_text(
        json.dumps(
            {
                "nres": args.nres,
                "nlow": args.nlow,
                "nhigh": args.nhigh,
                "k_crit": args.k_crit,
                "show_max": args.show_max,
                "runs": [
                    {
                        "alpha": float(item["alpha"]),
                        "k_crit": float(item["k_crit"]),
                        "target_velocity_slope": float(item["target_velocity_slope"]),
                        "fitted_fft_shell_slope": float(item["fitted_fft_shell_slope"]),
                        "init_seconds": item["init_seconds"],
                        "hydro_path": item["hydro_path"],
                        "diag_path": item["diag_path"],
                    }
                    for item in results
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Saved figure to {args.output.resolve()}")
    print(f"Saved summary to {args.summary.resolve()}")


if __name__ == "__main__":
    main()
