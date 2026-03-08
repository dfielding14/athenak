#!/usr/bin/env python3
"""Generate a 5-panel gallery of scalar_mixing velocity initializations.

The figure contains:
1. 2D projection-method slice
2. 2D stream-function slice
3. 3D projection-method midplane slice
4. 3D stream-function midplane slice
5. Radial velocity power spectra with the requested k^(-5/3) reference law

The 3D examples use aggressive importance sampling (`turb_k_crit = 2`) so the
gallery remains practical on a serial local build while still exercising the
`k_min = 2`, `k_max = 64`, `256^3` setup.
"""

from __future__ import annotations

import argparse
import glob
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from bin_convert_new import read_binary_as_athdf


KPLOT_MIN = 1
KPLOT_MAX = 128
KREF = 3
KDRIVE_MIN = 2
KDRIVE_MAX = 64
EXPO = 5.0 / 3.0


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _default_athena_path() -> Path:
    return _repo_root() / "build-codex-spectrum" / "src" / "athena"


def _case_specs() -> list[dict[str, object]]:
    return [
        {
            "name": "2D Projection",
            "basename": "scalar_mix_gallery_proj_2d",
            "input": "inputs/hydro/scalar_mixing_gallery_2d.athinput",
            "overrides": [],
            "slice_kind": "2d",
            "line_color": "#0b6e4f",
        },
        {
            "name": "2D Stream",
            "basename": "scalar_mix_gallery_stream_2d",
            "input": "inputs/hydro/scalar_mixing_gallery_2d.athinput",
            "overrides": ["problem/turb_use_stream_function=true"],
            "slice_kind": "2d",
            "line_color": "#bd632f",
        },
        {
            "name": "3D Projection",
            "basename": "scalar_mix_gallery_proj_3d",
            "input": "inputs/hydro/scalar_mixing_gallery_3d.athinput",
            "overrides": [],
            "slice_kind": "3d",
            "line_color": "#275dad",
        },
        {
            "name": "3D Stream",
            "basename": "scalar_mix_gallery_stream_3d",
            "input": "inputs/hydro/scalar_mixing_gallery_3d.athinput",
            "overrides": ["problem/turb_use_stream_function=true"],
            "slice_kind": "3d",
            "line_color": "#8b3fb3",
        },
    ]


def _run_case(athena: Path, case: dict[str, object], force: bool) -> Path:
    run_dir = athena.parent
    output_dir = run_dir / "bin"
    pattern = output_dir / f"{case['basename']}.hydro_w.*.bin"
    matches = sorted(glob.glob(str(pattern)))
    if matches and not force:
        return Path(matches[-1])

    for path in matches:
        Path(path).unlink()

    cmd = [
        str(athena),
        "-i",
        str(_repo_root() / str(case["input"])),
        f"job/basename={case['basename']}",
    ] + list(case["overrides"])
    subprocess.run(cmd, cwd=run_dir, check=True)

    matches = sorted(glob.glob(str(pattern)))
    if not matches:
        raise RuntimeError(f"No hydro_w output found for {case['name']}")
    return Path(matches[-1])


def _load_case(path: Path) -> dict[str, np.ndarray]:
    data = read_binary_as_athdf(str(path), quantities=["velx", "vely", "velz"],
                                dtype=np.float64)
    return {
        "velx": np.asarray(data["velx"], dtype=np.float64),
        "vely": np.asarray(data["vely"], dtype=np.float64),
        "velz": np.asarray(data["velz"], dtype=np.float64),
    }


def _midplane_slice(field: np.ndarray) -> np.ndarray:
    if field.shape[0] == 1:
        return field[0, :, :]
    return field[field.shape[0] // 2, :, :]


def _slice_payload(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    vx = _midplane_slice(data["velx"])
    vy = _midplane_slice(data["vely"])
    vz = _midplane_slice(data["velz"])
    speed = np.sqrt(vx*vx + vy*vy + vz*vz)
    return {"vx": vx, "vy": vy, "vz": vz, "speed": speed}


def _radial_spectrum(data: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    vx = data["velx"]
    vy = data["vely"]
    vz = data["velz"]

    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)
    energy = np.abs(uhat_x)**2 + np.abs(uhat_y)**2 + np.abs(uhat_z)**2

    nz, ny, nx = vx.shape
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny if ny > 1 else np.array([0.0])
    kz = np.fft.fftfreq(nz) * nz if nz > 1 else np.array([0.0])
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")
    shell = np.ceil(np.sqrt(kx3*kx3 + ky3*ky3 + kz3*kz3) - 1.0e-12).astype(np.int64)

    spectrum = np.bincount(shell.ravel(), weights=energy.ravel(), minlength=KPLOT_MAX + 2)
    kvals = np.arange(KPLOT_MIN, KPLOT_MAX + 1, dtype=np.float64)
    return kvals, spectrum[KPLOT_MIN:KPLOT_MAX + 1]


def _divergence_l2(data: dict[str, np.ndarray]) -> float:
    vx = data["velx"]
    vy = data["vely"]
    vz = data["velz"]

    uhat_x = np.fft.fftn(vx)
    uhat_y = np.fft.fftn(vy)
    uhat_z = np.fft.fftn(vz)

    nz, ny, nx = vx.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=1.0 / nx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=1.0 / ny) if ny > 1 else np.array([0.0])
    kz = 2.0 * np.pi * np.fft.fftfreq(nz, d=1.0 / nz) if nz > 1 else np.array([0.0])
    kz3, ky3, kx3 = np.meshgrid(kz, ky, kx, indexing="ij")

    div_hat = 1j * (kx3 * uhat_x + ky3 * uhat_y + kz3 * uhat_z)
    div = np.fft.ifftn(div_hat).real
    return float(np.sqrt(np.sum(div*div)))


def _plot_slice_panel(ax: plt.Axes, title: str, payload: dict[str, np.ndarray],
                      vmax: float, div_l2: float) -> None:
    image = ax.imshow(payload["speed"], origin="lower", cmap="magma",
                      extent=[0.0, 1.0, 0.0, 1.0], vmin=0.0, vmax=vmax)
    skip = max(1, payload["vx"].shape[0] // 16)
    xs = np.linspace(0.0, 1.0, payload["vx"].shape[1], endpoint=False) + 0.5/payload["vx"].shape[1]
    ys = np.linspace(0.0, 1.0, payload["vx"].shape[0], endpoint=False) + 0.5/payload["vx"].shape[0]
    grid_x, grid_y = np.meshgrid(xs[::skip], ys[::skip])
    ax.quiver(grid_x, grid_y, payload["vx"][::skip, ::skip], payload["vy"][::skip, ::skip],
              color="white", pivot="mid", scale=22.0, width=0.003, alpha=0.9)
    ax.set_title(f"{title}\nsqrt(sum div(v)^2) = {div_l2:.2e}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax._gallery_image = image  # type: ignore[attr-defined]


def _make_figure(cases: list[dict[str, object]], data_by_case: dict[str, dict[str, np.ndarray]],
                 spectra: dict[str, tuple[np.ndarray, np.ndarray]],
                 divergence_l2: dict[str, float], output: Path) -> None:
    payloads = {case["name"]: _slice_payload(data_by_case[str(case["name"])]) for case in cases}
    vmax = max(float(np.max(payload["speed"])) for payload in payloads.values())

    fig = plt.figure(figsize=(16, 10))
    grid = GridSpec(2, 3, figure=fig, width_ratios=[1.0, 1.0, 1.2], wspace=0.36, hspace=0.24)

    slice_axes = [
        fig.add_subplot(grid[0, 0]),
        fig.add_subplot(grid[0, 1]),
        fig.add_subplot(grid[1, 0]),
        fig.add_subplot(grid[1, 1]),
    ]
    spec_ax = fig.add_subplot(grid[:, 2])

    for ax, case in zip(slice_axes, cases):
        _plot_slice_panel(ax, str(case["name"]), payloads[str(case["name"])], vmax,
                          divergence_l2[str(case["name"])])

    colorbar = fig.colorbar(slice_axes[0]._gallery_image, ax=slice_axes, fraction=0.02, pad=0.02)
    colorbar.set_label(r"$|\mathbf{v}|$")
    colorbar.ax.yaxis.set_label_position("right")
    colorbar.ax.yaxis.tick_right()

    reference_k = np.arange(KDRIVE_MIN, KDRIVE_MAX + 1, dtype=np.float64)
    reference = (reference_k ** (-EXPO)) / (KREF ** (-EXPO))
    spec_ax.loglog(reference_k, reference, color="black", linestyle="--", linewidth=2.0,
                   label=r"reference $k^{-5/3}$")

    for case in cases:
        kvals, spectrum = spectra[str(case["name"])]
        ref_index = KREF - KPLOT_MIN
        if spectrum[ref_index] <= 0.0:
            raise RuntimeError(f"{case['name']} has non-positive E(k={KREF}); cannot normalize.")
        normalized = spectrum / spectrum[ref_index]
        spec_ax.loglog(kvals, normalized, linewidth=2.0, color=str(case["line_color"]),
                       label=str(case["name"]))

    spec_ax.axvline(KDRIVE_MIN, color="0.6", linestyle=":", linewidth=1.0)
    spec_ax.axvline(KDRIVE_MAX, color="0.6", linestyle=":", linewidth=1.0)
    spec_ax.set_xlabel("k")
    spec_ax.set_ylabel(r"$E(k) / E(3)$")
    spec_ax.set_title("Velocity Power Spectra")
    spec_ax.grid(True, which="both", alpha=0.25)
    spec_ax.legend(loc="lower left", frameon=False)
    spec_ax.set_xlim(KPLOT_MIN, KPLOT_MAX)
    spec_ax.set_ylim(bottom=1.0e-4)

    fig.suptitle("Scalar Mixing Velocity Initialization Gallery\n"
                 r"$256^2$ / $256^3$, $k_{\min}=2$, $k_{\max}=64$, target slope $-5/3$",
                 fontsize=16)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--athena", type=Path, default=_default_athena_path(),
                        help="Path to the AthenaK executable built with PROBLEM=scalar_mixing.")
    parser.add_argument("--output", type=Path,
                        default=_repo_root() / "docs" / "scalar_mixing_velocity_gallery.png",
                        help="Output figure path.")
    parser.add_argument("--force", action="store_true",
                        help="Re-run Athena even if matching hydro_w outputs already exist.")
    args = parser.parse_args()

    cases = _case_specs()
    data_by_case: dict[str, dict[str, np.ndarray]] = {}
    spectra: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    divergence_l2: dict[str, float] = {}

    for case in cases:
        output_path = _run_case(args.athena, case, args.force)
        case_data = _load_case(output_path)
        data_by_case[str(case["name"])] = case_data
        spectra[str(case["name"])] = _radial_spectrum(case_data)
        divergence_l2[str(case["name"])] = _divergence_l2(case_data)

    _make_figure(cases, data_by_case, spectra, divergence_l2, args.output)
    print(f"Saved gallery figure to {args.output}")


if __name__ == "__main__":
    main()
