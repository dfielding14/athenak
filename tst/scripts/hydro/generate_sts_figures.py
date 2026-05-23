"""Generate documentation figures for the RKL2 STS diffusion verification cases."""

from argparse import ArgumentParser
from pathlib import Path
import shutil
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "vis" / "python"))
import athena_read  # noqa: E402
import bin_convert  # noqa: E402


def run_case(executable, run_root, input_name, basename, *overrides):
    """Run one documentation case in an isolated output directory."""
    case_dir = run_root / basename
    case_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(executable),
        "-i",
        str(REPO_ROOT / "inputs" / "tests" / input_name),
        f"job/basename={basename}",
        *overrides,
    ]
    subprocess.run(command, cwd=case_dir, check=True, stdout=subprocess.DEVNULL)
    return case_dir


def uniform_scalar(binary_path):
    """Assemble a uniform two-dimensional scalar slice from binary MeshBlock data."""
    data = bin_convert.read_binary(str(binary_path))
    image = np.empty((data["Nx2"], data["Nx1"]))
    for block_id, block in enumerate(data["mb_data"]["s_00"]):
        logical_i, logical_j = data["mb_logical"][block_id][:2]
        nx1 = block.shape[-1]
        nx2 = block.shape[-2]
        i0 = int(logical_i) * nx1
        j0 = int(logical_j) * nx2
        image[j0:j0 + nx2, i0:i0 + nx1] = block[0, :, :]
    return image


def make_convergence_figure(executable, run_root, output_dir):
    """Measure scalar Fourier-mode errors and plot convergence with resolution."""
    resolutions = np.array([32, 64, 128])
    errors = []
    for resolution in resolutions:
        basename = f"figure_scalar_modes_{resolution}"
        case_dir = run_case(
            executable,
            run_root,
            "sts_scalar_modes.athinput",
            basename,
            f"mesh/nx1={resolution}",
            f"meshblock/nx1={min(resolution, 64)}",
        )
        errors.append(np.loadtxt(case_dir / f"{basename}-scalar-errors.dat")[2:4])
    errors = np.asarray(errors)

    fig, ax = plt.subplots(figsize=(5.7, 4.0), constrained_layout=True)
    ax.loglog(resolutions, errors[:, 0], "o-", label=r"$\kappa_0=0.005$")
    ax.loglog(resolutions, errors[:, 1], "s-", label=r"$\kappa_1=0.05$")
    reference = errors[0, 1] * (resolutions[0] / resolutions) ** 2
    ax.loglog(resolutions, reference, "k--", label=r"$\mathcal{O}(\Delta x^2)$")
    ax.set_xlabel("Number of cells")
    ax.set_ylabel(r"$L_1$ scalar error")
    ax.set_title("Independent scalar diffusivities with RKL2 STS")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(output_dir / "scalar_convergence.png", dpi=180)
    plt.close(fig)


def make_blob_figure(executable, run_root, output_dir):
    """Plot a two-dimensional advected and diffusing passive-scalar blob."""
    basename = "figure_scalar_blob"
    case_dir = run_case(
        executable,
        run_root,
        "sts_scalar_blob.athinput",
        basename,
    )
    initial = uniform_scalar(case_dir / "bin" / f"{basename}.hydro_w.00000.bin")
    final = uniform_scalar(case_dir / "bin" / f"{basename}.hydro_w.00001.bin")
    vmax = float(initial.max())

    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.55), constrained_layout=True)
    for ax, field, title in zip(
        axes,
        (initial, final),
        (r"$t=0$", r"$t=0.5$: advected and diffused"),
    ):
        image = ax.imshow(
            field,
            origin="lower",
            extent=(0.0, 1.0, 0.0, 1.0),
            cmap="magma",
            vmin=0.0,
            vmax=vmax,
        )
        ax.set_title(title)
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_aspect("equal")
    fig.colorbar(image, ax=axes, label="Scalar concentration", shrink=0.92)
    fig.savefig(output_dir / "scalar_blob_slice.png", dpi=180)
    plt.close(fig)


def make_thermal_front_figure(executable, run_root, output_dir):
    """Plot smoothing of a front under bounded power-law thermal conduction."""
    basename = "figure_thermal_front"
    case_dir = run_case(
        executable,
        run_root,
        "sts_thermal_front.athinput",
        basename,
    )
    initial = athena_read.tab(
        str(case_dir / "tab" / f"{basename}.hydro_w.00000.tab")
    )
    final = athena_read.tab(
        str(case_dir / "tab" / f"{basename}.hydro_w.00001.tab")
    )
    gm1 = 2.0 / 3.0
    initial_temp = gm1 * initial["eint"] / initial["dens"]
    final_temp = gm1 * final["eint"] / final["dens"]
    initial_order = np.argsort(initial["x1v"])
    final_order = np.argsort(final["x1v"])

    fig, ax = plt.subplots(figsize=(6.1, 3.75), constrained_layout=True)
    ax.plot(
        initial["x1v"][initial_order],
        initial_temp[initial_order],
        "--",
        label=r"$t=0$",
    )
    ax.plot(
        final["x1v"][final_order],
        final_temp[final_order],
        label=r"$t=0.04$",
    )
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$T$")
    ax.set_title(r"Bounded power-law conduction, $\kappa \propto T^{2.5}$")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(output_dir / "thermal_front_slice.png", dpi=180)
    plt.close(fig)


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--executable",
        type=Path,
        default=REPO_ROOT / "build-sts" / "src" / "athena",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "docs" / "source" / "_static" / "sts",
    )
    parser.add_argument(
        "--run-dir",
        type=Path,
        default=REPO_ROOT / "build-sts" / "sts-doc-runs",
    )
    args = parser.parse_args()
    executable = args.executable.resolve()
    if not executable.exists():
        parser.error(f"AthenaK executable does not exist: {executable}")

    output_dir = args.output_dir.resolve()
    run_root = args.run_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    shutil.rmtree(run_root, ignore_errors=True)
    run_root.mkdir(parents=True)

    make_convergence_figure(executable, run_root, output_dir)
    make_blob_figure(executable, run_root, output_dir)
    make_thermal_front_figure(executable, run_root, output_dir)

    for path in sorted(output_dir.glob("*.png")):
        print(path.relative_to(REPO_ROOT))


if __name__ == "__main__":
    main()
