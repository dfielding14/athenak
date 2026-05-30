#!/usr/bin/env python3
"""Generate self-gravity validation figures and the docs results page.

Run from the repository root:

  python docs/scripts/generate_selfgravity_validation.py \
      --build-cpu --with-mpi --build-mpi

The script intentionally reuses the regression-test input files and binary
output parser so the documentation figures track the test suite.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.ticker as mticker  # noqa: E402
import numpy as np  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
TST_DIR = REPO_ROOT / "tst"
BUILD_SRC = TST_DIR / "build" / "src"
STATIC_DIR = REPO_ROOT / "docs" / "source" / "_static" / "selfgravity"
RESULTS_PAGE = REPO_ROOT / "docs" / "source" / "modules" / "self_gravity_validation.md"

sys.path.insert(0, str(TST_DIR))
from test_suite.selfgravity import selfgravity_utils as sg  # noqa: E402


PROFILE_RE = re.compile(r"\[gravity-profile\]\s+(.*)")


def run_command(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command))
    result = subprocess.run(command, cwd=cwd, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"Command failed with exit code {result.returncode}: {' '.join(command)}"
        )


def build(extra_flags: list[str] | None = None) -> None:
    flags = extra_flags or []
    run_command(
        ["cmake", "-S", str(REPO_ROOT), "-B", str(TST_DIR / "build")] + flags, REPO_ROOT
    )
    run_command(
        ["cmake", "--build", str(TST_DIR / "build"), "-j", str(os.cpu_count() or 4)],
        REPO_ROOT,
    )


def ensure_build() -> None:
    athena = BUILD_SRC / "athena"
    if not athena.exists():
        raise FileNotFoundError(
            f"{athena} does not exist. Re-run with --build-cpu or build AthenaK first."
        )


def in_build_dir(func):
    def wrapper(*args, **kwargs):
        ensure_build()
        old_cwd = Path.cwd()
        os.chdir(BUILD_SRC)
        try:
            return func(*args, **kwargs)
        finally:
            os.chdir(old_cwd)

    return wrapper


@in_build_dir
def run_athena(*args, **kwargs):
    return sg.run_athena(*args, **kwargs)


@in_build_dir
def run_restart(*args, **kwargs):
    return sg.run_restart(*args, **kwargs)


@in_build_dir
def parse_latest_binary(basename: str, output_id: str):
    return sg.parse_binary_output(sg.latest_binary_output(basename, output_id))


@in_build_dir
def parse_first_binary(basename: str, output_id: str):
    return sg.parse_binary_output(sg.first_binary_output(basename, output_id))


@in_build_dir
def parse_latest_restart(basename: str):
    return sg.latest_restart_output(basename)


@in_build_dir
def cleanup(*basenames: str) -> None:
    sg.cleanup_outputs(*basenames)


@in_build_dir
def expect_failure(
    input_file: str | Path, extra_args: list[str] | None, expected: str
) -> str:
    return sg.expect_failure(input_file, extra_args=extra_args, expected_text=expected)


def centers_to_edges(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values)
    if values.size == 1:
        return np.array([values[0] - 0.5, values[0] + 0.5])
    edges = np.empty(values.size + 1)
    edges[1:-1] = 0.5 * (values[1:] + values[:-1])
    edges[0] = values[0] - 0.5 * (values[1] - values[0])
    edges[-1] = values[-1] + 0.5 * (values[-1] - values[-2])
    return edges


def add_image(ax, data, x, y, title, cmap="viridis", symmetric=False):
    finite = np.asarray(data)[np.isfinite(data)]
    vmin = vmax = None
    if symmetric and finite.size:
        vmax = np.max(np.abs(finite))
        vmin = -vmax
    mesh = ax.pcolormesh(
        centers_to_edges(x),
        centers_to_edges(y),
        data,
        shading="auto",
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    plt.colorbar(mesh, ax=ax, fraction=0.046, pad=0.04)


def save(fig, name: str) -> str:
    STATIC_DIR.mkdir(parents=True, exist_ok=True)
    path = STATIC_DIR / name
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return f"../_static/selfgravity/{name}"


def jeans_analytic(x1, x2, x3, amp: float) -> tuple[np.ndarray, float]:
    k_wave, four_pi_g, _ = sg.jeans_geometry(n_jeans=0.5)
    analytic = -four_pi_g * amp / (k_wave * k_wave)
    analytic *= np.sin(k_wave * sg.jeans_phase(x1, x2, x3))
    return analytic, k_wave


def parse_profile(output: str) -> dict[str, float]:
    for line in output.splitlines():
        match = PROFILE_RE.search(line)
        if not match:
            continue
        values: dict[str, float] = {}
        for token in match.group(1).split():
            if "=" not in token:
                continue
            key, raw = token.split("=", 1)
            if key == "stage":
                continue
            try:
                values[key] = float(raw)
            except ValueError:
                pass
        return values
    return {}


def uniform_sphere_potential(x, y, z, center, mass, radius, gconst):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    r = np.sqrt(dx * dx + dy * dy + dz * dz)
    outside = -gconst * mass / np.maximum(r, 1.0e-300)
    inside = -gconst * mass * (3.0 * radius * radius - r * r) / (2.0 * radius**3)
    return np.where(r < radius, inside, outside)


def point_mass_acceleration(x, y, z, center, mass, gconst):
    dx = x - center[0]
    dy = y - center[1]
    dz = z - center[2]
    r2 = dx * dx + dy * dy + dz * dz
    r3 = np.maximum(r2, 1.0e-300) ** 1.5
    return (-gconst * mass * dx / r3, -gconst * mass * dy / r3, -gconst * mass * dz / r3)


def periodic_jeans(metrics: dict[str, Any]) -> str:
    basename = "docs_selfgravity_periodic"
    try:
        run_athena(
            "inputs/tests/selfgravity.athinput", basename, expect_jeans_banner=True
        )
        output = parse_latest_binary(basename, "grav_phi")
        phi_grid, x, y, z = sg.uniform_grid_field(output, "grav_phi")
        zgrid, ygrid, xgrid = np.meshgrid(z, y, x, indexing="ij")
        analytic, _ = jeans_analytic(xgrid, ygrid, zgrid, 1.0e-3)
        numerical = phi_grid - np.mean(phi_grid)
        analytic -= np.mean(analytic)
        error = numerical - analytic
        rel_l2 = sg.relative_l2(error, analytic)
        initial = parse_first_binary(basename, "grav_phi")
        initial_phi, initial_x1, initial_x2, initial_x3 = sg.concatenate_field(
            initial, "grav_phi"
        )
        initial_analytic, _ = jeans_analytic(initial_x1, initial_x2, initial_x3, 1.0e-3)
        initial_phi -= np.mean(initial_phi)
        initial_analytic -= np.mean(initial_analytic)
        metrics["periodic_jeans"] = {
            "relative_l2": rel_l2,
            "initial_output_relative_l2": sg.relative_l2(
                initial_phi - initial_analytic, initial_analytic
            ),
        }

        mid = len(z) // 2
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        add_image(axes[0, 0], numerical[mid], x, y, "Numerical potential", "viridis")
        add_image(axes[0, 1], analytic[mid], x, y, "Analytic potential", "viridis")
        add_image(axes[1, 0], error[mid], x, y, "Potential error", "RdBu_r", True)
        axes[1, 1].scatter(analytic.ravel(), numerical.ravel(), s=2, alpha=0.35)
        lim = np.max(np.abs(analytic))
        axes[1, 1].plot([-lim, lim], [-lim, lim], color="black", linewidth=1)
        axes[1, 1].set_title(f"Pointwise comparison; rel. L2 = {rel_l2:.3e}")
        axes[1, 1].set_xlabel("analytic phi")
        axes[1, 1].set_ylabel("numerical phi")
        return save(fig, "periodic_jeans_potential.png")
    finally:
        cleanup(basename)


def hydro_mhd_equivalence(metrics: dict[str, Any]) -> str:
    hydro = "docs_selfgravity_equiv_hydro"
    mhd = "docs_selfgravity_equiv_mhd"
    try:
        run_athena("inputs/tests/selfgravity.athinput", hydro, expect_jeans_banner=True)
        run_athena("inputs/tests/selfgravity_mhd.athinput", mhd, expect_jeans_banner=True)
        h_grid, x, y, z = sg.uniform_grid_field(
            parse_latest_binary(hydro, "grav_phi"), "grav_phi"
        )
        m_grid, _, _, _ = sg.uniform_grid_field(
            parse_latest_binary(mhd, "grav_phi"), "grav_phi"
        )
        diff = h_grid - m_grid
        metrics["hydro_mhd_equivalence"] = {
            "max_abs_phi_difference": float(np.max(np.abs(diff)))
        }

        mid = len(z) // 2
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
        add_image(axes[0], h_grid[mid], x, y, "Hydro phi", "viridis")
        add_image(axes[1], m_grid[mid], x, y, "MHD zero-field phi", "viridis")
        add_image(axes[2], diff[mid], x, y, "Hydro - MHD", "RdBu_r", True)
        return save(fig, "hydro_mhd_equivalence.png")
    finally:
        cleanup(hydro, mhd)


def source_oracle(metrics: dict[str, Any]) -> str:
    on = "docs_selfgravity_source_on"
    off = "docs_selfgravity_source_off"
    try:
        run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            on,
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            off,
            extra_args=["hydro_srcterms/self_gravity=false"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        on_out = parse_latest_binary(on, "hydro_u")
        off_out = parse_latest_binary(off, "hydro_u")
        dens, x, y, z = sg.uniform_grid_field(on_out, "dens")
        zgrid, ygrid, xgrid = np.meshgrid(z, y, x, indexing="ij")
        k_wave, four_pi_g, direction = sg.jeans_geometry(n_jeans=0.5)
        accel = (four_pi_g * 1.0e-4 / k_wave) * np.cos(
            k_wave * sg.jeans_phase(xgrid, ygrid, zgrid)
        )
        time = on_out["metadata"]["time"]
        expected = [dens * accel * direction[n] * time for n in range(3)]
        rel = {}
        actual_fields = []
        for label, ref in zip(("mom1", "mom2", "mom3"), expected):
            on_mom, _, _, _ = sg.uniform_grid_field(on_out, label)
            off_mom, _, _, _ = sg.uniform_grid_field(off_out, label)
            actual = on_mom - off_mom
            actual_fields.append(actual)
            rel[label] = sg.relative_l2(actual - ref, ref)
        phi, _, _, _ = sg.uniform_grid_field(
            parse_latest_binary(on, "grav_phi"), "grav_phi"
        )
        dx = (x[1] - x[0], y[1] - y[0], z[1] - z[0])
        laplacian = (
            (np.roll(phi, -1, axis=2) - 2.0 * phi + np.roll(phi, 1, axis=2))
            / (dx[0] * dx[0])
            + (np.roll(phi, -1, axis=1) - 2.0 * phi + np.roll(phi, 1, axis=1))
            / (dx[1] * dx[1])
            + (np.roll(phi, -1, axis=0) - 2.0 * phi + np.roll(phi, 1, axis=0))
            / (dx[2] * dx[2])
        )
        poisson_source = four_pi_g * (dens - np.mean(dens))
        poisson_residual = laplacian - poisson_source
        rel["final_poisson_relative_l2"] = sg.relative_l2(
            poisson_residual, poisson_source
        )
        metrics["source_oracle"] = rel

        mid = len(z) // 2
        residual = actual_fields[0] - expected[0]
        fig, axes = plt.subplots(1, 4, figsize=(16, 3.8))
        add_image(
            axes[0], actual_fields[0][mid], x, y, "Actual mom1 increment", "viridis"
        )
        add_image(axes[1], expected[0][mid], x, y, "Analytic mom1 increment", "viridis")
        add_image(
            axes[2],
            residual[mid],
            x,
            y,
            f"Residual; rel. L2 = {rel['mom1']:.3e}",
            "RdBu_r",
            True,
        )
        add_image(
            axes[3],
            poisson_residual[mid],
            x,
            y,
            f"Final Poisson residual; rel. L2 = {rel['final_poisson_relative_l2']:.3e}",
            "RdBu_r",
            True,
        )
        return save(fig, "source_term_oracle.png")
    finally:
        cleanup(on, off)


def binary_multipole(metrics: dict[str, Any]) -> str:
    basename = "docs_selfgravity_binary"
    try:
        run_athena(
            "inputs/tests/selfgravity_binary.athinput", basename, max_final_defect=2.0e-6
        )
        output = parse_latest_binary(basename, "grav_phi")
        phi, x, y, z = sg.uniform_grid_field(output, "grav_phi")
        zgrid, ygrid, xgrid = np.meshgrid(z, y, x, indexing="ij")
        gconst = 1.0 / (4.0 * math.pi)
        radius = 0.10
        center1 = (0.15, 0.0, 0.0)
        center2 = (-0.15, 0.0, 0.0)
        analytic = uniform_sphere_potential(
            xgrid, ygrid, zgrid, center1, 1.0, radius, gconst
        ) + uniform_sphere_potential(xgrid, ygrid, zgrid, center2, 0.5, radius, gconst)
        r1 = np.sqrt((xgrid - center1[0]) ** 2 + ygrid * ygrid + zgrid * zgrid)
        r2 = np.sqrt((xgrid - center2[0]) ** 2 + ygrid * ygrid + zgrid * zgrid)
        far = (r1 > 1.5 * radius) & (r2 > 1.5 * radius)
        rel_l2 = sg.relative_l2(phi[far] - analytic[far], analytic[far])

        dphidz, dphidy, dphidx = np.gradient(phi, z, y, x, edge_order=2)
        ax1a, ax2a, ax3a = point_mass_acceleration(
            xgrid, ygrid, zgrid, center1, 1.0, gconst
        )
        bx1a, bx2a, bx3a = point_mass_acceleration(
            xgrid, ygrid, zgrid, center2, 0.5, gconst
        )
        ref = np.stack(((ax1a + bx1a), (ax2a + bx2a), (ax3a + bx3a)))
        num = np.stack((-dphidx, -dphidy, -dphidz))
        mask = (r1 > 2.0 * radius) & (r2 > 2.0 * radius)
        mask[:, :, :1] = False
        mask[:, :, -1:] = False
        mask[:, :1, :] = False
        mask[:, -1:, :] = False
        mask[:1, :, :] = False
        mask[-1:, :, :] = False
        force_cosine = float(
            np.sum(num[:, mask] * ref[:, mask])
            / (np.linalg.norm(num[:, mask]) * np.linalg.norm(ref[:, mask]))
        )
        metrics["binary_multipole"] = {
            "potential_relative_l2": rel_l2,
            "force_cosine": force_cosine,
        }

        mid = len(z) // 2
        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        add_image(axes[0, 0], phi[mid], x, y, "Numerical phi", "viridis")
        add_image(
            axes[0, 1], analytic[mid], x, y, "Uniform-sphere analytic phi", "viridis"
        )
        add_image(
            axes[1, 0], (phi - analytic)[mid], x, y, "Potential residual", "RdBu_r", True
        )
        dist = np.minimum(r1[far], r2[far])
        rel_err = np.abs(phi[far] - analytic[far]) / np.maximum(
            np.abs(analytic[far]), 1.0e-300
        )
        axes[1, 1].scatter(dist, rel_err, s=3, alpha=0.25)
        axes[1, 1].set_yscale("log")
        axes[1, 1].set_xlabel("distance to nearest sphere center")
        axes[1, 1].set_ylabel("|phi - phi_ref| / |phi_ref|")
        axes[1, 1].set_title(f"Potential error; force cosine = {force_cosine:.3f}")
        return save(fig, "binary_multipole_validation.png")
    finally:
        cleanup(basename)


def restart_equivalence(metrics: dict[str, Any]) -> str:
    straight = "docs_selfgravity_restart_straight"
    split = "docs_selfgravity_restart_split"
    try:
        run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            straight,
            extra_args=["time/nlim=2"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        run_athena(
            "inputs/tests/selfgravity_oracle.athinput",
            split,
            extra_args=["time/nlim=1"],
            expect_jeans_banner=True,
            max_final_defect=2.0e-8,
        )
        restart = parse_latest_restart(split)
        run_restart(restart, extra_args=["time/nlim=2"])
        straight_out = parse_latest_binary(straight, "hydro_u")
        split_out = parse_latest_binary(split, "hydro_u")
        diffs = {}
        diff_grids = {}
        for label in ("dens", "mom1", "mom2", "mom3"):
            a, x, y, z = sg.uniform_grid_field(straight_out, label)
            b, _, _, _ = sg.uniform_grid_field(split_out, label)
            diff_grids[label] = a - b
            diffs[label] = float(np.max(np.abs(diff_grids[label])))
        metrics["restart_equivalence"] = diffs

        mid = len(z) // 2
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
        add_image(
            axes[0], diff_grids["dens"][mid], x, y, "density difference", "RdBu_r", True
        )
        add_image(
            axes[1], diff_grids["mom1"][mid], x, y, "mom1 difference", "RdBu_r", True
        )
        axes[2].hist(diff_grids["dens"].ravel(), bins=40, alpha=0.65, label="dens")
        axes[2].hist(diff_grids["mom1"].ravel(), bins=40, alpha=0.65, label="mom1")
        axes[2].set_title("Straight run minus restart run")
        axes[2].set_xlabel("difference")
        axes[2].set_ylabel("cell count")
        axes[2].legend()
        return save(fig, "restart_equivalence.png")
    finally:
        cleanup(straight, split)


def controls_profile(metrics: dict[str, Any]) -> str:
    threshold = "docs_selfgravity_profile_threshold"
    fixed = "docs_selfgravity_profile_fixed"
    root_device = "docs_selfgravity_root_device"
    root_host = "docs_selfgravity_root_host"
    try:
        threshold_out = run_athena(
            "inputs/tests/selfgravity.athinput",
            threshold,
            extra_args=["gravity/profile=true"],
            expect_jeans_banner=True,
        )
        fixed_out = run_athena(
            "inputs/tests/selfgravity.athinput",
            fixed,
            extra_args=[
                "gravity/profile=true",
                "gravity/threshold=-1",
                "gravity/niteration=4",
            ],
            expect_jeans_banner=True,
            max_final_defect=None,
        )
        run_athena(
            "inputs/tests/selfgravity.athinput", root_device, expect_jeans_banner=True
        )
        run_athena(
            "inputs/tests/selfgravity.athinput",
            root_host,
            extra_args=["gravity/root_on_host=true"],
            expect_jeans_banner=True,
        )
        device_phi, _, _, _ = sg.uniform_grid_field(
            parse_latest_binary(root_device, "grav_phi"), "grav_phi"
        )
        host_phi, _, _, _ = sg.uniform_grid_field(
            parse_latest_binary(root_host, "grav_phi"), "grav_phi"
        )
        root_diff = float(np.max(np.abs(device_phi - host_phi)))
        threshold_defects = [float(v) for v in sg.DEFECT_RE.findall(threshold_out)]
        fixed_defects = [float(v) for v in sg.DEFECT_RE.findall(fixed_out)]
        profile = parse_profile(threshold_out)
        metrics["controls_profile"] = {
            "threshold_final_defect": threshold_defects[-1],
            "fixed_iteration_final_defect": fixed_defects[-1],
            "root_on_host_max_abs_phi_difference": root_diff,
            "profile_seconds": profile,
        }

        fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
        axes[0].plot(
            range(1, len(threshold_defects) + 1),
            threshold_defects,
            marker="o",
            label="threshold",
        )
        axes[0].plot(
            range(1, len(fixed_defects) + 1),
            fixed_defects,
            marker="o",
            label="fixed iteration",
        )
        axes[0].set_yscale("log")
        axes[0].set_xlabel("reported solve")
        axes[0].set_ylabel("final defect")
        axes[0].set_title("Convergence control")
        axes[0].legend()
        keys = [
            k
            for k in (
                "source_load",
                "setup",
                "root_transfer",
                "smooth",
                "boundary",
                "restrict_prolong",
                "result_retrieve",
            )
            if k in profile
        ]
        axes[1].barh(keys, [profile[k] for k in keys])
        axes[1].set_xlabel("seconds")
        axes[1].set_title("Profile phases")
        axes[1].xaxis.set_major_locator(mticker.MaxNLocator(4))
        axes[1].ticklabel_format(axis="x", style="sci", scilimits=(-2, 2))
        axes[2].bar(["root host/device"], [root_diff])
        axes[2].set_ylabel("max |delta phi|")
        axes[2].set_title("Root placement comparison")
        if root_diff > 0.0:
            axes[2].set_yscale("log")
        else:
            axes[2].set_ylim(0.0, 1.0)
            axes[2].text(0, 0.5, "exact match", ha="center", va="center")
        return save(fig, "controls_profile_root.png")
    finally:
        cleanup(threshold, fixed, root_device, root_host)


def be_static_refined(metrics: dict[str, Any]) -> str:
    basename = "docs_selfgravity_be_smr"
    try:
        out = run_athena(
            "inputs/tests/selfgravity_be_amr.athinput", basename, max_final_defect=2.0e-5
        )
        levels_match = re.search(r"Number of physical levels of refinement = (\d+)", out)
        levels = int(levels_match.group(1)) if levels_match else -1
        phi, x, y, z = sg.uniform_grid_field(
            parse_latest_binary(basename, "grav_phi"), "grav_phi"
        )
        dens, _, _, _ = sg.uniform_grid_field(
            parse_latest_binary(basename, "hydro_u"), "dens"
        )
        metrics["be_static_refined"] = {
            "physical_refinement_levels": levels,
            "phi_min": float(np.min(phi)),
            "phi_max": float(np.max(phi)),
            "density_min": float(np.min(dens)),
            "density_max": float(np.max(dens)),
        }

        mid = len(z) // 2
        fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
        add_image(axes[0], dens[mid], x, y, "BE density", "magma")
        add_image(axes[1], phi[mid], x, y, "BE potential", "viridis")
        axes[2].hist(dens.ravel(), bins=45, histtype="step", linewidth=1.5)
        axes[2].set_yscale("log")
        axes[2].set_xlabel("density")
        axes[2].set_ylabel("cell count")
        axes[2].set_title(f"Static refined hierarchy; levels = {levels}")
        return save(fig, "be_static_refined.png")
    finally:
        cleanup(basename)


def dynamic_amr_regrid(metrics: dict[str, Any]) -> str:
    basename = "docs_selfgravity_dynamic_amr"
    source = (REPO_ROOT / "inputs/tests/selfgravity.athinput").read_text()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "selfgravity_dynamic_amr.athinput"
        path.write_text(
            source.rstrip()
            + """

<mesh_refinement>
refinement = adaptive
max_nmb_per_rank = 128
num_levels = 2
ncycle_check = 1
refinement_interval = 1

<amr_criterion0>
method = min_max
variable = hydro_w_d
value_max = 0.5

<output2>
file_type = bin
variable = hydro_u
dt = 1.0
"""
        )
        try:
            out = run_athena(
                path,
                basename,
                extra_args=["mesh/nghost=4", "time/nlim=2", "gravity/threshold=1e-6"],
                max_final_defect=2.0e-6,
                expect_jeans_banner=True,
            )
            block_counts = re.findall(r"Current number of MeshBlocks = (\d+)", out)
            amr_counts = re.findall(
                r"(\d+) MeshBlocks created, (\d+) deleted by AMR", out
            )
            defects = [float(v) for v in sg.DEFECT_RE.findall(out)]
            phi, x, y, z = sg.uniform_grid_field(
                parse_latest_binary(basename, "grav_phi"), "grav_phi"
            )
            dens, _, _, _ = sg.uniform_grid_field(
                parse_latest_binary(basename, "hydro_u"), "dens"
            )
            metrics["dynamic_amr_regrid"] = {
                "final_meshblocks": int(block_counts[-1]) if block_counts else -1,
                "created_meshblocks": int(amr_counts[-1][0]) if amr_counts else -1,
                "deleted_meshblocks": int(amr_counts[-1][1]) if amr_counts else -1,
                "final_defect": defects[-1],
            }

            mid = len(z) // 2
            fig, axes = plt.subplots(1, 3, figsize=(12, 3.8))
            add_image(axes[0], dens[mid], x, y, "Dynamic AMR density", "magma")
            add_image(axes[1], phi[mid], x, y, "Dynamic AMR potential", "viridis")
            axes[2].hist(phi.ravel(), bins=45, histtype="step", linewidth=1.5)
            axes[2].set_xlabel("potential")
            axes[2].set_ylabel("cell count")
            final_meshblocks = metrics["dynamic_amr_regrid"]["final_meshblocks"]
            axes[2].set_title(f"Regrid smoke; final blocks = {final_meshblocks}")
            return save(fig, "dynamic_amr_regrid.png")
        finally:
            cleanup(basename)


def remove_block(text: str, block_name: str) -> str:
    lines = text.splitlines()
    out = []
    skipping = False
    for line in lines:
        stripped = line.strip()
        if stripped == f"<{block_name}>":
            skipping = True
            continue
        if skipping and stripped.startswith("<") and stripped.endswith(">"):
            skipping = False
        if not skipping:
            out.append(line)
    return "\n".join(out) + "\n"


def remove_parameters(text: str, block_name: str, parameters: set[str]) -> str:
    lines = text.splitlines()
    out = []
    in_block = False
    for line in lines:
        stripped = line.strip()
        if stripped == f"<{block_name}>":
            in_block = True
            out.append(line)
            continue
        if in_block and stripped.startswith("<") and stripped.endswith(">"):
            in_block = False
        key = stripped.split("=", 1)[0].strip()
        if in_block and key in parameters:
            continue
        out.append(line)
    return "\n".join(out) + "\n"


def failure_cases(metrics: dict[str, Any]) -> None:
    failures = [
        (
            "non-cubic MeshBlock",
            "inputs/tests/selfgravity.athinput",
            ["meshblock/nx3=8"],
            "logically cubic MeshBlocks",
        ),
        (
            "2D mesh",
            "inputs/tests/selfgravity.athinput",
            ["mesh/nx3=1", "meshblock/nx3=1"],
            "requires a 3D mesh",
        ),
        (
            "bad mg_bc",
            "inputs/tests/selfgravity.athinput",
            ["gravity/mg_bc=bogus"],
            "Unknown gravity/mg_bc",
        ),
        (
            "bad mporder",
            "inputs/tests/selfgravity_binary.athinput",
            ["gravity/mporder=3"],
            "gravity/mporder must be 2 or 4",
        ),
    ]
    table = []
    for label, input_file, args, expected in failures:
        expect_failure(input_file, args, expected)
        table.append({"case": label, "expected": expected})

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        source = (REPO_ROOT / "inputs/tests/selfgravity.athinput").read_text()
        no_gravity = tmp_path / "selfgravity_no_gravity.athinput"
        no_gravity.write_text(remove_block(source, "gravity"))
        expect_failure(no_gravity, None, "self_gravity source terms require")
        table.append(
            {
                "case": "missing gravity block",
                "expected": "self_gravity source terms require",
            }
        )

        no_controls = tmp_path / "selfgravity_no_controls.athinput"
        no_controls.write_text(
            remove_parameters(source, "gravity", {"threshold", "niteration"})
        )
        expect_failure(
            no_controls, None, "Set either gravity/threshold or gravity/niteration"
        )
        table.append(
            {
                "case": "missing convergence control",
                "expected": "Set either gravity/threshold or gravity/niteration",
            }
        )

        mhd_source = (REPO_ROOT / "inputs/tests/selfgravity_mhd.athinput").read_text()
        two_fluid = tmp_path / "selfgravity_two_fluid.athinput"
        two_fluid.write_text(
            mhd_source.rstrip()
            + """

<hydro>
eos = isothermal
reconstruct = plm
rsolver = llf
iso_sound_speed = 1.0

<ion-neutral>
drag_coeff = 1.0
"""
        )
        expect_failure(two_fluid, None, "Self-gravity currently supports one fluid")
        table.append(
            {
                "case": "two-fluid self-gravity",
                "expected": "Self-gravity currently supports one fluid",
            }
        )

    metrics["failure_cases"] = table


def mpi_summary(metrics: dict[str, Any]) -> str:
    mpi_metrics: dict[str, float] = {}
    for ranks in (2, 4):
        basename = f"docs_selfgravity_jeans_{ranks}rank"
        try:
            run_athena(
                "inputs/tests/selfgravity.athinput",
                basename,
                mpi_ranks=ranks,
                expect_jeans_banner=True,
            )
            output = parse_latest_binary(basename, "grav_phi")
            phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
            analytic, _ = jeans_analytic(x1, x2, x3, 1.0e-3)
            rel = sg.relative_l2((phi - np.mean(phi)) - analytic, analytic)
            mpi_metrics[f"{ranks}_rank_jeans_relative_l2"] = rel
        finally:
            cleanup(basename)

    device = "docs_selfgravity_mpi_root_device"
    host = "docs_selfgravity_mpi_root_host"
    try:
        run_athena(
            "inputs/tests/selfgravity.athinput",
            device,
            mpi_ranks=2,
            expect_jeans_banner=True,
        )
        run_athena(
            "inputs/tests/selfgravity.athinput",
            host,
            mpi_ranks=2,
            extra_args=["gravity/root_on_host=true"],
            expect_jeans_banner=True,
        )
        dphi, _, _, _ = sg.concatenate_field(
            parse_latest_binary(device, "grav_phi"), "grav_phi"
        )
        hphi, _, _, _ = sg.concatenate_field(
            parse_latest_binary(host, "grav_phi"), "grav_phi"
        )
        mpi_metrics["root_host_device_max_abs_phi_difference"] = float(
            np.max(np.abs(dphi - hphi))
        )
    finally:
        cleanup(device, host)

    basename = "docs_selfgravity_mpi_binary"
    try:
        run_athena(
            "inputs/tests/selfgravity_binary.athinput",
            basename,
            mpi_ranks=2,
            max_final_defect=2.0e-6,
        )
        output = parse_latest_binary(basename, "grav_phi")
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        gconst = 1.0 / (4.0 * math.pi)
        radius = 0.10
        r1 = np.sqrt((x1 - 0.15) ** 2 + x2 * x2 + x3 * x3)
        r2 = np.sqrt((x1 + 0.15) ** 2 + x2 * x2 + x3 * x3)
        analytic = -gconst / np.maximum(r1, 1.0e-300)
        analytic += -0.5 * gconst / np.maximum(r2, 1.0e-300)
        mask = (r1 > 2.0 * radius) & (r2 > 2.0 * radius)
        mpi_metrics["binary_multipole_relative_l2"] = sg.relative_l2(
            phi[mask] - analytic[mask], analytic[mask]
        )
    finally:
        cleanup(basename)

    metrics["mpi"] = mpi_metrics
    labels = list(mpi_metrics)
    values = [mpi_metrics[label] for label in labels]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    ax.barh(labels, values)
    ax.set_xscale("log")
    ax.set_xlabel("metric value")
    ax.set_title("MPI self-gravity validation metrics")
    return save(fig, "mpi_summary.png")


def gpu_summary(metrics: dict[str, Any]) -> str:
    gpu_metrics: dict[str, float] = {}
    basename = "docs_selfgravity_gpu_jeans"
    try:
        run_athena(
            "inputs/tests/selfgravity.athinput", basename, expect_jeans_banner=True
        )
        output = parse_latest_binary(basename, "grav_phi")
        phi, x1, x2, x3 = sg.concatenate_field(output, "grav_phi")
        analytic, _ = jeans_analytic(x1, x2, x3, 1.0e-3)
        gpu_metrics["periodic_jeans_relative_l2"] = sg.relative_l2(
            (phi - np.mean(phi)) - analytic, analytic
        )
    finally:
        cleanup(basename)

    device = "docs_selfgravity_gpu_root_device"
    host = "docs_selfgravity_gpu_root_host"
    try:
        run_athena("inputs/tests/selfgravity.athinput", device, expect_jeans_banner=True)
        run_athena(
            "inputs/tests/selfgravity.athinput",
            host,
            extra_args=["gravity/root_on_host=true"],
            expect_jeans_banner=True,
        )
        dphi, _, _, _ = sg.concatenate_field(
            parse_latest_binary(device, "grav_phi"), "grav_phi"
        )
        hphi, _, _, _ = sg.concatenate_field(
            parse_latest_binary(host, "grav_phi"), "grav_phi"
        )
        gpu_metrics["root_host_device_max_abs_phi_difference"] = float(
            np.max(np.abs(dphi - hphi))
        )
    finally:
        cleanup(device, host)

    metrics["gpu"] = gpu_metrics
    labels = list(gpu_metrics)
    values = [gpu_metrics[label] for label in labels]
    fig, ax = plt.subplots(figsize=(8, 3.5))
    ax.barh(labels, values)
    if all(value > 0.0 for value in values):
        ax.set_xscale("log")
    ax.set_xlabel("metric value")
    ax.set_title("GPU self-gravity validation metrics")
    return save(fig, "gpu_summary.png")


def format_metric(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4e}"
    if isinstance(value, int):
        return str(value)
    return str(value)


def metric_rows(metrics: dict[str, Any]) -> str:
    rows = ["| Group | Metric | Value |", "|---|---:|---:|"]
    for group, values in metrics.items():
        if group == "failure_cases" or not isinstance(values, dict):
            continue
        for key, value in values.items():
            if isinstance(value, dict):
                for subkey, subvalue in value.items():
                    rows.append(
                        f"| `{group}` | `{key}.{subkey}` | {format_metric(subvalue)} |"
                    )
            else:
                rows.append(f"| `{group}` | `{key}` | {format_metric(value)} |")
    return "\n".join(rows)


def failure_rows(metrics: dict[str, Any]) -> str:
    rows = ["| Failure case | Expected fatal diagnostic |", "|---|---|"]
    for item in metrics.get("failure_cases", []):
        rows.append(f"| {item['case']} | `{item['expected']}` |")
    return "\n".join(rows)


def write_results_page(
    figures: dict[str, str], metrics: dict[str, Any], with_mpi: bool, with_gpu: bool
) -> None:
    mpi_text = ""
    if with_mpi and "mpi_summary" in figures:
        mpi_text = f"""
## MPI Coverage

The MPI figure summarizes the 2-rank and 4-rank periodic Jeans tests, the
2-rank root host/device comparison, and the 2-rank isolated multipole reduction
test.

```{{figure}} {figures['mpi_summary']}
:alt: MPI self-gravity validation metrics
:width: 100%
```
"""
    else:
        mpi_text = """
## MPI Coverage

MPI figures were not generated in this artifact. Re-run the generator with
`--with-mpi --build-mpi` on a machine with MPI support.
"""

    gpu_text = """
## GPU Coverage

GPU smoke tests are defined in `tst/test_suite/selfgravity/test_selfgravity_gpu.py`.
This local artifact does not include CUDA results because no CUDA compiler/device
was available in the generation environment. On a CUDA node, run:

```bash
cd tst
python run_test_suite.py --gpu \\
    --test test_suite/selfgravity/test_selfgravity_gpu.py
```

Then regenerate this page with:

```bash
python docs/scripts/generate_selfgravity_validation.py \\
    --build-cpu --with-mpi --build-mpi --with-gpu --build-gpu
```
"""
    if with_gpu and "gpu_summary" in figures:
        gpu_text = f"""
## GPU Coverage

```{{figure}} {figures['gpu_summary']}
:alt: GPU self-gravity validation metrics
:width: 100%
```
"""

    coverage_rows = "\n".join(
        (
            "| Periodic Jeans solve | hydro/MHD convergence and analytic potential "
            "| potential slice, error slice, pointwise comparison |",
            "| Hydro/MHD equivalence | zero-field MHD equals hydro gravity "
            "| hydro/MHD slice and difference map |",
            "| Source coupling | momentum update against analytic Jeans force and "
            "final-density Poisson residual | momentum-increment and Poisson-residual "
            "slices |",
            "| Isolated boundaries | binary uniform-sphere multipole potential and "
            "force | potential slice, residual, radial error |",
            "| Restart | two cycles straight versus one cycle plus restart plus one "
            "cycle | density/momentum difference maps |",
            "| Controls and profiling | threshold mode, fixed iterations, root "
            "placement, profile phases | defect and timing plots |",
            "| Static refined hierarchy | BE collapse smoke on a static refined mesh "
            "| density and potential slices |",
            "| Dynamic AMR regrid | bounded experimental refine smoke "
            "| density and potential slices |",
            "| Invalid inputs | fatal diagnostics for bad user contracts "
            "| failure table |",
            "| MPI | 2-rank/4-rank Jeans, root placement, multipole reduction "
            "| MPI summary plot |",
            "| GPU | CUDA smoke hooks when available "
            "| GPU-node artifact pending in this local run |",
        )
    )

    text = f"""# Self-Gravity Validation Results

This page is a docs-facing companion to the self-gravity regression suite. It is
generated by `docs/scripts/generate_selfgravity_validation.py` from the same
input files and binary-output parser used by the tests.

## Test Coverage Matrix

| Test area | Regression coverage | Documentation artifact |
|---|---|---|
{coverage_rows}

## Quantitative Metrics

{metric_rows(metrics)}

## Periodic Jeans Potential

```{{figure}} {figures['periodic_jeans']}
:alt: Periodic Jeans self-gravity validation slices and pointwise comparison
:width: 100%
```

## Hydro/MHD Zero-Field Equivalence

```{{figure}} {figures['hydro_mhd']}
:alt: Hydro and zero-field MHD potential equivalence
:width: 100%
```

## Source-Term Oracle

```{{figure}} {figures['source_oracle']}
:alt: Source-term oracle comparison against analytic Jeans force
:width: 100%
```

## Isolated Binary Multipole

```{{figure}} {figures['binary_multipole']}
:alt: Isolated binary multipole potential and force validation
:width: 100%
```

## Restart Equivalence

```{{figure}} {figures['restart']}
:alt: Restart equivalence differences
:width: 100%
```

## Convergence Controls, Profiling, and Root Placement

```{{figure}} {figures['controls_profile']}
:alt: Self-gravity convergence controls and profile timing
:width: 100%
```

## Static Refined BE Smoke

```{{figure}} {figures['be_static_refined']}
:alt: Static refined BE self-gravity smoke slices
:width: 100%
```

## Dynamic AMR Regrid Smoke

```{{figure}} {figures['dynamic_amr_regrid']}
:alt: Dynamic AMR self-gravity regrid smoke slices
:width: 100%
```

## Failure-Mode Tests

{failure_rows(metrics)}

{mpi_text}

{gpu_text}
"""
    RESULTS_PAGE.write_text(text.rstrip() + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--build-cpu",
        action="store_true",
        help=(
            "configure and build a serial CPU AthenaK binary before generating "
            "CPU figures"
        ),
    )
    parser.add_argument(
        "--with-mpi", action="store_true", help="also generate MPI validation figures"
    )
    parser.add_argument(
        "--build-mpi",
        action="store_true",
        help="reconfigure and build with Athena_ENABLE_MPI=ON before MPI figures",
    )
    parser.add_argument(
        "--with-gpu",
        action="store_true",
        help="also generate CUDA-node validation figures",
    )
    parser.add_argument(
        "--build-gpu",
        action="store_true",
        help="reconfigure and build with Kokkos_ENABLE_CUDA=On before GPU figures",
    )
    args = parser.parse_args()

    if args.build_cpu:
        build(["-D", "Athena_ENABLE_MPI=OFF", "-D", "Kokkos_ENABLE_CUDA=OFF"])
    ensure_build()
    if STATIC_DIR.exists():
        shutil.rmtree(STATIC_DIR)
    STATIC_DIR.mkdir(parents=True)

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )

    metrics: dict[str, Any] = {}
    figures = {
        "periodic_jeans": periodic_jeans(metrics),
        "hydro_mhd": hydro_mhd_equivalence(metrics),
        "source_oracle": source_oracle(metrics),
        "binary_multipole": binary_multipole(metrics),
        "restart": restart_equivalence(metrics),
        "controls_profile": controls_profile(metrics),
        "be_static_refined": be_static_refined(metrics),
        "dynamic_amr_regrid": dynamic_amr_regrid(metrics),
    }
    failure_cases(metrics)

    if args.with_mpi:
        if args.build_mpi:
            build(["-D", "Athena_ENABLE_MPI=ON", "-D", "Kokkos_ENABLE_CUDA=OFF"])
        figures["mpi_summary"] = mpi_summary(metrics)

    if args.with_gpu:
        if args.build_gpu:
            build(["-D", "Athena_ENABLE_MPI=OFF", "-D", "Kokkos_ENABLE_CUDA=ON"])
        figures["gpu_summary"] = gpu_summary(metrics)

    write_results_page(figures, metrics, args.with_mpi, args.with_gpu)
    print(f"Wrote {RESULTS_PAGE}")
    print(f"Wrote figures under {STATIC_DIR}")


if __name__ == "__main__":
    main()
