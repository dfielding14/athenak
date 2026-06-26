#!/usr/bin/env python3
"""Build a strict Figure-8-style Q011 shock morphology figure."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import matplotlib as mpl
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

if __package__:
    from . import analyze_q011_section54_outputs as output_primitives
    from . import q011_section54_model as model
    from . import q011_section54_particles as particle_primitives
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    import analyze_q011_section54_outputs as output_primitives
    import q011_section54_model as model
    import q011_section54_particles as particle_primitives
    from pvtk_particles import ParticleVTKData, read_particle_vtk


SCHEMA_VERSION = 1
RECORD_TYPE = "q011_section54_figure8_metrics_v1"
MANIFEST_RECORD_TYPE = "q011_section54_figure8_manifest_v1"
MHD_FIELDS = ("dens", "velx", "vely", "velz", "eint", "bcc1", "bcc2", "bcc3")
PVTK_SCALARS = frozenset(
    {
        "gid",
        "ptag",
        "species",
        "cr_source",
        "macro_weight",
        "birth_time",
        "deltaf_f0",
        "deltaf_weight",
    }
)
PVTK_VECTORS = frozenset({"vel"})
SHOCK_SEARCH_OFFSETS = (-1200.0, 1200.0)
MAXIMUM_SHOCK_OFFSET = 600.0
DEFAULT_X_OFFSETS = (-1200.0, 1200.0)
_PVTK_HEADER = re.compile(
    rb"\A# vtk DataFile Version 2[.]0\n"
    rb"# AthenaK particle data at time= ([^ \n]+)  nranks= (0|[1-9][0-9]*)  "
    rb"cycle=(0|[1-9][0-9]*)  variables=([A-Za-z0-9_]+)\n"
)
_SAFE_STEM = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")


class Figure8Error(ValueError):
    """Reject incomplete, mismatched, or non-finite Figure 8 inputs."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise Figure8Error(message)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _artifact(path: Path) -> dict[str, object]:
    return {
        "path": str(path.resolve()),
        "sha256": _sha256(path),
        "byte_count": path.stat().st_size,
    }


def _canonical_json(value: Mapping[str, object]) -> str:
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def _finite_positive(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise Figure8Error(f"{label} must be a real scalar") from error
    _require(
        math.isfinite(result) and result > 0.0,
        f"{label} must be finite and positive",
    )
    return result


def _runtime_parameter(
    dataset: output_primitives.AthenaBinaryDataset, block: str, name: str
) -> float:
    try:
        value = dataset.input_parameters[block][name]
    except KeyError as error:
        raise Figure8Error(f"mhd_w_bcc header is missing {block}/{name}") from error
    return _finite_positive(value, f"{block}/{name}")


def _particle_header(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        prefix = stream.read(4096)
    match = _PVTK_HEADER.match(prefix)
    _require(match is not None, "prtcl_all PVTK execution header is malformed")
    assert match is not None
    try:
        time = float(match.group(1))
        nranks = int(match.group(2))
        cycle = int(match.group(3))
        variables = match.group(4).decode("ascii")
    except (UnicodeDecodeError, ValueError) as error:
        raise Figure8Error("prtcl_all PVTK execution header is invalid") from error
    _require(math.isfinite(time) and time >= 0.0, "prtcl_all time is invalid")
    _require(nranks > 0 and cycle >= 0, "prtcl_all rank or cycle metadata is invalid")
    _require(variables == "prtcl_all", "particle product is not prtcl_all")
    return {"time": time, "nranks": nranks, "cycle": cycle, "variables": variables}


def _compose_mhd_state(
    dataset: output_primitives.AthenaBinaryDataset,
) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray, np.ndarray]:
    _require(
        set(dataset.variable_names) == set(MHD_FIELDS),
        "mhd_w_bcc field inventory is incomplete or contains unexpected fields",
    )
    grids = {
        field: output_primitives.compose_leaf_field(dataset, field)
        for field in MHD_FIELDS
    }
    reference = grids["dens"]
    for field, grid in grids.items():
        _require(
            np.array_equal(grid.x1_faces, reference.x1_faces)
            and np.array_equal(grid.x2_faces, reference.x2_faces)
            and np.array_equal(grid.x3_faces, reference.x3_faces)
            and np.array_equal(grid.source_levels, reference.source_levels),
            f"mhd_w_bcc composite grid mismatch for {field}",
        )
        _require(grid.values.shape[0] == 1, "Figure 8 requires a 2D x1-x2 dataset")
    fields = {
        field: np.asarray(grid.values[0], dtype=np.float64)
        for field, grid in grids.items()
    }
    return fields, reference.x1_faces, reference.x2_faces, reference.source_levels[0]


def _validate_particles(data: ParticleVTKData) -> dict[str, np.ndarray]:
    _require(set(data.scalars) == PVTK_SCALARS, "prtcl_all scalar inventory drifted")
    _require(set(data.vectors) == PVTK_VECTORS, "prtcl_all vector inventory drifted")
    points = np.asarray(data.points, dtype=np.float64)
    velocity = np.asarray(data.vectors["vel"], dtype=np.float64)
    _require(
        points.ndim == 2 and points.shape[1:] == (3,) and velocity.shape == points.shape,
        "prtcl_all points and velocity must have shape (nparticle, 3)",
    )
    count = points.shape[0]
    _require(np.all(np.isfinite(points)) and np.all(np.isfinite(velocity)),
             "prtcl_all points or velocity are non-finite")
    _require(
        np.array_equal(points, points.astype(np.float32).astype(np.float64))
        and np.array_equal(velocity, velocity.astype(np.float32).astype(np.float64)),
        "prtcl_all positions and velocities must be exact decoded float32 values",
    )
    arrays = {"points": points, "velocity": velocity}
    for name, raw in data.scalars.items():
        values = np.asarray(raw)
        _require(values.shape == (count,), f"prtcl_all scalar {name} shape drifted")
        _require(np.all(np.isfinite(values)), f"prtcl_all scalar {name} is non-finite")
        arrays[name] = values
    _require(
        np.all(arrays["gid"] >= 0)
        and np.all(arrays["ptag"] >= 0)
        and np.unique(arrays["ptag"]).size == count,
        "prtcl_all gid/ptag identity is invalid",
    )
    _require(np.all(arrays["species"] == 0), "prtcl_all species inventory is not Q011")
    _require(
        np.all(arrays["macro_weight"] >= 0.0),
        "prtcl_all macro weights are negative",
    )
    return arrays


def analyze_figure8(
    mhd_path: Path,
    particle_path: Path,
    *,
    x_offsets: tuple[float, float] = DEFAULT_X_OFFSETS,
) -> dict[str, object]:
    """Read and reduce one exact common-cycle MHD/PVTK snapshot."""
    _require(x_offsets[0] < x_offsets[1], "near-shock x offsets must increase")
    try:
        dataset = output_primitives.read_athenak_binary(mhd_path)
    except (OSError, output_primitives.AnalysisError) as error:
        raise Figure8Error(f"unable to read mhd_w_bcc: {error}") from error
    header = _particle_header(particle_path)
    _require(dataset.cycle == header["cycle"], "mhd_w_bcc and prtcl_all cycles disagree")
    projected_time = float(format(float(header["time"]), ".6g"))
    _require(
        dataset.time == projected_time,
        "mhd_w_bcc time is not the six-significant-digit projection of PVTK time",
    )

    fields, x1_faces, x2_faces, source_levels = _compose_mhd_state(dataset)
    rho0 = _runtime_parameter(dataset, "problem", "ps_rho0")
    b0 = _runtime_parameter(dataset, "problem", "ps_b0")
    u0 = _runtime_parameter(dataset, "problem", "ps_u0")
    light_speed = _runtime_parameter(dataset, "particles", "pic_cr_light_speed")
    macro_mass = _runtime_parameter(dataset, "particles", "deposit_qscale")
    _require(
        u0 == model.UPSTREAM_SPEED_U0 and light_speed == model.LIGHT_SPEED,
        "mhd_w_bcc Q011 particle normalization drifted",
    )

    density = fields["dens"] / rho0
    bmag = np.sqrt(fields["bcc1"] ** 2 + fields["bcc2"] ** 2 + fields["bcc3"] ** 2)
    bmag /= b0
    _require(np.all(density > 0.0), "normalized gas density must be positive")
    _require(np.any(bmag > 0.0), "normalized magnetic magnitude has no positive cells")
    x1_centers = 0.5 * (x1_faces[:-1] + x1_faces[1:])
    density_profile = np.mean(density, axis=0)
    ideal_x = model.x_ideal(header["time"])
    try:
        front = output_primitives.detect_shock_front(
            x1_centers,
            density_profile,
            search_window=(
                ideal_x + SHOCK_SEARCH_OFFSETS[0],
                ideal_x + SHOCK_SEARCH_OFFSETS[1],
            ),
            gradient_sign="negative",
        )
    except output_primitives.AnalysisError as error:
        raise Figure8Error(f"shock-front detection failed: {error}") from error
    _require(
        abs(front.x1 - ideal_x) <= MAXIMUM_SHOCK_OFFSET,
        "detected shock front exceeds the Q011 ideal-surface offset bound",
    )
    x_limits = (front.x1 + x_offsets[0], front.x1 + x_offsets[1])
    _require(
        x_limits[0] >= x1_faces[0] and x_limits[1] <= x1_faces[-1],
        "near-shock x extent lies outside the retained MHD domain",
    )
    x_indices = np.flatnonzero(
        (x1_centers >= x_limits[0]) & (x1_centers <= x_limits[1])
    )
    _require(x_indices.size >= 2, "near-shock x extent selects fewer than two cells")
    _require(np.all(np.diff(x_indices) == 1), "near-shock x selection is not contiguous")
    first_x = int(x_indices[0])
    last_x = int(x_indices[-1]) + 1
    panel_x_faces = x1_faces[first_x : last_x + 1]

    try:
        particle_data = read_particle_vtk(particle_path)
    except (OSError, ValueError) as error:
        raise Figure8Error(f"unable to read prtcl_all: {error}") from error
    arrays = _validate_particles(particle_data)
    source_admitted = arrays["cr_source"] == 1
    birth_admitted = source_admitted & (arrays["birth_time"] >= 45.0)
    selected = birth_admitted & (arrays["macro_weight"] > 0.0)
    _require(np.any(selected), "prtcl_all has no selected positive-weight Q011 particles")
    try:
        chi = particle_primitives.reconstruct_chi_from_pvtk_velocity(arrays["velocity"])
    except particle_primitives.ParticleReducerError as error:
        raise Figure8Error(f"particle chi reconstruction failed: {error}") from error
    _require(np.all(chi[selected] > 0.0), "selected particle chi must be positive")
    momentum_squared = chi * u0**2
    specific_energy = momentum_squared / (
        np.sqrt(1.0 + momentum_squared / light_speed**2) + 1.0
    )
    energy_weight = macro_mass * arrays["macro_weight"] * specific_energy
    _require(
        np.all(np.isfinite(energy_weight)) and np.all(energy_weight[selected] > 0.0),
        "selected particle energy weights must be finite and positive",
    )

    log10_chi = np.log10(chi[selected])
    log10_chi_edges = np.log10(np.asarray(particle_primitives.CHI_BIN_EDGES))
    phase_energy, _, _ = np.histogram2d(
        arrays["points"][selected, 0],
        log10_chi,
        bins=(panel_x_faces, log10_chi_edges),
        weights=energy_weight[selected],
    )
    phase_counts, _, _ = np.histogram2d(
        arrays["points"][selected, 0],
        log10_chi,
        bins=(panel_x_faces, log10_chi_edges),
    )
    total_selected_energy = float(np.sum(energy_weight[selected]))
    phase_energy_total = float(np.sum(phase_energy))
    bin_measure = np.diff(panel_x_faces)[:, None] * np.diff(log10_chi_edges)[None, :]
    phase_density = phase_energy / total_selected_energy / bin_measure
    _require(
        np.any(phase_energy > 0.0) and np.all(np.isfinite(phase_density)),
        "energy-weighted phase-space histogram is empty or non-finite",
    )

    wrong_source = arrays["cr_source"] != 1
    early_birth = source_admitted & (arrays["birth_time"] < 45.0)
    nonpositive_weight = birth_admitted & (arrays["macro_weight"] <= 0.0)
    metrics: dict[str, object] = {
        "record_type": RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "qualification_effect": "derived_figure_only_no_claim_closure",
        "snapshot": {
            "cycle": int(dataset.cycle),
            "mhd_time": float(dataset.time),
            "particle_time": float(header["time"]),
            "particle_nranks": int(header["nranks"]),
            "mesh_time_projection": projected_time,
        },
        "normalization": {
            "rho0": rho0,
            "b0": b0,
            "u0": u0,
            "particle_light_speed": light_speed,
            "particle_macro_mass": macro_mass,
            "chi_formula": "gamma(v)^2 * |v|^2 / u0^2",
            "specific_kinetic_energy_formula": (
                "p_over_m_squared / (sqrt(1 + p_over_m_squared / C_squared) + 1)"
            ),
            "phase_space_weight": (
                "particle_macro_mass * macro_weight * "
                "relativistic_specific_kinetic_energy"
            ),
        },
        "shock": {
            "detector": "unique_strongest_negative_density_gradient",
            "ideal_surface_x1": ideal_x,
            "detected_front_x1": front.x1,
            "density_gradient": front.density_gradient,
            "offset_from_ideal_surface": front.x1 - ideal_x,
            "near_shock_x_offsets": list(x_offsets),
            "near_shock_x_limits": list(x_limits),
        },
        "mesh": {
            "near_shock_cell_shape_y_x": [int(density.shape[0]), int(x_indices.size)],
            "source_refinement_levels": sorted(
                int(value) for value in np.unique(source_levels[:, first_x:last_x])
            ),
            "density_over_rho0_min": float(np.min(density[:, first_x:last_x])),
            "density_over_rho0_max": float(np.max(density[:, first_x:last_x])),
            "bmag_over_b0_min": float(np.min(bmag[:, first_x:last_x])),
            "bmag_over_b0_max": float(np.max(bmag[:, first_x:last_x])),
        },
        "particle_filter": {
            "all_particle_count": int(arrays["points"].shape[0]),
            "rejected_wrong_source_count": int(np.count_nonzero(wrong_source)),
            "rejected_early_birth_count": int(np.count_nonzero(early_birth)),
            "rejected_nonpositive_weight_count": int(
                np.count_nonzero(nonpositive_weight)
            ),
            "selected_particle_count": int(np.count_nonzero(selected)),
            "selection": ["cr_source == 1", "birth_time >= 45", "macro_weight > 0"],
        },
        "phase_space": {
            "x1_edges": panel_x_faces.tolist(),
            "log10_chi_edges": log10_chi_edges.tolist(),
            "counts": phase_counts.astype(np.int64).tolist(),
            "energy_weighted_histogram": phase_energy.tolist(),
            "energy_probability_density_per_x1_per_log10_chi": phase_density.tolist(),
            "selected_total_energy_weight": total_selected_energy,
            "histogrammed_energy_weight": phase_energy_total,
            "histogrammed_energy_fraction": phase_energy_total / total_selected_energy,
            "histogrammed_particle_count": int(np.sum(phase_counts)),
            "selected_chi_min": float(np.min(chi[selected])),
            "selected_chi_max": float(np.max(chi[selected])),
        },
    }
    return {
        "metrics": metrics,
        "density": density[:, first_x:last_x],
        "bmag": bmag[:, first_x:last_x],
        "x1_faces": panel_x_faces,
        "x2_faces": x2_faces,
        "log10_chi_edges": log10_chi_edges,
        "phase_density": phase_density,
    }


def _style() -> None:
    plt.style.use("default")
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "font.family": "serif",
            "font.serif": ["DejaVu Serif", "Times New Roman", "Times"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _positive_log_norm(values: np.ndarray) -> colors.LogNorm:
    positive = values[values > 0.0]
    _require(positive.size > 0, "logarithmic panel contains no positive values")
    lower = float(np.percentile(positive, 1.0))
    upper = float(np.max(positive))
    if lower == upper:
        lower = max(np.nextafter(upper, 0.0), upper * 0.5)
    return colors.LogNorm(vmin=lower, vmax=upper)


def render_figure8(
    reduction: Mapping[str, object], output_paths: Sequence[Path], dpi: int
) -> None:
    """Render one already validated Figure 8 reduction."""
    _require(dpi > 0, "dpi must be positive")
    _style()
    density = np.asarray(reduction["density"])
    bmag = np.asarray(reduction["bmag"])
    x1_faces = np.asarray(reduction["x1_faces"])
    x2_faces = np.asarray(reduction["x2_faces"])
    log10_chi_edges = np.asarray(reduction["log10_chi_edges"])
    phase_density = np.asarray(reduction["phase_density"])
    metrics = reduction["metrics"]
    shock = metrics["shock"]
    snapshot = metrics["snapshot"]

    fig, axes = plt.subplots(
        3, 1, figsize=(11.2, 7.4), sharex=True, constrained_layout=True
    )
    panels = (
        (density, "viridis", None, r"$\rho/\rho_0$"),
        (bmag, "magma", _positive_log_norm(bmag), r"$|B|/B_0$"),
    )
    for index, (values, cmap, norm, label) in enumerate(panels):
        mesh = axes[index].pcolormesh(
            x1_faces, x2_faces, values, shading="flat", cmap=cmap, norm=norm
        )
        colorbar = fig.colorbar(mesh, ax=axes[index], pad=0.012)
        colorbar.set_label(label)
        axes[index].set_ylabel(r"$x_2/(c/\omega_{pi})$")
        axes[index].text(
            0.012, 0.90, f"({chr(ord('a') + index)})", transform=axes[index].transAxes,
            color="white", fontweight="bold", fontsize=11,
        )

    phase_plot = np.where(phase_density.T > 0.0, phase_density.T, np.nan)
    phase = axes[2].pcolormesh(
        x1_faces,
        log10_chi_edges,
        phase_plot,
        shading="flat",
        cmap="cividis",
        norm=_positive_log_norm(phase_density),
    )
    colorbar = fig.colorbar(phase, ax=axes[2], pad=0.012)
    colorbar.set_label(r"energy-weighted probability density")
    axes[2].set_ylabel(r"$\log_{10}\chi$")
    axes[2].set_xlabel(r"$x_1/(c/\omega_{pi})$")
    axes[2].text(
        0.012, 0.90, "(c)", transform=axes[2].transAxes,
        color="white", fontweight="bold", fontsize=11,
    )
    for axis in axes:
        axis.axvline(shock["ideal_surface_x1"], color="#f0f0f0", lw=0.8, ls=":")
        axis.axvline(shock["detected_front_x1"], color="#ffffff", lw=0.9, ls="--")
        axis.set_xlim(x1_faces[0], x1_faces[-1])
    fig.suptitle(
        rf"Q011 Section 5.4: $\Omega_0 t={snapshot['particle_time']:.6g}$, "
        rf"cycle {snapshot['cycle']}"
    )
    for path in output_paths:
        if path.suffix.lower() == ".png":
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
        else:
            fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def make_figure8(
    mhd_path: Path,
    particle_path: Path,
    output_dir: Path,
    *,
    stem: str = "q011_section54_figure8_v1",
    x_offsets: tuple[float, float] = DEFAULT_X_OFFSETS,
    dpi: int = 300,
) -> list[Path]:
    """Analyze exact inputs, render PNG/PDF, and write metrics plus manifest."""
    _require(_SAFE_STEM.fullmatch(stem) is not None, "output stem is unsafe")
    mhd_path = mhd_path.resolve(strict=True)
    particle_path = particle_path.resolve(strict=True)
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    reduction = analyze_figure8(mhd_path, particle_path, x_offsets=x_offsets)
    png_path = output_dir / f"{stem}.png"
    pdf_path = output_dir / f"{stem}.pdf"
    metrics_path = output_dir / f"{stem}_metrics.json"
    manifest_path = output_dir / f"{stem}_manifest.json"
    render_figure8(reduction, (png_path, pdf_path), dpi)
    metrics_path.write_text(_canonical_json(reduction["metrics"]), encoding="utf-8")
    manifest = {
        "record_type": MANIFEST_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "qualification_effect": "derived_figure_only_no_claim_closure",
        "generator": _artifact(Path(__file__)),
        "inputs": [_artifact(mhd_path), _artifact(particle_path)],
        "outputs": [_artifact(path) for path in (png_path, pdf_path, metrics_path)],
        "snapshot": reduction["metrics"]["snapshot"],
        "analysis_parameters": {
            "near_shock_x_offsets": list(x_offsets),
            "dpi": dpi,
            "particle_filter": [
                "cr_source == 1",
                "birth_time >= 45",
                "macro_weight > 0",
            ],
            "cross_cycle_substitution": False,
        },
    }
    manifest_path.write_text(_canonical_json(manifest), encoding="utf-8")
    return [png_path, pdf_path, metrics_path, manifest_path]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mhd", type=Path, required=True, help="Exact mhd_w_bcc Athena binary"
    )
    parser.add_argument(
        "--particles", type=Path, required=True, help="Exact prtcl_all PVTK"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stem", default="q011_section54_figure8_v1")
    parser.add_argument("--x-offset-min", type=float, default=DEFAULT_X_OFFSETS[0])
    parser.add_argument("--x-offset-max", type=float, default=DEFAULT_X_OFFSETS[1])
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    outputs = make_figure8(
        args.mhd,
        args.particles,
        args.output_dir,
        stem=args.stem,
        x_offsets=(args.x_offset_min, args.x_offset_max),
        dpi=args.dpi,
    )
    print(_canonical_json({"outputs": [str(path) for path in outputs]}), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
