#!/usr/bin/env python3
"""Render a strict time-ordered Q011 Figure-8-style PNG sequence."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import tempfile
from typing import Mapping, Sequence

import matplotlib as mpl

mpl.use("Agg", force=True)

from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np

if __package__:
    from . import make_q011_section54_figure8_v1 as figure8
    from . import pvtk_particles
else:
    import make_q011_section54_figure8_v1 as figure8
    import pvtk_particles


SCHEMA_VERSION = 1
MANIFEST_RECORD_TYPE = "q011_section54_evolution_manifest_v1"
DEFAULT_STEM = "q011_section54_evolution_v1"
DEFAULT_X_OFFSETS = (-800.0, 1200.0)
FRAME_DIGITS = 4
_SAFE_STEM = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}")
_MHD_NAME = re.compile(
    r"(?P<basename>[^/]+)[.]mhd_w_bcc[.](?P<index>[0-9]{5,})[.]bin"
)
_PARTICLE_NAME = re.compile(
    r"(?P<basename>[^/]+)[.]prtcl_all[.](?P<index>[0-9]{5,})[.]part[.]vtk"
)


class EvolutionError(ValueError):
    """Reject an incomplete, ambiguous, or inconsistent evolution sequence."""


@dataclass(frozen=True)
class ProductPair:
    """One exact same-output-index MHD and particle product pair."""

    output_index: int
    index_token: str
    basename: str
    mhd_path: Path
    particle_path: Path


@dataclass(frozen=True)
class FrameReduction:
    """The bounded data retained from one strict Figure 8 analysis."""

    pair: ProductPair
    cycle: int
    mhd_time: float
    particle_time: float
    particle_nranks: int
    mesh_time_projection: float
    ideal_surface_x1: float
    detected_front_x1: float
    density_gradient: float
    offset_from_ideal_surface: float
    selected_particle_count: int
    histogrammed_energy_fraction: float
    density: np.ndarray
    bmag: np.ndarray
    x1_faces: np.ndarray
    x2_faces: np.ndarray
    log10_chi_edges: np.ndarray
    phase_density: np.ndarray


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EvolutionError(message)


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


def _validate_options(
    stem: str, x_offsets: tuple[float, float], dpi: int
) -> None:
    _require(_SAFE_STEM.fullmatch(stem) is not None, "output stem is unsafe")
    _require(
        all(math.isfinite(value) for value in x_offsets)
        and x_offsets[0] < x_offsets[1],
        "shock-relative x offsets must be finite and increasing",
    )
    _require(type(dpi) is int and dpi > 0, "dpi must be a positive integer")


def _scan_products(
    directory: Path, marker: str, pattern: re.Pattern[str], label: str
) -> dict[int, tuple[str, str, Path]]:
    _require(directory.is_dir(), f"run root is missing the {directory.name} directory")
    products: dict[int, tuple[str, str, Path]] = {}
    try:
        entries = sorted(directory.iterdir(), key=lambda path: path.name)
    except OSError as error:
        raise EvolutionError(f"unable to inventory {directory}: {error}") from error
    for path in entries:
        if marker not in path.name:
            continue
        _require(path.is_file(), f"{label} product is not a regular file: {path}")
        match = pattern.fullmatch(path.name)
        _require(match is not None, f"malformed {label} product name: {path.name}")
        assert match is not None
        index_token = match.group("index")
        output_index = int(index_token)
        _require(
            output_index not in products,
            f"duplicate {label} product for output index {output_index}",
        )
        try:
            resolved = path.resolve(strict=True)
        except OSError as error:
            raise EvolutionError(
                f"unable to resolve {label} product {path}: {error}"
            ) from error
        products[output_index] = (match.group("basename"), index_token, resolved)
    _require(products, f"run root contains no {label} products")
    return products


def discover_product_pairs(run_root: Path) -> tuple[ProductPair, ...]:
    """Discover one unambiguous, gap-free MHD/PVTK pair per output index."""
    try:
        root = run_root.resolve(strict=True)
    except OSError as error:
        raise EvolutionError(f"unable to resolve run root {run_root}: {error}") from error
    _require(root.is_dir(), "run root must be a directory")
    mhd = _scan_products(root / "bin", ".mhd_w_bcc.", _MHD_NAME, "mhd_w_bcc")
    particles = _scan_products(
        root / "pvtk", ".prtcl_all.", _PARTICLE_NAME, "prtcl_all"
    )
    mhd_indices = set(mhd)
    particle_indices = set(particles)
    _require(
        mhd_indices == particle_indices,
        "mhd_w_bcc/prtcl_all output-index inventories disagree: "
        f"missing prtcl_all={sorted(mhd_indices - particle_indices)}, "
        f"missing mhd_w_bcc={sorted(particle_indices - mhd_indices)}",
    )
    ordered_indices = sorted(mhd_indices)
    _require(ordered_indices[0] == 0, "output inventory is missing index 00000")
    missing = sorted(set(range(ordered_indices[-1] + 1)) - mhd_indices)
    _require(not missing, f"output inventory has index gaps: {missing}")

    pairs = []
    basenames = set()
    for output_index in ordered_indices:
        mhd_basename, mhd_token, mhd_path = mhd[output_index]
        particle_basename, particle_token, particle_path = particles[output_index]
        _require(
            mhd_token == particle_token,
            f"output index {output_index} uses mismatched index tokens",
        )
        _require(
            mhd_basename == particle_basename,
            f"output index {output_index} uses mismatched product basenames",
        )
        basenames.add(mhd_basename)
        pairs.append(
            ProductPair(
                output_index=output_index,
                index_token=mhd_token,
                basename=mhd_basename,
                mhd_path=mhd_path,
                particle_path=particle_path,
            )
        )
    _require(len(basenames) == 1, "product basename changes within the run")
    return tuple(pairs)


def _validate_initial_pair(pair: ProductPair) -> dict[str, object]:
    """Validate the canonical t=0 pair without asking Figure 8 for particles."""
    _require(pair.output_index == 0, "initial exclusion must use output index 00000")
    try:
        header = figure8._particle_header(pair.particle_path)
        dataset = figure8.output_primitives.read_athenak_binary(pair.mhd_path)
    except (
        OSError,
        figure8.Figure8Error,
        figure8.output_primitives.AnalysisError,
    ) as error:
        raise EvolutionError(
            f"unable to validate index 00000 metadata: {error}"
        ) from error
    particle_time = float(header["time"])
    particle_cycle = int(header["cycle"])
    projected_time = float(format(particle_time, ".6g"))
    _require(
        particle_time == 0.0 and particle_cycle == 0,
        "output index 00000 is not the canonical cycle-zero t=0 snapshot",
    )
    _require(
        dataset.cycle == particle_cycle and dataset.time == projected_time,
        "cycle-zero mhd_w_bcc and prtcl_all metadata disagree",
    )
    try:
        particle_data = figure8.read_particle_vtk(pair.particle_path)
        arrays = figure8._validate_particles(particle_data)
    except (OSError, ValueError, figure8.Figure8Error) as error:
        raise EvolutionError(
            f"unable to validate index 00000 particles: {error}"
        ) from error
    particle_count = int(arrays["points"].shape[0])
    _require(
        particle_count == 0,
        "cycle-zero t=0 prtcl_all product is not a no-particle initial frame",
    )
    return {
        "output_index": pair.output_index,
        "index_token": pair.index_token,
        "reason": "cycle_zero_t0_no_particle_initial_frame",
        "cycle": particle_cycle,
        "mhd_time": float(dataset.time),
        "particle_time": particle_time,
        "particle_count": particle_count,
    }


def _mapping(value: object, label: str) -> Mapping[str, object]:
    _require(isinstance(value, Mapping), f"strict analyzer returned invalid {label}")
    assert isinstance(value, Mapping)
    return value


def _array(value: object, label: str) -> np.ndarray:
    # Detach panel slices from the strict analyzer's full-domain backing arrays.
    array = np.array(value, dtype=np.float64, copy=True)
    _require(np.all(np.isfinite(array)), f"strict analyzer returned non-finite {label}")
    return array


def _compact_reduction(
    pair: ProductPair, reduction: Mapping[str, object]
) -> FrameReduction:
    """Drop bulky JSON histogram tables while retaining rendered reductions."""
    metrics = _mapping(reduction.get("metrics"), "metrics")
    snapshot = _mapping(metrics.get("snapshot"), "snapshot metrics")
    shock = _mapping(metrics.get("shock"), "shock metrics")
    particle_filter = _mapping(
        metrics.get("particle_filter"), "particle-filter metrics"
    )
    phase_space = _mapping(metrics.get("phase_space"), "phase-space metrics")

    density = _array(reduction.get("density"), "density")
    bmag = _array(reduction.get("bmag"), "magnetic magnitude")
    x1_faces = _array(reduction.get("x1_faces"), "x1 faces")
    x2_faces = _array(reduction.get("x2_faces"), "x2 faces")
    log10_chi_edges = _array(
        reduction.get("log10_chi_edges"), "log10 chi edges"
    )
    phase_density = _array(reduction.get("phase_density"), "phase density")
    _require(
        density.ndim == 2 and bmag.shape == density.shape,
        "strict analyzer returned incompatible MHD panel shapes",
    )
    ny, nx = density.shape
    _require(
        x1_faces.shape == (nx + 1,) and x2_faces.shape == (ny + 1,),
        "strict analyzer returned incompatible MHD face coordinates",
    )
    _require(
        phase_density.ndim == 2
        and phase_density.shape[0] == nx
        and log10_chi_edges.shape == (phase_density.shape[1] + 1,),
        "strict analyzer returned incompatible phase-space panel shapes",
    )
    _require(
        np.all(np.diff(x1_faces) > 0.0)
        and np.all(np.diff(x2_faces) > 0.0)
        and np.all(np.diff(log10_chi_edges) > 0.0),
        "strict analyzer returned non-increasing panel coordinates",
    )
    _require(np.all(density > 0.0), "strict analyzer returned non-positive density")
    _require(
        np.all(bmag >= 0.0) and np.any(bmag > 0.0),
        "strict analyzer returned invalid magnetic magnitude",
    )
    _require(
        np.all(phase_density >= 0.0) and np.any(phase_density > 0.0),
        "strict analyzer returned invalid phase density",
    )

    frame = FrameReduction(
        pair=pair,
        cycle=int(snapshot["cycle"]),
        mhd_time=float(snapshot["mhd_time"]),
        particle_time=float(snapshot["particle_time"]),
        particle_nranks=int(snapshot["particle_nranks"]),
        mesh_time_projection=float(snapshot["mesh_time_projection"]),
        ideal_surface_x1=float(shock["ideal_surface_x1"]),
        detected_front_x1=float(shock["detected_front_x1"]),
        density_gradient=float(shock["density_gradient"]),
        offset_from_ideal_surface=float(shock["offset_from_ideal_surface"]),
        selected_particle_count=int(particle_filter["selected_particle_count"]),
        histogrammed_energy_fraction=float(
            phase_space["histogrammed_energy_fraction"]
        ),
        density=density,
        bmag=bmag,
        x1_faces=x1_faces,
        x2_faces=x2_faces,
        log10_chi_edges=log10_chi_edges,
        phase_density=phase_density,
    )
    scalar_values = (
        frame.mhd_time,
        frame.particle_time,
        frame.mesh_time_projection,
        frame.ideal_surface_x1,
        frame.detected_front_x1,
        frame.density_gradient,
        frame.offset_from_ideal_surface,
        frame.histogrammed_energy_fraction,
    )
    _require(
        all(math.isfinite(value) for value in scalar_values),
        "strict analyzer returned non-finite frame metrics",
    )
    _require(
        frame.cycle > 0
        and frame.particle_time > 0.0
        and frame.particle_nranks > 0
        and frame.selected_particle_count > 0,
        "strict analyzer returned a non-renderable frame",
    )
    return frame


def analyze_sequence(
    pairs: Sequence[ProductPair],
    *,
    x_offsets: tuple[float, float] = DEFAULT_X_OFFSETS,
) -> tuple[tuple[FrameReduction, ...], dict[str, object]]:
    """Strictly analyze every nonzero pair and order reductions by physical time."""
    _require(bool(pairs), "product-pair inventory is empty")
    excluded = _validate_initial_pair(pairs[0])
    frames = []
    for pair in pairs[1:]:
        try:
            reduction = figure8.analyze_figure8(
                pair.mhd_path, pair.particle_path, x_offsets=x_offsets
            )
        except (OSError, ValueError, figure8.Figure8Error) as error:
            raise EvolutionError(
                f"strict Figure 8 analysis failed for output index "
                f"{pair.index_token}: {error}"
            ) from error
        frames.append(_compact_reduction(pair, reduction))
    _require(frames, "run contains no nonzero particle frames")
    frames.sort(key=lambda frame: frame.particle_time)
    times = [frame.particle_time for frame in frames]
    cycles = [frame.cycle for frame in frames]
    _require(
        all(right > left for left, right in zip(times, times[1:])),
        "rendered particle times are duplicate or non-increasing",
    )
    _require(
        all(right > left for left, right in zip(cycles, cycles[1:])),
        "rendered cycles do not increase with physical time",
    )
    return tuple(frames), excluded


def _positive_log_limits(arrays: Sequence[np.ndarray], label: str) -> tuple[float, float]:
    chunks = [values[values > 0.0] for values in arrays]
    _require(all(chunk.size > 0 for chunk in chunks), f"{label} has an empty frame")
    positive = np.concatenate(chunks)
    lower = float(np.percentile(positive, 1.0))
    upper = float(np.max(positive))
    if lower == upper:
        lower = max(float(np.nextafter(upper, 0.0)), upper * 0.5)
    _require(
        math.isfinite(lower) and math.isfinite(upper) and 0.0 < lower < upper,
        f"unable to construct common logarithmic {label} limits",
    )
    return lower, upper


def common_color_normalization(
    frames: Sequence[FrameReduction],
) -> dict[str, dict[str, object]]:
    """Compute one display normalization per panel across all frames."""
    _require(bool(frames), "cannot normalize an empty frame sequence")
    density_min = min(float(np.min(frame.density)) for frame in frames)
    density_max = max(float(np.max(frame.density)) for frame in frames)
    if density_min == density_max:
        padding = max(abs(density_min) * 0.01, 1.0e-12)
        density_min -= padding
        density_max += padding
    bmag_min, bmag_max = _positive_log_limits(
        [frame.bmag for frame in frames], "magnetic magnitude"
    )
    phase_min, phase_max = _positive_log_limits(
        [frame.phase_density for frame in frames], "phase density"
    )
    return {
        "density_over_rho0": {
            "scale": "linear_global_min_to_max",
            "vmin": density_min,
            "vmax": density_max,
        },
        "bmag_over_b0": {
            "scale": "log_pooled_positive_1st_percentile_to_global_max",
            "vmin": bmag_min,
            "vmax": bmag_max,
        },
        "phase_space": {
            "scale": "log_pooled_positive_1st_percentile_to_global_max",
            "vmin": phase_min,
            "vmax": phase_max,
        },
    }


def _norm(
    normalization: Mapping[str, Mapping[str, object]], name: str
) -> colors.Normalize:
    values = normalization[name]
    vmin = float(values["vmin"])
    vmax = float(values["vmax"])
    if str(values["scale"]).startswith("log_"):
        return colors.LogNorm(vmin=vmin, vmax=vmax)
    return colors.Normalize(vmin=vmin, vmax=vmax)


def render_evolution_frame(
    frame: FrameReduction,
    output_path: Path,
    normalization: Mapping[str, Mapping[str, object]],
    *,
    x_offsets: tuple[float, float] = DEFAULT_X_OFFSETS,
    dpi: int = 300,
) -> None:
    """Render one reduced frame using sequence-wide color limits."""
    _require(output_path.suffix.lower() == ".png", "evolution frames must be PNGs")
    _validate_options(DEFAULT_STEM, x_offsets, dpi)
    figure8._style()
    relative_x_faces = frame.x1_faces - frame.detected_front_x1
    ideal_relative_x = frame.ideal_surface_x1 - frame.detected_front_x1
    fig, axes = plt.subplots(
        3, 1, figsize=(11.2, 7.4), sharex=True, constrained_layout=True
    )
    try:
        panels = (
            (
                frame.density,
                "viridis",
                _norm(normalization, "density_over_rho0"),
                r"$\rho/\rho_0$",
            ),
            (
                frame.bmag,
                "magma",
                _norm(normalization, "bmag_over_b0"),
                r"$|B|/B_0$",
            ),
        )
        for index, (values, cmap, norm, label) in enumerate(panels):
            mesh = axes[index].pcolormesh(
                relative_x_faces,
                frame.x2_faces,
                values,
                shading="flat",
                cmap=cmap,
                norm=norm,
            )
            colorbar = fig.colorbar(mesh, ax=axes[index], pad=0.012)
            colorbar.set_label(label)
            axes[index].set_ylabel(r"$x_2/(c/\omega_{pi})$")
            axes[index].text(
                0.012,
                0.90,
                f"({chr(ord('a') + index)})",
                transform=axes[index].transAxes,
                color="white",
                fontweight="bold",
                fontsize=11,
            )

        phase_plot = np.where(frame.phase_density.T > 0.0, frame.phase_density.T, np.nan)
        phase = axes[2].pcolormesh(
            relative_x_faces,
            frame.log10_chi_edges,
            phase_plot,
            shading="flat",
            cmap="cividis",
            norm=_norm(normalization, "phase_space"),
        )
        colorbar = fig.colorbar(phase, ax=axes[2], pad=0.012)
        colorbar.set_label(r"energy-weighted probability density")
        axes[2].set_ylabel(r"$\log_{10}\chi$")
        axes[2].set_xlabel(r"$(x_1-x_{\rm sh})/(c/\omega_{pi})$")
        axes[2].text(
            0.012,
            0.90,
            "(c)",
            transform=axes[2].transAxes,
            color="white",
            fontweight="bold",
            fontsize=11,
        )
        for axis in axes:
            axis.axvline(ideal_relative_x, color="#f0f0f0", lw=0.8, ls=":")
            axis.axvline(0.0, color="#ffffff", lw=0.9, ls="--")
            axis.set_xlim(*x_offsets)
        fig.suptitle(
            rf"Q011 Section 5.4: $\Omega_0 t={frame.particle_time:.6g}$, "
            rf"cycle {frame.cycle}"
        )
        fig.savefig(output_path, dpi=dpi, bbox_inches="tight")
    finally:
        plt.close(fig)


def _frame_metrics(
    sequence_index: int,
    frame: FrameReduction,
    output: Mapping[str, object],
) -> dict[str, object]:
    return {
        "sequence_index": sequence_index,
        "output_index": frame.pair.output_index,
        "input_index_token": frame.pair.index_token,
        "output": dict(output),
        "cycle": frame.cycle,
        "time": {
            "mhd": frame.mhd_time,
            "particle": frame.particle_time,
            "mesh_time_projection": frame.mesh_time_projection,
        },
        "front": {
            "ideal_surface_x1": frame.ideal_surface_x1,
            "detected_front_x1": frame.detected_front_x1,
            "density_gradient": frame.density_gradient,
            "offset_from_ideal_surface": frame.offset_from_ideal_surface,
        },
        "particle_metrics": {
            "selected_particle_count": frame.selected_particle_count,
            "histogrammed_energy_fraction": frame.histogrammed_energy_fraction,
        },
    }


def _source_dependencies() -> list[dict[str, object]]:
    paths = (
        Path(figure8.__file__),
        Path(figure8.output_primitives.__file__),
        Path(figure8.model.__file__),
        Path(figure8.particle_primitives.__file__),
        Path(pvtk_particles.__file__),
    )
    return [_artifact(path) for path in paths]


def _reject_stale_frames(
    output_dir: Path, stem: str, expected_names: set[str]
) -> None:
    prefix = f"{stem}_frame_"
    existing = {
        path.name
        for path in output_dir.iterdir()
        if path.name.startswith(prefix) and path.suffix.lower() == ".png"
    }
    stale = sorted(existing - expected_names)
    _require(not stale, f"output directory contains stale sequence frames: {stale}")


def make_evolution_sequence(
    run_root: Path,
    output_dir: Path,
    *,
    stem: str = DEFAULT_STEM,
    x_offsets: tuple[float, float] = DEFAULT_X_OFFSETS,
    dpi: int = 300,
) -> list[Path]:
    """Analyze a run root, render PNG frames, and publish a hashed manifest."""
    _validate_options(stem, x_offsets, dpi)
    pairs = discover_product_pairs(run_root)
    input_records = []
    for pair in pairs:
        input_records.append(
            {
                "output_index": pair.output_index,
                "index_token": pair.index_token,
                "mhd_w_bcc": _artifact(pair.mhd_path),
                "prtcl_all": _artifact(pair.particle_path),
            }
        )
    frames, excluded = analyze_sequence(pairs, x_offsets=x_offsets)
    _require(
        len(frames) <= 10**FRAME_DIGITS,
        f"frame count exceeds {FRAME_DIGITS}-digit sequence namespace",
    )
    normalization = common_color_normalization(frames)

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    frame_paths = [
        output_dir / f"{stem}_frame_{index:0{FRAME_DIGITS}d}.png"
        for index in range(len(frames))
    ]
    _reject_stale_frames(output_dir, stem, {path.name for path in frame_paths})
    with tempfile.TemporaryDirectory(
        prefix=f".{stem}-frames-", dir=output_dir
    ) as staging_name:
        staging_dir = Path(staging_name)
        staged_paths = [staging_dir / path.name for path in frame_paths]
        for frame, staged_path in zip(frames, staged_paths):
            render_evolution_frame(
                frame,
                staged_path,
                normalization,
                x_offsets=x_offsets,
                dpi=dpi,
            )
            _require(
                staged_path.is_file() and staged_path.stat().st_size > 0,
                f"renderer did not create a nonempty frame: {staged_path.name}",
            )
        for staged_path, output_path in zip(staged_paths, frame_paths):
            staged_path.replace(output_path)

    output_artifacts = [_artifact(path) for path in frame_paths]
    frame_records = [
        _frame_metrics(index, frame, output_artifacts[index])
        for index, frame in enumerate(frames)
    ]
    rendered_indices = {
        frame.pair.output_index: index for index, frame in enumerate(frames)
    }
    for record in input_records:
        output_index = int(record["output_index"])
        if output_index == int(excluded["output_index"]):
            record["disposition"] = "excluded_t0_no_particle_initial_frame"
        else:
            record["disposition"] = "rendered"
            record["sequence_index"] = rendered_indices[output_index]

    manifest = {
        "record_type": MANIFEST_RECORD_TYPE,
        "schema_version": SCHEMA_VERSION,
        "qualification_effect": "derived_sequence_only_no_claim_closure",
        "generator": _artifact(Path(__file__)),
        "source_dependencies": _source_dependencies(),
        "run_root": str(pairs[0].mhd_path.parent.parent),
        "inputs": input_records,
        "outputs": output_artifacts,
        "excluded_frames": [excluded],
        "frames": frame_records,
        "analysis_parameters": {
            "near_shock_x_offsets": list(x_offsets),
            "shock_relative_x_axis": True,
            "common_color_normalization": normalization,
            "frame_order": "strictly_increasing_particle_time",
            "frame_name_pattern": f"{stem}_frame_%0{FRAME_DIGITS}d.png",
            "dpi": dpi,
            "particle_filter": [
                "cr_source == 1",
                "birth_time >= 45",
                "macro_weight > 0",
            ],
            "cross_cycle_substitution": False,
        },
    }
    manifest_path = output_dir / f"{stem}_manifest.json"
    with tempfile.TemporaryDirectory(
        prefix=f".{stem}-manifest-", dir=output_dir
    ) as staging_name:
        staged_manifest = Path(staging_name) / manifest_path.name
        staged_manifest.write_text(_canonical_json(manifest), encoding="utf-8")
        staged_manifest.replace(manifest_path)
    return [*frame_paths, manifest_path]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--stem", default=DEFAULT_STEM)
    parser.add_argument("--x-offset-min", type=float, default=DEFAULT_X_OFFSETS[0])
    parser.add_argument("--x-offset-max", type=float, default=DEFAULT_X_OFFSETS[1])
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args(argv)
    outputs = make_evolution_sequence(
        args.run_root,
        args.output_dir,
        stem=args.stem,
        x_offsets=(args.x_offset_min, args.x_offset_max),
        dpi=args.dpi,
    )
    print(_canonical_json({"outputs": [str(path) for path in outputs]}), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
