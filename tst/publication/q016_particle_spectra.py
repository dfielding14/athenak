"""Bounded species/source/birth-cohort particle spectra from AthenaK PVTK output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

if __package__:
    from .pvtk_particles import ParticleVTKData, read_particle_vtk
else:
    from pvtk_particles import ParticleVTKData, read_particle_vtk

_DELTA_F_DEFINITIONS = {
    "full_f": (
        "Histogram weight is macro_weight. This is the sampled full-f particle "
        "spectrum; deltaf_f0 and deltaf_weight are not applied."
    ),
    "delta_f_perturbation": (
        "Histogram weight is macro_weight * deltaf_weight. The result is the "
        "signed sampled perturbation spectrum only; the analytic background "
        "represented by deltaf_f0 is intentionally not added."
    ),
}

_SOURCE_LABELS = {
    0: "initial",
    1: "shock_injected",
}


def _as_edges(values: Iterable[float], name: str) -> np.ndarray:
    edges = np.asarray(list(values), dtype=np.float64)
    if edges.ndim != 1 or edges.size < 2:
        raise ValueError(f"{name} must contain at least two edges")
    if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0.0):
        raise ValueError(f"{name} must be finite and strictly increasing")
    return edges


def _scalar(data: ParticleVTKData, name: str) -> np.ndarray:
    if name not in data.scalars:
        raise ValueError(f"Particle VTK scalar '{name}' is required")
    values = np.asarray(data.scalars[name])
    if values.ndim != 1 or values.size != data.points.shape[0]:
        raise ValueError(f"Particle VTK scalar '{name}' has the wrong shape")
    return values


def _weights(data: ParticleVTKData, delta_f_semantics: str) -> np.ndarray:
    if delta_f_semantics not in _DELTA_F_DEFINITIONS:
        raise ValueError(f"Unsupported delta-f semantics: {delta_f_semantics}")
    macro_weight = _scalar(data, "macro_weight").astype(np.float64)
    if not np.all(np.isfinite(macro_weight)) or np.any(macro_weight < 0.0):
        raise ValueError("macro_weight must contain finite non-negative values")
    if delta_f_semantics == "full_f":
        return macro_weight

    deltaf_f0 = _scalar(data, "deltaf_f0").astype(np.float64)
    deltaf_weight = _scalar(data, "deltaf_weight").astype(np.float64)
    if not np.all(np.isfinite(deltaf_f0)) or np.any(deltaf_f0 <= 0.0):
        raise ValueError("deltaf_f0 must contain finite positive values")
    if not np.all(np.isfinite(deltaf_weight)):
        raise ValueError("deltaf_weight must contain finite values")
    return macro_weight * deltaf_weight


def _histogram_record(
    axis: np.ndarray,
    weights: np.ndarray,
    mask: np.ndarray,
    edges: np.ndarray,
) -> dict[str, object]:
    selected_axis = axis[mask]
    selected_weights = weights[mask]
    weighted_sum, _ = np.histogram(selected_axis, bins=edges, weights=selected_weights)
    counts, _ = np.histogram(selected_axis, bins=edges)
    underflow = selected_axis < edges[0]
    overflow = selected_axis > edges[-1]
    return {
        "particle_count_total": int(np.count_nonzero(mask)),
        "particle_count_in_bins": int(np.sum(counts)),
        "weighted_sum_in_bins": weighted_sum.astype(float).tolist(),
        "weighted_density_per_bin_width": (weighted_sum / np.diff(edges)).tolist(),
        "weighted_sum_underflow": float(np.sum(selected_weights[underflow])),
        "weighted_sum_overflow": float(np.sum(selected_weights[overflow])),
    }


def build_weighted_species_cohort_spectra(
    data: ParticleVTKData,
    bin_edges: Iterable[float],
    birth_time_edges: Iterable[float],
    *,
    delta_f_semantics: str = "full_f",
) -> dict[str, object]:
    """Return bounded physical-speed spectra resolved by species/source/birth cohort."""
    edges = _as_edges(bin_edges, "bin_edges")
    birth_edges = _as_edges(birth_time_edges, "birth_time_edges")
    velocity = np.asarray(data.vectors.get("vel"), dtype=np.float64)
    if velocity.shape != (data.points.shape[0], 3):
        raise ValueError("Particle VTK vector 'vel' is required with shape (nparticle, 3)")
    if not np.all(np.isfinite(velocity)):
        raise ValueError("Particle VTK velocity contains non-finite values")

    tracking_id = _scalar(data, "ptag").astype(np.int64)
    species = _scalar(data, "species").astype(np.int64)
    source = _scalar(data, "cr_source").astype(np.int64)
    birth_time = _scalar(data, "birth_time").astype(np.float64)
    if np.unique(tracking_id).size != tracking_id.size:
        raise ValueError("ptag must be unique within one particle VTK snapshot")
    if np.any(species < 0) or np.any(source < 0):
        raise ValueError("species and cr_source must contain non-negative values")
    if not np.all(np.isfinite(birth_time)):
        raise ValueError("birth_time must contain finite values")
    if np.any(birth_time < birth_edges[0]) or np.any(birth_time > birth_edges[-1]):
        raise ValueError("birth_time_edges must bound every particle birth_time")

    weights = _weights(data, delta_f_semantics)
    speed = np.linalg.norm(velocity, axis=1)
    all_mask = np.ones(speed.size, dtype=bool)
    groups = []
    for sp in np.unique(species):
        for src in np.unique(source[species == sp]):
            for cohort in range(birth_edges.size - 1):
                mask = (species == sp) & (source == src)
                mask &= birth_time >= birth_edges[cohort]
                if cohort + 1 == birth_edges.size - 1:
                    mask &= birth_time <= birth_edges[cohort + 1]
                else:
                    mask &= birth_time < birth_edges[cohort + 1]
                if not np.any(mask):
                    continue
                groups.append({
                    "species": int(sp),
                    "cr_source": int(src),
                    "cr_source_label": _SOURCE_LABELS.get(int(src), "unknown"),
                    "birth_cohort": int(cohort),
                    "birth_time_min": float(birth_edges[cohort]),
                    "birth_time_max": float(birth_edges[cohort + 1]),
                    "spectrum": _histogram_record(speed, weights, mask, edges),
                })

    return {
        "schema_version": 1,
        "quantity": "physical_speed",
        "bin_edges": edges.tolist(),
        "birth_time_edges": birth_edges.tolist(),
        "normalization": "weighted_sum_in_bins / diff(bin_edges)",
        "delta_f_semantics": delta_f_semantics,
        "delta_f_definition": _DELTA_F_DEFINITIONS[delta_f_semantics],
        "metadata_fields": {
            "tracking_id": "ptag",
            "species": "species",
            "source": "cr_source",
            "birth_time": "birth_time",
            "macro_weight": "macro_weight",
            "delta_f_background": "deltaf_f0",
            "delta_f_weight": "deltaf_weight",
        },
        "all_particles": _histogram_record(speed, weights, all_mask, edges),
        "groups": groups,
    }


def write_weighted_species_cohort_spectra(
    data: ParticleVTKData,
    output_path: str | Path,
    bin_edges: Iterable[float],
    birth_time_edges: Iterable[float],
    *,
    delta_f_semantics: str = "full_f",
) -> dict[str, object]:
    """Write one JSON spectrum record and return the emitted payload."""
    payload = build_weighted_species_cohort_spectra(
        data,
        bin_edges,
        birth_time_edges,
        delta_f_semantics=delta_f_semantics,
    )
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return payload


def _parse_edges(value: str) -> list[float]:
    return [float(item) for item in value.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("pvtk")
    parser.add_argument("output_json")
    parser.add_argument("--bin-edges", required=True, type=_parse_edges)
    parser.add_argument("--birth-time-edges", required=True, type=_parse_edges)
    parser.add_argument(
        "--delta-f-semantics",
        choices=sorted(_DELTA_F_DEFINITIONS),
        default="full_f",
    )
    args = parser.parse_args()
    write_weighted_species_cohort_spectra(
        read_particle_vtk(args.pvtk),
        args.output_json,
        args.bin_edges,
        args.birth_time_edges,
        delta_f_semantics=args.delta_f_semantics,
    )


if __name__ == "__main__":
    main()
