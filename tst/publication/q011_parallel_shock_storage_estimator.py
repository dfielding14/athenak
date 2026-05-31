#!/usr/bin/env python3
"""Deterministic raw-PVTK storage estimator for the Q-011 Section 5.4 deck.

This sidecar is an analytical planning tool, not runtime evidence.  It assumes
that shock-injected particles do not escape, that the startup cohort with birth
time below the configured cutoff is removed, and that PVTK snapshots are
retained at t=0 and every configured PVTK cadence through tlim.  It projects
only the binary particle-array payload written by ``vtk_prtcl.cpp``.  ASCII VTK
headers, MHD ``bin`` outputs, restart files, logs, filesystem allocation,
replication, and campaign safety margin are deliberately excluded.

Particle creation is modeled as a continuous swept-mass analytical rate.  The
runtime creates integral particles with a timestep-carried mass reservoir, so
individual runtime snapshots can differ by sub-particle quantization and
timestep-edge effects.  The full-run payload estimate is suitable for storage
planning under the documented no-escape assumption.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
from pathlib import Path
import re
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DECK = (
    REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper.athinput"
)
DEFAULT_GRID_COUNT = 3
DEFAULT_SEED_COUNT = 8

_BLOCK_RE = re.compile(r"<([A-Za-z_][A-Za-z0-9_]*)>")
_PARAMETER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
_INTEGER_RE = re.compile(r"[+-]?[0-9]+")
_OUTPUT_BLOCK_RE = re.compile(r"output[0-9]+")

PVTK_LAYOUT = {
    "position_float_components": 3,
    "integer_scalar_fields": ["gid", "ptag", "species", "cr_source"],
    "float_scalar_fields": [
        "macro_weight",
        "birth_time",
        "deltaf_f0",
        "deltaf_weight",
    ],
    "velocity_float_components": 3,
    "bytes_per_scalar": 4,
}


class EstimatorError(ValueError):
    """Raised when a deck cannot be estimated without making new assumptions."""


def parse_athinput(path: Path) -> dict[str, dict[str, str]]:
    """Parse Athena input blocks strictly and reject malformed or duplicate data."""
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise EstimatorError(f"{path}: cannot read deck as UTF-8 text") from exc

    blocks: dict[str, dict[str, str]] = {}
    current: str | None = None
    for lineno, raw_line in enumerate(lines, 1):
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<"):
            match = _BLOCK_RE.fullmatch(line)
            if match is None:
                raise EstimatorError(f"{path}:{lineno}: malformed block header")
            current = match.group(1)
            if current in blocks:
                raise EstimatorError(f"{path}:{lineno}: duplicate block {current}")
            blocks[current] = {}
            continue
        if current is None or line.count("=") != 1:
            raise EstimatorError(f"{path}:{lineno}: malformed parameter line")
        name, value = (item.strip() for item in line.split("=", 1))
        if _PARAMETER_RE.fullmatch(name) is None or not value:
            raise EstimatorError(f"{path}:{lineno}: malformed parameter")
        if name in blocks[current]:
            raise EstimatorError(f"{path}:{lineno}: duplicate {current}/{name}")
        blocks[current][name] = value
    if not blocks:
        raise EstimatorError(f"{path}: deck contains no input blocks")
    return blocks


def _value(
    blocks: dict[str, dict[str, str]], block: str, parameter: str
) -> str:
    try:
        return blocks[block][parameter]
    except KeyError as exc:
        raise EstimatorError(f"missing required parameter {block}/{parameter}") from exc


def _fraction(
    blocks: dict[str, dict[str, str]], block: str, parameter: str
) -> Fraction:
    value = _value(blocks, block, parameter)
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise EstimatorError(
            f"{block}/{parameter}: expected a finite decimal number"
        ) from exc


def _integer(
    blocks: dict[str, dict[str, str]], block: str, parameter: str
) -> int:
    value = _value(blocks, block, parameter)
    if _INTEGER_RE.fullmatch(value) is None:
        raise EstimatorError(f"{block}/{parameter}: expected an integer")
    return int(value)


def _boolean(
    blocks: dict[str, dict[str, str]], block: str, parameter: str
) -> bool:
    value = _value(blocks, block, parameter)
    if value == "true":
        return True
    if value == "false":
        return False
    raise EstimatorError(f"{block}/{parameter}: expected true or false")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise EstimatorError(message)


def _json_quantity(value: Fraction) -> int | dict[str, int | float]:
    if value.denominator == 1:
        return value.numerator
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "approx": float(value),
    }


def _positive_count(label: str, value: int) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise EstimatorError(f"{label} must be a positive integer")
    return value


def pvtk_bytes_per_particle() -> int:
    """Return the raw binary array payload emitted per cosmic-ray particle."""
    scalar_count = (
        PVTK_LAYOUT["position_float_components"]
        + len(PVTK_LAYOUT["integer_scalar_fields"])
        + len(PVTK_LAYOUT["float_scalar_fields"])
        + PVTK_LAYOUT["velocity_float_components"]
    )
    return scalar_count * PVTK_LAYOUT["bytes_per_scalar"]


def _find_particle_output(
    blocks: dict[str, dict[str, str]]
) -> tuple[str, Fraction]:
    matches = []
    for block, parameters in blocks.items():
        if _OUTPUT_BLOCK_RE.fullmatch(block) is None:
            continue
        if parameters.get("file_type") == "pvtk":
            matches.append(block)
    _require(len(matches) == 1, "deck must contain exactly one PVTK output block")
    block = matches[0]
    _require(
        _value(blocks, block, "variable") == "prtcl_all",
        f"{block}/variable must be prtcl_all",
    )
    cadence = _fraction(blocks, block, "dt")
    _require(cadence > 0, f"{block}/dt must be positive")
    return block, cadence


def _snapshot_times(tlim: Fraction, cadence: Fraction) -> list[Fraction]:
    _require(tlim > 0, "time/tlim must be positive")
    ratio = tlim / cadence
    _require(
        ratio.denominator == 1,
        "PVTK cadence must divide time/tlim exactly for a closed projection",
    )
    return [index * cadence for index in range(ratio.numerator + 1)]


def _retained_duration(
    snapshot_time: Fraction,
    injection_start: Fraction,
    injection_stop: Fraction,
    removal_cutoff: Fraction,
) -> Fraction:
    retained_start = max(injection_start, removal_cutoff)
    retained_stop = min(snapshot_time, injection_stop)
    return max(Fraction(0), retained_stop - retained_start)


def estimate_storage(
    deck_path: Path = DEFAULT_DECK,
    *,
    grid_count: int = DEFAULT_GRID_COUNT,
    seed_count: int = DEFAULT_SEED_COUNT,
) -> dict[str, Any]:
    """Project raw PVTK payload for one run and a grid-by-seed campaign."""
    grid_count = _positive_count("grid_count", grid_count)
    seed_count = _positive_count("seed_count", seed_count)
    blocks = parse_athinput(deck_path)

    _require(
        _value(blocks, "problem", "pgen_name") == "pic_parallel_shock",
        "problem/pgen_name must be pic_parallel_shock",
    )
    _require(
        _value(blocks, "particles", "particle_type") == "cosmic_ray",
        "particles/particle_type must be cosmic_ray",
    )
    _require(
        _fraction(blocks, "particles", "ppc") == 0,
        "particles/ppc must be zero; initial particles are not modeled",
    )
    _require(
        _boolean(blocks, "particles", "pic_enable_2d3v"),
        "particles/pic_enable_2d3v must be true",
    )
    _require(
        _integer(blocks, "mesh", "nx3") == 1,
        "mesh/nx3 must be one for the Section 5.4 2D3V estimator",
    )
    _require(
        _value(blocks, "problem", "ps_shock_speed_model") == "ideal_surface",
        "problem/ps_shock_speed_model must be ideal_surface",
    )
    _require(
        _boolean(blocks, "problem", "ps_enable_injection"),
        "problem/ps_enable_injection must be true",
    )
    _require(
        not _boolean(blocks, "problem", "ps_enable_frame_tracking"),
        "problem/ps_enable_frame_tracking must be false",
    )

    tlim = _fraction(blocks, "time", "tlim")
    _, cadence = _find_particle_output(blocks)
    snapshots = _snapshot_times(tlim, cadence)
    gamma = _fraction(blocks, "mhd", "gamma")
    rho0 = _fraction(blocks, "problem", "ps_rho0")
    u0 = _fraction(blocks, "problem", "ps_u0")
    eta = _fraction(blocks, "problem", "ps_eta")
    injection_start = _fraction(blocks, "problem", "ps_inject_t_start")
    injection_stop = _fraction(blocks, "problem", "ps_inject_t_stop")
    removal_cutoff = _fraction(
        blocks, "problem", "ps_remove_birth_time_before"
    )
    x1min = _fraction(blocks, "mesh", "x1min")
    x1max = _fraction(blocks, "mesh", "x1max")
    x2min = _fraction(blocks, "mesh", "x2min")
    x2max = _fraction(blocks, "mesh", "x2max")
    x3min = _fraction(blocks, "mesh", "x3min")
    x3max = _fraction(blocks, "mesh", "x3max")
    inject_species = _integer(blocks, "problem", "ps_inject_species")
    species_block = f"species{inject_species}"
    species_mass = _fraction(blocks, species_block, "mass")
    deposit_qscale = _fraction(blocks, "particles", "deposit_qscale")

    _require(gamma > 1, "mhd/gamma must be greater than one")
    _require(rho0 > 0, "problem/ps_rho0 must be positive")
    _require(u0 > 0, "problem/ps_u0 must be positive")
    _require(eta >= 0, "problem/ps_eta must be non-negative")
    _require(species_mass > 0, f"{species_block}/mass must be positive")
    _require(deposit_qscale > 0, "particles/deposit_qscale must be positive")
    _require(injection_start >= 0, "problem/ps_inject_t_start must be non-negative")
    _require(injection_stop >= injection_start, "injection stop precedes start")
    _require(removal_cutoff >= injection_start, "startup removal cutoff precedes injection")
    _require(removal_cutoff <= tlim, "startup removal cutoff exceeds time/tlim")
    _require(x1max > x1min, "mesh x1 extent must be positive")
    _require(x2max > x2min, "mesh x2 extent must be positive")
    _require(x3max > x3min, "mesh x3 extent must be positive")

    shock_speed = (gamma - 1) * u0 / 2
    shock_position_at_tlim = x1min + shock_speed * tlim
    _require(
        shock_position_at_tlim <= x1max,
        "ideal shock surface leaves the mesh before time/tlim",
    )
    transverse_area = (x2max - x2min) * (x3max - x3min)
    macro_particle_mass = deposit_qscale * species_mass
    sweep_speed = u0 + shock_speed
    injected_mass_rate = eta * rho0 * sweep_speed * transverse_area
    particle_rate = injected_mass_rate / macro_particle_mass
    bytes_per_particle = pvtk_bytes_per_particle()

    snapshot_reports = []
    payload_per_run = Fraction(0)
    for snapshot in snapshots:
        duration = _retained_duration(
            snapshot,
            injection_start,
            injection_stop,
            removal_cutoff,
        )
        particles = particle_rate * duration
        payload = particles * bytes_per_particle
        payload_per_run += payload
        snapshot_reports.append(
            {
                "time": _json_quantity(snapshot),
                "retained_particle_count_analytical": _json_quantity(particles),
                "raw_pvtk_payload_bytes_before_overhead": _json_quantity(payload),
            }
        )

    final_particles = particle_rate * _retained_duration(
        tlim,
        injection_start,
        injection_stop,
        removal_cutoff,
    )
    campaign_runs = grid_count * seed_count
    campaign_payload = payload_per_run * campaign_runs
    return {
        "schema_version": 1,
        "record_type": "q011_parallel_shock_storage_estimate",
        "qualification_effect": "none_analytical_storage_planning_only",
        "deck": str(deck_path.resolve()),
        "method": "continuous_swept_mass_no_escape_analytical_projection",
        "assumptions": {
            "particle_retention": (
                "no particle escape; all shock-injected particles with birth_time "
                "at or above the configured cutoff remain retained"
            ),
            "cadence": (
                "retain PVTK at t=0 and every configured PVTK dt through time/tlim"
            ),
            "quantization": (
                "continuous swept-mass rate; runtime timestep-edge and carried "
                "mass-reservoir sub-particle quantization are excluded"
            ),
            "payload_scope": (
                "raw PVTK binary particle arrays only; excludes ASCII headers, "
                "MHD bin outputs, restart files, logs, filesystem allocation, "
                "replication, and safety margin"
            ),
        },
        "parsed_inputs": {
            "tlim": _json_quantity(tlim),
            "pvtk_dt": _json_quantity(cadence),
            "startup_removal_birth_time_before": _json_quantity(removal_cutoff),
            "transverse_area": _json_quantity(transverse_area),
            "ideal_surface_speed": _json_quantity(shock_speed),
            "upstream_relative_sweep_speed": _json_quantity(sweep_speed),
            "macro_particle_mass": _json_quantity(macro_particle_mass),
            "injected_mass_rate": _json_quantity(injected_mass_rate),
            "injected_particle_rate": _json_quantity(particle_rate),
        },
        "pvtk_layout": {
            **PVTK_LAYOUT,
            "bytes_per_particle": bytes_per_particle,
        },
        "per_run": {
            "pvtk_snapshot_count": len(snapshot_reports),
            "snapshots": snapshot_reports,
            "retained_particle_count_at_tlim_analytical": _json_quantity(
                final_particles
            ),
            "raw_pvtk_payload_bytes_before_overhead": _json_quantity(
                payload_per_run
            ),
            "raw_pvtk_payload_gb_decimal_before_overhead": (
                float(payload_per_run) / 1.0e9
            ),
        },
        "campaign_parameters": {
            "grid_count": grid_count,
            "seed_count": seed_count,
            "run_count": campaign_runs,
        },
        "campaign": {
            "raw_pvtk_payload_bytes_before_overhead": _json_quantity(
                campaign_payload
            ),
            "raw_pvtk_payload_tb_decimal_before_overhead": (
                float(campaign_payload) / 1.0e12
            ),
        },
    }


def _parse_positive_count(value: str) -> int:
    if _INTEGER_RE.fullmatch(value) is None:
        raise argparse.ArgumentTypeError("expected a positive integer")
    count = int(value)
    if count <= 0:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return count


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, default=DEFAULT_DECK)
    parser.add_argument("--grid-count", type=_parse_positive_count, default=DEFAULT_GRID_COUNT)
    parser.add_argument("--seed-count", type=_parse_positive_count, default=DEFAULT_SEED_COUNT)
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    report = estimate_storage(
        args.deck,
        grid_count=args.grid_count,
        seed_count=args.seed_count,
    )
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
