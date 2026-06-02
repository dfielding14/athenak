#!/usr/bin/env python3
"""Deterministic storage-review estimator for the Q-011 Section 5.4 deck.

This sidecar is an analytical planning tool, not runtime evidence.  It assumes
that shock-injected particles do not escape, that the startup cohort with birth
time below the configured cutoff is removed, and that PVTK snapshots are
retained at t=0 and every configured PVTK cadence through tlim.  It projects
the binary particle-array payload written by ``vtk_prtcl.cpp`` exactly as the
predecessor sidecar did, then layers a conservative full-campaign storage
planning envelope around that retained calculation.

Particle creation is modeled as a continuous swept-mass analytical rate.  The
runtime creates integral particles with a timestep-carried mass reservoir, so
individual runtime snapshots can differ by sub-particle quantization and
timestep-edge effects.  The full-run payload estimate is suitable for storage
planning under the documented no-escape assumption.  Mesh binary products are
projected from the deck and the binary-output layout.  Restart, log/telemetry,
filesystem-allocation, replication, and safety-margin values are explicit
planning allowances, not manufactured runtime measurements.
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
    REPO_ROOT / "inputs/publication/pic_parallel_shock_section54_paper_vl2_tsc.athinput"
)
DEFAULT_GRID_COUNT = 3
DEFAULT_SEED_COUNT = 8
DEFAULT_LOG_TELEMETRY_ALLOWANCE_BYTES_PER_RUN = 64 * 1024**2
DEFAULT_FILESYSTEM_ALLOCATION_OVERHEAD_FRACTION = Fraction(1, 10)
DEFAULT_REPLICATION_COPY_COUNT = 1
DEFAULT_SAFETY_MARGIN_FRACTION = Fraction(1, 4)

MESH_BIN_FLOAT_BYTES = 4
MESH_BIN_INT32_BYTES = 4
MESH_BIN_LOCATION_REAL_BYTES_PLANNING_ALLOWANCE = 8
MESH_BIN_BLOCK_METADATA_BYTES = (
    10 * MESH_BIN_INT32_BYTES + 6 * MESH_BIN_LOCATION_REAL_BYTES_PLANNING_ALLOWANCE
)
MESH_BIN_PRODUCTS = {
    "rho": "mhd_w_d",
    "bmag": "mhd_bmag",
    "prtcl_jx": "prtcl_jx",
    "j2": "mhd_j2",
}

RESTART_CR_REAL_FIELD_COUNT = 26
RESTART_CR_INTEGER_FIELD_COUNT = 4
RESTART_REAL_BYTES_PLANNING_ALLOWANCE = 8
RESTART_INTEGER_BYTES_PLANNING_ALLOWANCE = 4
RESTART_PARTICLE_BYTES_PLANNING_ALLOWANCE = (
    RESTART_CR_REAL_FIELD_COUNT * RESTART_REAL_BYTES_PLANNING_ALLOWANCE
    + RESTART_CR_INTEGER_FIELD_COUNT * RESTART_INTEGER_BYTES_PLANNING_ALLOWANCE
)
RESTART_MESH_BYTES_PER_CELL_PLANNING_ALLOWANCE = 256
RESTART_FIXED_PUBLICATION_BYTES_PER_CHECKPOINT_ALLOWANCE = 16 * 1024**2

ORION_ONLY_REPLICATION_POLICY = {
    "copy_count": DEFAULT_REPLICATION_COPY_COUNT,
    "selected_destination": "/lustre/orion/ast207/proj-shared/dfielding/PIC",
    "role": "user_selected_sole_bulk_evidence_root",
    "durability_risk": (
        "Orion-only retention is user-directed and does not provide an "
        "institutional or approved off-site durable archive."
    ),
}

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


def _quantity_fraction(value: int | dict[str, int | float]) -> Fraction:
    if isinstance(value, int):
        return Fraction(value)
    return Fraction(value["numerator"], value["denominator"])


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


def _find_mesh_bin_outputs(
    blocks: dict[str, dict[str, str]], tlim: Fraction
) -> list[dict[str, Any]]:
    outputs = []
    remaining = dict(MESH_BIN_PRODUCTS)
    for block, parameters in blocks.items():
        if _OUTPUT_BLOCK_RE.fullmatch(block) is None:
            continue
        if parameters.get("file_type") != "bin":
            continue
        output_id = parameters.get("id")
        if output_id not in remaining:
            raise EstimatorError(f"{block}/id is not a required Section 5.4 bin product")
        expected_variable = remaining.pop(output_id)
        _require(
            _value(blocks, block, "variable") == expected_variable,
            f"{block}/variable must be {expected_variable}",
        )
        _require(
            _value(blocks, block, "ghost_zones") == "false",
            f"{block}/ghost_zones must be false",
        )
        cadence = _fraction(blocks, block, "dt")
        _require(cadence > 0, f"{block}/dt must be positive")
        snapshots = _snapshot_times(tlim, cadence)
        outputs.append(
            {
                "block": block,
                "id": output_id,
                "variable": expected_variable,
                "cadence": cadence,
                "snapshot_count": len(snapshots),
            }
        )
    _require(
        not remaining,
        "deck must contain exactly the required Section 5.4 bin products",
    )
    _require(
        len(outputs) == len(MESH_BIN_PRODUCTS),
        "deck must contain exactly the required Section 5.4 bin products",
    )
    return sorted(outputs, key=lambda item: item["id"])


def _find_restart_output(
    blocks: dict[str, dict[str, str]], tlim: Fraction
) -> tuple[str, Fraction, list[Fraction]]:
    matches = []
    for block, parameters in blocks.items():
        if _OUTPUT_BLOCK_RE.fullmatch(block) is None:
            continue
        if parameters.get("file_type") == "rst":
            matches.append(block)
    _require(len(matches) == 1, "deck must contain exactly one restart output block")
    block = matches[0]
    cadence = _fraction(blocks, block, "dt")
    _require(cadence > 0, f"{block}/dt must be positive")
    ratio = tlim / cadence
    _require(
        ratio.denominator == 1,
        "restart cadence must divide time/tlim exactly for a closed projection",
    )
    return block, cadence, [
        index * cadence for index in range(1, ratio.numerator + 1)
    ]


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


def _mesh_variant_plans(
    blocks: dict[str, dict[str, str]], grid_count: int
) -> list[dict[str, Any]]:
    root_cells = (
        _integer(blocks, "mesh", "nx1"),
        _integer(blocks, "mesh", "nx2"),
        _integer(blocks, "mesh", "nx3"),
    )
    meshblock_cells = (
        _integer(blocks, "meshblock", "nx1"),
        _integer(blocks, "meshblock", "nx2"),
        _integer(blocks, "meshblock", "nx3"),
    )
    _require(all(value > 0 for value in root_cells), "mesh cell counts must be positive")
    _require(
        all(value > 0 for value in meshblock_cells),
        "meshblock cell counts must be positive",
    )
    _require(
        all(root % block == 0 for root, block in zip(root_cells, meshblock_cells)),
        "root mesh cell counts must divide evenly into MeshBlocks",
    )
    _require(
        _value(blocks, "mesh_refinement", "refinement") == "adaptive",
        "mesh_refinement/refinement must be adaptive",
    )
    num_levels = _integer(blocks, "mesh_refinement", "num_levels")
    _require(num_levels == 3, "mesh_refinement/num_levels must be three")
    finest_multiplier = 2 ** (num_levels - 1)
    finest_cells = (
        root_cells[0] * finest_multiplier,
        root_cells[1] * finest_multiplier,
        root_cells[2],
    )

    canonical = [
        (
            "coarse_uniform_dx12",
            root_cells,
            "uniform_root_grid_exact",
        ),
        (
            "three_level_amr_root_dx12_finest_dx3",
            finest_cells,
            "full_domain_finest_level_upper_bound_not_runtime_amr_occupancy",
        ),
        (
            "fine_uniform_dx3",
            finest_cells,
            "uniform_fine_grid_exact_from_frozen_variant_override",
        ),
    ]
    if grid_count != len(canonical):
        canonical = [
            (
                f"noncanonical_grid_{index + 1}_finest_equivalent",
                finest_cells,
                "noncanonical_grid_count_conservative_full_domain_finest_level_upper_bound",
            )
            for index in range(grid_count)
        ]

    plans = []
    for variant, cells_by_axis, method in canonical:
        _require(
            all(cells % block == 0 for cells, block in zip(cells_by_axis, meshblock_cells)),
            f"{variant}: projected cell counts must divide evenly into MeshBlocks",
        )
        meshblocks_by_axis = tuple(
            cells // block for cells, block in zip(cells_by_axis, meshblock_cells)
        )
        plans.append(
            {
                "variant": variant,
                "projection_method": method,
                "cells_by_axis": list(cells_by_axis),
                "cell_count": cells_by_axis[0] * cells_by_axis[1] * cells_by_axis[2],
                "meshblocks_by_axis": list(meshblocks_by_axis),
                "meshblock_count": (
                    meshblocks_by_axis[0]
                    * meshblocks_by_axis[1]
                    * meshblocks_by_axis[2]
                ),
            }
        )
    return plans


def _mesh_bin_payload_for_variant(
    variant: dict[str, Any], outputs: list[dict[str, Any]]
) -> dict[str, Any]:
    product_reports = []
    total = 0
    for output in outputs:
        bytes_per_snapshot = variant["meshblock_count"] * (
            MESH_BIN_BLOCK_METADATA_BYTES + variant["cell_count"] // variant["meshblock_count"]
            * MESH_BIN_FLOAT_BYTES
        )
        payload = bytes_per_snapshot * output["snapshot_count"]
        total += payload
        product_reports.append(
            {
                "id": output["id"],
                "variable": output["variable"],
                "cadence": _json_quantity(output["cadence"]),
                "snapshot_count": output["snapshot_count"],
                "binary_payload_bytes_before_ascii_headers": payload,
            }
        )
    return {
        "products": product_reports,
        "binary_payload_bytes_before_ascii_headers": total,
    }


def _restart_payload_for_variant(
    variant: dict[str, Any],
    checkpoint_times: list[Fraction],
    *,
    particle_rate: Fraction,
    injection_start: Fraction,
    injection_stop: Fraction,
    removal_cutoff: Fraction,
) -> dict[str, Any]:
    checkpoints = []
    total = Fraction(0)
    mesh_allowance = (
        variant["cell_count"] * RESTART_MESH_BYTES_PER_CELL_PLANNING_ALLOWANCE
    )
    for checkpoint in checkpoint_times:
        particles = particle_rate * _retained_duration(
            checkpoint,
            injection_start,
            injection_stop,
            removal_cutoff,
        )
        particle_allowance = particles * RESTART_PARTICLE_BYTES_PLANNING_ALLOWANCE
        payload = (
            particle_allowance
            + mesh_allowance
            + RESTART_FIXED_PUBLICATION_BYTES_PER_CHECKPOINT_ALLOWANCE
        )
        total += payload
        checkpoints.append(
            {
                "time": _json_quantity(checkpoint),
                "retained_particle_count_analytical": _json_quantity(particles),
                "particle_layout_bytes_allowance": _json_quantity(particle_allowance),
                "mesh_state_bytes_allowance": mesh_allowance,
                "fixed_publication_metadata_bytes_allowance": (
                    RESTART_FIXED_PUBLICATION_BYTES_PER_CHECKPOINT_ALLOWANCE
                ),
                "restart_payload_bytes_planning_allowance": _json_quantity(payload),
            }
        )
    return {
        "checkpoint_count": len(checkpoints),
        "checkpoints": checkpoints,
        "restart_payload_bytes_planning_allowance": _json_quantity(total),
    }


def estimate_storage(
    deck_path: Path = DEFAULT_DECK,
    *,
    grid_count: int = DEFAULT_GRID_COUNT,
    seed_count: int = DEFAULT_SEED_COUNT,
) -> dict[str, Any]:
    """Project retained raw PVTK bytes and a conservative campaign envelope."""
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
    mesh_bin_outputs = _find_mesh_bin_outputs(blocks, tlim)
    _, restart_cadence, restart_checkpoints = _find_restart_output(blocks, tlim)
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
    variant_plans = _mesh_variant_plans(blocks, grid_count)
    variant_envelopes = []
    logical_campaign_bytes = Fraction(0)
    for variant in variant_plans:
        mesh_bins = _mesh_bin_payload_for_variant(variant, mesh_bin_outputs)
        restarts = _restart_payload_for_variant(
            variant,
            restart_checkpoints,
            particle_rate=particle_rate,
            injection_start=injection_start,
            injection_stop=injection_stop,
            removal_cutoff=removal_cutoff,
        )
        logical_bytes = (
            payload_per_run
            + mesh_bins["binary_payload_bytes_before_ascii_headers"]
            + _quantity_fraction(restarts["restart_payload_bytes_planning_allowance"])
            + DEFAULT_LOG_TELEMETRY_ALLOWANCE_BYTES_PER_RUN
        )
        logical_campaign_bytes += logical_bytes * seed_count
        variant_envelopes.append(
            {
                **variant,
                "mesh_bin_products": mesh_bins,
                "restart_payload": restarts,
                "raw_pvtk_payload_bytes_before_overhead": _json_quantity(payload_per_run),
                "log_telemetry_bytes_planning_allowance": (
                    DEFAULT_LOG_TELEMETRY_ALLOWANCE_BYTES_PER_RUN
                ),
                "logical_bytes_per_run_before_filesystem_overhead": _json_quantity(
                    logical_bytes
                ),
            }
        )
    filesystem_overhead = (
        logical_campaign_bytes * DEFAULT_FILESYSTEM_ALLOCATION_OVERHEAD_FRACTION
    )
    filesystem_enveloped = logical_campaign_bytes + filesystem_overhead
    replicated = filesystem_enveloped * DEFAULT_REPLICATION_COPY_COUNT
    safety_margin = replicated * DEFAULT_SAFETY_MARGIN_FRACTION
    reservation_envelope = replicated + safety_margin
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
                "the predecessor raw-PVTK binary-array calculation is retained "
                "unchanged; the separate planning envelope adds mesh-bin binary "
                "payload, restart planning allowances, log/telemetry allowance, "
                "filesystem-allocation overhead, Orion-only replication policy, "
                "and safety margin"
            ),
        },
        "parsed_inputs": {
            "tlim": _json_quantity(tlim),
            "pvtk_dt": _json_quantity(cadence),
            "restart_dt": _json_quantity(restart_cadence),
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
        "planning_envelope": {
            "evidence_boundary": (
                "deterministic conservative storage planning only; allowances are "
                "not observed runtime bytes, a filesystem quota reservation, "
                "campaign authorization, or qualification evidence"
            ),
            "mesh_bin_policy": {
                "required_products": dict(MESH_BIN_PRODUCTS),
                "float_bytes_per_cell": MESH_BIN_FLOAT_BYTES,
                "meshblock_metadata_bytes": MESH_BIN_BLOCK_METADATA_BYTES,
                "ascii_headers": (
                    "ASCII headers are not modeled byte-exactly; covered conservatively by the "
                    "filesystem-allocation overhead allowance"
                ),
            },
            "restart_policy": {
                "checkpoint_times": [
                    _json_quantity(checkpoint) for checkpoint in restart_checkpoints
                ],
                "particle_layout": {
                    "real_field_count": RESTART_CR_REAL_FIELD_COUNT,
                    "integer_field_count": RESTART_CR_INTEGER_FIELD_COUNT,
                    "real_bytes_planning_allowance": RESTART_REAL_BYTES_PLANNING_ALLOWANCE,
                    "integer_bytes_planning_allowance": (
                        RESTART_INTEGER_BYTES_PLANNING_ALLOWANCE
                    ),
                    "bytes_per_particle_planning_allowance": (
                        RESTART_PARTICLE_BYTES_PLANNING_ALLOWANCE
                    ),
                    "classification": (
                        "source-layout-derived allowance, not measured checkpoint bytes"
                    ),
                },
                "mesh_bytes_per_cell_planning_allowance": (
                    RESTART_MESH_BYTES_PER_CELL_PLANNING_ALLOWANCE
                ),
                "fixed_publication_bytes_per_checkpoint_allowance": (
                    RESTART_FIXED_PUBLICATION_BYTES_PER_CHECKPOINT_ALLOWANCE
                ),
                "classification": (
                    "restart planning allowance only, not measured runtime bytes; "
                    "runtime AMR occupancy, rank layout, serialized checkpoint bytes, "
                    "manifests, and completion markers must still be measured and retained"
                ),
            },
            "log_telemetry_policy": {
                "bytes_per_run_allowance": DEFAULT_LOG_TELEMETRY_ALLOWANCE_BYTES_PER_RUN,
                "classification": "planning allowance, not observed stdout or telemetry bytes",
            },
            "filesystem_allocation_policy": {
                "overhead_fraction": _json_quantity(
                    DEFAULT_FILESYSTEM_ALLOCATION_OVERHEAD_FRACTION
                ),
                "classification": (
                    "conservative planning allowance for headers, allocation "
                    "granularity, inventories, manifests, markers, and small "
                    "derived records; not a measured Orion filesystem overhead"
                ),
            },
            "replication_policy": dict(ORION_ONLY_REPLICATION_POLICY),
            "safety_margin_policy": {
                "fraction": _json_quantity(DEFAULT_SAFETY_MARGIN_FRACTION),
                "classification": "planning allowance, not measured variance",
            },
            "variants": variant_envelopes,
            "campaign": {
                "logical_bytes_before_filesystem_overhead": _json_quantity(
                    logical_campaign_bytes
                ),
                "filesystem_allocation_overhead_bytes_allowance": _json_quantity(
                    filesystem_overhead
                ),
                "bytes_after_filesystem_overhead": _json_quantity(filesystem_enveloped),
                "replication_copy_count": DEFAULT_REPLICATION_COPY_COUNT,
                "bytes_after_replication_policy": _json_quantity(replicated),
                "safety_margin_bytes_allowance": _json_quantity(safety_margin),
                "reservation_envelope_bytes": _json_quantity(reservation_envelope),
                "reservation_envelope_tb_decimal": float(reservation_envelope) / 1.0e12,
            },
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
