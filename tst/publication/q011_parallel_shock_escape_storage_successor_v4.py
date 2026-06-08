#!/usr/bin/env python3
"""Q011 storage successor including compact raw particle-escape evidence."""

from __future__ import annotations

import argparse
import copy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from typing import Any

if __package__:
    from . import q011_parallel_shock_storage_estimator as predecessor
else:
    import q011_parallel_shock_storage_estimator as predecessor


DEFAULT_DECK = (
    predecessor.REPO_ROOT
    / "inputs/publication/"
    "pic_parallel_shock_section54_production_science_successor_v1_vl2_tsc.athinput"
)
ESCAPE_EVENT_BYTES = 176
DEFAULT_MAXIMUM_RANKS_PER_RUN = 64
DEFAULT_MAXIMUM_SEGMENTS_PER_RUN = 16
ESCAPE_STREAM_FRAMING_BYTES_PER_RANK_SEGMENT_ALLOWANCE = 4096
PRODUCTION_MESH_BIN_COMPONENTS = {
    "mhd_w_bcc": ("mhd_w_bcc", 8),
    "prtcl_rho": ("prtcl_rho", 1),
    "prtcl_jx": ("prtcl_jx", 1),
    "prtcl_jy": ("prtcl_jy", 1),
    "prtcl_jz": ("prtcl_jz", 1),
    "mhd_j2": ("mhd_j2", 1),
}
PHYSICAL_BLOCKS = (
    "mesh",
    "meshblock",
    "mesh_refinement",
    "time",
    "mhd",
    "coord",
    "particles",
    "species0",
    "problem",
)


class EscapeStorageError(ValueError):
    """Reject an escape-storage projection requiring unstated assumptions."""


def _positive_integer(value: int, *, label: str) -> int:
    if type(value) is not int or value <= 0:
        raise EscapeStorageError(f"{label} must be a positive integer")
    return value


def _ceil(value: Fraction) -> int:
    return -(-value.numerator // value.denominator)


def _canonical_sha256(value: object) -> str:
    payload = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _production_mesh_outputs(
    blocks: dict[str, dict[str, str]], tlim: Fraction
) -> list[dict[str, Any]]:
    observed = {}
    hst_count = 0
    for block, parameters in blocks.items():
        if predecessor._OUTPUT_BLOCK_RE.fullmatch(block) is None:
            continue
        file_type = parameters.get("file_type")
        if file_type == "bin":
            output_id = parameters.get("id")
            if output_id not in PRODUCTION_MESH_BIN_COMPONENTS or output_id in observed:
                raise EscapeStorageError(f"{block}: unexpected production bin output")
            variable, components = PRODUCTION_MESH_BIN_COMPONENTS[output_id]
            if (
                parameters.get("variable") != variable
                or parameters.get("ghost_zones") != "false"
            ):
                raise EscapeStorageError(f"{block}: production bin semantics drifted")
            cadence = predecessor._fraction(blocks, block, "dt")
            observed[output_id] = {
                "id": output_id,
                "variable": variable,
                "components": components,
                "cadence": cadence,
                "snapshot_count": len(predecessor._snapshot_times(tlim, cadence)),
            }
        elif file_type == "hst":
            hst_count += 1
    if set(observed) != set(PRODUCTION_MESH_BIN_COMPONENTS) or hst_count != 1:
        raise EscapeStorageError("production output inventory is incomplete or drifted")
    predecessor._find_particle_output(blocks)
    predecessor._find_restart_output(blocks, tlim)
    return [observed[key] for key in sorted(observed)]


def _production_mesh_payload(
    variant: dict[str, Any], outputs: list[dict[str, Any]]
) -> dict[str, Any]:
    products = []
    total = 0
    cells_per_block = variant["cell_count"] // variant["meshblock_count"]
    for output in outputs:
        bytes_per_snapshot = variant["meshblock_count"] * (
            predecessor.MESH_BIN_BLOCK_METADATA_BYTES
            + cells_per_block
            * predecessor.MESH_BIN_FLOAT_BYTES
            * output["components"]
        )
        payload = bytes_per_snapshot * output["snapshot_count"]
        total += payload
        products.append(
            {
                "id": output["id"],
                "variable": output["variable"],
                "component_count": output["components"],
                "cadence": predecessor._json_quantity(output["cadence"]),
                "snapshot_count": output["snapshot_count"],
                "binary_payload_bytes_before_ascii_headers": payload,
            }
        )
    return {
        "products": products,
        "binary_payload_bytes_before_ascii_headers": total,
    }


def estimate_storage_with_escape_evidence(
    deck_path: Path = DEFAULT_DECK,
    *,
    grid_count: int = predecessor.DEFAULT_GRID_COUNT,
    seed_count: int = predecessor.DEFAULT_SEED_COUNT,
    maximum_ranks_per_run: int = DEFAULT_MAXIMUM_RANKS_PER_RUN,
    maximum_segments_per_run: int = DEFAULT_MAXIMUM_SEGMENTS_PER_RUN,
) -> dict[str, Any]:
    """Extend the historical envelope with a fail-safe all-particles-escape bound."""
    maximum_ranks_per_run = _positive_integer(
        maximum_ranks_per_run, label="maximum_ranks_per_run"
    )
    maximum_segments_per_run = _positive_integer(
        maximum_segments_per_run, label="maximum_segments_per_run"
    )
    try:
        base = predecessor.estimate_storage(
            predecessor.DEFAULT_DECK,
            grid_count=grid_count,
            seed_count=seed_count,
        )
        blocks = predecessor.parse_athinput(deck_path)
        predecessor_blocks = predecessor.parse_athinput(predecessor.DEFAULT_DECK)
        if any(
            blocks.get(name) != predecessor_blocks.get(name)
            for name in PHYSICAL_BLOCKS
        ):
            raise EscapeStorageError(
                "production deck physical controls differ from historical estimator deck"
            )
        tlim = predecessor._fraction(blocks, "time", "tlim")
        production_outputs = _production_mesh_outputs(blocks, tlim)
        start = predecessor._fraction(blocks, "problem", "ps_inject_t_start")
        stop = predecessor._fraction(blocks, "problem", "ps_inject_t_stop")
        particle_rate = predecessor._quantity_fraction(
            base["parsed_inputs"]["injected_particle_rate"]
        )
    except predecessor.EstimatorError as error:
        raise EscapeStorageError(str(error)) from error
    duration = max(Fraction(0), min(tlim, stop) - start)
    maximum_event_count = _ceil(particle_rate * duration)
    event_payload = maximum_event_count * ESCAPE_EVENT_BYTES
    framing = (
        maximum_ranks_per_run
        * maximum_segments_per_run
        * ESCAPE_STREAM_FRAMING_BYTES_PER_RANK_SEGMENT_ALLOWANCE
    )
    per_run = event_payload + framing
    run_count = int(base["campaign_parameters"]["run_count"])
    campaign_escape = per_run * run_count

    result = copy.deepcopy(base)
    result["schema_version"] = 2
    result["record_type"] = "q011_parallel_shock_escape_storage_estimate_successor_v4"
    result["deck"] = str(deck_path.resolve())
    result["method"] = (
        "historical_continuous_swept_mass_projection_plus_all_injected_particles_"
        "escape_binary_evidence_bound"
    )
    result["predecessor_estimate_sha256"] = _canonical_sha256(base)
    result["assumptions"]["escape_evidence"] = (
        "storage is reserved for one 176-byte raw event for every analytically "
        "injected particle, even though only particles physically destroyed at "
        "outer x1 are written; fixed stream framing is separately bounded for "
        "the configured maximum ranks and restart segments"
    )
    result["escape_evidence"] = {
        "binary_schema": "Q011ESC1/Q011END1",
        "event_bytes": ESCAPE_EVENT_BYTES,
        "maximum_injected_particle_count_per_run": maximum_event_count,
        "maximum_event_payload_bytes_per_run": event_payload,
        "maximum_ranks_per_run": maximum_ranks_per_run,
        "maximum_segments_per_run": maximum_segments_per_run,
        "framing_bytes_per_rank_segment_allowance": (
            ESCAPE_STREAM_FRAMING_BYTES_PER_RANK_SEGMENT_ALLOWANCE
        ),
        "maximum_framing_bytes_per_run": framing,
        "maximum_escape_evidence_bytes_per_run": per_run,
        "maximum_escape_evidence_gb_decimal_per_run": per_run / 1.0e9,
        "maximum_escape_evidence_bytes_campaign": campaign_escape,
        "maximum_escape_evidence_tb_decimal_campaign": campaign_escape / 1.0e12,
        "classification": (
            "source-schema-derived conservative allowance, not observed escape count "
            "or measured filesystem allocation"
        ),
    }

    logical = Fraction(0)
    production_mesh_increment = 0
    for variant in result["planning_envelope"]["variants"]:
        old_mesh = variant["mesh_bin_products"][
            "binary_payload_bytes_before_ascii_headers"
        ]
        production_mesh = _production_mesh_payload(variant, production_outputs)
        historical_logical = predecessor._quantity_fraction(
            variant["logical_bytes_per_run_before_filesystem_overhead"]
        )
        variant_logical = (
            historical_logical
            - old_mesh
            + production_mesh["binary_payload_bytes_before_ascii_headers"]
            + per_run
        )
        production_mesh_increment += (
            production_mesh["binary_payload_bytes_before_ascii_headers"] - old_mesh
        ) * seed_count
        variant["historical_mesh_bin_products"] = variant["mesh_bin_products"]
        variant["mesh_bin_products"] = production_mesh
        variant["escape_evidence_bytes_per_run_allowance"] = per_run
        variant["logical_bytes_per_run_before_filesystem_overhead"] = (
            predecessor._json_quantity(variant_logical)
        )
        logical += variant_logical * seed_count

    base_campaign = predecessor._quantity_fraction(
        base["planning_envelope"]["campaign"][
            "logical_bytes_before_filesystem_overhead"
        ]
    )
    filesystem_overhead = (
        logical * predecessor.DEFAULT_FILESYSTEM_ALLOCATION_OVERHEAD_FRACTION
    )
    filesystem_enveloped = logical + filesystem_overhead
    replicated = filesystem_enveloped * predecessor.DEFAULT_REPLICATION_COPY_COUNT
    safety_margin = replicated * predecessor.DEFAULT_SAFETY_MARGIN_FRACTION
    reservation = replicated + safety_margin
    result["planning_envelope"]["campaign"] = {
        "predecessor_logical_bytes_before_escape_evidence": (
            predecessor._json_quantity(base_campaign)
        ),
        "production_mesh_bin_increment_bytes": production_mesh_increment,
        "escape_evidence_bytes_allowance": campaign_escape,
        "logical_bytes_before_filesystem_overhead": predecessor._json_quantity(
            logical
        ),
        "filesystem_allocation_overhead_bytes_allowance": (
            predecessor._json_quantity(filesystem_overhead)
        ),
        "bytes_after_filesystem_overhead": predecessor._json_quantity(
            filesystem_enveloped
        ),
        "replication_copy_count": predecessor.DEFAULT_REPLICATION_COPY_COUNT,
        "bytes_after_replication_policy": predecessor._json_quantity(replicated),
        "safety_margin_bytes_allowance": predecessor._json_quantity(safety_margin),
        "reservation_envelope_bytes": predecessor._json_quantity(reservation),
        "reservation_envelope_tb_decimal": float(reservation) / 1.0e12,
    }
    return result


def _positive_argument(value: str) -> int:
    try:
        parsed = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("expected a positive integer") from error
    if parsed <= 0:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return parsed


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deck", type=Path, default=DEFAULT_DECK)
    parser.add_argument("--grid-count", type=_positive_argument, default=3)
    parser.add_argument("--seed-count", type=_positive_argument, default=8)
    parser.add_argument(
        "--maximum-ranks-per-run",
        type=_positive_argument,
        default=DEFAULT_MAXIMUM_RANKS_PER_RUN,
    )
    parser.add_argument(
        "--maximum-segments-per-run",
        type=_positive_argument,
        default=DEFAULT_MAXIMUM_SEGMENTS_PER_RUN,
    )
    args = parser.parse_args()
    print(
        json.dumps(
            estimate_storage_with_escape_evidence(
                args.deck,
                grid_count=args.grid_count,
                seed_count=args.seed_count,
                maximum_ranks_per_run=args.maximum_ranks_per_run,
                maximum_segments_per_run=args.maximum_segments_per_run,
            ),
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
