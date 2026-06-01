#!/usr/bin/env python3
"""Standalone Sun & Bai paper_smooth receiver-resolution TSC oracle.

The refinement interface is at x=0.  Fine receivers live to its left and
coarse receivers live to its right.  Each raw contribution uses the receiving
cell's resolution.  The cross-interface totals are intentionally not
renormalized.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from fractions import Fraction
import json
from typing import Iterable, Sequence


@dataclass(frozen=True)
class Receiver:
    """One cell center and the resolution used for its raw TSC weight."""

    level: str
    center: Fraction
    dx: Fraction


FINE_RECEIVER_CENTERS = (
    Fraction("-1.75"),
    Fraction("-1.25"),
    Fraction("-0.75"),
    Fraction("-0.25"),
)
COARSE_RECEIVER_CENTERS = (
    Fraction("0.5"),
    Fraction("1.5"),
    Fraction("2.5"),
)
FINE_DX = Fraction("0.5")
COARSE_DX = Fraction("1.0")

RECEIVERS = tuple(
    Receiver("fine", center, FINE_DX) for center in FINE_RECEIVER_CENTERS
) + tuple(
    Receiver("coarse", center, COARSE_DX)
    for center in COARSE_RECEIVER_CENTERS
)

PARTICLES = {
    "a": Fraction("-1.45"),
    "b": Fraction("-0.55"),
    "c": Fraction("-0.25"),
    "d": Fraction("0.15"),
    "e": Fraction("0.8"),
    "f": Fraction("1.5"),
}

# Explicit frozen rows are ordered like RECEIVERS.  Keeping these independent
# of the evaluator makes the script useful as an implementation oracle.
EXPECTED_RAW_WEIGHTS = {
    "a": (
        Fraction(81, 200),
        Fraction(59, 100),
        Fraction(1, 200),
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(0),
    ),
    "b": (
        Fraction(0),
        Fraction(1, 200),
        Fraction(59, 100),
        Fraction(81, 200),
        Fraction(81, 800),
        Fraction(0),
        Fraction(0),
    ),
    "c": (
        Fraction(0),
        Fraction(0),
        Fraction(1, 8),
        Fraction(3, 4),
        Fraction(9, 32),
        Fraction(0),
        Fraction(0),
    ),
    "d": (
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(49, 200),
        Fraction(251, 400),
        Fraction(9, 800),
        Fraction(0),
    ),
    "e": (
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(33, 50),
        Fraction(8, 25),
        Fraction(0),
    ),
    "f": (
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(0),
        Fraction(1, 8),
        Fraction(3, 4),
        Fraction(1, 8),
    ),
}

EXPECTED_RAW_TOTALS = {
    "a": Fraction(1),
    "b": Fraction(881, 800),
    "c": Fraction(37, 32),
    "d": Fraction(707, 800),
    "e": Fraction(49, 50),
    "f": Fraction(1),
}

NON_UNIT_TOTAL_PARTICLES = ("b", "c", "d", "e")
EXPECTED_TENSOR_PRODUCT_WEIGHT = Fraction(234171, 4000000)


class OracleError(RuntimeError):
    """Raised when a frozen paper_smooth oracle invariant drifts."""


def _fraction(value: Fraction | float | int | str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, float):
        return Fraction(str(value))
    return Fraction(value)


def raw_tsc_weight(
    particle_coordinate: Fraction | float | int | str,
    receiver_center: Fraction | float | int | str,
    receiver_dx: Fraction | float | int | str,
) -> Fraction:
    """Return one raw 1D TSC weight using the receiver's resolution."""

    coordinate = _fraction(particle_coordinate)
    center = _fraction(receiver_center)
    dx = _fraction(receiver_dx)
    if dx <= 0:
        raise ValueError("receiver_dx must be positive")

    distance = abs(coordinate - center) / dx
    if distance < Fraction(1, 2):
        return Fraction(3, 4) - distance * distance
    if distance < Fraction(3, 2):
        return Fraction(1, 2) * (Fraction(3, 2) - distance) ** 2
    return Fraction(0)


def receiver_raw_weights(
    particle_coordinate: Fraction | float | int | str,
    receivers: Iterable[Receiver] = RECEIVERS,
) -> tuple[Fraction, ...]:
    """Return unnormalized weights for the supplied refinement-interface cells."""

    return tuple(
        raw_tsc_weight(particle_coordinate, receiver.center, receiver.dx)
        for receiver in receivers
    )


def tensor_product_weight(
    particle_coordinates: Sequence[Fraction | float | int | str],
    receiver_centers: Sequence[Fraction | float | int | str],
    receiver_dx: Sequence[Fraction | float | int | str],
) -> Fraction:
    """Return the multidimensional raw TSC tensor-product weight."""

    if not (
        len(particle_coordinates) == len(receiver_centers) == len(receiver_dx)
    ):
        raise ValueError("tensor-product coordinates and dx must have equal length")
    weight = Fraction(1)
    for coordinate, center, dx in zip(
        particle_coordinates, receiver_centers, receiver_dx
    ):
        weight *= raw_tsc_weight(coordinate, center, dx)
    return weight


def _fraction_record(value: Fraction) -> dict[str, object]:
    return {
        "fraction": str(value),
        "float": float(value),
    }


def _validate_support_bounds() -> None:
    for dx in (FINE_DX, COARSE_DX):
        if raw_tsc_weight(0, 0, dx) != Fraction(3, 4):
            raise OracleError("TSC center weight drifted")
        if raw_tsc_weight(0, dx / 2, dx) != Fraction(1, 2):
            raise OracleError("TSC inner piece boundary drifted")
        if raw_tsc_weight(0, 3 * dx / 2, dx) != 0:
            raise OracleError("TSC support must exclude distance=1.5*dx")
        if raw_tsc_weight(0, 2 * dx, dx) != 0:
            raise OracleError("TSC support must exclude distance>1.5*dx")


def _validate_tensor_product() -> None:
    particle = (PARTICLES["b"], Fraction("0.2"), Fraction("-0.4"))
    center = (Fraction("-0.25"), Fraction(0), Fraction(0))
    dx = (FINE_DX, FINE_DX, FINE_DX)
    if tensor_product_weight(particle, center, dx) != EXPECTED_TENSOR_PRODUCT_WEIGHT:
        raise OracleError("multidimensional TSC weight is not a tensor product")
    outside = (PARTICLES["b"], Fraction("0.2"), Fraction("0.75"))
    if tensor_product_weight(outside, center, dx) != 0:
        raise OracleError("tensor-product support must vanish outside any axis")


def validate_oracle() -> None:
    """Validate the frozen coefficient grid and concise shape invariants."""

    for label, coordinate in PARTICLES.items():
        measured = receiver_raw_weights(coordinate)
        expected = EXPECTED_RAW_WEIGHTS[label]
        if measured != expected:
            raise OracleError(
                f"raw receiver-resolution weights drifted for particle {label}: "
                f"measured={measured!r} expected={expected!r}"
            )
        total = sum(measured, Fraction(0))
        if total != EXPECTED_RAW_TOTALS[label]:
            raise OracleError(
                f"raw unnormalized total drifted for particle {label}: "
                f"measured={total} expected={EXPECTED_RAW_TOTALS[label]}"
            )

    if tuple(
        label for label in PARTICLES if EXPECTED_RAW_TOTALS[label] != 1
    ) != NON_UNIT_TOTAL_PARTICLES:
        raise OracleError("frozen non-unit raw-total cases drifted")

    _validate_support_bounds()
    _validate_tensor_product()


def build_report() -> dict[str, object]:
    """Return the checked standalone oracle report."""

    validate_oracle()
    return {
        "oracle": "sun_bai_paper_smooth_receiver_resolution_raw_tsc",
        "status": "pass",
        "policy": {
            "interface_x": 0.0,
            "resolution_rule": "evaluate_each_weight_at_receiver_resolution",
            "renormalize_cross_interface_totals": False,
            "shape": "raw_tsc_tensor_product",
            "support_bound_receiver_dx": "abs(x_receiver-x_particle) < 1.5*dx",
        },
        "receivers": [
            {
                "level": receiver.level,
                "center": float(receiver.center),
                "dx": float(receiver.dx),
            }
            for receiver in RECEIVERS
        ],
        "particles": {
            label: {
                "position": float(PARTICLES[label]),
                "raw_weights": [
                    _fraction_record(weight)
                    for weight in receiver_raw_weights(PARTICLES[label])
                ],
                "raw_total": _fraction_record(EXPECTED_RAW_TOTALS[label]),
                "non_unit_raw_total": label in NON_UNIT_TOTAL_PARTICLES,
            }
            for label in PARTICLES
        },
        "checks": {
            "frozen_receiver_weight_grid": "pass",
            "raw_non_unit_totals_preserved": "pass",
            "support_bound": "pass",
            "tensor_product": "pass",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--indent",
        type=int,
        default=2,
        help="JSON indentation level (default: 2)",
    )
    args = parser.parse_args()
    print(json.dumps(build_report(), indent=args.indent, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
