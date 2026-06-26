#!/usr/bin/env python3
"""Measure particle transport and exact cohort retention in a coupled shock."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
from typing import Any

import numpy as np

from tst.publication.pvtk_particles import read_particle_vtk
from tst.publication.q011_section54_particles import reconstruct_chi_from_pvtk_velocity


UPSTREAM_SPEED = 30.0
INJECTION_SPEED_RATIO = 3.16227766017
INJECTION_MOMENTUM = UPSTREAM_SPEED * INJECTION_SPEED_RATIO
_HEADER_RE = re.compile(
    rb"time=\s*([^\s]+)\s+nranks=\s*([0-9]+)\s+cycle=([0-9]+)"
)
_INDEX_RE = re.compile(r"[.]([0-9]+)[.]part[.]vtk$")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _header(path: Path) -> tuple[float, int, int]:
    match = _HEADER_RE.search(path.read_bytes()[:512])
    if match is None:
        raise RuntimeError(f"missing AthenaK particle header: {path}")
    return float(match.group(1)), int(match.group(3)), int(match.group(2))


def _percentiles(values: np.ndarray) -> dict[str, float | None]:
    if values.size == 0:
        return {name: None for name in ("p50", "p90", "p99", "pmax")}
    return {
        "p50": float(np.percentile(values, 50.0)),
        "p90": float(np.percentile(values, 90.0)),
        "p99": float(np.percentile(values, 99.0)),
        "pmax": float(np.max(values)),
    }


def analyze(run_root: Path, basename: str) -> dict[str, Any]:
    paths = sorted(
        (run_root / "pvtk").glob(f"{basename}.prtcl_all.*.part.vtk"),
        key=lambda path: int(_INDEX_RE.search(path.name).group(1)),
    )
    if len(paths) < 3:
        raise RuntimeError(f"shock quicklook needs at least 3 particle snapshots, found {len(paths)}")
    snapshots: list[dict[str, Any]] = []
    particle_states: list[dict[str, np.ndarray]] = []
    for path in paths:
        time, cycle, ranks = _header(path)
        data = read_particle_vtk(path)
        required_scalars = {"cr_source", "birth_time", "macro_weight", "ptag"}
        if not required_scalars.issubset(data.scalars) or "vel" not in data.vectors:
            raise RuntimeError(f"particle field inventory is incomplete: {path}")
        mask = np.asarray(data.scalars["cr_source"]) == 1
        velocity = np.asarray(data.vectors["vel"], dtype=np.float64)[mask]
        chi = reconstruct_chi_from_pvtk_velocity(velocity)
        momentum = UPSTREAM_SPEED * np.sqrt(chi)
        statistics = _percentiles(momentum)
        particle_states.append(
            {
                "ptag": np.asarray(data.scalars["ptag"], dtype=np.int64)[mask],
                "momentum": momentum,
            }
        )
        snapshots.append(
            {
                "path": str(path.relative_to(run_root)),
                "sha256": _sha256(path),
                "time": time,
                "cycle": cycle,
                "rank_count": ranks,
                "particle_count": int(mask.sum()),
                "oldest_birth_time": (
                    float(np.min(data.scalars["birth_time"][mask])) if mask.any() else None
                ),
                **statistics,
                "simulation_frame_p99_over_injection": (
                    statistics["p99"] / INJECTION_MOMENTUM
                    if statistics["p99"] is not None
                    else None
                ),
            }
        )
    nonempty = [index for index, row in enumerate(snapshots) if row["particle_count"] > 0]
    if not nonempty:
        raise RuntimeError("shock run produced no injected particles")
    initial_index = nonempty[0]
    initial_state = particle_states[initial_index]
    initial_order = np.argsort(initial_state["ptag"])
    initial_tags = initial_state["ptag"][initial_order]
    initial_momentum = initial_state["momentum"][initial_order]
    cohort_history = []
    maximum_matched_gain = 0.0
    for index in nonempty:
        state = particle_states[index]
        common, first, current = np.intersect1d(
            initial_tags, state["ptag"], return_indices=True
        )
        delta = 100.0 * (
            state["momentum"][current] / initial_momentum[first] - 1.0
        )
        maximum_matched_gain = max(maximum_matched_gain, float(np.max(delta)))
        cohort_history.append(
            {
                "time": snapshots[index]["time"],
                "cycle": snapshots[index]["cycle"],
                "matched_count": int(common.size),
                "surviving_fraction": float(common.size / initial_tags.size),
                "momentum_change_percentiles": {
                    f"p{quantile}": float(np.percentile(delta, quantile))
                    for quantile in (10, 50, 90, 99)
                },
            }
        )
    transport_ready = bool(
        len(snapshots) >= 3
        and initial_tags.size > 0
        and cohort_history[-1]["matched_count"] > 0
        and snapshots[-1]["time"] > snapshots[initial_index]["time"]
    )
    return {
        "record_type": "q011_shock_transport_engineering_quicklook_v1",
        "run_root": str(run_root),
        "basename": basename,
        "snapshot_count": len(snapshots),
        "injection_momentum": INJECTION_MOMENTUM,
        "initial_cohort_count": int(initial_tags.size),
        "final_cohort_surviving_fraction": cohort_history[-1]["surviving_fraction"],
        "maximum_matched_momentum_gain_percent": maximum_matched_gain,
        "transport_ready": transport_ready,
        "snapshots": snapshots,
        "cohort_history": cohort_history,
        "claim_scope": (
            "compact_uniform_coupled_shock_injection_and_transport_engineering_evidence;"
            "no_diffusive_shock_acceleration_or_section54_qualification_claim"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_root", type=Path)
    parser.add_argument("basename")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--require-transport", action="store_true")
    args = parser.parse_args()
    report = analyze(args.run_root.resolve(strict=True), args.basename)
    payload = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")
    return 0 if (not args.require_transport or report["transport_ready"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
