#!/usr/bin/env python3
"""Run R15 without turning its finite-limiter hard-bound diagnostic fatal."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import cgl_lf_stage_i_fast as fast


CASE_ID = "R15"
NODES = 24
RUNS_RELATIVE = Path(
    "runs/mks24-stage-i-fast-relaxed/E03-forcing-policy/24-node"
)
STRICT_OVERRIDE = "mhd/cgl_lf_strict_admissibility=false"


class RelaxedR15Error(RuntimeError):
    """The finite-limiter diagnostic run cannot be launched or continued."""


def configure() -> None:
    fast.FAST_RUNS_RELATIVE = RUNS_RELATIVE
    fast.CASE_NODES[CASE_ID] = NODES
    fast.CASE_SEEDS.pop(CASE_ID, None)
    fast.__file__ = str(Path(__file__).resolve())


def patch_segment(segment: Path) -> None:
    script = segment / "manifest/run.sbatch"
    text = script.read_text(encoding="utf-8")
    needle = " time/tlim=10.0\n"
    if text.count(needle) != 1:
        raise RelaxedR15Error(f"unexpected batch command shape: {script}")
    script.write_text(
        text.replace(needle, f" time/tlim=10.0 {STRICT_OVERRIDE}\n"),
        encoding="utf-8",
    )
    manifest_path = fast.segment_manifest(segment)
    manifest = fast.load_json(manifest_path)
    manifest.update({
        "variant": "finite_limiter_hard_bound_diagnostic_nonfatal",
        "command_line_overrides": [STRICT_OVERRIDE],
        "claim_scope": (
            "R15 finite-rate limiter scan; every hard-bound event remains "
            "retained, but hard-bound occupancy alone is not fatal"
        ),
    })
    fast.write_json(manifest_path, manifest)


def prepare_segment(
    root: Path,
    sequence: int,
    start_time: float,
    restart: Path | None,
    restart_sha: str | None,
) -> Path:
    configure()
    segment = fast.prepare_segment(
        root,
        fast.campaign_cases()[CASE_ID],
        sequence,
        start_time,
        restart,
        restart_sha,
    )
    patch_segment(segment)
    return segment


def launch(root: Path, *, submit: bool) -> Path:
    configure()
    segments = fast.segment_directories(root, CASE_ID)
    if segments:
        raise RelaxedR15Error(f"relaxed R15 lineage already exists: {segments[-1]}")
    segment = prepare_segment(root, 0, 0.0, None, None)
    if submit:
        fast.submit_segment(segment)
    else:
        print(segment)
    return segment


def continue_segment(segment: Path, *, submit: bool) -> Path | None:
    configure()
    manifest = fast.load_json(fast.segment_manifest(segment))
    if str(manifest["case_id"]) != CASE_ID:
        raise RelaxedR15Error(f"expected {CASE_ID} segment: {segment}")
    restart = fast.terminal_product_group(
        segment / "output" / "rst",
        ".rst",
        int(manifest["ranks"]),
    )
    restart_time = float(restart["physical_time"])
    start_time = float(manifest["start_time"])
    target_time = float(manifest["target_time"])
    if not math.isfinite(restart_time) or restart_time <= start_time:
        raise RelaxedR15Error(f"restart did not advance beyond segment start: {segment}")
    if restart_time >= target_time - 1.0e-12:
        print(f"complete {CASE_ID} relaxed diagnostic: t={restart_time}")
        return None
    next_segment = prepare_segment(
        Path(str(manifest["root"])),
        int(manifest["sequence"]) + 1,
        restart_time,
        Path(str(restart["rank_zero"])),
        str(restart["sha256"]),
    )
    if submit:
        fast.submit_segment(next_segment)
    else:
        print(next_segment)
    return next_segment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=fast.DEFAULT_ROOT)
    commands = parser.add_subparsers(dest="command", required=True)
    launch_command = commands.add_parser("launch", help="launch relaxed R15")
    launch_command.add_argument("--submit", action="store_true")
    advance = commands.add_parser("advance", help="continue relaxed R15")
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")
    args = parser.parse_args()

    try:
        if args.command == "launch":
            launch(args.root.resolve(), submit=args.submit)
        else:
            continue_segment(args.segment.resolve(), submit=args.submit)
        return 0
    except (
        RelaxedR15Error,
        fast.FastRunError,
        OSError,
        ValueError,
        KeyError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
