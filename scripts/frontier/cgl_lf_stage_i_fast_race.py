#!/usr/bin/env python3
"""Launch and continue independent high-node races for direct Stage I cases."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import cgl_lf_stage_i_fast as fast


class RaceError(RuntimeError):
    """A race segment cannot be launched or continued."""


def configure(root: Path, runs_relative: Path, case_id: str, nodes: int) -> None:
    if runs_relative.is_absolute() or ".." in runs_relative.parts:
        raise RaceError("--runs-relative must remain beneath the campaign root")
    if nodes < 1:
        raise RaceError("--nodes must be positive")
    fast.FAST_RUNS_RELATIVE = runs_relative
    fast.CASE_NODES[case_id] = nodes
    fast.__file__ = str(Path(__file__).resolve())
    if not (root / runs_relative).resolve().is_relative_to(root.resolve()):
        raise RaceError("race run directory escapes the campaign root")


def configure_from_segment(segment: Path) -> dict[str, object]:
    manifest = fast.load_json(fast.segment_manifest(segment))
    root = Path(str(manifest["root"])).resolve()
    case_id = str(manifest["case_id"])
    nodes = int(manifest["nodes"])
    case_parent = segment.resolve().parent
    try:
        runs_relative = case_parent.parent.relative_to(root)
    except ValueError as error:
        raise RaceError("race segment is outside its declared campaign root") from error
    configure(root, runs_relative, case_id, nodes)
    return manifest


def continue_segment(segment: Path, *, submit: bool) -> Path | None:
    manifest = configure_from_segment(segment)
    restart = fast.terminal_product_group(
        segment / "output" / "rst",
        ".rst",
        int(manifest["ranks"]),
    )
    restart_time = float(restart["physical_time"])
    start_time = float(manifest["start_time"])
    target_time = float(manifest["target_time"])
    if not math.isfinite(restart_time) or restart_time <= start_time:
        raise RaceError(f"restart did not advance beyond segment start: {segment}")
    if restart_time >= target_time - 1.0e-12:
        print(f"complete {manifest['case_id']} race: t={restart_time}")
        return None

    next_segment = fast.prepare_segment(
        Path(str(manifest["root"])),
        fast.campaign_cases()[str(manifest["case_id"])],
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

    launch = commands.add_parser("launch", help="launch an independent race")
    launch.add_argument("case_id", choices=sorted(fast.CASE_NODES))
    launch.add_argument("--nodes", required=True, type=int)
    launch.add_argument(
        "--runs-relative",
        required=True,
        type=Path,
        help="race run root relative to --root, above the RXX directory",
    )
    launch.add_argument("--submit", action="store_true")

    advance = commands.add_parser("advance", help="continue a completed race segment")
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")
    args = parser.parse_args()

    try:
        root = args.root.resolve()
        if args.command == "launch":
            configure(root, args.runs_relative, args.case_id, args.nodes)
            fast.launch_case(root, args.case_id, submit=args.submit)
        else:
            continue_segment(args.segment.resolve(), submit=args.submit)
        return 0
    except (RaceError, fast.FastRunError, OSError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
