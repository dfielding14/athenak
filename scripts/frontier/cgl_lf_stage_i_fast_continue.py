#!/usr/bin/env python3
"""Continue direct CGL-LF Stage I runs from their latest complete restart."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import sys

import cgl_lf_stage_i_fast as fast


class ContinueError(RuntimeError):
    """A direct continuation cannot be prepared."""


def continue_segment(segment: Path, *, submit: bool) -> Path | None:
    manifest = fast.load_json(fast.segment_manifest(segment))
    restart = fast.terminal_product_group(
        segment / "output" / "rst",
        ".rst",
        int(manifest["ranks"]),
    )
    restart_time = float(restart["physical_time"])
    start_time = float(manifest["start_time"])
    target_time = float(manifest["target_time"])
    if not math.isfinite(restart_time) or restart_time <= start_time:
        raise ContinueError(f"restart did not advance beyond segment start: {segment}")
    if restart_time >= target_time - 1.0e-12:
        print(f"complete {manifest['case_id']}: t={restart_time}")
        return None

    # Make the generated batch script pin and invoke this lean continuation path.
    fast.__file__ = str(Path(__file__).resolve())
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
    advance = commands.add_parser("advance", help="continue from a completed segment")
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")
    args = parser.parse_args()

    try:
        if args.command == "advance":
            continue_segment(args.segment.resolve(), submit=args.submit)
        return 0
    except (
        ContinueError,
        fast.FastRunError,
        OSError,
        ValueError,
        KeyError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
