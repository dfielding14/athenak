#!/usr/bin/env python3
"""Advance R14/R15 with short native Athena checkpoint segments.

The finite-limiter diagnostics have shown reproducible turbulence-driver
failures only after long uninterrupted process lifetimes.  This adapter keeps
the corrected campaign's executable, inputs, ranks, physics overrides, and
restart lineage unchanged while requesting a clean Athena stop every 30
minutes.  Each successful stop writes a complete restart and self-submits the
next exact continuation.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys


SCRIPT_PATH = Path(__file__).resolve()
CORRECTED_LAUNCHER = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_corrected.py")
BOUNDED_CASES = ("R14", "R15")
REQUESTED_WALLTIME = "00:40:00"
ATHENA_WALLTIME = "00:30:00"
SEGMENTATION_POLICY = "native_wallclock_checkpoint_30m_v1"


def _load_corrected_launcher():
    name = "cgl_lf_stage_i_fast_corrected_bounded_base"
    spec = importlib.util.spec_from_file_location(name, CORRECTED_LAUNCHER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import corrected launcher: {CORRECTED_LAUNCHER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


corrected = _load_corrected_launcher()
_CORRECTED_PREPARE_SEGMENT = corrected.prepare_segment


class BoundedContinuationError(RuntimeError):
    """A bounded-continuation operational contract failed."""


def require_bounded_case(case_id: str) -> None:
    if case_id not in BOUNDED_CASES:
        raise BoundedContinuationError(
            f"bounded continuation supports only {', '.join(BOUNDED_CASES)}: {case_id}"
        )


def configure(
    root: Path, node_overrides: dict[str, int] | None = None
) -> Path:
    corrected.SCRIPT_PATH = SCRIPT_PATH
    corrected.prepare_segment = prepare_segment
    root = corrected.configure(root, node_overrides)
    corrected.fast.__file__ = str(SCRIPT_PATH)
    corrected.fast.prepare_segment = prepare_segment
    corrected.fast.REQUESTED_WALLTIME = REQUESTED_WALLTIME
    corrected.fast.ATHENA_WALLTIME = ATHENA_WALLTIME
    return root


def prepare_segment(
    root: Path,
    case: dict[str, object],
    sequence: int,
    start_time: float,
    restart: Path | None,
    restart_sha: str | None,
) -> Path:
    case_id = str(case["id"])
    require_bounded_case(case_id)
    if sequence < 1:
        raise BoundedContinuationError(
            f"{case_id} bounded operation requires an authenticated continuation"
        )
    segment = _CORRECTED_PREPARE_SEGMENT(
        root,
        case,
        sequence,
        start_time,
        restart,
        restart_sha,
    )
    manifest_path = corrected.fast.segment_manifest(segment)
    manifest = corrected.fast.load_json(manifest_path)
    manifest.update(
        {
            "runtime_segmentation_policy": SEGMENTATION_POLICY,
            "requested_walltime": REQUESTED_WALLTIME,
            "athena_walltime": ATHENA_WALLTIME,
            "runtime_segmentation_changes_physics": False,
        }
    )
    corrected.fast.write_json(manifest_path, manifest)
    validate_bounded_segment(segment)
    return segment


def validate_bounded_segment(segment: Path) -> dict[str, object]:
    manifest = corrected.validate_segment(segment)
    case_id = str(manifest["case_id"])
    require_bounded_case(case_id)
    exact = {
        "fast_script": str(SCRIPT_PATH),
        "runtime_segmentation_policy": SEGMENTATION_POLICY,
        "requested_walltime": REQUESTED_WALLTIME,
        "athena_walltime": ATHENA_WALLTIME,
        "runtime_segmentation_changes_physics": False,
        "variant": corrected.FINITE_LIMITER_VARIANT,
        "command_line_overrides": [corrected.FINITE_LIMITER_OVERRIDE],
    }
    for key, expected in exact.items():
        if manifest.get(key) != expected:
            raise BoundedContinuationError(
                f"bounded segment {key} differs: {segment}"
            )
    script = (segment / "manifest/run.sbatch").read_text(encoding="utf-8")
    required_fragments = (
        f"#SBATCH -t {REQUESTED_WALLTIME}",
        f"-t {ATHENA_WALLTIME}",
        str(SCRIPT_PATH),
    )
    for fragment in required_fragments:
        if fragment not in script:
            raise BoundedContinuationError(
                f"bounded batch script lacks {fragment!r}: {segment}"
            )
    return manifest


def launch_case(root: Path, case_id: str, *, submit: bool) -> Path | None:
    require_bounded_case(case_id)
    if not corrected.fast.segment_directories(root, case_id):
        raise BoundedContinuationError(
            f"{case_id} has no corrected lineage to continue"
        )
    return corrected.launch_case(root, case_id, submit=submit)


def next_from_segment(segment: Path, *, submit: bool) -> Path | None:
    manifest = corrected.validate_segment(segment)
    require_bounded_case(str(manifest["case_id"]))
    return corrected.next_from_segment(segment, submit=submit)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=corrected.CAMPAIGN_ROOT,
        help="fixed corrected-production campaign root",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    commands.add_parser("preflight", help="authenticate bounded production inputs")

    launch = commands.add_parser(
        "launch", help="continue and optionally submit R14/R15"
    )
    launch.add_argument("cases", nargs="+", choices=BOUNDED_CASES)
    launch.add_argument(
        "--nodes",
        action="append",
        default=[],
        metavar="RXX=N",
    )
    launch.add_argument("--submit", action="store_true")

    advance = commands.add_parser(
        "advance", help="analyze and continue a bounded segment"
    )
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")

    analyze = commands.add_parser("analyze", help="analyze one R14/R15 segment")
    analyze.add_argument("--segment", required=True, type=Path)

    status = commands.add_parser("status", help="show R14/R15 states")
    status.add_argument("cases", nargs="*", default=list(BOUNDED_CASES))

    args = parser.parse_args()
    try:
        root = configure(args.root)
        if args.command == "preflight":
            print(
                json.dumps(
                    corrected.validate_provenance(root, BOUNDED_CASES),
                    indent=2,
                    sort_keys=True,
                )
            )
        elif args.command == "launch":
            cases = list(args.cases)
            root = configure(
                root, corrected.parse_node_overrides(args.nodes, cases)
            )
            corrected.validate_provenance(root, cases)
            for case_id in cases:
                launch_case(root, case_id, submit=args.submit)
        elif args.command == "advance":
            next_from_segment(args.segment.resolve(), submit=args.submit)
        elif args.command == "analyze":
            manifest = corrected.validate_segment(args.segment.resolve())
            require_bounded_case(str(manifest["case_id"]))
            print(
                json.dumps(
                    corrected.analyze_segment(args.segment.resolve()),
                    indent=2,
                    sort_keys=True,
                )
            )
        elif args.command == "status":
            cases = list(args.cases)
            for case_id in cases:
                require_bounded_case(case_id)
            corrected.fast.command_status(root, cases)
        return 0
    except (
        BoundedContinuationError,
        corrected.CorrectedFastError,
        corrected.fast.FastRunError,
        OSError,
        subprocess.SubprocessError,
        ValueError,
        KeyError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


configure(corrected.CAMPAIGN_ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
