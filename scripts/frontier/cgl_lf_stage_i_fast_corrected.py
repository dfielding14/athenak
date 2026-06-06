#!/usr/bin/env python3
"""Launch corrected-production active CGL Stage I cases directly on Frontier.

This is an isolated production adapter over ``cgl_lf_stage_i_fast.py``.  It
pins the ppar^2 fast-discriminant source and executable, forces every initial
case lineage to start from t=0, and permits continuations only from this
corrected campaign's own output.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
from typing import Iterable


SCRIPT_PATH = Path(__file__).resolve()
BASE_LAUNCHER = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast.py")


def _load_base_launcher():
    name = "cgl_lf_stage_i_fast_corrected_base"
    spec = importlib.util.spec_from_file_location(name, BASE_LAUNCHER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import direct-fast base launcher: {BASE_LAUNCHER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fast = _load_base_launcher()

CAMPAIGN_ID = "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
CAMPAIGN_ROOT = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
)
FROZEN_SOURCE = Path(
    "/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-corrected-0c406312f"
)
SOURCE_REVISION = "0c406312fa35d5e1c7041d80333b0ea24b0127ae"
MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
ATHENA = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/build/"
    "frontier-hip-0c406312fa35-cpe25.09-cce20-rocm6.4.2/src/athena"
)
ATHENA_SHA256 = "0f032379c4d829fc86353cc734a716f6ef6d6eda8951a5b79a3a97ad671e9b4b"
IDENTITY = CAMPAIGN_ROOT / "campaign-identity.json"
IDENTITY_SHA256 = "fad99d9651e6da155bb5dfae1b81a8675cce232fb17e31a36e9be36d62810643"
RUNS_RELATIVE = Path("runs/mks24-stage-i-fast/E03-forcing-policy")
LOGS_RELATIVE = Path("logs/slurm-fast")

ACTIVE_CASES = (
    "R02",
    "R03",
    "R04",
    "R05",
    "R10",
    "R11",
    "R12",
    "R13",
    "R14",
    "R15",
    "R16",
    "R17",
)
FINITE_LIMITER_DIAGNOSTIC_CASES = frozenset({"R14", "R15"})
FINITE_LIMITER_OVERRIDE = "mhd/cgl_lf_strict_admissibility=false"
FINITE_LIMITER_VARIANT = "finite_limiter_hard_bound_diagnostic_nonfatal"

# These defaults use the successful high-node campaign profiles while keeping
# rank counts below each input's meshblock count.
DEFAULT_CASE_NODES = {
    "R02": 24,
    "R03": 24,
    "R04": 24,
    "R05": 24,
    "R10": 4,
    "R11": 24,
    "R12": 24,
    "R13": 4,
    "R14": 24,
    "R15": 24,
    "R16": 3,
    "R17": 128,
}
MAX_USEFUL_NODES = {
    **{case_id: 27 for case_id in ACTIVE_CASES if case_id not in {"R16", "R17"}},
    "R16": 3,
    "R17": 216,
}
CASE_INPUT_SHA256 = {
    "R02": "c310509aa1638418bab7117427e8cca06e5a651c60e15212e4763baad5101206",
    "R03": "997f449abb3c2e4d1de50509efffa127c6230e3ecb2632612200010b8fa6c0b0",
    "R04": "571eea2ccec5d069ccb1b49d132ba4c5b8bdac7ce15ee8d8a04c92327c453137",
    "R05": "6527a2072105d515287701c904c3af2cebb07e10ab0b45c3fb41f29d3467622c",
    "R10": "93d7ed9b0846019f59bce34feeb9b22b098afffdfc2dbaa59aaf5913655cd705",
    "R11": "ab0dba23f80ea2d6c173ecbb724dbcb332b13d3b1751b8776ad05313f9ccaf6c",
    "R12": "98ddea4b4f7fec18cc40abdbf5f7c8ba5b583a91f84f4dae00e1411f23d42e7c",
    "R13": "a190129ee34c46a4fe83a392a19120a58b9a9a4064742ac58511441370c59c35",
    "R14": "2b8d5837f8a7f3070ca2ef56f8b44e53d048839918185a8eef155cb084736578",
    "R15": "9a698a60bef4c558ccee4635d69c3acbf3fea478401bd633d09bbc0526d943d8",
    "R16": "c0ac4b54248e8f8dfb0f5fd34c0cfb4414b5330529cbf2836961c5277af3f2d1",
    "R17": "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04",
}

_BASE_PREPARE_SEGMENT = fast.prepare_segment
_BASE_BATCH_SCRIPT_TEXT = fast.batch_script_text


class CorrectedFastError(RuntimeError):
    """A corrected-production provenance or isolation contract failed."""


def require_campaign_root(root: Path) -> Path:
    actual = root.expanduser().resolve()
    expected = CAMPAIGN_ROOT.expanduser().resolve()
    if actual != expected:
        raise CorrectedFastError(
            f"campaign root must be exactly {expected}; received {actual}"
        )
    return actual


def require_beneath(path: Path, parent: Path, label: str) -> Path:
    resolved = path.expanduser().resolve()
    try:
        relative = resolved.relative_to(parent.expanduser().resolve())
    except ValueError as error:
        raise CorrectedFastError(f"{label} escapes corrected campaign root: {resolved}") from error
    if not relative.parts:
        raise CorrectedFastError(f"{label} must be below corrected campaign root: {resolved}")
    return resolved


def git_read(arguments: list[str]) -> str:
    completed = subprocess.run(
        ["/usr/bin/git", *arguments],
        check=True,
        text=True,
        capture_output=True,
    )
    return completed.stdout.strip()


def require_source_revision() -> None:
    if not FROZEN_SOURCE.is_dir():
        raise CorrectedFastError(f"missing corrected frozen source: {FROZEN_SOURCE}")
    head = git_read(["-C", str(FROZEN_SOURCE), "rev-parse", "--verify", "HEAD"])
    if head != SOURCE_REVISION:
        raise CorrectedFastError(
            f"corrected frozen source revision mismatch: {head} != {SOURCE_REVISION}"
        )
    status = git_read(
        ["-C", str(FROZEN_SOURCE), "status", "--porcelain", "--untracked-files=all"]
    )
    if status:
        raise CorrectedFastError(f"corrected frozen source is dirty: {FROZEN_SOURCE}")


def validate_node_profile(case_id: str, nodes: int) -> None:
    if case_id not in DEFAULT_CASE_NODES:
        raise CorrectedFastError(f"unsupported corrected active case: {case_id}")
    maximum = MAX_USEFUL_NODES[case_id]
    if nodes < 1 or nodes > maximum:
        raise CorrectedFastError(
            f"{case_id} nodes must be in [1, {maximum}] to keep ranks useful"
        )


def parse_node_overrides(values: Iterable[str], selected_cases: Iterable[str]) -> dict[str, int]:
    selected = set(selected_cases)
    result: dict[str, int] = {}
    for value in values:
        match = re.fullmatch(r"(R\d{2})=([1-9]\d*)", value)
        if match is None:
            raise CorrectedFastError(f"invalid --nodes override: {value}; expected RXX=N")
        case_id, raw_nodes = match.groups()
        if case_id not in selected:
            raise CorrectedFastError(f"--nodes override is not selected for launch: {case_id}")
        nodes = int(raw_nodes)
        validate_node_profile(case_id, nodes)
        result[case_id] = nodes
    return result


def configure(root: Path, node_overrides: dict[str, int] | None = None) -> Path:
    root = require_campaign_root(root)
    nodes = dict(DEFAULT_CASE_NODES)
    for case_id, value in (node_overrides or {}).items():
        validate_node_profile(case_id, value)
        nodes[case_id] = value

    fast.DEFAULT_ROOT = root
    fast.FROZEN_SOURCE = FROZEN_SOURCE
    fast.MATRIX_RELATIVE = MATRIX_RELATIVE
    fast.MATRIX_SHA256 = MATRIX_SHA256
    fast.ATHENA = ATHENA
    fast.ATHENA_SHA256 = ATHENA_SHA256
    fast.FAST_RUNS_RELATIVE = RUNS_RELATIVE
    fast.LOGS_RELATIVE = LOGS_RELATIVE
    fast.CASE_NODES = nodes
    fast.CASE_INPUT_SHA256 = dict(CASE_INPUT_SHA256)
    fast.CASE_SEEDS = {}
    fast.__file__ = str(SCRIPT_PATH)
    fast.prepare_segment = prepare_segment
    fast.batch_script_text = batch_script_text
    return root


def validate_common_provenance(root: Path) -> None:
    require_campaign_root(root)
    require_source_revision()
    fast.require_sha(IDENTITY, IDENTITY_SHA256, "corrected campaign identity")
    fast.require_sha(FROZEN_SOURCE / MATRIX_RELATIVE, MATRIX_SHA256, "Stage I matrix")
    fast.require_sha(ATHENA, ATHENA_SHA256, "corrected executable")


def validate_provenance(root: Path, cases: Iterable[str] = ACTIVE_CASES) -> dict[str, object]:
    root = require_campaign_root(root)
    validate_common_provenance(root)
    matrix_cases = fast.campaign_cases()
    selected = fast.expand_cases(cases)
    for case_id in selected:
        case = matrix_cases[case_id]
        fast.require_sha(
            FROZEN_SOURCE / str(case["input"]),
            CASE_INPUT_SHA256[case_id],
            f"{case_id} corrected frozen input",
        )
    return {
        "campaign_id": CAMPAIGN_ID,
        "campaign_root": str(root),
        "source": str(FROZEN_SOURCE),
        "source_revision": SOURCE_REVISION,
        "campaign_identity": str(IDENTITY),
        "campaign_identity_sha256": IDENTITY_SHA256,
        "matrix_sha256": MATRIX_SHA256,
        "executable": str(ATHENA),
        "executable_sha256": ATHENA_SHA256,
        "fresh_initial_time": 0.0,
        "active_cases": selected,
        "case_nodes": {case_id: fast.CASE_NODES[case_id] for case_id in selected},
    }


def batch_script_text(manifest: dict[str, object]) -> str:
    # The retained base renderer reads module globals. Install corrected values
    # here as well as in configure() so direct rendering remains deterministic.
    fast.ATHENA = ATHENA
    fast.ATHENA_SHA256 = ATHENA_SHA256
    fast.__file__ = str(SCRIPT_PATH)
    text = _BASE_BATCH_SCRIPT_TEXT(manifest)
    text = text.replace("#SBATCH -J cglf_", "#SBATCH -J cglc_", 1)

    athena_assignment = f"ATHENA={shlex.quote(str(ATHENA))}\n"
    if text.count(athena_assignment) != 1:
        raise CorrectedFastError("unexpected direct-fast ATHENA assignment")
    text = text.replace(
        athena_assignment,
        athena_assignment
        + f"FROZEN_SOURCE={shlex.quote(str(FROZEN_SOURCE))}\n"
        + f"MATRIX={shlex.quote(str(FROZEN_SOURCE / MATRIX_RELATIVE))}\n"
        + f"IDENTITY={shlex.quote(str(IDENTITY))}\n"
        + f"CAMPAIGN_ROOT={shlex.quote(str(CAMPAIGN_ROOT.resolve()))}\n",
        1,
    )

    executable_check = f'require_sha {ATHENA_SHA256} "$ATHENA" executable\n'
    if text.count(executable_check) != 1:
        raise CorrectedFastError("unexpected direct-fast executable check")
    corrected_checks = f"""require_sha {ATHENA_SHA256} "$ATHENA" executable
require_sha {IDENTITY_SHA256} "$IDENTITY" corrected_campaign_identity
require_sha {MATRIX_SHA256} "$MATRIX" stage_i_matrix
source_revision="$(/usr/bin/git -C "$FROZEN_SOURCE" rev-parse --verify HEAD)"
test "$source_revision" = {SOURCE_REVISION} || {{
  echo "corrected frozen source revision mismatch: $source_revision" >&2
  exit 1
}}
test -z "$(/usr/bin/git -C "$FROZEN_SOURCE" status --porcelain --untracked-files=all)" || {{
  echo "corrected frozen source is dirty: $FROZEN_SOURCE" >&2
  exit 1
}}
case "$RUN_DIR" in
  "$CAMPAIGN_ROOT"/runs/*) ;;
  *) echo "run directory escapes corrected campaign root: $RUN_DIR" >&2; exit 1 ;;
esac
case "$OUT_DIR" in
  "$RUN_DIR"/output) ;;
  *) echo "output directory escapes corrected run directory: $OUT_DIR" >&2; exit 1 ;;
esac
"""
    return text.replace(executable_check, corrected_checks, 1)


def prepare_segment(
    root: Path,
    case: dict[str, object],
    sequence: int,
    start_time: float,
    restart: Path | None,
    restart_sha: str | None,
) -> Path:
    root = require_campaign_root(root)
    case_id = str(case["id"])
    if case_id not in ACTIVE_CASES:
        raise CorrectedFastError(f"only active corrected-production cases are supported: {case_id}")
    validate_node_profile(case_id, int(fast.CASE_NODES[case_id]))
    validate_common_provenance(root)

    if sequence == 0:
        if start_time != 0.0 or restart is not None or restart_sha is not None:
            raise CorrectedFastError(
                f"{case_id} initial corrected-production segment must start fresh from t=0"
            )
        launch_origin = "fresh_t0_corrected_production"
    else:
        if sequence < 1 or start_time <= 0.0 or restart is None or restart_sha is None:
            raise CorrectedFastError(f"{case_id} continuation metadata is incomplete")
        require_beneath(
            restart,
            root / RUNS_RELATIVE / case_id,
            f"{case_id} continuation restart",
        )
        launch_origin = "corrected_campaign_continuation"

    segment = _BASE_PREPARE_SEGMENT(
        root,
        case,
        sequence,
        start_time,
        restart,
        restart_sha,
    )
    script = segment / "manifest/run.sbatch"
    command_line_overrides: list[str] = []
    strict_admissibility = True
    variant = "corrected_production_strict"
    continuation_policy = "strict_base_scientific_sanity"
    if case_id in FINITE_LIMITER_DIAGNOSTIC_CASES:
        text = script.read_text(encoding="utf-8")
        needle = f" time/tlim={fast.TARGET_TIME}\n"
        if text.count(needle) != 1:
            raise CorrectedFastError(f"unexpected corrected diagnostic batch command: {script}")
        script.write_text(
            text.replace(needle, f" time/tlim={fast.TARGET_TIME} {FINITE_LIMITER_OVERRIDE}\n"),
            encoding="utf-8",
        )
        command_line_overrides = [FINITE_LIMITER_OVERRIDE]
        strict_admissibility = False
        variant = FINITE_LIMITER_VARIANT
        continuation_policy = "finite_progress_complete_terminal_products"

    manifest_path = fast.segment_manifest(segment)
    manifest = fast.load_json(manifest_path)
    manifest.update(
        {
            "campaign_id": CAMPAIGN_ID,
            "campaign_root": str(root),
            "source": str(FROZEN_SOURCE),
            "source_revision": SOURCE_REVISION,
            "campaign_identity": str(IDENTITY),
            "campaign_identity_sha256": IDENTITY_SHA256,
            "launch_origin": launch_origin,
            "fresh_lineage_root": sequence == 0,
            "legacy_restart_permitted": False,
            "variant": variant,
            "strict_admissibility": strict_admissibility,
            "command_line_overrides": command_line_overrides,
            "continuation_policy": continuation_policy,
        }
    )
    fast.write_json(manifest_path, manifest)
    return segment


def validate_segment(segment: Path) -> dict[str, object]:
    root = require_campaign_root(CAMPAIGN_ROOT)
    segment = require_beneath(segment, root / RUNS_RELATIVE, "corrected segment")
    manifest = fast.load_json(fast.segment_manifest(segment))
    case_id = str(manifest.get("case_id"))
    if case_id not in ACTIVE_CASES:
        raise CorrectedFastError(f"segment is not an active corrected case: {case_id}")

    exact = {
        "root": str(root),
        "campaign_root": str(root),
        "campaign_id": CAMPAIGN_ID,
        "source": str(FROZEN_SOURCE),
        "source_revision": SOURCE_REVISION,
        "campaign_identity": str(IDENTITY),
        "campaign_identity_sha256": IDENTITY_SHA256,
        "matrix_sha256": MATRIX_SHA256,
        "executable": str(ATHENA),
        "executable_sha256": ATHENA_SHA256,
        "input_sha256": CASE_INPUT_SHA256[case_id],
        "legacy_restart_permitted": False,
    }
    for key, expected in exact.items():
        if manifest.get(key) != expected:
            raise CorrectedFastError(f"corrected segment {key} differs: {segment}")

    script_text = (segment / "manifest/run.sbatch").read_text(encoding="utf-8")
    if case_id in FINITE_LIMITER_DIAGNOSTIC_CASES:
        diagnostic_exact = {
            "variant": FINITE_LIMITER_VARIANT,
            "strict_admissibility": False,
            "command_line_overrides": [FINITE_LIMITER_OVERRIDE],
            "continuation_policy": "finite_progress_complete_terminal_products",
        }
        for key, expected in diagnostic_exact.items():
            if manifest.get(key) != expected:
                raise CorrectedFastError(f"corrected diagnostic segment {key} differs: {segment}")
        if script_text.count(FINITE_LIMITER_OVERRIDE) != 1:
            raise CorrectedFastError(f"corrected diagnostic command-line override differs: {segment}")
    else:
        strict_exact = {
            "variant": "corrected_production_strict",
            "strict_admissibility": True,
            "command_line_overrides": [],
            "continuation_policy": "strict_base_scientific_sanity",
        }
        for key, expected in strict_exact.items():
            if manifest.get(key) != expected:
                raise CorrectedFastError(f"corrected strict segment {key} differs: {segment}")
        if FINITE_LIMITER_OVERRIDE in script_text:
            raise CorrectedFastError(f"strict corrected case has diagnostic override: {segment}")

    nodes = int(manifest["nodes"])
    validate_node_profile(case_id, nodes)
    if int(manifest["ranks"]) != nodes * fast.RANKS_PER_NODE:
        raise CorrectedFastError(f"corrected segment rank count differs: {segment}")

    sequence = int(manifest["sequence"])
    restart = manifest.get("restart")
    if sequence == 0:
        if (
            float(manifest["start_time"]) != 0.0
            or restart is not None
            or manifest.get("restart_sha256") is not None
            or manifest.get("launch_origin") != "fresh_t0_corrected_production"
            or manifest.get("fresh_lineage_root") is not True
        ):
            raise CorrectedFastError(f"corrected initial segment is not fresh from t=0: {segment}")
    else:
        if not isinstance(restart, str) or not restart:
            raise CorrectedFastError(f"corrected continuation lacks a restart: {segment}")
        require_beneath(
            Path(restart),
            root / RUNS_RELATIVE / case_id,
            f"{case_id} continuation restart",
        )
        if (
            manifest.get("launch_origin") != "corrected_campaign_continuation"
            or manifest.get("fresh_lineage_root") is not False
        ):
            raise CorrectedFastError(f"corrected continuation lineage differs: {segment}")
    return manifest


def launch_case(root: Path, case_id: str, *, submit: bool) -> Path | None:
    root = require_campaign_root(root)
    validate_common_provenance(root)
    cases = fast.campaign_cases()
    segments = fast.segment_directories(root, case_id)
    if segments:
        latest = segments[-1]
        manifest = validate_segment(latest)
        job_id = manifest.get("job_id")
        if job_id is None:
            if submit:
                fast.submit_segment(latest)
            else:
                print(latest)
            return latest
        if fast.job_state(str(job_id)) in {
            "PENDING",
            "RUNNING",
            "CONFIGURING",
            "COMPLETING",
        }:
            print(f"active {case_id}: {job_id}")
            return None
        return next_from_segment(latest, submit=submit)

    segment = prepare_segment(root, cases[case_id], 0, 0.0, None, None)
    if submit:
        fast.submit_segment(segment)
    else:
        print(segment)
    return segment


def next_from_segment(segment: Path, *, submit: bool) -> Path | None:
    manifest = validate_segment(segment)
    fast.CASE_NODES[str(manifest["case_id"])] = int(manifest["nodes"])
    analysis = analyze_segment(segment)
    if not analysis["passed"]:
        raise CorrectedFastError(f"corrected continuation gate failed: {segment}")
    if analysis["complete"]:
        print(f"complete {manifest['case_id']}: t={analysis['final_time']}")
        return None
    restart = analysis["terminal_restart"]
    if not isinstance(restart, dict):
        raise CorrectedFastError(f"missing corrected terminal restart: {segment}")
    next_segment = prepare_segment(
        CAMPAIGN_ROOT,
        fast.campaign_cases()[str(manifest["case_id"])],
        int(manifest["sequence"]) + 1,
        float(analysis["final_time"]),
        Path(str(restart["rank_zero"])),
        str(restart["sha256"]),
    )
    if submit:
        fast.submit_segment(next_segment)
    else:
        print(next_segment)
    return next_segment


def analyze_segment(segment: Path) -> dict[str, object]:
    manifest = validate_segment(segment)
    result = fast.analyze_segment(segment.resolve(), save=False)
    if str(manifest["case_id"]) in FINITE_LIMITER_DIAGNOSTIC_CASES:
        restart = result.get("terminal_restart")
        snapshot = result.get("terminal_snapshot")
        ranks = int(manifest["ranks"])
        final_time = float(result["final_time"])
        start_time = float(result["start_time"])
        finite_progress = math.isfinite(final_time) and final_time > start_time
        restart_complete = (
            isinstance(restart, dict)
            and int(restart.get("rank_count", -1)) == ranks
            and math.isclose(
                float(restart.get("physical_time", math.nan)),
                final_time,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        )
        snapshot_complete = (
            isinstance(snapshot, dict)
            and int(snapshot.get("rank_count", -1)) == ranks
        )
        result["base_strict_passed"] = result["passed"]
        result["continuation_policy"] = "finite_progress_complete_terminal_products"
        result["continuation_gate"] = {
            "finite_progress": finite_progress,
            "complete_terminal_restart": restart_complete,
            "complete_terminal_snapshot": snapshot_complete,
            "hard_bound_zero_required": False,
        }
        result["passed"] = finite_progress and restart_complete and snapshot_complete
    fast.write_json(segment / "manifest/fast_analysis.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        default=CAMPAIGN_ROOT,
        help="fixed corrected-production campaign root; any other root is rejected",
    )
    commands = parser.add_subparsers(dest="command", required=True)

    commands.add_parser("preflight", help="authenticate corrected source, executable, and inputs")

    launch = commands.add_parser("launch", help="prepare and optionally submit corrected active cases")
    launch.add_argument("cases", nargs="*", default=list(ACTIVE_CASES))
    launch.add_argument(
        "--nodes",
        action="append",
        default=[],
        metavar="RXX=N",
        help="override a selected case up to its useful meshblock-derived node ceiling",
    )
    launch.add_argument("--submit", action="store_true")

    advance = commands.add_parser("advance", help="analyze and continue a corrected segment")
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")

    analyze = commands.add_parser("analyze", help="analyze one corrected segment")
    analyze.add_argument("--segment", required=True, type=Path)

    status = commands.add_parser("status", help="show corrected-production case states")
    status.add_argument("cases", nargs="*", default=list(ACTIVE_CASES))

    args = parser.parse_args()
    try:
        root = configure(args.root)
        if args.command == "preflight":
            print(json.dumps(validate_provenance(root), indent=2, sort_keys=True))
        elif args.command == "launch":
            cases = fast.expand_cases(args.cases)
            root = configure(root, parse_node_overrides(args.nodes, cases))
            validate_provenance(root, cases)
            for case_id in cases:
                launch_case(root, case_id, submit=args.submit)
        elif args.command == "advance":
            next_from_segment(args.segment.resolve(), submit=args.submit)
        elif args.command == "analyze":
            print(json.dumps(analyze_segment(args.segment.resolve()), indent=2, sort_keys=True))
        elif args.command == "status":
            fast.command_status(root, fast.expand_cases(args.cases))
        return 0
    except (
        CorrectedFastError,
        fast.FastRunError,
        OSError,
        subprocess.SubprocessError,
        ValueError,
        KeyError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


configure(CAMPAIGN_ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
