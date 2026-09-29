#!/usr/bin/env python3
"""Build deterministic Stage I CT inventories from accepted whole-case bundles.

The builder accepts no free-form segment or restart paths.  It authenticates
the exact whole-case bundle, selects the preregistered retained states from its
accepted segment lineage, and emits the schema consumed by
``cgl_lf_stage_i_scientific_acceptance.py audit-ct-divb``.  Output is a
non-authorizing immutable candidate; this utility never publishes or mutates
canonical campaign state.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path
import sys


ACCEPTANCE_PATH = Path(__file__).with_name("cgl_lf_stage_i_scientific_acceptance.py")
ACCEPTANCE_SPEC = importlib.util.spec_from_file_location(
    "_cgl_lf_stage_i_ct_inventory_acceptance", ACCEPTANCE_PATH
)
if ACCEPTANCE_SPEC is None or ACCEPTANCE_SPEC.loader is None:
    raise RuntimeError(f"cannot load scientific acceptance utility: {ACCEPTANCE_PATH}")
acceptance = importlib.util.module_from_spec(ACCEPTANCE_SPEC)
sys.modules[ACCEPTANCE_SPEC.name] = acceptance
ACCEPTANCE_SPEC.loader.exec_module(acceptance)


SCHEMA_VERSION = 2
RECORD_TYPE = "stage-i-restart-ct-inventory"


class CtInventoryError(ValueError):
    """Raised when an accepted bundle cannot produce an exact CT inventory."""


def write_candidate(path: Path, inventory: dict[str, object]) -> None:
    """Write one immutable candidate while forbidding canonical-root mutation."""

    absolute = path.expanduser().absolute()
    try:
        absolute.resolve(strict=False).relative_to(
            acceptance.CANONICAL_CAMPAIGN_ROOT.resolve(strict=False)
        )
    except ValueError:
        pass
    else:
        raise CtInventoryError("candidate output beneath canonical root is forbidden")
    acceptance.write_candidate(absolute, inventory)


def normalized_rank_files(
    terminal: dict[str, object],
    expected_rank_count: int,
) -> list[dict[str, object]]:
    """Authenticate and normalize one exact terminal rank-local restart set."""

    if terminal.get("storage") != "per_rank":
        raise CtInventoryError("terminal restart storage must be per_rank")
    declared = acceptance.require_list(
        terminal.get("rank_files"), "terminal restart rank files"
    )
    if len(declared) != expected_rank_count:
        raise CtInventoryError("terminal restart rank count differs from allocation")
    ranks: list[dict[str, object]] = []
    paths: set[str] = set()
    digests: set[str] = set()
    for value in declared:
        observed = acceptance.verify_declared_binding(value, "terminal restart rank file")
        rank = acceptance.rank_identity(Path(str(observed["path"])))
        path = str(observed["path"])
        digest = str(observed["sha256"])
        if path in paths or digest in digests:
            raise CtInventoryError("terminal restart contains duplicated path or copied rank bytes")
        paths.add(path)
        digests.add(digest)
        ranks.append({**observed, "rank": rank})
    ranks.sort(key=lambda value: int(value["rank"]))
    if [int(value["rank"]) for value in ranks] != list(range(expected_rank_count)):
        raise CtInventoryError("terminal restart canonical ranks are not exact and contiguous")
    rank_zero = ranks[0]
    terminal_observed = acceptance.verify_declared_binding(
        terminal, "terminal restart representative"
    )
    if terminal_observed != {
        key: rank_zero[key] for key in ("path", "size_bytes", "sha256")
    }:
        raise CtInventoryError("terminal restart representative differs from rank zero")
    return ranks


def accepted_segment_record(
    policy: dict[str, object],
    case_id: str,
    expected_case: dict[str, object],
    segment_path: Path,
) -> tuple[dict[str, object], dict[str, object], float, str]:
    """Authenticate one accepted segment named by the whole-case bundle."""

    segment, binding = acceptance.load_json(segment_path, "accepted segment manifest")
    accounting = acceptance.require_dict(segment.get("accounting"), "segment accounting")
    command = acceptance.require_dict(segment.get("command"), "segment command")
    inspection = acceptance.require_dict(
        segment.get("scientific_inspection"), "segment scientific inspection"
    )
    final_time = acceptance.require_finite(
        inspection.get("final_time"), "segment final time"
    )
    executable = acceptance.require_sha256(
        command.get("executable_sha256"), "segment executable SHA-256"
    )
    if (
        accounting.get("result") != "accepted"
        or accounting.get("case_id") != case_id
        or accounting.get("case_name") != expected_case.get("name")
        or accounting.get("executable_sha256") != executable
        or inspection.get("schema_version") != 3
        or inspection.get("accepted") is not True
        or inspection.get("case_id") != case_id
        or acceptance.resolve_bound_path(inspection.get("manifest")).resolve(strict=True)
        != Path(str(binding["path"]))
        or command.get("matrix_sha256")
        != policy["verified_sources"]["stage_i_manifest"]["sha256"]
    ):
        raise CtInventoryError("whole-case bundle segment is not an accepted CT source")
    return segment, binding, final_time, executable


def state_from_segment(
    segment: dict[str, object],
    segment_binding: dict[str, object],
    required_time: float,
) -> dict[str, object]:
    """Return one exact required retained state from an accepted segment."""

    inspection = acceptance.require_dict(
        segment.get("scientific_inspection"), "required state scientific inspection"
    )
    accounting = acceptance.require_dict(segment.get("accounting"), "required state accounting")
    final_time = acceptance.require_finite(
        inspection.get("final_time"), "required state final time"
    )
    if (
        final_time != required_time
        or acceptance.require_finite(
            inspection.get("required_time"), "required state required time"
        )
        != required_time
        or acceptance.require_finite(
            inspection.get("terminal_restart_time"), "required state terminal restart time"
        )
        != required_time
        or inspection.get("segment") != accounting.get("segment")
    ):
        raise CtInventoryError("accepted segment does not terminate at the exact required state")
    checks = acceptance.require_dict(inspection.get("checks"), "required state checks")
    terminal = acceptance.require_dict(
        inspection.get("terminal_restart"), "required state terminal restart"
    )
    if (
        checks.get("required_time_reached") is not True
        or checks.get("restart_retained") is not True
        or checks.get("terminal_restart_physical_time_matches_final") is not True
        or required_time
        not in [
            acceptance.require_finite(value, "required state restart time")
            for value in acceptance.require_list(
                inspection.get("restart_times"), "required state restart times"
            )
        ]
        or terminal
        not in acceptance.require_list(
            inspection.get("restarts"), "required state retained restarts"
        )
    ):
        raise CtInventoryError("accepted segment lacks the exact retained terminal restart")
    allocation = acceptance.require_dict(segment.get("allocation"), "required state allocation")
    expected_rank_count = acceptance.require_int(
        allocation.get("nodes"), "required state nodes", minimum=1
    ) * acceptance.require_int(
        allocation.get("ranks_per_node"), "required state ranks_per_node", minimum=1
    )
    return {
        "time": required_time,
        "accepted_segment_manifest": segment_binding,
        "rank_files": normalized_rank_files(terminal, expected_rank_count),
    }


def build_ct_inventory(
    policy: dict[str, object],
    case_id: str,
    bundle_path: Path,
    expected_bundle_sha256: str,
    authority_mode: str,
) -> dict[str, object]:
    """Build one deterministic schema-2 CT inventory from an accepted bundle only."""

    case_id = acceptance.require_case_id(case_id)
    expected_bundle_sha256 = acceptance.require_sha256(
        expected_bundle_sha256, "expected bundle SHA-256"
    )
    if authority_mode not in ("canonical", "offline"):
        raise CtInventoryError("authority_mode must be canonical or offline")
    bundle, bundle_binding = acceptance.load_json(bundle_path, "accepted whole-case bundle")
    if bundle_binding["sha256"] != expected_bundle_sha256:
        raise CtInventoryError("accepted whole-case bundle SHA-256 differs")
    if authority_mode == "canonical" and not acceptance.path_is_exact(
        Path(str(bundle_binding["path"])), acceptance.canonical_bundle_path(case_id)
    ):
        raise CtInventoryError("canonical authority requires the exact canonical bundle path")
    expected_case = acceptance.case_manifest_record(policy, case_id)
    selected_cases = acceptance.require_list(bundle.get("cases"), "accepted bundle cases")
    if (
        bundle.get("workflow") != "paper-mks24-stage-i-production"
        or bundle.get("status") != "accepted_for_analysis"
        or bundle.get("production_case_id") != case_id
        or acceptance.require_finite(
            bundle.get("required_final_time"), "bundle required final time"
        )
        != 10.0
        or acceptance.require_finite(
            bundle.get("accepted_final_time"), "bundle accepted final time"
        )
        != 10.0
        or len(selected_cases) != 1
    ):
        raise CtInventoryError("accepted whole-case bundle identity differs")
    selected_case = acceptance.require_dict(selected_cases[0], "accepted bundle case")
    if (
        selected_case.get("name") != expected_case.get("name")
        or selected_case.get("input") != expected_case.get("input")
        or selected_case.get("status") != "passed"
    ):
        raise CtInventoryError("accepted whole-case bundle selected case differs")

    segment_paths = [
        acceptance.resolve_bound_path(value).resolve(strict=True)
        for value in acceptance.require_list(
            bundle.get("production_segment_manifests"), "production segment manifests"
        )
    ]
    if not segment_paths or len(segment_paths) != len(set(segment_paths)):
        raise CtInventoryError("accepted whole-case bundle segment inventory is empty or duplicated")
    segments: list[tuple[dict[str, object], dict[str, object], float, str]] = []
    for path in segment_paths:
        segments.append(accepted_segment_record(policy, case_id, expected_case, path))
    times = [record[2] for record in segments]
    if times != sorted(times) or len(times) != len(set(times)) or times[-1] != 10.0:
        raise CtInventoryError("accepted whole-case bundle segment times are not exact and ordered")
    executables = {record[3] for record in segments}
    if len(executables) != 1:
        raise CtInventoryError("accepted whole-case bundle executable lineage differs")

    ct_policy = acceptance.require_dict(
        policy["criteria"].get("ct_divb_policy"), "CT-divB policy"
    )
    required_times = [
        acceptance.require_finite(value, "required CT state time")
        for value in acceptance.require_list(
            ct_policy.get("required_state_times"), "required CT state times"
        )
    ]
    states: list[dict[str, object]] = []
    all_paths: set[str] = set()
    all_digests: set[str] = set()
    for required_time in required_times:
        matches = [record for record in segments if record[2] == required_time]
        if len(matches) != 1:
            raise CtInventoryError(
                f"accepted whole-case bundle must contain exactly one t={required_time:g} state"
            )
        state = state_from_segment(matches[0][0], matches[0][1], required_time)
        for rank in acceptance.require_list(state["rank_files"], "state rank files"):
            record = acceptance.require_dict(rank, "state rank file")
            path = str(record["path"])
            digest = str(record["sha256"])
            if path in all_paths or digest in all_digests:
                raise CtInventoryError("inventory contains duplicated path or copied rank bytes")
            all_paths.add(path)
            all_digests.add(digest)
        states.append(state)
    return {
        "schema_version": SCHEMA_VERSION,
        "record_type": RECORD_TYPE,
        "case_id": case_id,
        "coverage_complete": True,
        "required_state_times": required_times,
        "executable_sha256": next(iter(executables)),
        "accepted_bundle_manifest": bundle_binding,
        "states": states,
    }


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line interface."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--criteria", type=Path, default=acceptance.DEFAULT_CRITERIA)
    parser.add_argument(
        "--criteria-review", type=Path, default=acceptance.DEFAULT_CRITERIA_REVIEW
    )
    parser.add_argument("--case-id", required=True)
    parser.add_argument("--bundle-manifest", type=Path, required=True)
    parser.add_argument("--expected-bundle-sha256", required=True)
    parser.add_argument("--authority-mode", choices=("canonical", "offline"), required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    """Build one immutable non-authorizing CT inventory candidate."""

    args = build_parser().parse_args(argv)
    try:
        policy = acceptance.load_validated_policy(args.criteria, args.criteria_review)
        inventory = build_ct_inventory(
            policy,
            args.case_id,
            args.bundle_manifest,
            args.expected_bundle_sha256,
            args.authority_mode,
        )
        write_candidate(args.output, inventory)
    except (CtInventoryError, acceptance.AcceptanceError, OSError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    sys.stdout.write(json.dumps({
        "case_id": inventory["case_id"],
        "output": str(args.output.expanduser().absolute()),
        "output_sha256": acceptance.sha256_bytes(acceptance.stable_json(inventory)),
        "required_state_times": inventory["required_state_times"],
        "state_count": len(inventory["states"]),
    }, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
