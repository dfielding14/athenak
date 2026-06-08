#!/usr/bin/env python3
"""Produce final reviewed science from the corrected/composite Stage I report.

This adapter fills the narrow gap between corrected downstream analysis and
the existing direct-fast acceptance/science tools.  It admits exactly the
corrected active R02--R05 and R10--R17 trajectories plus authenticated legacy
R06--R09 passive controls.  It delegates all scientific evaluation to the
reviewed fast-acceptance and fast-science utilities.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
CORRECTED_DOWNSTREAM_TOOL = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_fast_corrected_downstream.py"
)
CORRECTED_REPORT_TOOL = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_fast_corrected_report.py"
)
FAST_ACCEPTANCE_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_acceptance.py")
FAST_SCIENCE_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_science.py")
REVIEWED_ACCEPTANCE_TOOL = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_scientific_acceptance.py"
)
DEFAULT_CRITERIA = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.json"
)
DEFAULT_CRITERIA_REVIEW = REPO_ROOT / (
    "inputs/cgl_lf_paper/mks24_stage_i_scientific_acceptance_criteria.review.json"
)
DEFAULT_CAMPAIGN = Path(
    "/lustre/orion/ast207/proj-shared/dfielding/CGL/campaigns/"
    "mks24-stage-i-eos-fastdisc-ppar2-corrected-v1"
)
DEFAULT_IDENTITY = DEFAULT_CAMPAIGN / "campaign-identity.json"
DEFAULT_INVENTORY = DEFAULT_CAMPAIGN / (
    "analysis/mks24-stage-i-composite/E03-forcing-policy/R02-R17-final/inventory.json"
)
DEFAULT_ACCEPTANCE = DEFAULT_INVENTORY.parent.with_name(
    f"{DEFAULT_INVENTORY.parent.name}-acceptance"
)
DEFAULT_SCIENCE = DEFAULT_INVENTORY.parent.with_name(
    f"{DEFAULT_INVENTORY.parent.name}-science"
)

SCHEMA = "athenak-cgl-corrected-composite-reviewed-science"
SCHEMA_VERSION = 1
RECORD_NAME = "corrected-composite-science.json"
ACTIVE_EVIDENCE_CLASS = "corrected_active_production"
PASSIVE_EVIDENCE_CLASS = "authenticated_legacy_passive_control"
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
PASSIVE_CASES = ("R06", "R07", "R08", "R09")
ALL_CASES = tuple(f"R{number:02d}" for number in range(2, 18))
VALID_RESULTS = {"pass", "fail", "inconclusive"}
DEPENDENCY_PINS = (
    ("corrected_downstream_sha256", "corrected downstream", CORRECTED_DOWNSTREAM_TOOL),
    ("corrected_report_sha256", "corrected report", CORRECTED_REPORT_TOOL),
    ("fast_acceptance_sha256", "fast acceptance", FAST_ACCEPTANCE_TOOL),
    ("fast_science_sha256", "fast science", FAST_SCIENCE_TOOL),
    ("reviewed_acceptance_sha256", "reviewed acceptance", REVIEWED_ACCEPTANCE_TOOL),
)
EXPECTED_SCIENCE_TABLES = tuple(
    f"tables/{stem}.{suffix}"
    for stem in ("cases", "contrasts", "gates", "resolution", "mks24")
    for suffix in ("csv", "md")
)


class CorrectedScienceError(RuntimeError):
    """The corrected/composite reviewed-science contract failed."""


def load_module(name: str, path: Path):
    """Import one local workflow dependency by exact path."""

    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CorrectedScienceError(f"cannot import workflow dependency: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one file."""

    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha256(value: object, label: str) -> str:
    """Return one exact lowercase SHA-256 digest."""

    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise CorrectedScienceError(f"{label} must be a lowercase SHA-256 digest")
    return value


def verify_dependency_pins(args: argparse.Namespace) -> dict[str, dict[str, object]]:
    """Verify every directly imported or invoked workflow dependency."""

    verified: dict[str, dict[str, object]] = {}
    for argument, label, path in DEPENDENCY_PINS:
        expected = require_sha256(getattr(args, argument), f"{label} pin")
        current = binding(path)
        if current["sha256"] != expected:
            raise CorrectedScienceError(
                f"{label} SHA-256 differs: expected={expected} "
                f"actual={current['sha256']} path={current['path']}"
            )
        verified[argument] = current
    return verified


def binding(path: Path) -> dict[str, object]:
    """Bind one current regular file."""

    resolved = path.expanduser().absolute().resolve(strict=True)
    if not resolved.is_file():
        raise CorrectedScienceError(f"expected a regular file: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise CorrectedScienceError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise CorrectedScienceError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise CorrectedScienceError(f"{label} must be a nonempty string")
    return value


def exact_scope_records(
    reviewed: object, policy: dict[str, object]
) -> tuple[dict[str, object], dict[str, object]]:
    """Return exact validated current-review and total-intervention scopes."""

    limitation = policy.get("current_science_scope_limitation")
    if (
        not isinstance(limitation, dict)
        or limitation != getattr(reviewed, "CURRENT_SCIENCE_SCOPE_LIMITATION", None)
    ):
        raise CorrectedScienceError("current science scope limitation differs")
    try:
        intervention = policy["criteria"]["family_gates"]["active_passive"][
            "intervention_scope"
        ]
    except (KeyError, TypeError) as error:
        raise CorrectedScienceError(
            "active/passive intervention scope is unavailable"
        ) from error
    if (
        not isinstance(intervention, dict)
        or intervention != getattr(reviewed, "ACTIVE_PASSIVE_INTERVENTION_SCOPE", None)
    ):
        raise CorrectedScienceError("active/passive intervention scope differs")
    return limitation, intervention


def require_exact_scope_records(
    value: dict[str, object],
    limitation: dict[str, object],
    intervention: dict[str, object],
    label: str,
) -> None:
    """Require one retained product to carry both exact science scopes."""

    if value.get("current_science_scope_limitation") != limitation:
        raise CorrectedScienceError(f"{label} current science scope limitation differs")
    if value.get("active_passive_intervention_scope") != intervention:
        raise CorrectedScienceError(f"{label} active/passive intervention scope differs")


def validated_scope_records(
    criteria: Path, criteria_review: Path
) -> tuple[dict[str, object], dict[str, object]]:
    """Load and return the exact scopes for independent downstream validation."""

    reviewed = load_module("_cgl_corrected_science_reviewed", REVIEWED_ACCEPTANCE_TOOL)
    try:
        policy = reviewed.load_validated_policy(criteria, criteria_review)
    except Exception as error:
        raise CorrectedScienceError(
            f"cannot load reviewed scientific policy: {error}"
        ) from error
    return exact_scope_records(reviewed, policy)


def verify_binding(value: object, label: str) -> dict[str, object]:
    """Verify and normalize one path/SHA/size binding."""

    declared = require_dict(value, f"{label} binding")
    current = binding(Path(require_text(declared.get("path"), f"{label} path")))
    if (
        declared.get("sha256") != current["sha256"]
        or declared.get("size_bytes") != current["size_bytes"]
    ):
        raise CorrectedScienceError(f"{label} differs from its declared binding")
    return current


def require_same_binding(value: object, expected: object, label: str) -> None:
    if verify_binding(value, label) != verify_binding(expected, f"expected {label}"):
        raise CorrectedScienceError(f"{label} differs from the selected artifact")


def load_json(path: Path, label: str) -> dict[str, object]:
    """Load one finite JSON object."""

    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CorrectedScienceError(f"cannot load {label}: {path}") from error
    return require_dict(value, label)


def stable_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def write_json(path: Path, value: object) -> None:
    """Atomically write deterministic JSON."""

    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as stream:
            staged = Path(stream.name)
            stream.write(stable_json(value))
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def path_contains(parent: Path, child: Path) -> bool:
    return child == parent or parent in child.parents


def validate_output_layout(
    context: dict[str, object], acceptance_output: Path, science_output: Path
) -> tuple[Path, Path]:
    """Require separate final acceptance and science output trees."""

    report = Path(str(context["inventory_output"])).resolve(strict=True)
    acceptance = acceptance_output.expanduser().absolute().resolve(strict=False)
    science = science_output.expanduser().absolute().resolve(strict=False)
    for output, label in ((acceptance, "acceptance"), (science, "science")):
        if path_contains(report, output) or path_contains(output, report):
            raise CorrectedScienceError(
                f"{label} output must be separate from the composite report"
            )
    if path_contains(acceptance, science) or path_contains(science, acceptance):
        raise CorrectedScienceError(
            "acceptance and science output trees must be separate"
        )
    return acceptance, science


def validate_corrected_context(
    identity: Path, inventory: Path, inventory_sha256: str
) -> dict[str, object]:
    """Authenticate the exact corrected-active/legacy-passive composite split."""

    downstream = load_module(
        "_cgl_corrected_science_downstream", CORRECTED_DOWNSTREAM_TOOL
    )
    try:
        context = downstream.validate_inventory(identity, inventory, inventory_sha256)
    except Exception as error:
        raise CorrectedScienceError(
            f"corrected/composite inventory authentication failed: {error}"
        ) from error
    if context.get("inventory_kind") != "corrected-composite":
        raise CorrectedScienceError(
            "final reviewed science requires an authenticated corrected composite inventory"
        )
    if context.get("active_cases") != list(ACTIVE_CASES):
        raise CorrectedScienceError("corrected active case classification differs")
    if context.get("passive_cases") != list(PASSIVE_CASES):
        raise CorrectedScienceError(
            "only authenticated legacy R06-R09 passive controls may be admitted"
        )
    if context.get("selected_cases") != list(ALL_CASES):
        raise CorrectedScienceError("final reviewed science requires exact R02-R17 coverage")

    record = require_dict(context.get("inventory"), "composite inventory")
    if record.get("record_type") != "cgl_lf_stage_i_corrected_composite_report":
        raise CorrectedScienceError("composite inventory record type differs")
    expected_classes = {
        ACTIVE_EVIDENCE_CLASS: list(ACTIVE_CASES),
        PASSIVE_EVIDENCE_CLASS: list(PASSIVE_CASES),
    }
    if record.get("evidence_classes") != expected_classes:
        raise CorrectedScienceError("composite evidence-class inventory differs")
    compatibility = require_dict(
        record.get("passive_compatibility"), "passive compatibility"
    )
    contract_id = require_text(
        compatibility.get("contract_id"), "passive compatibility contract ID"
    )
    if (
        compatibility.get("result") != "pass"
        or compatibility.get("scope") != list(PASSIVE_CASES)
    ):
        raise CorrectedScienceError("legacy passive compatibility is not an exact pass")

    cases = require_dict(record.get("cases"), "composite cases")
    if set(cases) != set(ALL_CASES):
        raise CorrectedScienceError("composite case coverage differs from R02-R17")
    for case_id in ALL_CASES:
        case = require_dict(cases[case_id], f"{case_id} composite case")
        authority = require_dict(
            case.get("execution_authority"), f"{case_id} execution authority"
        )
        expected = (
            ACTIVE_EVIDENCE_CLASS if case_id in ACTIVE_CASES else PASSIVE_EVIDENCE_CLASS
        )
        if case.get("evidence_class") != expected or authority.get(
            "evidence_class"
        ) != expected:
            raise CorrectedScienceError(
                f"{case_id} evidence class is not authorized for final science"
            )
        compatibility_id = authority.get("compatibility_contract_id")
        if case_id in PASSIVE_CASES and compatibility_id != contract_id:
            raise CorrectedScienceError(
                f"{case_id} lacks the authenticated passive compatibility contract"
            )
        if case_id in ACTIVE_CASES and compatibility_id is not None:
            raise CorrectedScienceError(
                f"{case_id} corrected active evidence was relabeled as passive-compatible"
            )
    return context


def workflow_commands(
    context: dict[str, object],
    acceptance_output: Path,
    science_output: Path,
    criteria: Path,
    criteria_review: Path,
    python: Path,
) -> tuple[list[str], list[str]]:
    """Build exact existing-tool invocations for all R02-R17 cases."""

    inventory = Path(
        require_text(
            require_dict(context["inventory_binding"], "inventory binding").get("path"),
            "inventory path",
        )
    )
    cases = ",".join(ALL_CASES)
    acceptance = [
        str(python),
        str(FAST_ACCEPTANCE_TOOL),
        "--inventory",
        str(inventory),
        "--output",
        str(acceptance_output),
        "--criteria",
        str(criteria),
        "--criteria-review",
        str(criteria_review),
        "--inventory-only",
        "--cases",
        cases,
    ]
    science = [
        str(python),
        str(FAST_SCIENCE_TOOL),
        "--inventory",
        str(inventory),
        "--acceptance",
        str(acceptance_output),
        "--output",
        str(science_output),
        "--criteria",
        str(criteria),
        "--criteria-review",
        str(criteria_review),
        "--cases",
        cases,
    ]
    return acceptance, science


def run_checked(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise CorrectedScienceError(
            f"reviewed-science stage failed ({completed.returncode}): "
            f"{' '.join(command)}\n{detail}"
        )
    if completed.stdout.strip():
        print(completed.stdout.strip())


def validate_products(
    context: dict[str, object],
    acceptance_output: Path,
    science_output: Path,
    criteria: Path,
    criteria_review: Path,
) -> dict[str, object]:
    """Revalidate final existing-tool outputs and build their binding record."""

    acceptance_root = acceptance_output.resolve(strict=True)
    science_root = science_output.resolve(strict=True)
    inventory = verify_binding(context["inventory_binding"], "composite inventory")
    inventory_path = Path(str(inventory["path"]))
    fast_science = load_module("_cgl_corrected_science_aggregator", FAST_SCIENCE_TOOL)
    reviewed = load_module("_cgl_corrected_science_reviewed", REVIEWED_ACCEPTANCE_TOOL)
    try:
        policy = reviewed.load_validated_policy(criteria, criteria_review)
        limitation, intervention = exact_scope_records(reviewed, policy)
        fast_science.validate_inventory(inventory_path, policy, list(ALL_CASES))
        summaries, campaign, acceptance_provenance, _ = (
            fast_science.validate_acceptance_root(
                acceptance_root, inventory_path, policy
            )
        )
    except Exception as error:
        raise CorrectedScienceError(
            f"acceptance/science provenance validation failed: {error}"
        ) from error
    if set(summaries) != set(ALL_CASES):
        raise CorrectedScienceError("acceptance evidence does not cover exact R02-R17")
    if campaign.get("result") not in VALID_RESULTS:
        raise CorrectedScienceError("acceptance campaign result is malformed")
    require_exact_scope_records(
        campaign, limitation, intervention, "acceptance campaign"
    )

    acceptance_summary = load_json(
        acceptance_root / "summary.json", "acceptance summary"
    )
    if (
        acceptance_summary.get("record_type")
        != "cgl-lf-stage-i-direct-fast-acceptance-summary"
        or acceptance_summary.get("selected_cases") != list(ALL_CASES)
        or set(
            require_dict(
                acceptance_summary.get("case_results"),
                "acceptance summary case results",
            )
        )
        != set(ALL_CASES)
    ):
        raise CorrectedScienceError("acceptance summary coverage differs from R02-R17")
    require_exact_scope_records(
        acceptance_summary, limitation, intervention, "acceptance summary"
    )
    acceptance_outputs = require_dict(
        acceptance_provenance.get("outputs"), "acceptance provenance outputs"
    )
    require_same_binding(
        acceptance_outputs.get("summary"),
        binding(acceptance_root / "summary.json"),
        "acceptance summary",
    )
    require_same_binding(
        acceptance_outputs.get("campaign_evidence"),
        binding(acceptance_root / "campaign_evidence.json"),
        "acceptance campaign evidence",
    )
    case_acceptance = require_dict(
        acceptance_outputs.get("case_acceptance"), "case acceptance outputs"
    )
    if set(case_acceptance) != set(ALL_CASES):
        raise CorrectedScienceError("case acceptance output coverage differs from R02-R17")
    for case_id, value in case_acceptance.items():
        require_same_binding(
            value,
            binding(acceptance_root / "cases" / case_id / "case_acceptance.json"),
            f"{case_id} acceptance output",
        )

    science_record = load_json(science_root / "science.json", "reviewed science")
    try:
        reviewed.verify_evidence_digest(science_record, "reviewed science")
    except Exception as error:
        raise CorrectedScienceError(f"reviewed science evidence is forged: {error}") from error
    if (
        science_record.get("record_type")
        != "cgl-lf-stage-i-direct-fast-reviewed-science-comparisons"
        or science_record.get("selected_cases") != list(ALL_CASES)
        or science_record.get("result") not in VALID_RESULTS
    ):
        raise CorrectedScienceError("reviewed science identity or coverage differs")
    require_exact_scope_records(
        science_record, limitation, intervention, "direct reviewed science"
    )
    dispositions = require_dict(
        science_record.get("case_dispositions"), "reviewed science case dispositions"
    )
    if set(dispositions) != set(ALL_CASES):
        raise CorrectedScienceError("reviewed science dispositions do not cover R02-R17")
    for case_id, value in dispositions.items():
        disposition = require_dict(value, f"{case_id} science disposition")
        if (
            disposition.get("inventory_status") != "complete"
            or disposition.get("diagnostics_available") is not True
            or disposition.get("snapshot_products_authenticated") is not True
        ):
            raise CorrectedScienceError(
                f"{case_id} lacks complete authenticated downstream science products"
            )

    science_provenance = load_json(
        science_root / "provenance.json", "reviewed science provenance"
    )
    if (
        science_provenance.get("record_type")
        != "cgl-lf-stage-i-direct-fast-reviewed-science-provenance"
        or science_provenance.get("inputs") != science_record.get("provenance")
    ):
        raise CorrectedScienceError("reviewed science provenance identity differs")
    provenance_inputs = require_dict(
        science_provenance.get("inputs"), "reviewed science provenance inputs"
    )
    require_same_binding(
        provenance_inputs.get("inventory"), inventory, "reviewed science inventory"
    )
    require_same_binding(
        provenance_inputs.get("acceptance_provenance"),
        binding(acceptance_root / "provenance.json"),
        "reviewed science acceptance provenance",
    )
    require_same_binding(
        provenance_inputs.get("acceptance_campaign_evidence"),
        binding(acceptance_root / "campaign_evidence.json"),
        "reviewed science acceptance campaign",
    )
    require_same_binding(
        provenance_inputs.get("aggregator"),
        binding(FAST_SCIENCE_TOOL),
        "reviewed science aggregator",
    )
    for key in ("case_acceptance", "case_lineages", "case_diagnostics"):
        values = require_dict(provenance_inputs.get(key), f"reviewed science {key}")
        if set(values) != set(ALL_CASES):
            raise CorrectedScienceError(
                f"reviewed science {key} coverage differs from R02-R17"
            )
        for case_id, value in values.items():
            verify_binding(value, f"{case_id} reviewed science {key}")

    science_outputs = require_dict(
        science_provenance.get("outputs"), "reviewed science provenance outputs"
    )
    require_same_binding(
        science_outputs.get("science"),
        binding(science_root / "science.json"),
        "reviewed science record",
    )
    tables = require_list(science_outputs.get("tables"), "reviewed science tables")
    table_bindings = {
        str(Path(str(verify_binding(value, "reviewed science table")["path"])).relative_to(
            science_root
        )): verify_binding(value, "reviewed science table")
        for value in tables
    }
    if set(table_bindings) != set(EXPECTED_SCIENCE_TABLES):
        raise CorrectedScienceError("reviewed science table inventory differs")

    compatibility = require_dict(
        require_dict(context["inventory"], "composite inventory").get(
            "passive_compatibility"
        ),
        "passive compatibility",
    )
    reviewed_case_outputs = require_dict(
        acceptance_outputs.get("reviewed_case_evidence"),
        "reviewed case evidence outputs",
    )
    if not set(reviewed_case_outputs) <= set(ALL_CASES):
        raise CorrectedScienceError("reviewed case evidence contains unsupported cases")
    for case_id, value in reviewed_case_outputs.items():
        verify_binding(value, f"{case_id} reviewed case evidence")

    return {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "record_type": "cgl-lf-stage-i-corrected-composite-reviewed-science",
        "status": "complete",
        "result": science_record["result"],
        "campaign_kind": "corrected-composite",
        "active_passive_intervention_scope": intervention,
        "current_science_scope_limitation": limitation,
        "campaign_identity": context["identity_binding"],
        "inventory": inventory,
        "case_classification": {
            "corrected_active": list(ACTIVE_CASES),
            "authenticated_legacy_passive": list(PASSIVE_CASES),
            "selected": list(ALL_CASES),
        },
        "passive_compatibility": {
            "contract_id": compatibility["contract_id"],
            "result": compatibility["result"],
            "scope": compatibility["scope"],
        },
        "tools": {
            "adapter": binding(SCRIPT_PATH),
            "corrected_downstream": binding(CORRECTED_DOWNSTREAM_TOOL),
            "corrected_report": binding(CORRECTED_REPORT_TOOL),
            "fast_acceptance": binding(FAST_ACCEPTANCE_TOOL),
            "fast_science": binding(FAST_SCIENCE_TOOL),
            "reviewed_acceptance": binding(REVIEWED_ACCEPTANCE_TOOL),
        },
        "acceptance": {
            "root": str(acceptance_root),
            "provenance": binding(acceptance_root / "provenance.json"),
            "summary": binding(acceptance_root / "summary.json"),
            "campaign_evidence": binding(acceptance_root / "campaign_evidence.json"),
            "case_acceptance": {
                case_id: verify_binding(value, f"{case_id} case acceptance")
                for case_id, value in sorted(case_acceptance.items())
            },
            "reviewed_case_evidence": {
                case_id: verify_binding(value, f"{case_id} reviewed case evidence")
                for case_id, value in sorted(reviewed_case_outputs.items())
            },
        },
        "science": {
            "root": str(science_root),
            "provenance": binding(science_root / "provenance.json"),
            "record": binding(science_root / "science.json"),
            "tables": dict(sorted(table_bindings.items())),
        },
    }


def run_workflow(args: argparse.Namespace) -> Path:
    """Run existing acceptance/science tools and write the corrected binding record."""

    verify_dependency_pins(args)
    context = validate_corrected_context(
        args.identity, args.inventory, args.inventory_sha256
    )
    acceptance, science = validate_output_layout(
        context, args.acceptance_output, args.science_output
    )
    python = args.python.expanduser().absolute().resolve(strict=True)
    criteria = args.criteria.expanduser().absolute().resolve(strict=True)
    review = args.criteria_review.expanduser().absolute().resolve(strict=True)
    commands = workflow_commands(
        context, acceptance, science, criteria, review, python
    )
    run_checked(commands[0])
    verify_dependency_pins(args)
    validate_corrected_context(args.identity, args.inventory, args.inventory_sha256)
    run_checked(commands[1])
    verify_dependency_pins(args)
    context = validate_corrected_context(
        args.identity, args.inventory, args.inventory_sha256
    )
    record = validate_products(context, acceptance, science, criteria, review)
    path = science / RECORD_NAME
    write_json(path, record)
    if load_json(path, "corrected reviewed-science record") != validate_products(
        context, acceptance, science, criteria, review
    ):
        raise CorrectedScienceError("corrected reviewed-science record differs")
    verify_dependency_pins(args)
    print(f"wrote corrected/composite reviewed-science record: {path}")
    return path


def validate_workflow(args: argparse.Namespace) -> Path:
    """Revalidate the corrected/composite binding record and all source evidence."""

    verify_dependency_pins(args)
    context = validate_corrected_context(
        args.identity, args.inventory, args.inventory_sha256
    )
    acceptance, science = validate_output_layout(
        context, args.acceptance_output, args.science_output
    )
    expected = validate_products(
        context,
        acceptance,
        science,
        args.criteria.expanduser().absolute().resolve(strict=True),
        args.criteria_review.expanduser().absolute().resolve(strict=True),
    )
    path = science / RECORD_NAME
    if load_json(path, "corrected reviewed-science record") != expected:
        raise CorrectedScienceError("corrected reviewed-science record differs")
    verify_dependency_pins(args)
    print(f"validated corrected/composite reviewed-science record: {path}")
    return path


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--inventory-sha256", required=True)
    parser.add_argument("--acceptance-output", type=Path, default=DEFAULT_ACCEPTANCE)
    parser.add_argument("--science-output", type=Path, default=DEFAULT_SCIENCE)
    parser.add_argument("--criteria", type=Path, default=DEFAULT_CRITERIA)
    parser.add_argument("--criteria-review", type=Path, default=DEFAULT_CRITERIA_REVIEW)
    parser.add_argument("--corrected-downstream-sha256", required=True)
    parser.add_argument("--corrected-report-sha256", required=True)
    parser.add_argument("--fast-acceptance-sha256", required=True)
    parser.add_argument("--fast-science-sha256", required=True)
    parser.add_argument("--reviewed-acceptance-sha256", required=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    run = commands.add_parser("run", help="produce final acceptance and science records")
    add_common_arguments(run)
    run.add_argument("--python", type=Path, default=Path(sys.executable))
    validate = commands.add_parser("validate", help="revalidate final science records")
    add_common_arguments(validate)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "run":
            run_workflow(args)
        else:
            validate_workflow(args)
    except (
        CorrectedScienceError,
        OSError,
        ValueError,
        TypeError,
        KeyError,
        subprocess.SubprocessError,
    ) as error:
        print(f"corrected science error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
