#!/usr/bin/env python3
"""Orchestrate corrected-production Stage I downstream work.

The workflow accepts either a corrected-active-only assembled inventory or a
composite inventory that additionally reuses passive R06--R09 evidence.  Every
active case is reauthenticated against the corrected campaign identity and its
selected segment manifests before any downstream job is prepared.  Legacy
active evidence is therefore rejected even when it appears in an otherwise
well-formed composite inventory.

``prepare`` creates the all-snapshot literature-correct hyperbolicity and
analysis attempts through the existing launchers and creates CT/publication
Slurm jobs.  It never submits jobs.  ``submit`` is an explicit later action.
``retry-publication`` validates the latest successful upstream attempts before
preparing and submitting a fresh publication attempt bound to their job IDs.
``complete`` writes immutable corrected-production completion records only
after all products and their exact inventory, executable, and formula bindings
have been revalidated.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shlex
import stat
import subprocess
import sys
import tempfile
from typing import Callable, Iterable


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[2]
IDENTITY_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_corrected_identity.py")
HYPER_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_hyperbolicity.py")
ANALYSIS_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_analyze.py")
CT_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_ct_audit.py")
PUBLICATION_TOOL = SCRIPT_PATH.with_name("cgl_lf_stage_i_fast_publication.py")
CORRECTED_REPORT_TOOL = SCRIPT_PATH.with_name(
    "cgl_lf_stage_i_fast_corrected_report.py"
)

SCHEMA = "athenak-cgl-corrected-downstream-workflow"
SCHEMA_VERSION = 1
COMPLETION_SCHEMA = "athenak-cgl-corrected-downstream-completion"
POINTER_SCHEMA = "athenak-cgl-corrected-downstream-pointer"
FORMULA_ID = "literature-correct"
EXECUTABLE_FORMULA_ID = "literature-correct"
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
ALL_CASES = tuple(sorted((*ACTIVE_CASES, *PASSIVE_CASES)))
SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")
JOB_ID_PATTERN = re.compile(r"[1-9]\d*(?:;[A-Za-z0-9_.-]+)?")
SUBMITTED_JOB_ID_PATTERN = re.compile(r"[1-9]\d*")
ATTEMPT_PATTERN = re.compile(r"attempt-([0-9]{3})")


class CorrectedDownstreamError(RuntimeError):
    """A corrected-production downstream contract failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def stable_json(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    ).encode("utf-8")


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as stream:
            staged = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(staged, path)
        staged = None
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def write_json(path: Path, value: object) -> None:
    atomic_write(path, stable_json(value))


def write_immutable_json(path: Path, value: object) -> None:
    payload = stable_json(value)
    if path.exists():
        if not path.is_file() or path.read_bytes() != payload:
            raise CorrectedDownstreamError(
                f"refusing to replace differing retained artifact: {path}"
            )
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    staged: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=path.parent,
            prefix=f".{path.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            staged = Path(stream.name)
            stream.write(payload)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(staged, path)
    except FileExistsError as error:
        raise CorrectedDownstreamError(
            f"retained artifact appeared concurrently: {path}"
        ) from error
    finally:
        if staged is not None:
            staged.unlink(missing_ok=True)


def require_dict(value: object, label: str) -> dict[str, object]:
    if not isinstance(value, dict):
        raise CorrectedDownstreamError(f"{label} must be an object")
    return value


def require_list(value: object, label: str) -> list[object]:
    if not isinstance(value, list):
        raise CorrectedDownstreamError(f"{label} must be a list")
    return value


def require_text(value: object, label: str) -> str:
    if not isinstance(value, str) or not value:
        raise CorrectedDownstreamError(f"{label} must be a nonempty string")
    return value


def require_sha256(value: object, label: str) -> str:
    if not isinstance(value, str) or SHA256_PATTERN.fullmatch(value) is None:
        raise CorrectedDownstreamError(f"{label} must be a lowercase SHA-256")
    return value


def require_job_id(value: object, label: str) -> str:
    if not isinstance(value, str) or SUBMITTED_JOB_ID_PATTERN.fullmatch(value) is None:
        raise CorrectedDownstreamError(f"{label} must be a submitted Slurm job ID")
    return value


def require_job_ids(value: object, label: str) -> list[str]:
    job_ids = [
        require_job_id(item, f"{label} item")
        for item in require_list(value, label)
    ]
    if len(set(job_ids)) != len(job_ids):
        raise CorrectedDownstreamError(f"{label} contains duplicate job IDs")
    if job_ids != sorted(job_ids):
        raise CorrectedDownstreamError(f"{label} must be sorted")
    return job_ids


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.resolve(strict=False).relative_to(parent.resolve(strict=False))
        return True
    except ValueError:
        return False


def artifact_binding(path: Path) -> dict[str, object]:
    try:
        resolved = path.expanduser().resolve(strict=True)
    except OSError as error:
        raise CorrectedDownstreamError(f"artifact does not resolve: {path}") from error
    before = resolved.stat()
    if not stat.S_ISREG(before.st_mode):
        raise CorrectedDownstreamError(f"artifact is not a regular file: {resolved}")
    digest = sha256_file(resolved)
    after = resolved.stat()
    if (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    ) != (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    ):
        raise CorrectedDownstreamError(f"artifact changed while hashing: {resolved}")
    return {
        "path": str(resolved),
        "size_bytes": after.st_size,
        "sha256": digest,
    }


def verify_binding(value: object, label: str) -> dict[str, object]:
    declared = require_dict(value, f"{label} binding")
    path = Path(require_text(declared.get("path"), f"{label} path"))
    current = artifact_binding(path)
    if declared.get("sha256") != current["sha256"]:
        raise CorrectedDownstreamError(f"{label} SHA-256 differs")
    if declared.get("size_bytes") != current["size_bytes"]:
        raise CorrectedDownstreamError(f"{label} size differs")
    return current


def binding_identity(value: object, label: str) -> tuple[str, int, str]:
    binding = require_dict(value, f"{label} binding")
    path = str(Path(require_text(binding.get("path"), f"{label} path")).resolve())
    size = binding.get("size_bytes")
    if not isinstance(size, int) or isinstance(size, bool) or size < 0:
        raise CorrectedDownstreamError(f"{label} size_bytes is malformed")
    digest = require_sha256(binding.get("sha256"), f"{label} SHA-256")
    return path, size, digest


def same_binding(first: object, second: object, label: str) -> bool:
    return binding_identity(first, f"{label} first") == binding_identity(
        second, f"{label} second"
    )


def same_binding_lists(first: object, second: object, label: str) -> bool:
    left = require_list(first, f"{label} first list")
    right = require_list(second, f"{label} second list")
    return len(left) == len(right) and all(
        same_binding(a, b, f"{label} item {index}")
        for index, (a, b) in enumerate(zip(left, right))
    )


def load_bound_json(
    path: Path, label: str, expected_sha256: str | None = None
) -> tuple[dict[str, object], dict[str, object]]:
    binding = artifact_binding(path)
    if expected_sha256 is not None and binding["sha256"] != require_sha256(
        expected_sha256, f"{label} expected SHA-256"
    ):
        raise CorrectedDownstreamError(f"{label} SHA-256 differs from authority")
    try:
        value = json.loads(Path(str(binding["path"])).read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise CorrectedDownstreamError(f"{label} is not valid UTF-8 JSON") from error
    return require_dict(value, label), binding


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise CorrectedDownstreamError(f"cannot import workflow dependency: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def identity_validation(identity_path: Path) -> dict[str, object]:
    identity = load_module("_cgl_corrected_downstream_identity", IDENTITY_TOOL)
    try:
        return identity.validate_identity(identity_path)
    except identity.IdentityError as error:
        raise CorrectedDownstreamError(str(error)) from error


def passive_delta(case: dict[str, object], case_id: str) -> bool:
    model = require_dict(case.get("model_choices"), f"{case_id} model choices")
    value = require_text(model.get("passive_delta"), f"{case_id} passive_delta").lower()
    if value not in {"true", "false"}:
        raise CorrectedDownstreamError(f"{case_id} passive_delta is not boolean")
    return value == "true"


def selected_active_roots(identity: dict[str, object]) -> dict[str, Path]:
    root = Path(require_text(identity.get("campaign_root"), "campaign root"))
    selected: dict[str, Path] = {}
    for relative_value in require_list(
        identity.get("selected_run_paths"), "identity selected run paths"
    ):
        relative = Path(require_text(relative_value, "selected run path"))
        case_id = relative.name
        if case_id not in ACTIVE_CASES:
            raise CorrectedDownstreamError(
                f"identity selects a non-active corrected run path: {relative}"
            )
        if case_id in selected:
            raise CorrectedDownstreamError(
                f"identity selects duplicate corrected run roots for {case_id}"
            )
        selected[case_id] = (root / relative).resolve(strict=False)
    if set(selected) != set(ACTIVE_CASES):
        raise CorrectedDownstreamError(
            "identity must select exactly corrected active R02-R05 and R10-R17"
        )
    return selected


def active_manifest_bindings(
    case_id: str,
    case: dict[str, object],
    identity: dict[str, object],
    identity_path: Path,
    identity_sha256: str,
    selected_root: Path,
) -> list[dict[str, object]]:
    artifacts = require_dict(identity.get("artifacts"), "identity artifacts")
    executable = require_dict(artifacts.get("executable"), "identity executable")
    matrix = require_dict(artifacts.get("matrix"), "identity matrix")
    exact = {
        "campaign_id": identity.get("campaign_id"),
        "campaign_root": identity.get("campaign_root"),
        "root": identity.get("campaign_root"),
        "source_revision": identity.get("source_revision"),
        "campaign_identity": str(identity_path),
        "campaign_identity_sha256": identity_sha256,
        "matrix_sha256": matrix.get("sha256"),
        "executable": executable.get("path"),
        "executable_sha256": executable.get("sha256"),
        "legacy_restart_permitted": False,
        "case_id": case_id,
    }
    identities = require_dict(
        case.get("lineage_identities"), f"{case_id} lineage identities"
    )
    if identities.get("matrix_sha256") != [matrix.get("sha256")]:
        raise CorrectedDownstreamError(f"{case_id} matrix lineage is not corrected")
    if identities.get("executable_sha256") != [executable.get("sha256")]:
        raise CorrectedDownstreamError(f"{case_id} executable lineage is not corrected")

    lineage = require_list(case.get("lineage"), f"{case_id} selected lineage")
    if not lineage:
        raise CorrectedDownstreamError(f"{case_id} selected lineage is empty")
    bindings: list[dict[str, object]] = []
    sequences: list[int] = []
    for order, value in enumerate(lineage):
        segment = require_dict(value, f"{case_id} lineage segment {order}")
        if segment.get("kind") != "fast" or segment.get("order") != order:
            raise CorrectedDownstreamError(
                f"{case_id} active lineage contains legacy or malformed segment {order}"
            )
        segment_dir = Path(
            require_text(segment.get("segment_dir"), f"{case_id} segment directory")
        ).resolve(strict=False)
        if not is_relative_to(segment_dir, selected_root):
            raise CorrectedDownstreamError(
                f"{case_id} active segment escapes corrected selected run root"
            )
        binding = verify_binding(
            segment.get("manifest"), f"{case_id} segment {order} manifest"
        )
        manifest, _ = load_bound_json(
            Path(str(binding["path"])), f"{case_id} segment {order} manifest"
        )
        for key, expected in exact.items():
            if manifest.get(key) != expected:
                raise CorrectedDownstreamError(
                    f"{case_id} active segment {order} {key} is not corrected"
                )
        sequence = manifest.get("sequence")
        if not isinstance(sequence, int) or isinstance(sequence, bool):
            raise CorrectedDownstreamError(f"{case_id} segment sequence is malformed")
        sequences.append(sequence)
        if order == 0:
            if (
                sequence != 0
                or float(manifest.get("start_time", -1.0)) != 0.0
                or manifest.get("restart") is not None
                or manifest.get("launch_origin") != "fresh_t0_corrected_production"
                or manifest.get("fresh_lineage_root") is not True
            ):
                raise CorrectedDownstreamError(
                    f"{case_id} corrected active lineage does not start fresh at t=0"
                )
        else:
            restart = manifest.get("restart")
            if (
                sequence != sequences[-2] + 1
                or not isinstance(restart, str)
                or not is_relative_to(Path(restart), selected_root)
                or manifest.get("launch_origin")
                != "corrected_campaign_continuation"
                or manifest.get("fresh_lineage_root") is not False
            ):
                raise CorrectedDownstreamError(
                    f"{case_id} corrected continuation lineage is invalid"
                )
        bindings.append(binding)
    return bindings


def validate_inventory(
    identity_path: Path,
    inventory_path: Path,
    expected_inventory_sha256: str,
) -> dict[str, object]:
    identity_result = identity_validation(identity_path)
    identity, identity_binding = load_bound_json(
        Path(str(identity_result["identity"]["path"])),
        "corrected campaign identity",
        str(identity_result["identity"]["sha256"]),
    )
    identity_artifacts: dict[str, dict[str, object]] = {}
    for name, value in require_dict(
        identity_result.get("artifacts"), "identity artifacts"
    ).items():
        declared = require_dict(value, f"identity {name} artifact")
        current = artifact_binding(Path(require_text(
            declared.get("path"), f"identity {name} artifact path"
        )))
        if current["sha256"] != declared.get("sha256"):
            raise CorrectedDownstreamError(
                f"identity {name} artifact SHA-256 differs"
            )
        identity_artifacts[name] = current
    inventory, inventory_binding = load_bound_json(
        inventory_path, "corrected/composite inventory", expected_inventory_sha256
    )
    if inventory.get("schema_version") != 1:
        raise CorrectedDownstreamError("assembled inventory schema differs")
    output = Path(
        require_text(inventory.get("output"), "assembled inventory output")
    ).resolve(strict=True)
    if Path(str(inventory_binding["path"])) != output / "inventory.json":
        raise CorrectedDownstreamError(
            "assembled inventory is outside its declared output root"
        )
    campaign_root = Path(str(identity_result["campaign_root"])).resolve()
    if is_relative_to(output, campaign_root / "runs"):
        raise CorrectedDownstreamError("assembled inventory is beneath corrected runs/")

    matrix = verify_binding(inventory.get("matrix"), "assembled matrix")
    identity_matrix = require_dict(identity_artifacts.get("matrix"), "identity matrix")
    if matrix["sha256"] != identity_matrix.get("sha256"):
        raise CorrectedDownstreamError("assembled inventory matrix is not identity-bound")
    reporter = verify_binding(inventory.get("adapter"), "assembled reporter")
    record_type = inventory.get("record_type")
    composite_validation: dict[str, object] | None = None
    composite_adapter: dict[str, object] | None = None
    if record_type == "cgl_lf_stage_i_corrected_composite_report":
        if not CORRECTED_REPORT_TOOL.is_file():
            raise CorrectedDownstreamError(
                "corrected composite inventory validator is unavailable"
            )
        composite_adapter = verify_binding(
            inventory.get("composite_adapter"), "corrected composite adapter"
        )
        if Path(str(composite_adapter["path"])) != CORRECTED_REPORT_TOOL.resolve():
            raise CorrectedDownstreamError(
                "corrected composite adapter path differs from reviewed tool"
            )
        composite = load_module(
            "_cgl_corrected_downstream_composite_validator",
            CORRECTED_REPORT_TOOL,
        )
        try:
            composite_validation = composite.validate_composite_inventory(
                composite.DEFAULT_CONFIG, output
            )
        except composite.CompositeReportError as error:
            raise CorrectedDownstreamError(str(error)) from error
    elif record_type is not None:
        raise CorrectedDownstreamError(
            f"unsupported corrected inventory record_type: {record_type}"
        )
    cases = require_dict(inventory.get("cases"), "assembled cases")
    unsupported = set(cases) - set(ALL_CASES)
    if unsupported:
        raise CorrectedDownstreamError(
            f"assembled inventory contains unsupported cases: {sorted(unsupported)}"
        )
    missing_active = set(ACTIVE_CASES) - set(cases)
    if missing_active:
        raise CorrectedDownstreamError(
            f"assembled inventory lacks corrected active cases: {sorted(missing_active)}"
        )
    present_passive = set(cases) & set(PASSIVE_CASES)
    if present_passive and composite_validation is None:
        raise CorrectedDownstreamError(
            "passive controls require an authenticated corrected composite inventory"
        )
    if composite_validation is not None and set(cases) != set(ALL_CASES):
        raise CorrectedDownstreamError(
            "corrected composite inventory must contain exact R02-R17 coverage"
        )

    selected_roots = selected_active_roots(identity)
    active_bindings: dict[str, list[dict[str, object]]] = {}
    passive_executables: dict[str, list[str]] = {}
    selected_cases = sorted(cases)
    for case_id in selected_cases:
        case = require_dict(cases[case_id], f"{case_id} assembled case")
        lineage_path = output / "cases" / case_id / "lineage.json"
        lineage, _ = load_bound_json(lineage_path, f"{case_id} lineage")
        if lineage != case:
            raise CorrectedDownstreamError(
                f"{case_id} inventory record differs from lineage.json"
            )
        if case.get("case_id") != case_id or case.get("status") != "complete":
            raise CorrectedDownstreamError(f"{case_id} is not a complete assembled case")
        if case.get("errors") != []:
            raise CorrectedDownstreamError(f"{case_id} retains assembly errors")
        is_passive = passive_delta(case, case_id)
        if case_id in ACTIVE_CASES:
            if is_passive:
                raise CorrectedDownstreamError(
                    f"{case_id} corrected active case is marked passive"
                )
            active_bindings[case_id] = active_manifest_bindings(
                case_id,
                case,
                identity,
                Path(str(identity_binding["path"])),
                str(identity_binding["sha256"]),
                selected_roots[case_id],
            )
        else:
            if not is_passive:
                raise CorrectedDownstreamError(
                    f"{case_id} reused passive case is marked active"
                )
            identities = require_dict(
                case.get("lineage_identities"), f"{case_id} lineage identities"
            )
            executables = require_list(
                identities.get("executable_sha256"),
                f"{case_id} passive executable identities",
            )
            passive_executables[case_id] = [
                require_sha256(value, f"{case_id} passive executable SHA-256")
                for value in executables
            ]

    return {
        "identity": identity,
        "identity_validation": identity_result,
        "identity_artifacts": identity_artifacts,
        "identity_binding": identity_binding,
        "inventory": inventory,
        "inventory_binding": inventory_binding,
        "inventory_output": output,
        "matrix": matrix,
        "reporter": reporter,
        "inventory_kind": (
            "corrected-composite" if composite_validation is not None
            else "corrected-active-only"
        ),
        "composite_adapter": composite_adapter,
        "composite_validation": composite_validation,
        "cases": cases,
        "selected_cases": selected_cases,
        "active_cases": list(ACTIVE_CASES),
        "passive_cases": sorted(set(selected_cases) & set(PASSIVE_CASES)),
        "active_manifest_bindings": active_bindings,
        "passive_executable_sha256": passive_executables,
    }


def tool_bindings(context: dict[str, object]) -> dict[str, dict[str, object]]:
    tools = {
        "orchestrator": SCRIPT_PATH,
        "identity": IDENTITY_TOOL,
        "hyperbolicity": HYPER_TOOL,
        "analysis": ANALYSIS_TOOL,
        "ct": CT_TOOL,
        "publication": PUBLICATION_TOOL,
    }
    bindings = {name: artifact_binding(path) for name, path in tools.items()}
    if context.get("inventory_kind") == "corrected-composite":
        bindings["corrected_report"] = artifact_binding(CORRECTED_REPORT_TOOL)
    audit = require_dict(context["identity_artifacts"], "identity artifacts").get(
        "audit"
    )
    bindings["audit"] = verify_binding(audit, "identity audit")
    hyper = load_module("_cgl_corrected_downstream_hyper_preflight", HYPER_TOOL)
    if Path(hyper.DEFAULT_AUDIT).resolve() != Path(str(bindings["audit"]["path"])):
        raise CorrectedDownstreamError(
            "hyperbolicity launcher default audit differs from corrected identity"
        )
    return bindings


def exact_formula_binding(context: dict[str, object]) -> dict[str, object]:
    artifacts = require_dict(context["identity_artifacts"], "identity artifacts")
    return {
        "audit_formula_id": FORMULA_ID,
        "executable_formula_id": EXECUTABLE_FORMULA_ID,
        "compatibility": "compatible",
        "audit_script": artifacts["audit"],
        "corrected_eos": artifacts["eos"],
    }


def require_workflow_root(path: Path, context: dict[str, object]) -> Path:
    root = path.expanduser().resolve(strict=False)
    campaign = Path(
        require_text(
            require_dict(context["identity_validation"], "identity validation").get(
                "campaign_root"
            ),
            "campaign root",
        )
    )
    output = Path(str(context["inventory_output"]))
    if is_relative_to(root, campaign / "runs"):
        raise CorrectedDownstreamError("workflow root must be outside corrected runs/")
    if is_relative_to(root, output) or is_relative_to(output, root):
        raise CorrectedDownstreamError(
            "workflow root and assembled inventory output must be separate"
        )
    return root


def workflow_commands(
    context: dict[str, object],
    workflow_root: Path,
    python: Path,
    account: str,
    partition: str,
    walltime: str,
    cpus_per_task: int,
) -> dict[str, object]:
    inventory = Path(str(context["inventory_binding"]["path"]))
    inventory_sha = str(context["inventory_binding"]["sha256"])
    analysis = Path(str(context["inventory_output"]))
    selected = [str(value) for value in context["selected_cases"]]
    active = [str(value) for value in context["active_cases"]]
    downstream_root = analysis / "corrected-downstream"
    hyper_jobs = downstream_root / "hyperbolicity"
    analysis_jobs = downstream_root / "analysis"
    ct_output = workflow_root / "ct"
    publication_output = workflow_root / "publication"
    hyper = [
        str(python),
        str(HYPER_TOOL),
        str(inventory),
        "--inventory-sha256",
        inventory_sha,
        "--jobs-dir",
        str(hyper_jobs),
        "--formula",
        FORMULA_ID,
        "launch",
        *active,
        "--exploratory",
        "--snapshot-policy",
        "all",
        "--account",
        account,
        "--partition",
        partition,
        "--walltime",
        walltime,
        "--cpus-per-task",
        str(cpus_per_task),
        "--python",
        str(python),
    ]
    analysis_command = [
        str(python),
        str(ANALYSIS_TOOL),
        str(inventory),
        "--jobs-dir",
        str(analysis_jobs),
        "launch",
        *selected,
        "--account",
        account,
        "--partition",
        partition,
        "--walltime",
        walltime,
        "--cpus-per-task",
        str(cpus_per_task),
        "--python",
        str(python),
    ]
    ct = [
        str(python),
        str(CT_TOOL),
        "--inventory",
        str(inventory),
        "--expected-inventory-sha256",
        inventory_sha,
        "--output",
        str(ct_output),
        "--cases",
        ",".join(selected),
        "--snapshot-policy",
        "all",
    ]
    publication = [
        str(python),
        str(PUBLICATION_TOOL),
        str(analysis),
        "--output",
        str(publication_output),
        "--acceptance",
        str(ct_output),
    ]
    return {
        "hyperbolicity": {"jobs_dir": str(hyper_jobs), "command": hyper},
        "analysis": {"jobs_dir": str(analysis_jobs), "command": analysis_command},
        "ct": {
            "job_dir": str(workflow_root / "jobs/ct/attempt-000"),
            "output": str(ct_output),
            "command": ct,
        },
        "publication": {
            "job_dir": str(workflow_root / "jobs/publication/attempt-000"),
            "output": str(publication_output),
            "command": publication,
        },
    }


def run_checked(command: list[str]) -> None:
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise CorrectedDownstreamError(
            f"downstream preparation command failed ({completed.returncode}): "
            f"{shlex.join(command)}\n{detail}"
        )
    if completed.stdout.strip():
        print(completed.stdout.strip())


def attempt_directories(jobs: Path, case_id: str) -> list[Path]:
    return sorted(
        path
        for path in (jobs / case_id).glob("attempt-[0-9][0-9][0-9]")
        if path.is_dir() and (path / "manifest.json").is_file()
    )


def matching_attempt(
    jobs: Path,
    case_id: str,
    predicate: Callable[[dict[str, object]], bool],
) -> Path:
    for attempt in reversed(attempt_directories(jobs, case_id)):
        manifest, _ = load_bound_json(attempt / "manifest.json", f"{case_id} attempt")
        if predicate(manifest):
            return attempt
    raise CorrectedDownstreamError(f"{case_id} lacks a matching prepared attempt")


def prepared_attempts(
    commands: dict[str, object], context: dict[str, object]
) -> dict[str, dict[str, str]]:
    inventory_path = str(context["inventory_binding"]["path"])
    inventory_sha = str(context["inventory_binding"]["sha256"])
    hyper_jobs = Path(str(require_dict(commands["hyperbolicity"], "hyper stage")["jobs_dir"]))
    analysis_jobs = Path(str(require_dict(commands["analysis"], "analysis stage")["jobs_dir"]))
    result: dict[str, dict[str, str]] = {"hyperbolicity": {}, "analysis": {}}
    for case_id in context["active_cases"]:
        attempt = matching_attempt(
            hyper_jobs,
            str(case_id),
            lambda manifest: (
                manifest.get("record_type")
                == "cgl_lf_stage_i_direct_fast_hyperbolicity_job"
                and require_dict(manifest.get("inventory"), "hyper inventory").get("path")
                == inventory_path
                and require_dict(manifest.get("inventory"), "hyper inventory").get("sha256")
                == inventory_sha
                and manifest.get("formula_id") == FORMULA_ID
                and manifest.get("snapshot_policy") == "all"
                and require_dict(
                    manifest.get("snapshot_coverage"), "hyper snapshot coverage"
                ).get("all_complete_retained_snapshots_selected")
                is True
            ),
        )
        result["hyperbolicity"][str(case_id)] = str(attempt)
    for case_id in context["selected_cases"]:
        attempt = matching_attempt(
            analysis_jobs,
            str(case_id),
            lambda manifest: (
                manifest.get("record_type")
                == "cgl_lf_stage_i_direct_fast_analysis_job"
                and require_dict(
                    manifest.get("inventory"), "analysis inventory"
                ).get("path")
                == inventory_path
                and require_dict(
                    manifest.get("inventory"), "analysis inventory"
                ).get("sha256")
                == inventory_sha
            ),
        )
        result["analysis"][str(case_id)] = str(attempt)
    return result


def formula_executable_evidence(context: dict[str, object]) -> dict[str, object]:
    artifacts = require_dict(context["identity_artifacts"], "identity artifacts")
    return {
        "executable_formula_id": EXECUTABLE_FORMULA_ID,
        "formula_executable_compatibility": "compatible",
        "formula_executable_compatibility_reason": (
            "the corrected campaign identity binds the ppar^2 EOS source and "
            "the exact executable used by every selected active lineage"
        ),
        "formula_executable_binding": {
            "campaign_identity": context["identity_binding"],
            "corrected_executable": artifacts["executable"],
            "corrected_eos": artifacts["eos"],
            "binder": artifact_binding(SCRIPT_PATH),
        },
    }


def validate_formula_executable_evidence(
    manifest: dict[str, object],
    result: dict[str, object],
    context: dict[str, object],
    case_id: str,
) -> None:
    expected = formula_executable_evidence(context)
    if manifest.get("formula_executable_evidence") != expected:
        raise CorrectedDownstreamError(
            f"{case_id} hyperbolicity manifest lacks exact formula/executable binding"
        )
    provenance = require_dict(result.get("provenance"), f"{case_id} result provenance")
    for key, value in expected.items():
        if provenance.get(key) != value:
            raise CorrectedDownstreamError(
                f"{case_id} result {key} differs from corrected formula binding"
            )


def enrich_hyper_result(identity_path: Path, manifest_path: Path) -> Path:
    manifest, _ = load_bound_json(manifest_path, "hyperbolicity manifest")
    inventory = verify_binding(manifest.get("inventory"), "hyperbolicity inventory")
    context = validate_inventory(
        identity_path,
        Path(str(inventory["path"])),
        str(inventory["sha256"]),
    )
    case_id = require_text(manifest.get("case_id"), "hyperbolicity case ID")
    if case_id not in ACTIVE_CASES:
        raise CorrectedDownstreamError(
            f"formula/executable binding is forbidden for non-active case {case_id}"
        )
    if (
        manifest.get("formula_id") != FORMULA_ID
        or manifest.get("snapshot_policy") != "all"
        or not same_binding_lists(
            manifest.get("segment_manifests"),
            context["active_manifest_bindings"][case_id],
            f"{case_id} corrected segment manifests",
        )
    ):
        raise CorrectedDownstreamError(
            f"{case_id} hyperbolicity manifest is not corrected-production exact"
        )
    hyper = load_module("_cgl_corrected_downstream_hyper_binder", HYPER_TOOL)
    attempt = manifest_path.resolve().parent
    validation_error = hyper.validate_result(attempt, manifest)
    if validation_error is not None:
        raise CorrectedDownstreamError(
            f"{case_id} result cannot be formula-bound: {validation_error}"
        )
    result_record = require_dict(manifest.get("result"), "hyperbolicity result")
    result_path = Path(require_text(result_record.get("path"), "result path"))
    result_sha_path = Path(
        require_text(result_record.get("sha256_path"), "result SHA-256 path")
    )
    result, _ = load_bound_json(result_path, f"{case_id} hyperbolicity result")
    provenance = require_dict(result.get("provenance"), f"{case_id} result provenance")
    evidence = formula_executable_evidence(context)
    for key, value in evidence.items():
        existing = provenance.get(key)
        if existing is not None and existing != value:
            raise CorrectedDownstreamError(
                f"{case_id} result retains conflicting {key}"
            )
        provenance[key] = value
    result["provenance"] = provenance
    write_json(result_path, result)
    atomic_write(
        result_sha_path,
        f"{sha256_file(result_path)}  {result_path.name}\n".encode("ascii"),
    )
    retained, _ = load_bound_json(result_path, f"{case_id} enriched result")
    validate_formula_executable_evidence(manifest, retained, context, case_id)
    validation_error = hyper.validate_result(attempt, manifest)
    if validation_error is not None:
        raise CorrectedDownstreamError(
            f"{case_id} enriched result is invalid: {validation_error}"
        )
    return result_path


def patch_hyper_attempt(attempt: Path, context: dict[str, object]) -> None:
    manifest_path = attempt / "manifest.json"
    manifest, _ = load_bound_json(manifest_path, "hyperbolicity manifest")
    case_id = require_text(manifest.get("case_id"), "hyperbolicity case ID")
    evidence = formula_executable_evidence(context)
    existing = manifest.get("formula_executable_evidence")
    if existing is not None and existing != evidence:
        raise CorrectedDownstreamError(
            f"{case_id} hyperbolicity manifest has conflicting formula binding"
        )
    if existing is None:
        if manifest.get("job_id") is not None or (attempt / "exit_code.txt").exists():
            raise CorrectedDownstreamError(
                f"{case_id} hyperbolicity attempt cannot be patched after execution"
            )
        manifest["formula_executable_evidence"] = evidence
        write_json(manifest_path, manifest)

    script = attempt / "run.sbatch"
    text = script.read_text(encoding="utf-8")
    marker = "# corrected formula/executable binding"
    if marker in text:
        return
    if manifest.get("job_id") is not None or (attempt / "exit_code.txt").exists():
        raise CorrectedDownstreamError(
            f"{case_id} hyperbolicity batch script cannot be patched after execution"
        )
    exit_path = attempt / "exit_code.txt"
    needle = (
        f"printf '%s\\n' \"$status\" > {shlex.quote(str(exit_path))}\n"
        "exit \"$status\"\n"
    )
    if text.count(needle) != 1:
        raise CorrectedDownstreamError(
            f"{case_id} hyperbolicity batch-script exit block differs"
        )
    binder_command = [
        require_text(manifest.get("python"), "hyperbolicity Python"),
        str(SCRIPT_PATH),
        "bind-hyper-result",
        "--identity",
        str(context["identity_binding"]["path"]),
        "--manifest",
        str(manifest_path),
    ]
    replacement = "\n".join([
        marker,
        'if [[ "$status" -eq 0 ]]; then',
        "  set +e",
        f"  {shlex.join(binder_command)}",
        "  status=$?",
        "  set -e",
        "fi",
        f"printf '%s\\n' \"$status\" > {shlex.quote(str(exit_path))}",
        'exit "$status"',
        "",
    ])
    atomic_write(script, text.replace(needle, replacement, 1).encode("utf-8"))
    script.chmod(0o755)


def patch_hyper_attempts(
    attempts: dict[str, dict[str, str]], context: dict[str, object]
) -> None:
    for attempt in require_dict(
        attempts.get("hyperbolicity"), "prepared hyperbolicity attempts"
    ).values():
        patch_hyper_attempt(Path(require_text(attempt, "hyperbolicity attempt")), context)


def generic_job_manifest(
    stage: str,
    job_dir: Path,
    command: list[str],
    context: dict[str, object],
    tools: dict[str, dict[str, object]],
    account: str,
    partition: str,
    walltime: str,
    cpus_per_task: int,
) -> dict[str, object]:
    executable = require_dict(context["identity_artifacts"], "identity artifacts")[
        "executable"
    ]
    return {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_corrected_downstream_job",
        "stage": stage,
        "job_dir": str(job_dir),
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "corrected_executable": executable,
        "formula_id": FORMULA_ID,
        "executable_formula_id": EXECUTABLE_FORMULA_ID,
        "formula_binding": exact_formula_binding(context),
        "tools": tools,
        "command": command,
        "account": account,
        "partition": partition,
        "walltime": walltime,
        "nodes": 1,
        "ntasks": 1,
        "cpus_per_task": cpus_per_task,
        "job_name": f"cglcd_{stage}",
        "slurm_log": str(job_dir / "%x.%j.log"),
        "job_id": None,
    }


def generic_batch_script(manifest: dict[str, object]) -> str:
    command = [
        require_text(value, "job command argument")
        for value in require_list(manifest.get("command"), "job command")
    ]
    job_dir = Path(require_text(manifest.get("job_dir"), "job directory"))
    checks: list[str] = []
    for label, binding_value in (
        ("campaign_identity", manifest["campaign_identity"]),
        ("inventory", manifest["inventory"]),
        ("corrected_executable", manifest["corrected_executable"]),
    ):
        binding = require_dict(binding_value, f"{label} binding")
        checks.append(
            f"require_sha {binding['sha256']} {shlex.quote(str(binding['path']))} {label}"
        )
    return "\n".join([
        "#!/bin/bash",
        f"#SBATCH --account={manifest['account']}",
        f"#SBATCH --partition={manifest['partition']}",
        f"#SBATCH --time={manifest['walltime']}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={manifest['cpus_per_task']}",
        f"#SBATCH --job-name={manifest['job_name']}",
        f"#SBATCH --output={manifest['slurm_log']}",
        "",
        "set -uo pipefail",
        "require_sha() {",
        "  expected=$1",
        "  path=$2",
        "  label=$3",
        "  observed=$(/usr/bin/sha256sum \"$path\" | /usr/bin/awk '{print $1}')",
        "  test \"$observed\" = \"$expected\" || {",
        "    echo \"$label SHA-256 differs: $observed\" >&2",
        "    exit 1",
        "  }",
        "}",
        *checks,
        f"cd {shlex.quote(str(REPO_ROOT))}",
        f"export OMP_NUM_THREADS={int(manifest['cpus_per_task'])}",
        "set +e",
        (
            "srun --nodes=1 --ntasks=1 "
            f"--cpus-per-task={int(manifest['cpus_per_task'])} "
            + shlex.join(command)
        ),
        "status=$?",
        f"printf '%s\\n' \"$status\" > {shlex.quote(str(job_dir / 'exit_code.txt'))}",
        "exit \"$status\"",
        "",
    ])


def prepare_generic_job(manifest: dict[str, object]) -> None:
    job_dir = Path(str(manifest["job_dir"]))
    job_dir.mkdir(parents=True, exist_ok=True)
    write_immutable_json(job_dir / "manifest.json", manifest)
    script = job_dir / "run.sbatch"
    payload = generic_batch_script(manifest).encode("utf-8")
    if script.exists() and script.read_bytes() != payload:
        raise CorrectedDownstreamError(f"refusing to replace differing job script: {script}")
    if not script.exists():
        atomic_write(script, payload)
        script.chmod(0o755)


def prepare_workflow(args: argparse.Namespace) -> Path:
    context = validate_inventory(args.identity, args.inventory, args.inventory_sha256)
    workflow_root = require_workflow_root(args.workflow_root, context)
    python = args.python.expanduser().resolve(strict=True)
    tools = tool_bindings(context)
    commands = workflow_commands(
        context,
        workflow_root,
        python,
        args.account,
        args.partition,
        args.walltime,
        args.cpus_per_task,
    )
    run_checked(list(require_dict(commands["hyperbolicity"], "hyper stage")["command"]))
    run_checked(list(require_dict(commands["analysis"], "analysis stage")["command"]))
    attempts = prepared_attempts(commands, context)
    patch_hyper_attempts(attempts, context)
    for stage in ("ct", "publication"):
        spec = require_dict(commands[stage], f"{stage} stage")
        prepare_generic_job(generic_job_manifest(
            stage,
            Path(str(spec["job_dir"])),
            list(spec["command"]),
            context,
            tools,
            args.account,
            args.partition,
            args.walltime,
            args.cpus_per_task,
        ))
    executable = require_dict(context["identity_artifacts"], "identity artifacts")[
        "executable"
    ]
    workflow = {
        "schema": SCHEMA,
        "schema_version": SCHEMA_VERSION,
        "campaign_kind": "corrected-production",
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "inventory_output": str(context["inventory_output"]),
        "workflow_root": str(workflow_root),
        "corrected_executable": executable,
        "formula_binding": exact_formula_binding(context),
        "inventory_kind": context["inventory_kind"],
        "case_classification": {
            "corrected_active": context["active_cases"],
            "reused_passive": context["passive_cases"],
            "selected": context["selected_cases"],
        },
        "active_execution_manifests": context["active_manifest_bindings"],
        "passive_executable_sha256": context["passive_executable_sha256"],
        "tools": tools,
        "commands": commands,
        "submission_policy": (
            "prepare never submits; submit is an explicit separate command; "
            "publication depends on hyperbolicity, analysis, and CT"
        ),
    }
    workflow_path = workflow_root / "workflow.json"
    write_immutable_json(workflow_path, workflow)
    print(f"prepared corrected downstream workflow: {workflow_path}")
    return workflow_path


def load_workflow(workflow_root: Path) -> tuple[dict[str, object], dict[str, object]]:
    root = workflow_root.expanduser().resolve(strict=True)
    workflow, binding = load_bound_json(root / "workflow.json", "corrected workflow")
    if (
        workflow.get("schema") != SCHEMA
        or workflow.get("schema_version") != SCHEMA_VERSION
        or workflow.get("campaign_kind") != "corrected-production"
        or workflow.get("workflow_root") != str(root)
    ):
        raise CorrectedDownstreamError("corrected workflow identity differs")
    return workflow, binding


def workflow_context(workflow: dict[str, object]) -> dict[str, object]:
    identity = verify_binding(workflow.get("campaign_identity"), "workflow identity")
    inventory = verify_binding(workflow.get("inventory"), "workflow inventory")
    context = validate_inventory(
        Path(str(identity["path"])),
        Path(str(inventory["path"])),
        str(inventory["sha256"]),
    )
    executable = verify_binding(
        workflow.get("corrected_executable"), "workflow corrected executable"
    )
    identity_executable = require_dict(
        context["identity_artifacts"], "identity artifacts"
    )["executable"]
    if executable != identity_executable:
        raise CorrectedDownstreamError(
            "workflow corrected executable differs from campaign identity"
        )
    formula = require_dict(workflow.get("formula_binding"), "workflow formula binding")
    if formula != exact_formula_binding(context):
        raise CorrectedDownstreamError("workflow formula/executable binding differs")
    if workflow.get("inventory_kind") != context["inventory_kind"]:
        raise CorrectedDownstreamError("workflow inventory kind differs")
    for binding in require_dict(workflow.get("tools"), "workflow tools").values():
        verify_binding(binding, "workflow tool")
    expected = {
        "corrected_active": context["active_cases"],
        "reused_passive": context["passive_cases"],
        "selected": context["selected_cases"],
    }
    if workflow.get("case_classification") != expected:
        raise CorrectedDownstreamError("workflow case classification differs")
    if workflow.get("active_execution_manifests") != context["active_manifest_bindings"]:
        raise CorrectedDownstreamError("workflow active execution bindings differ")
    return context


def submit_manifest_job(job_dir: Path, dependency_ids: Iterable[str] = ()) -> str:
    manifest_path = job_dir / "manifest.json"
    manifest, _ = load_bound_json(manifest_path, "job manifest")
    requested_dependencies = [
        require_job_id(dependency, "job dependency")
        for dependency in dependency_ids
    ]
    if len(set(requested_dependencies)) != len(requested_dependencies):
        raise CorrectedDownstreamError("job submission contains duplicate dependencies")
    dependencies = sorted(requested_dependencies)
    declared_dependencies = manifest.get("dependency_job_ids")
    if declared_dependencies is not None:
        if require_job_ids(
            declared_dependencies, "job manifest dependency IDs"
        ) != dependencies:
            raise CorrectedDownstreamError(
                f"job manifest dependencies differ from submission request: {job_dir}"
            )
    job_id = manifest.get("job_id")
    if isinstance(job_id, str):
        return require_job_id(job_id, "job manifest job ID")
    if job_id is not None:
        raise CorrectedDownstreamError(f"job manifest job ID is malformed: {job_dir}")
    if (job_dir / "exit_code.txt").is_file():
        raise CorrectedDownstreamError(f"cannot submit completed job again: {job_dir}")
    if dependencies and declared_dependencies is None:
        manifest["dependency_job_ids"] = dependencies
        write_json(manifest_path, manifest)
    command = ["/usr/bin/sbatch", "--parsable"]
    if dependencies:
        command.append(f"--dependency=afterok:{':'.join(dependencies)}")
    command.append(str(job_dir / "run.sbatch"))
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    response = completed.stdout.strip()
    if JOB_ID_PATTERN.fullmatch(response) is None:
        raise CorrectedDownstreamError(f"unexpected sbatch response: {response!r}")
    job_id = response.split(";", 1)[0]
    manifest["job_id"] = job_id
    manifest["submitted_utc"] = utc_now()
    write_json(manifest_path, manifest)
    return job_id


def submit_workflow(workflow_root: Path) -> Path:
    workflow, workflow_binding = load_workflow(workflow_root)
    context = workflow_context(workflow)
    attempts = prepared_attempts(
        require_dict(workflow.get("commands"), "workflow commands"), context
    )
    patch_hyper_attempts(attempts, context)
    upstream_ids: list[str] = []
    jobs: dict[str, object] = {"hyperbolicity": {}, "analysis": {}}
    for stage in ("hyperbolicity", "analysis"):
        values = require_dict(attempts.get(stage), f"{stage} attempts")
        for case_id, attempt in sorted(values.items()):
            job_id = submit_manifest_job(Path(require_text(attempt, f"{case_id} attempt")))
            require_dict(jobs[stage], f"{stage} submitted jobs")[case_id] = job_id
            upstream_ids.append(job_id)
    commands = require_dict(workflow.get("commands"), "workflow commands")
    ct_dir = Path(require_text(require_dict(commands["ct"], "ct stage")["job_dir"], "CT job"))
    ct_id = submit_manifest_job(ct_dir)
    upstream_ids.append(ct_id)
    if len(set(upstream_ids)) != len(upstream_ids):
        raise CorrectedDownstreamError("submitted upstream jobs contain duplicate IDs")
    publication_dependencies = sorted(upstream_ids)
    publication_dir = Path(require_text(
        require_dict(commands["publication"], "publication stage")["job_dir"],
        "publication job",
    ))
    publication_id = submit_manifest_job(publication_dir, publication_dependencies)
    jobs["ct"] = ct_id
    jobs["publication"] = publication_id
    submission = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_corrected_downstream_submission",
        "workflow": workflow_binding,
        "jobs": jobs,
        "publication_dependency": publication_dependencies,
    }
    path = Path(str(workflow["workflow_root"])) / "submission.json"
    write_immutable_json(path, submission)
    print(f"submitted corrected downstream workflow: {path}")
    return path


def publication_stage(workflow: dict[str, object]) -> dict[str, object]:
    return require_dict(
        require_dict(workflow.get("commands"), "workflow commands").get("publication"),
        "publication stage",
    )


def publication_attempt_directories(workflow: dict[str, object]) -> list[Path]:
    initial = Path(
        require_text(publication_stage(workflow).get("job_dir"), "publication job")
    )
    if initial.name != "attempt-000":
        raise CorrectedDownstreamError(
            "initial publication job directory must be attempt-000"
        )
    return sorted(
        path
        for path in initial.parent.glob("attempt-[0-9][0-9][0-9]")
        if path.is_dir() and ATTEMPT_PATTERN.fullmatch(path.name) is not None
    )


def publication_contract_manifest(
    workflow: dict[str, object],
    context: dict[str, object],
    template: dict[str, object],
    attempt: Path,
) -> dict[str, object]:
    stage = publication_stage(workflow)
    return generic_job_manifest(
        "publication",
        attempt,
        [
            require_text(value, "publication command argument")
            for value in require_list(stage.get("command"), "publication command")
        ],
        context,
        require_dict(workflow.get("tools"), "workflow tools"),
        require_text(template.get("account"), "publication account"),
        require_text(template.get("partition"), "publication partition"),
        require_text(template.get("walltime"), "publication walltime"),
        int(template["cpus_per_task"]),
    )


def publication_template_manifest(
    workflow: dict[str, object], context: dict[str, object]
) -> dict[str, object]:
    stage = publication_stage(workflow)
    initial = Path(require_text(stage.get("job_dir"), "publication job"))
    manifest, _ = load_bound_json(initial / "manifest.json", "publication template")
    cpus_per_task = manifest.get("cpus_per_task")
    if (
        not isinstance(cpus_per_task, int)
        or isinstance(cpus_per_task, bool)
        or cpus_per_task < 1
    ):
        raise CorrectedDownstreamError(
            "publication template cpus_per_task must be positive"
        )
    expected = publication_contract_manifest(workflow, context, manifest, initial)
    if any(
        manifest.get(key) != value
        for key, value in expected.items()
        if key != "job_id"
    ):
        raise CorrectedDownstreamError(
            "initial publication job differs from the corrected workflow contract"
        )
    if manifest.get("job_id") is not None:
        require_job_id(manifest.get("job_id"), "publication template job ID")
    return manifest


def publication_attempt_matches(
    workflow: dict[str, object],
    context: dict[str, object],
    template: dict[str, object],
    attempt: Path,
    manifest: dict[str, object],
) -> bool:
    expected = publication_contract_manifest(workflow, context, template, attempt)
    return all(
        manifest.get(key) == value
        for key, value in expected.items()
        if key != "job_id"
    )


def publication_attempt_dependency_ids(manifest: dict[str, object]) -> list[str]:
    return require_job_ids(
        manifest.get("dependency_job_ids"), "publication dependency job IDs"
    )


def bound_job_id(
    binding_value: object, expected_path: Path | None, label: str
) -> str:
    declared = require_dict(binding_value, f"{label} binding")
    path = Path(require_text(declared.get("path"), f"{label} path"))
    if expected_path is not None and path.resolve() != expected_path.resolve():
        raise CorrectedDownstreamError(f"{label} manifest path differs")
    manifest, current = load_bound_json(path, label)
    if not same_binding(declared, current, label):
        raise CorrectedDownstreamError(f"{label} binding differs")
    return require_job_id(manifest.get("job_id"), f"{label} job ID")


def successful_upstream_job_ids(
    hyper: dict[str, object],
    analysis: dict[str, object],
    ct: dict[str, object],
) -> list[str]:
    job_ids: list[str] = []
    for stage, records in (("hyperbolicity", hyper), ("analysis", analysis)):
        for case_id, value in sorted(records.items()):
            record = require_dict(value, f"{case_id} {stage} completion")
            attempt = Path(
                require_text(record.get("attempt"), f"{case_id} {stage} attempt")
            )
            job_ids.append(
                bound_job_id(
                    record.get("manifest"),
                    attempt / "manifest.json",
                    f"{case_id} {stage} manifest",
                )
            )
    job_ids.append(bound_job_id(ct.get("job_manifest"), None, "CT job manifest"))
    if len(set(job_ids)) != len(job_ids):
        raise CorrectedDownstreamError(
            "successful upstream attempts contain duplicate job IDs"
        )
    return sorted(job_ids)


def matching_publication_attempt(
    workflow: dict[str, object],
    context: dict[str, object],
    dependencies: list[str],
) -> tuple[Path, dict[str, object], dict[str, object]]:
    template = publication_template_manifest(workflow, context)
    for attempt in reversed(publication_attempt_directories(workflow)):
        manifest_path = attempt / "manifest.json"
        if not manifest_path.is_file():
            continue
        manifest, binding = load_bound_json(
            manifest_path, "publication attempt manifest"
        )
        if not publication_attempt_matches(
            workflow, context, template, attempt, manifest
        ):
            continue
        if publication_attempt_dependency_ids(manifest) != dependencies:
            continue
        require_job_id(manifest.get("job_id"), "publication attempt job ID")
        if not (attempt / "exit_code.txt").is_file():
            raise CorrectedDownstreamError(
                f"matching publication attempt is still pending: {attempt}"
            )
        if job_exit_code(attempt, "publication") == 0:
            return attempt, manifest, binding
    raise CorrectedDownstreamError(
        "publication lacks a successful attempt matching current upstream jobs"
    )


def job_exit_code(job_dir: Path, label: str) -> int:
    path = job_dir / "exit_code.txt"
    try:
        return int(path.read_text(encoding="utf-8").strip())
    except (OSError, ValueError) as error:
        raise CorrectedDownstreamError(f"{label} lacks a valid exit code") from error


def zero_exit(job_dir: Path, label: str) -> None:
    value = job_exit_code(job_dir, label)
    if value != 0:
        raise CorrectedDownstreamError(f"{label} failed with exit code {value}")


def completed_hyperbolicity(
    workflow: dict[str, object], context: dict[str, object]
) -> dict[str, object]:
    hyper = load_module("_cgl_corrected_downstream_hyper_complete", HYPER_TOOL)
    jobs = Path(str(require_dict(
        require_dict(workflow["commands"], "workflow commands")["hyperbolicity"],
        "hyper stage",
    )["jobs_dir"]))
    inventory = context["inventory_binding"]
    records: dict[str, object] = {}
    for case_id in context["active_cases"]:
        attempt = matching_attempt(
            jobs,
            str(case_id),
            lambda manifest: (
                manifest.get("formula_id") == FORMULA_ID
                and manifest.get("snapshot_policy") == "all"
                and same_binding(
                    manifest.get("inventory"), inventory, "hyper inventory"
                )
            ),
        )
        manifest, manifest_binding = load_bound_json(
            attempt / "manifest.json", f"{case_id} hyperbolicity manifest"
        )
        zero_exit(attempt, f"{case_id} hyperbolicity audit")
        error = hyper.validate_result(attempt, manifest)
        if error is not None:
            raise CorrectedDownstreamError(
                f"{case_id} hyperbolicity result is invalid: {error}"
            )
        if not same_binding_lists(
            manifest.get("segment_manifests"),
            context["active_manifest_bindings"][case_id],
            f"{case_id} corrected segment manifests",
        ):
            raise CorrectedDownstreamError(
                f"{case_id} hyperbolicity audit is not bound to corrected execution"
            )
        result_path = Path(str(require_dict(manifest["result"], "hyper result")["path"]))
        result_binding = artifact_binding(result_path)
        result, _ = load_bound_json(result_path, f"{case_id} hyperbolicity result")
        validate_formula_executable_evidence(manifest, result, context, str(case_id))
        records[str(case_id)] = {
            "attempt": str(attempt),
            "manifest": manifest_binding,
            "result": result_binding,
            "formula_id": FORMULA_ID,
            "executable_sha256": require_dict(
                context["identity_artifacts"], "identity artifacts"
            )["executable"]["sha256"],
            "formula_executable_compatibility": "compatible",
            "disposition": hyper.audit_disposition(attempt, manifest),
        }
    return records


def completed_analysis(
    workflow: dict[str, object], context: dict[str, object]
) -> dict[str, object]:
    analysis = load_module("_cgl_corrected_downstream_analysis_complete", ANALYSIS_TOOL)
    jobs = Path(str(require_dict(
        require_dict(workflow["commands"], "workflow commands")["analysis"],
        "analysis stage",
    )["jobs_dir"]))
    inventory = context["inventory_binding"]
    output = Path(str(context["inventory_output"]))
    records: dict[str, object] = {}
    for case_id in context["selected_cases"]:
        attempt = matching_attempt(
            jobs,
            str(case_id),
            lambda manifest: same_binding(
                manifest.get("inventory"), inventory, "analysis inventory"
            ),
        )
        manifest, manifest_binding = load_bound_json(
            attempt / "manifest.json", f"{case_id} analysis manifest"
        )
        zero_exit(attempt, f"{case_id} analysis")
        errors = analysis.analysis_completion_errors(attempt, manifest)
        if errors:
            raise CorrectedDownstreamError(
                f"{case_id} analysis is incomplete: {'; '.join(errors)}"
            )
        records[str(case_id)] = {
            "attempt": str(attempt),
            "manifest": manifest_binding,
            "diagnostics": artifact_binding(
                output / "cases" / str(case_id) / "diagnostics.json"
            ),
        }
    return records


def completed_ct(
    workflow: dict[str, object], context: dict[str, object]
) -> dict[str, object]:
    stage = require_dict(
        require_dict(workflow["commands"], "workflow commands")["ct"], "CT stage"
    )
    job_dir = Path(str(stage["job_dir"]))
    zero_exit(job_dir, "CT audit")
    audit_path = Path(str(stage["output"])) / "ct_audit.json"
    audit, binding = load_bound_json(audit_path, "CT audit")
    if not same_binding(audit.get("inventory"), context["inventory_binding"], "CT inventory"):
        raise CorrectedDownstreamError("CT audit inventory binding differs")
    selection = require_dict(audit.get("selection"), "CT selection")
    if (
        selection.get("cases") != context["selected_cases"]
        or selection.get("snapshot_policy") != "all"
        or audit.get("result") == "authentication_failed"
    ):
        raise CorrectedDownstreamError("CT audit coverage or authentication differs")
    cases = require_dict(audit.get("cases"), "CT cases")
    if set(cases) != set(context["selected_cases"]) or any(
        require_dict(value, "CT case").get("provenance_authenticated") is not True
        for value in cases.values()
    ):
        raise CorrectedDownstreamError("CT audit lacks authenticated selected coverage")
    return {"job_manifest": artifact_binding(job_dir / "manifest.json"), "audit": binding}


def retry_publication(workflow_root: Path) -> Path:
    workflow, _ = load_workflow(workflow_root)
    context = workflow_context(workflow)
    hyper = completed_hyperbolicity(workflow, context)
    analysis = completed_analysis(workflow, context)
    ct = completed_ct(workflow, context)
    dependencies = successful_upstream_job_ids(hyper, analysis, ct)
    template = publication_template_manifest(workflow, context)

    for attempt in reversed(publication_attempt_directories(workflow)):
        manifest_path = attempt / "manifest.json"
        if not manifest_path.is_file():
            continue
        manifest, _ = load_bound_json(manifest_path, "publication attempt manifest")
        if not publication_attempt_matches(
            workflow, context, template, attempt, manifest
        ):
            continue
        if publication_attempt_dependency_ids(manifest) != dependencies:
            continue
        job_id = manifest.get("job_id")
        if job_id is None:
            if (attempt / "exit_code.txt").exists():
                raise CorrectedDownstreamError(
                    f"unsubmitted publication attempt has an exit code: {attempt}"
                )
            submit_manifest_job(attempt, dependencies)
            print(f"submitted prepared publication retry: {attempt}")
            return attempt
        require_job_id(job_id, "publication retry job ID")
        if not (attempt / "exit_code.txt").is_file():
            raise CorrectedDownstreamError(
                f"matching publication attempt is still pending: {attempt}"
            )
        if job_exit_code(attempt, "publication") != 0:
            break
        raise CorrectedDownstreamError(
            f"matching publication attempt already succeeded: {attempt}"
        )

    attempts = publication_attempt_directories(workflow)
    last_index = max(
        int(ATTEMPT_PATTERN.fullmatch(attempt.name).group(1))
        for attempt in attempts
    )
    if last_index >= 999:
        raise CorrectedDownstreamError("publication attempt namespace is exhausted")
    attempt = attempts[0].parent / f"attempt-{last_index + 1:03d}"
    try:
        attempt.mkdir()
    except FileExistsError as error:
        raise CorrectedDownstreamError(
            f"publication retry attempt appeared concurrently: {attempt}"
        ) from error
    manifest = publication_contract_manifest(workflow, context, template, attempt)
    manifest["dependency_job_ids"] = dependencies
    prepare_generic_job(manifest)
    submit_manifest_job(attempt, dependencies)
    print(f"prepared and submitted fresh publication retry: {attempt}")
    return attempt


def completed_publication(
    workflow: dict[str, object],
    context: dict[str, object],
    hyper: dict[str, object],
    analysis: dict[str, object],
    ct: dict[str, object],
) -> dict[str, object]:
    dependencies = successful_upstream_job_ids(hyper, analysis, ct)
    job_dir, _, job_binding = matching_publication_attempt(
        workflow, context, dependencies
    )
    stage = publication_stage(workflow)
    manifest_path = Path(str(stage["output"])) / "manifest.json"
    manifest, binding = load_bound_json(manifest_path, "publication manifest")
    if (
        manifest.get("record_type") != "cgl_lf_stage_i_fast_publication_products"
        or manifest.get("analysis_output") != str(context["inventory_output"])
    ):
        raise CorrectedDownstreamError("publication manifest identity differs")
    source_paths = {
        require_text(require_dict(value, "publication source").get("path"), "source path")
        for value in require_list(manifest.get("sources"), "publication sources")
    }
    required_sources = {str(ct["audit"]["path"])}
    required_sources.update(
        str(record[key]["path"])
        for record in hyper.values()
        for key in ("manifest", "result")
    )
    required_sources.update(
        str(record["diagnostics"]["path"]) for record in analysis.values()
    )
    if not required_sources <= source_paths:
        raise CorrectedDownstreamError(
            "publication products omit corrected downstream evidence"
        )
    for value in require_list(manifest.get("products"), "publication products"):
        verify_binding(value, "publication product")
    return {
        "attempt": str(job_dir),
        "job_manifest": job_binding,
        "dependency_job_ids": dependencies,
        "manifest": binding,
    }


def completion_record(
    workflow: dict[str, object],
    workflow_binding: dict[str, object],
    context: dict[str, object],
) -> dict[str, object]:
    hyper = completed_hyperbolicity(workflow, context)
    analysis = completed_analysis(workflow, context)
    ct = completed_ct(workflow, context)
    publication = completed_publication(workflow, context, hyper, analysis, ct)
    executable = require_dict(context["identity_artifacts"], "identity artifacts")[
        "executable"
    ]
    return {
        "schema": COMPLETION_SCHEMA,
        "schema_version": 1,
        "status": "complete",
        "campaign_kind": "corrected-production",
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "workflow": workflow_binding,
        "case_classification": workflow["case_classification"],
        "corrected_executable": executable,
        "formula_binding": {
            **exact_formula_binding(context),
            "active_cases": context["active_cases"],
        },
        "hyperbolicity": hyper,
        "analysis": analysis,
        "ct": ct,
        "publication": publication,
    }


def complete_workflow(workflow_root: Path) -> Path:
    workflow, workflow_binding = load_workflow(workflow_root)
    context = workflow_context(workflow)
    record = completion_record(workflow, workflow_binding, context)
    root = Path(str(workflow["workflow_root"]))
    inventory_sha = str(context["inventory_binding"]["sha256"])
    record_path = root / "completion/records" / f"{inventory_sha}.json"
    write_immutable_json(record_path, record)
    pointer = {
        "schema": POINTER_SCHEMA,
        "schema_version": 1,
        "status": "complete",
        "campaign_kind": "corrected-production",
        "campaign_identity": context["identity_binding"],
        "inventory": context["inventory_binding"],
        "corrected_executable": record["corrected_executable"],
        "formula_binding": record["formula_binding"],
        "completion_record": artifact_binding(record_path),
    }
    pointer_path = root / "completion/corrected-production-complete.json"
    write_immutable_json(pointer_path, pointer)
    print(f"wrote corrected-only completion pointer: {pointer_path}")
    return pointer_path


def validate_completion(workflow_root: Path) -> Path:
    workflow, workflow_binding = load_workflow(workflow_root)
    context = workflow_context(workflow)
    root = Path(str(workflow["workflow_root"]))
    pointer_path = root / "completion/corrected-production-complete.json"
    pointer, _ = load_bound_json(pointer_path, "completion pointer")
    if (
        pointer.get("schema") != POINTER_SCHEMA
        or pointer.get("status") != "complete"
        or pointer.get("campaign_identity") != context["identity_binding"]
        or pointer.get("inventory") != context["inventory_binding"]
    ):
        raise CorrectedDownstreamError("completion pointer identity differs")
    record_binding = verify_binding(pointer.get("completion_record"), "completion record")
    record, _ = load_bound_json(Path(str(record_binding["path"])), "completion record")
    expected = completion_record(workflow, workflow_binding, context)
    if record != expected:
        raise CorrectedDownstreamError("completion record differs from current evidence")
    print(f"validated corrected-only completion pointer: {pointer_path}")
    return pointer_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subcommands = parser.add_subparsers(dest="command", required=True)

    prepare = subcommands.add_parser(
        "prepare", help="prepare corrected-only downstream jobs without submitting"
    )
    prepare.add_argument("--identity", type=Path, required=True)
    prepare.add_argument("--inventory", type=Path, required=True)
    prepare.add_argument("--inventory-sha256", required=True)
    prepare.add_argument("--workflow-root", type=Path, required=True)
    prepare.add_argument("--python", type=Path, default=Path(sys.executable))
    prepare.add_argument("--account", default="ast207")
    prepare.add_argument("--partition", default="batch")
    prepare.add_argument("--walltime", default="02:00:00")
    prepare.add_argument("--cpus-per-task", type=int, default=56)

    for name, help_text in (
        ("submit", "explicitly submit a previously prepared workflow"),
        (
            "retry-publication",
            "validate current upstream outputs and submit a fresh publication attempt",
        ),
        ("complete", "validate all outputs and write corrected-only pointers"),
        ("validate", "revalidate an existing corrected-only completion pointer"),
    ):
        command = subcommands.add_parser(name, help=help_text)
        command.add_argument("--workflow-root", type=Path, required=True)
    bind = subcommands.add_parser("bind-hyper-result", help=argparse.SUPPRESS)
    bind.add_argument("--identity", type=Path, required=True)
    bind.add_argument("--manifest", type=Path, required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.command == "prepare":
            if args.cpus_per_task < 1:
                raise CorrectedDownstreamError("--cpus-per-task must be positive")
            prepare_workflow(args)
        elif args.command == "submit":
            submit_workflow(args.workflow_root)
        elif args.command == "retry-publication":
            retry_publication(args.workflow_root)
        elif args.command == "complete":
            complete_workflow(args.workflow_root)
        elif args.command == "validate":
            validate_completion(args.workflow_root)
        else:
            enrich_hyper_result(args.identity, args.manifest)
    except (
        CorrectedDownstreamError,
        OSError,
        subprocess.CalledProcessError,
        ValueError,
        TypeError,
        KeyError,
    ) as error:
        print(f"corrected downstream error: {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
