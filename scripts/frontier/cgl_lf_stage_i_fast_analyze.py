#!/usr/bin/env python3
"""Submit independent direct-fast per-case analysis jobs to Slurm.

This is a lean analysis launcher, not a campaign controller.  It reads an
assembled direct-fast inventory, writes job records outside simulation roots,
and runs one ``cgl_lf_stage_i_fast_report.py analyze-case`` process per job.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
from typing import Iterable


DEFAULT_ACCOUNT = "ast207"
DEFAULT_PARTITION = "batch"
DEFAULT_WALLTIME = "02:00:00"
DEFAULT_CPUS_PER_TASK = 56
DEFAULT_SNAPSHOT_WORKERS = 4
DEFAULT_SNAPSHOT_MEMORY_BUDGET_GIB = 384.0
REPORTER = Path(__file__).resolve().with_name("cgl_lf_stage_i_fast_report.py")
CORRECTED_REPORTER = Path(__file__).resolve().with_name(
    "cgl_lf_stage_i_fast_corrected_report.py"
)
CASE_ID = re.compile(r"R\d{2}")
JOB_ID = re.compile(r"[1-9]\d*(?:;[A-Za-z0-9_.-]+)?")
FAILED_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "TIMEOUT",
}
FAILED_ANALYSIS_OUTPUT = "FAILED_ANALYSIS_OUTPUT"


class AnalysisLaunchError(RuntimeError):
    """A direct-fast analysis launch precondition failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise AnalysisLaunchError(f"cannot read JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise AnalysisLaunchError(f"expected a JSON object: {path}")
    return value


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        delete=False,
    ) as stream:
        stream.write(payload)
        temporary = Path(stream.name)
    temporary.replace(path)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact_binding(path: Path) -> dict[str, object]:
    resolved = path.resolve()
    stat = resolved.stat()
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "size_bytes": stat.st_size,
    }


def binding_match_errors(
    observed: object, expected: object, label: str
) -> list[str]:
    if not isinstance(observed, dict) or not isinstance(expected, dict):
        return [f"{label} binding is malformed"]
    errors: list[str] = []
    for key in ("path", "sha256", "size_bytes"):
        if key not in observed or key not in expected:
            errors.append(f"{label} binding lacks {key}")
            continue
        left = observed[key]
        right = expected[key]
        if key == "path" and isinstance(left, str) and isinstance(right, str):
            left = str(Path(left).expanduser().resolve())
            right = str(Path(right).expanduser().resolve())
        if left != right:
            errors.append(f"{label} binding differs in {key}")
    return errors


def current_binding_errors(
    binding: object, label: str, expected_path: Path | None = None
) -> list[str]:
    if not isinstance(binding, dict) or not isinstance(binding.get("path"), str):
        return [f"{label} binding is malformed"]
    path = Path(str(binding["path"])).expanduser().resolve()
    errors = []
    if expected_path is not None and path != expected_path.expanduser().resolve():
        errors.append(f"{label} path differs from the expected artifact")
    try:
        current = artifact_binding(path)
    except OSError as error:
        errors.append(f"cannot bind current {label}: {error}")
        return errors
    errors.extend(binding_match_errors(binding, current, label))
    return errors


def load_corrected_reporter() -> object:
    """Load the exact corrected reporter that owns terminal authorization."""

    name = "_cgl_lf_stage_i_fast_analyze_corrected_reporter"
    existing = sys.modules.get(name)
    if existing is not None:
        return existing
    spec = importlib.util.spec_from_file_location(name, CORRECTED_REPORTER)
    if spec is None or spec.loader is None:
        raise AnalysisLaunchError(f"cannot import corrected reporter: {CORRECTED_REPORTER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def authenticated_terminal_cases(
    inventory_path: Path, inventory: dict[str, object]
) -> set[str]:
    """Revalidate a corrected composite before admitting terminal partials."""

    cases = inventory.get("cases")
    if not isinstance(cases, dict) or not any(
        isinstance(record, dict) and record.get("status") == "failed_partial"
        for record in cases.values()
    ):
        return set()
    if inventory.get("record_type") != "cgl_lf_stage_i_corrected_composite_report":
        return set()
    try:
        corrected = load_corrected_reporter()
        expected_adapter = artifact_binding(CORRECTED_REPORTER)
        adapter_errors = binding_match_errors(
            inventory.get("composite_adapter"),
            expected_adapter,
            "corrected composite adapter",
        )
        if adapter_errors:
            raise AnalysisLaunchError("; ".join(adapter_errors))
        validation = corrected.validate_composite_inventory(
            corrected.DEFAULT_CONFIG, inventory_path.parent
        )
        validation_errors = binding_match_errors(
            validation.get("inventory"),
            artifact_binding(inventory_path),
            "corrected composite validation",
        )
        if validation_errors:
            raise AnalysisLaunchError("; ".join(validation_errors))
        dispositions = validation.get("terminal_dispositions")
        if not isinstance(dispositions, dict):
            raise AnalysisLaunchError("corrected terminal dispositions are malformed")
        return set(dispositions)
    except AnalysisLaunchError:
        raise
    except Exception as error:
        raise AnalysisLaunchError(
            f"corrected composite terminal validation failed: {error}"
        ) from error


def analysis_completion_errors(
    attempt_dir: Path, manifest: dict[str, object]
) -> list[str]:
    """Return reasons a zero-exit analysis attempt is not complete."""

    errors: list[str] = []
    case_id = manifest.get("case_id")
    analysis_value = manifest.get("analysis_output")
    attempt_value = manifest.get("attempt_dir")
    if not isinstance(case_id, str) or not CASE_ID.fullmatch(case_id):
        return ["job manifest case ID is malformed"]
    if not isinstance(analysis_value, str):
        return ["job manifest analysis output is malformed"]
    if (
        not isinstance(attempt_value, str)
        or Path(attempt_value).expanduser().resolve() != attempt_dir.resolve()
    ):
        errors.append("job manifest attempt directory differs from the attempt")
    analysis = Path(analysis_value).expanduser().resolve()
    inventory_path = analysis / "inventory.json"
    errors.extend(
        current_binding_errors(
            manifest.get("inventory"), "job inventory", inventory_path
        )
    )
    errors.extend(
        current_binding_errors(manifest.get("reporter"), "job reporter", REPORTER)
    )

    inventory: dict[str, object] | None = None
    try:
        inventory = load_json(inventory_path)
    except AnalysisLaunchError as error:
        errors.append(str(error))

    diagnostics_path = analysis / "cases" / case_id / "diagnostics.json"
    try:
        diagnostics = load_json(diagnostics_path)
    except AnalysisLaunchError as error:
        errors.append(str(error))
        return errors
    if diagnostics.get("case_id") != case_id:
        errors.append("case diagnostics identity differs from the job")
    prepared_value = manifest.get("prepared_utc")
    analyzed_value = diagnostics.get("analyzed_utc")
    try:
        prepared = datetime.fromisoformat(str(prepared_value))
        analyzed = datetime.fromisoformat(str(analyzed_value))
        if prepared.tzinfo is None or analyzed.tzinfo is None:
            raise ValueError("timestamps must include time zones")
        if analyzed < prepared:
            errors.append("case diagnostics predate the analysis attempt")
    except ValueError:
        errors.append("case diagnostics or job preparation timestamp is malformed")
    inventory_case: dict[str, object] | None = None
    if inventory is not None:
        cases = inventory.get("cases")
        value = cases.get(case_id) if isinstance(cases, dict) else None
        inventory_case = value if isinstance(value, dict) else None
    terminal_cases: set[str] = set()
    if inventory is not None:
        try:
            terminal_cases = authenticated_terminal_cases(inventory_path, inventory)
        except AnalysisLaunchError as error:
            errors.append(str(error))
    terminal_failure = case_id in terminal_cases
    expected_assembly = "failed_partial" if terminal_failure else "complete"
    expected_analysis = "partial" if terminal_failure else "complete"
    if diagnostics.get("assembly_status") != expected_assembly:
        errors.append(
            f"case diagnostics assembly_status is not {expected_assembly}"
        )
    if diagnostics.get("analysis_status") != expected_analysis:
        errors.append(
            f"case diagnostics analysis_status is not {expected_analysis}"
        )
    if diagnostics.get("snapshot_analysis_status") != "complete":
        errors.append("case diagnostics snapshot_analysis_status is not complete")
    if diagnostics.get("analysis_errors") != []:
        errors.append(
            "case diagnostics retain or lack a valid empty analysis_errors list"
        )

    provenance = diagnostics.get("provenance")
    if not isinstance(provenance, dict):
        errors.append("case diagnostics provenance is malformed")
        return errors
    errors.extend(
        binding_match_errors(
            provenance.get("adapter"),
            manifest.get("reporter"),
            "diagnostics reporter",
        )
    )

    lineage_path = analysis / "cases" / case_id / "lineage.json"
    errors.extend(
        current_binding_errors(
            provenance.get("lineage"), "diagnostics lineage", lineage_path
        )
    )
    if inventory is not None:
        if not isinstance(inventory_case, dict):
            errors.append("job inventory lacks the expected case record")
        else:
            if inventory_case.get("status") != expected_assembly:
                errors.append(
                    f"job inventory case status is not {expected_assembly}"
                )
            if diagnostics.get("case_name") != inventory_case.get("case_name"):
                errors.append("case diagnostics name differs from the job inventory")
            try:
                lineage = load_json(lineage_path)
            except AnalysisLaunchError as error:
                errors.append(str(error))
            else:
                if lineage != inventory_case:
                    errors.append("diagnostics lineage differs from the job inventory")
    return errors


def is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def inventory_context(
    analysis_argument: Path, jobs_argument: Path | None
) -> tuple[Path, Path, Path, dict[str, object], set[str]]:
    supplied = analysis_argument.expanduser().resolve()
    inventory_path = supplied if supplied.is_file() else supplied / "inventory.json"
    analysis = inventory_path.parent
    inventory = load_json(inventory_path)
    cases = inventory.get("cases")
    if not isinstance(cases, dict):
        raise AnalysisLaunchError(f"inventory has no case mapping: {inventory_path}")
    declared_output = inventory.get("output")
    if not isinstance(declared_output, str):
        raise AnalysisLaunchError(f"inventory has no declared output: {inventory_path}")
    if Path(declared_output).resolve() != analysis:
        raise AnalysisLaunchError(
            f"inventory output does not match its directory: {inventory_path}"
        )
    root_value = inventory.get("root")
    if not isinstance(root_value, str):
        raise AnalysisLaunchError(f"inventory has no campaign root: {inventory_path}")
    simulation_root = Path(root_value).resolve() / "runs"
    jobs = (
        jobs_argument.expanduser().resolve()
        if jobs_argument is not None
        else analysis / "slurm-analysis"
    )
    if is_relative_to(analysis, simulation_root) or is_relative_to(
        jobs, simulation_root
    ):
        raise AnalysisLaunchError("analysis jobs and output must be outside runs/")
    return (
        analysis,
        inventory_path,
        jobs,
        inventory,
        authenticated_terminal_cases(inventory_path, inventory),
    )


def expand_case_tokens(tokens: Iterable[str], available: set[str]) -> list[str]:
    selected: set[str] = set()
    for raw_token in tokens:
        for token in raw_token.split(","):
            token = token.strip().upper()
            if not token:
                continue
            match = re.fullmatch(r"(R\d{2})-(R\d{2})", token)
            if match:
                start = int(match.group(1)[1:])
                end = int(match.group(2)[1:])
                if start > end:
                    raise AnalysisLaunchError(f"descending case range: {token}")
                selected.update(f"R{number:02d}" for number in range(start, end + 1))
            elif CASE_ID.fullmatch(token):
                selected.add(token)
            else:
                raise AnalysisLaunchError(f"invalid case selector: {token}")
    if not selected:
        selected = set(available)
    unknown = selected - available
    if unknown:
        raise AnalysisLaunchError(
            f"cases absent from inventory: {', '.join(sorted(unknown))}"
        )
    return sorted(selected)


def eligible_cases(
    inventory: dict[str, object],
    selectors: Iterable[str],
    include_all: bool,
    terminal_cases: set[str] | None = None,
) -> list[str]:
    records = inventory["cases"]
    assert isinstance(records, dict)
    selected = expand_case_tokens(selectors, set(records))
    if include_all:
        return selected
    eligible = []
    for case_id in selected:
        record = records[case_id]
        if (
            isinstance(record, dict)
            and record.get("status") == "complete"
        ) or case_id in (terminal_cases or set()):
            eligible.append(case_id)
        else:
            print(f"skip {case_id}: assembled lineage is not complete")
    return eligible


def analysis_arguments(args: argparse.Namespace) -> list[str]:
    if args.snapshot_workers < 1:
        raise AnalysisLaunchError("--snapshot-workers must be positive")
    if (
        not math.isfinite(args.snapshot_memory_budget_gib)
        or args.snapshot_memory_budget_gib <= 0.0
    ):
        raise AnalysisLaunchError(
            "--snapshot-memory-budget-gib must be positive and finite"
        )
    options: list[str] = []
    if args.skip_snapshots:
        options.append("--skip-snapshots")
    options.extend(["--snapshot-time-start", str(args.snapshot_time_start)])
    options.extend(["--snapshot-time-end", str(args.snapshot_time_end)])
    options.extend(["--snapshot-workers", str(args.snapshot_workers)])
    options.extend([
        "--snapshot-memory-budget-gib",
        str(args.snapshot_memory_budget_gib),
    ])
    options.extend(["--pdf-bins", str(args.pdf_bins)])
    options.extend(["--alignment-shells", args.alignment_shells])
    if args.eddy_samples is not None:
        options.extend(["--eddy-samples", str(args.eddy_samples)])
    options.extend(["--eddy-bins", str(args.eddy_bins)])
    options.extend(["--eddy-seed", str(args.eddy_seed)])
    return options


def validate_retry_analysis_overrides(args: argparse.Namespace) -> None:
    if args.snapshot_workers is not None and args.snapshot_workers < 1:
        raise AnalysisLaunchError("--snapshot-workers must be positive")
    if args.snapshot_memory_budget_gib is not None and (
        not math.isfinite(args.snapshot_memory_budget_gib)
        or args.snapshot_memory_budget_gib <= 0.0
    ):
        raise AnalysisLaunchError(
            "--snapshot-memory-budget-gib must be positive and finite"
        )


def retry_report_options(
    inherited: list[str], args: argparse.Namespace
) -> list[str]:
    validate_retry_analysis_overrides(args)
    updated = list(inherited)
    overrides = (
        ("--snapshot-workers", args.snapshot_workers),
        ("--snapshot-memory-budget-gib", args.snapshot_memory_budget_gib),
    )
    for option, value in overrides:
        if value is None:
            continue
        positions = [
            index for index, token in enumerate(updated) if token == option
        ]
        if len(positions) != 1 or positions[0] + 1 >= len(updated):
            raise AnalysisLaunchError(
                f"inherited report options must contain exactly one {option} value"
            )
        updated[positions[0] + 1] = str(value)
    return updated


def initial_retry_report_options(args: argparse.Namespace) -> list[str]:
    validate_retry_analysis_overrides(args)
    options: list[str] = []
    if args.snapshot_workers is not None:
        options.extend(["--snapshot-workers", str(args.snapshot_workers)])
    if args.snapshot_memory_budget_gib is not None:
        options.extend([
            "--snapshot-memory-budget-gib",
            str(args.snapshot_memory_budget_gib),
        ])
    return options


def attempt_directories(jobs: Path, case_id: str) -> list[Path]:
    parent = jobs / case_id
    if not parent.is_dir():
        return []
    return sorted(
        path
        for path in parent.glob("attempt-[0-9][0-9][0-9]")
        if path.is_dir() and (path / "manifest.json").is_file()
    )


def batch_script(manifest: dict[str, object]) -> str:
    command = manifest["command"]
    if not isinstance(command, list) or not all(
        isinstance(item, str) for item in command
    ):
        raise AnalysisLaunchError("internal error: malformed analysis command")
    attempt_dir = Path(str(manifest["attempt_dir"]))
    lines = [
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
        f"cd {shlex.quote(str(manifest['repo']))}",
        f"export OMP_NUM_THREADS={int(manifest['cpus_per_task'])}",
        "set +e",
        (
            "srun --nodes=1 --ntasks=1 "
            f"--cpus-per-task={int(manifest['cpus_per_task'])} "
            + shlex.join(command)
        ),
        "status=$?",
        (
            f"printf '%s\\n' \"$status\" > "
            f"{shlex.quote(str(attempt_dir / 'exit_code.txt'))}"
        ),
        "exit \"$status\"",
        "",
    ]
    return "\n".join(lines)


def prepare_attempt(
    *,
    analysis: Path,
    inventory_path: Path,
    jobs: Path,
    case_id: str,
    report_options: list[str],
    account: str,
    partition: str,
    walltime: str,
    cpus_per_task: int,
    python: Path,
    retry_of: str | None = None,
) -> Path:
    if cpus_per_task < 1:
        raise AnalysisLaunchError("--cpus-per-task must be positive")
    existing = attempt_directories(jobs, case_id)
    number = len(existing)
    attempt_dir = jobs / case_id / f"attempt-{number:03d}"
    if attempt_dir.exists():
        raise AnalysisLaunchError(f"attempt path already exists: {attempt_dir}")
    attempt_dir.mkdir(parents=True)
    log_dir = jobs / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[2]
    command = [
        str(python),
        str(REPORTER),
        "--output",
        str(analysis),
        "analyze-case",
        case_id,
        *report_options,
    ]
    manifest: dict[str, object] = {
        "schema_version": 1,
        "record_type": "cgl_lf_stage_i_direct_fast_analysis_job",
        "prepared_utc": utc_now(),
        "case_id": case_id,
        "attempt": number,
        "attempt_dir": str(attempt_dir),
        "analysis_output": str(analysis),
        "inventory": artifact_binding(inventory_path),
        "repo": str(repo),
        "python": str(python),
        "reporter": artifact_binding(REPORTER),
        "launcher": artifact_binding(Path(__file__)),
        "report_options": report_options,
        "command": command,
        "account": account,
        "partition": partition,
        "walltime": walltime,
        "nodes": 1,
        "ntasks": 1,
        "cpus_per_task": cpus_per_task,
        "job_name": f"cglfa_{case_id}_a{number:03d}",
        "slurm_log": str(log_dir / f"{case_id}.attempt-{number:03d}.%j.log"),
        "retry_of": retry_of,
        "job_id": None,
    }
    write_json(attempt_dir / "manifest.json", manifest)
    script = attempt_dir / "run.sbatch"
    script.write_text(batch_script(manifest), encoding="utf-8")
    script.chmod(0o755)
    return attempt_dir


def submit_attempt(attempt_dir: Path) -> str:
    manifest_path = attempt_dir / "manifest.json"
    manifest = load_json(manifest_path)
    if manifest.get("job_id") is not None:
        raise AnalysisLaunchError(f"attempt already submitted: {attempt_dir}")
    completed = subprocess.run(
        ["/usr/bin/sbatch", "--parsable", str(attempt_dir / "run.sbatch")],
        check=True,
        text=True,
        capture_output=True,
    )
    response = completed.stdout.strip()
    if not JOB_ID.fullmatch(response):
        raise AnalysisLaunchError(f"unexpected sbatch response: {response!r}")
    job_id = response.split(";", 1)[0]
    manifest["job_id"] = job_id
    manifest["submitted_utc"] = utc_now()
    write_json(manifest_path, manifest)
    return job_id


def update_prepared_attempt(
    attempt_dir: Path,
    manifest: dict[str, object],
    account: str | None,
    partition: str | None,
    walltime: str | None,
    cpus_per_task: int | None,
    report_options: list[str] | None,
) -> None:
    if (
        account is None
        and partition is None
        and walltime is None
        and cpus_per_task is None
        and report_options is None
    ):
        return
    if account is not None:
        manifest["account"] = account
    if partition is not None:
        manifest["partition"] = partition
    if walltime is not None:
        manifest["walltime"] = walltime
    if cpus_per_task is not None:
        if cpus_per_task < 1:
            raise AnalysisLaunchError("--cpus-per-task must be positive")
        manifest["cpus_per_task"] = cpus_per_task
    if report_options is not None:
        reporter = manifest.get("reporter")
        if not isinstance(reporter, dict) or not isinstance(
            reporter.get("path"), str
        ):
            raise AnalysisLaunchError("prepared manifest reporter is malformed")
        required_strings = {
            key: manifest.get(key)
            for key in ("python", "analysis_output", "case_id")
        }
        if not all(isinstance(value, str) for value in required_strings.values()):
            raise AnalysisLaunchError("prepared manifest command fields are malformed")
        manifest["report_options"] = report_options
        manifest["command"] = [
            str(required_strings["python"]),
            str(reporter["path"]),
            "--output",
            str(required_strings["analysis_output"]),
            "analyze-case",
            str(required_strings["case_id"]),
            *report_options,
        ]
        manifest["launcher"] = artifact_binding(Path(__file__))
    write_json(attempt_dir / "manifest.json", manifest)
    script = attempt_dir / "run.sbatch"
    script.write_text(batch_script(manifest), encoding="utf-8")
    script.chmod(0o755)


def normalized_state(value: str) -> str:
    return value.strip().split("+", 1)[0].split()[0].upper() if value.strip() else ""


def attempt_state(attempt_dir: Path, manifest: dict[str, object]) -> str:
    exit_path = attempt_dir / "exit_code.txt"
    if exit_path.is_file():
        try:
            exit_code = int(exit_path.read_text(encoding="utf-8").strip())
        except (OSError, ValueError):
            return "INVALID_EXIT_RECORD"
        if exit_code != 0:
            return f"FAILED_EXIT_{exit_code}"
        return (
            FAILED_ANALYSIS_OUTPUT
            if analysis_completion_errors(attempt_dir, manifest)
            else "COMPLETED"
        )
    job_id = manifest.get("job_id")
    if not isinstance(job_id, str):
        return "PREPARED"
    queued = subprocess.run(
        ["/usr/bin/squeue", "-h", "-j", job_id, "-o", "%T"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if queued:
        return normalized_state(queued.splitlines()[0])
    accounted = subprocess.run(
        ["/usr/bin/sacct", "-X", "-n", "-P", "-j", job_id, "-o", "State"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if accounted:
        return normalized_state(accounted.splitlines()[0].split("|", 1)[0])
    return "UNKNOWN"


def print_attempt(case_id: str, attempt_dir: Path) -> str:
    manifest = load_json(attempt_dir / "manifest.json")
    state = attempt_state(attempt_dir, manifest)
    job_id = str(manifest.get("job_id") or "-")
    print(
        f"{case_id}\t{attempt_dir.name}\t{job_id}\t{state}\t"
        f"{manifest.get('slurm_log')}"
    )
    return state


def launch(args: argparse.Namespace) -> int:
    analysis, inventory_path, jobs, inventory, terminal_cases = inventory_context(
        args.analysis, args.jobs_dir
    )
    cases = eligible_cases(inventory, args.cases, args.all, terminal_cases)
    report_options = analysis_arguments(args)
    python = args.python.expanduser().absolute()
    for case_id in cases:
        existing = attempt_directories(jobs, case_id)
        if existing:
            print(f"skip {case_id}: attempt exists; use retry")
            continue
        attempt = prepare_attempt(
            analysis=analysis,
            inventory_path=inventory_path,
            jobs=jobs,
            case_id=case_id,
            report_options=report_options,
            account=args.account,
            partition=args.partition,
            walltime=args.walltime,
            cpus_per_task=args.cpus_per_task,
            python=python,
        )
        if args.submit:
            print(f"submitted {case_id}: {submit_attempt(attempt)}")
        else:
            print(f"prepared {case_id}: {attempt}")
    return 0


def status(args: argparse.Namespace) -> int:
    _, _, jobs, inventory, _ = inventory_context(args.analysis, args.jobs_dir)
    records = inventory["cases"]
    assert isinstance(records, dict)
    cases = expand_case_tokens(args.cases, set(records))
    print("CASE\tATTEMPT\tJOB_ID\tSTATE\tLOG")
    for case_id in cases:
        attempts = attempt_directories(jobs, case_id)
        selected = attempts if args.all_attempts else attempts[-1:]
        if not selected:
            print(f"{case_id}\t-\t-\tNOT_PREPARED\t-")
        for attempt in selected:
            print_attempt(case_id, attempt)
    return 0


def retry(args: argparse.Namespace) -> int:
    analysis, inventory_path, jobs, inventory, terminal_cases = inventory_context(
        args.analysis, args.jobs_dir
    )
    cases = eligible_cases(inventory, args.cases, args.all, terminal_cases)
    python = args.python.expanduser().absolute()
    validate_retry_analysis_overrides(args)
    for case_id in cases:
        attempts = attempt_directories(jobs, case_id)
        if not attempts:
            attempt = prepare_attempt(
                analysis=analysis,
                inventory_path=inventory_path,
                jobs=jobs,
                case_id=case_id,
                report_options=initial_retry_report_options(args),
                account=(
                    DEFAULT_ACCOUNT if args.account is None else args.account
                ),
                partition=(
                    DEFAULT_PARTITION if args.partition is None else args.partition
                ),
                walltime=(
                    DEFAULT_WALLTIME if args.walltime is None else args.walltime
                ),
                cpus_per_task=(
                    DEFAULT_CPUS_PER_TASK
                    if args.cpus_per_task is None
                    else args.cpus_per_task
                ),
                python=python,
            )
        else:
            latest = attempts[-1]
            manifest = load_json(latest / "manifest.json")
            state = attempt_state(latest, manifest)
            if state == "PREPARED":
                options = manifest.get("report_options")
                if not isinstance(options, list) or not all(
                    isinstance(item, str) for item in options
                ):
                    raise AnalysisLaunchError(
                        f"malformed report options in {latest / 'manifest.json'}"
                    )
                overridden_options = retry_report_options(options, args)
                attempt = latest
                update_prepared_attempt(
                    attempt,
                    manifest,
                    args.account,
                    args.partition,
                    args.walltime,
                    args.cpus_per_task,
                    (
                        overridden_options
                        if overridden_options != options
                        else None
                    ),
                )
            elif (
                state.startswith("FAILED_EXIT_")
                or state == FAILED_ANALYSIS_OUTPUT
                or state in FAILED_STATES
            ):
                options = manifest.get("report_options")
                if not isinstance(options, list) or not all(
                    isinstance(item, str) for item in options
                ):
                    raise AnalysisLaunchError(
                        f"malformed report options in {latest / 'manifest.json'}"
                    )
                options = retry_report_options(options, args)
                attempt = prepare_attempt(
                    analysis=analysis,
                    inventory_path=inventory_path,
                    jobs=jobs,
                    case_id=case_id,
                    report_options=options,
                    account=(
                        args.account
                        if args.account is not None
                        else str(manifest.get("account", DEFAULT_ACCOUNT))
                    ),
                    partition=(
                        args.partition
                        if args.partition is not None
                        else str(manifest.get("partition", DEFAULT_PARTITION))
                    ),
                    walltime=(
                        args.walltime
                        if args.walltime is not None
                        else str(manifest.get("walltime", DEFAULT_WALLTIME))
                    ),
                    cpus_per_task=(
                        args.cpus_per_task
                        if args.cpus_per_task is not None
                        else int(
                            manifest.get(
                                "cpus_per_task", DEFAULT_CPUS_PER_TASK
                            )
                        )
                    ),
                    python=Path(str(manifest.get("python", python))),
                    retry_of=str(latest),
                )
            else:
                print(f"skip {case_id}: latest attempt is {state}")
                continue
        if args.submit:
            print(f"submitted {case_id}: {submit_attempt(attempt)}")
        else:
            print(f"retry-ready {case_id}: {attempt}")
    return 0


def add_selection_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("cases", nargs="*", help="case IDs or ranges; default: all")
    parser.add_argument(
        "--all",
        action="store_true",
        help="include partial/in-progress lineages instead of complete cases only",
    )


def add_slurm_options(
    parser: argparse.ArgumentParser, *, retry_overrides: bool = False
) -> None:
    parser.add_argument(
        "--account", default=None if retry_overrides else DEFAULT_ACCOUNT
    )
    parser.add_argument(
        "--partition", default=None if retry_overrides else DEFAULT_PARTITION
    )
    parser.add_argument(
        "--walltime", default=None if retry_overrides else DEFAULT_WALLTIME
    )
    parser.add_argument(
        "--cpus-per-task",
        type=int,
        default=None if retry_overrides else DEFAULT_CPUS_PER_TASK,
    )
    parser.add_argument(
        "--python",
        type=Path,
        default=Path(sys.executable),
        help="Python used inside analysis jobs; default: current Python",
    )
    parser.add_argument(
        "--submit",
        action="store_true",
        help="submit with sbatch; otherwise only prepare job records",
    )


def add_analysis_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--skip-snapshots", action="store_true")
    parser.add_argument("--snapshot-time-start", type=float, default=8.0)
    parser.add_argument("--snapshot-time-end", type=float, default=10.0)
    parser.add_argument(
        "--snapshot-workers",
        type=int,
        default=DEFAULT_SNAPSHOT_WORKERS,
        help="maximum independent snapshot-analysis processes inside each case job",
    )
    parser.add_argument(
        "--snapshot-memory-budget-gib",
        type=float,
        default=DEFAULT_SNAPSHOT_MEMORY_BUDGET_GIB,
        help="aggregate coordinator-plus-worker memory budget per case job",
    )
    parser.add_argument("--pdf-bins", type=int, default=64)
    parser.add_argument(
        "--alignment-shells", default="2,4,6,8,12,16,24,32,64,128"
    )
    parser.add_argument("--eddy-samples", type=int)
    parser.add_argument("--eddy-bins", type=int, default=24)
    parser.add_argument("--eddy-seed", type=int, default=731)


def add_retry_analysis_overrides(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--snapshot-workers",
        type=int,
        default=None,
        help="override inherited maximum snapshot-analysis processes",
    )
    parser.add_argument(
        "--snapshot-memory-budget-gib",
        type=float,
        default=None,
        help="override inherited aggregate snapshot-analysis memory budget",
    )


def parser() -> argparse.ArgumentParser:
    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "analysis",
        type=Path,
        help="direct-fast analysis output directory or its inventory.json",
    )
    command.add_argument(
        "--jobs-dir",
        type=Path,
        help="job-record directory; default: ANALYSIS/slurm-analysis",
    )
    subcommands = command.add_subparsers(dest="command", required=True)

    launch_command = subcommands.add_parser(
        "launch", help="prepare or submit one job per eligible case"
    )
    add_selection_options(launch_command)
    add_slurm_options(launch_command)
    add_analysis_options(launch_command)

    status_command = subcommands.add_parser(
        "status", help="show latest prepared/submitted job states"
    )
    status_command.add_argument("cases", nargs="*", help="case IDs or ranges")
    status_command.add_argument(
        "--all-attempts", action="store_true", help="show every attempt"
    )

    retry_command = subcommands.add_parser(
        "retry", help="prepare or submit missing and terminally failed jobs"
    )
    add_selection_options(retry_command)
    add_slurm_options(retry_command, retry_overrides=True)
    add_retry_analysis_overrides(retry_command)
    return command


def main() -> int:
    args = parser().parse_args()
    try:
        if args.command == "launch":
            return launch(args)
        if args.command == "status":
            return status(args)
        return retry(args)
    except (
        AnalysisLaunchError,
        OSError,
        subprocess.CalledProcessError,
        ValueError,
    ) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
