#!/usr/bin/env python3
"""Run and summarize reproducible CGL Landau-fluid validation workflows."""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


ROOT_DIR = Path(__file__).resolve().parents[1]
SCHEMA_VERSION = 1
HISTORY_LABEL = re.compile(r"\[\d+\]=(\S+)")
FULL_INPUTS = (
    "inputs/unit_tests/cgl_lf_quant_parallel.athinput",
    "inputs/unit_tests/cgl_lf_quant_parallel_collisional.athinput",
    "inputs/unit_tests/cgl_lf_quant_perp.athinput",
    "inputs/unit_tests/cgl_lf_quant_perp_collisional.athinput",
    "inputs/unit_tests/cgl_lf_quant_grad_b.athinput",
    "inputs/unit_tests/cgl_lf_flux_limiter.athinput",
    "inputs/unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput",
    "inputs/unit_tests/cgl_lf_limiter_mirror.athinput",
    "inputs/unit_tests/cgl_lf_limiter_firehose.athinput",
    "inputs/unit_tests/cgl_lf_field_wave.athinput",
    "inputs/unit_tests/cgl_lf_paper_oblique_wave.athinput",
    "inputs/unit_tests/cgl_pure_paper_oblique_wave.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_alfven.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_slow.athinput",
    "inputs/unit_tests/cgl_pure_paper_eigen_fast.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_alfven.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_slow.athinput",
    "inputs/unit_tests/cgl_lf_paper_eigen_fast.athinput",
)
UNSAFE_LF_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)


@dataclass(frozen=True)
class CaseSpec:
    """Description of one AthenaK run within a workflow."""

    name: str
    input_path: str
    overrides: tuple[str, ...] = ()
    lf_active: bool = True
    validation_output: bool = False
    comparison_group: str | None = None
    comparison_role: str | None = None
    amr: bool = False


@dataclass
class RunPaths:
    """Paths owned by one workflow result bundle."""

    root: Path
    logs: Path = field(init=False)
    history: Path = field(init=False)
    data: Path = field(init=False)
    figures: Path = field(init=False)
    inputs: Path = field(init=False)
    cases: Path = field(init=False)

    def __post_init__(self) -> None:
        self.logs = self.root / "logs"
        self.history = self.root / "history"
        self.data = self.root / "data"
        self.figures = self.root / "figures"
        self.inputs = self.root / "inputs"
        self.cases = self.root / "cases"

    def create(self) -> None:
        for path in (
            self.root,
            self.logs,
            self.history,
            self.data,
            self.figures,
            self.inputs,
            self.cases,
        ):
            path.mkdir(parents=True, exist_ok=True)


def display_path(path: Path, base: Path) -> str:
    """Return a readable relative path when path is inside base."""

    try:
        return str(path.relative_to(base))
    except ValueError:
        return str(path)


def resolve_from_root(value: str | Path) -> Path:
    """Resolve relative command-line paths from the repository root."""

    path = Path(value)
    return path if path.is_absolute() else ROOT_DIR / path


def timestamp() -> str:
    """UTC timestamp appropriate for run-directory names and manifests."""

    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def parse_history(path: Path) -> dict[str, list[float]]:
    """Parse an AthenaK history file without importing analysis packages."""

    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#  [1]="):
                labels = HISTORY_LABEL.findall(line)
            elif line and not line.startswith("#") and line.strip():
                rows.append([float(value) for value in line.split()])
    if not labels or not rows:
        raise ValueError(f"history output is incomplete: {path}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def parse_table(path: Path) -> dict[str, list[float]]:
    """Parse an AthenaK formatted-table output."""

    labels: list[str] = []
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("# gid"):
                labels = line[1:].split()
            elif line and not line.startswith("#") and line.strip():
                rows.append([float(value) for value in line.split()])
    if not labels or not rows:
        raise ValueError(f"table output is incomplete: {path}")
    return {
        label: [row[index] for row in rows]
        for index, label in enumerate(labels)
    }


def final_value(history: dict[str, list[float]], key: str) -> float:
    """Return the last value of a named history column."""

    if key not in history:
        raise ValueError(f"required history column is missing: {key}")
    return history[key][-1]


def lf_diagnostics(path: Path) -> dict[str, object]:
    """Extract LF safety and limiter occupancy diagnostics from history."""

    history = parse_history(path)
    nstage = final_value(history, "lf_nstage")
    counters = {
        key: final_value(history, key)
        for key in (
            "lf_nstage",
            "lf_dfloor",
            "lf_pfloor",
            "lf_nonfin",
            "lf_nonpos",
            "lf_mirror",
            "lf_firehs",
            "lf_hardbd",
        )
    }
    clean = all(counters[key] == 0.0 for key in UNSAFE_LF_COLUMNS)
    denominator = max(nstage, 1.0)
    return {
        "clean": clean,
        "counters": counters,
        "mirror_occupancy": counters["lf_mirror"] / denominator,
        "firehose_occupancy": counters["lf_firehs"] / denominator,
    }


def append_history_output(source_text: str) -> str:
    """Ensure a staged input emits history diagnostics."""

    has_hst = re.search(r"(?ms)^<output\d+>.*?^\s*file_type\s*=\s*hst\b",
                        source_text)
    if has_hst:
        return source_text
    output_ids = [
        int(value)
        for value in re.findall(r"(?m)^<output(\d+)>", source_text)
    ]
    output_id = max(output_ids, default=0) + 1
    return (
        source_text.rstrip()
        + "\n\n"
        + f"<output{output_id}>\n"
        + "file_type = hst\n"
        + "data_format = %24.16e\n"
        + "dcycle = 1\n"
    )


def stage_input(spec: CaseSpec, paths: RunPaths) -> Path:
    """Create the execution input retained with a result bundle."""

    source = ROOT_DIR / spec.input_path
    staged = paths.inputs / f"{spec.name}.athinput"
    staged.write_text(
        append_history_output(source.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    return staged


def workflow_cases(workflow: str) -> list[CaseSpec]:
    """Return the concrete cases for an executable workflow."""

    if workflow == "quick":
        return [
            CaseSpec(
                "cgl_lf_quick",
                "inputs/unit_tests/cgl_lf_quant_parallel.athinput",
                validation_output=True,
            )
        ]
    if workflow == "amr":
        return [
            CaseSpec(
                "cgl_lf_amr",
                "inputs/tests/cgl_lf_amr_2d.athinput",
                amr=True,
            )
        ]
    if workflow == "compare":
        source = "inputs/tests/cgl_lf_decay.athinput"
        return [
            CaseSpec(
                "compare_collisionless_sts",
                source,
                ("time/sts_max_dt_ratio=1.0", "time/nlim=-1"),
                comparison_group="collisionless",
                comparison_role="sts",
            ),
            CaseSpec(
                "compare_collisionless_explicit",
                source,
                (
                    "mhd/cgl_heat_flux_integrator=explicit",
                    "time/sts_integrator=none",
                    "time/nlim=-1",
                ),
                comparison_group="collisionless",
                comparison_role="explicit",
            ),
            CaseSpec(
                "compare_collision_sts",
                source,
                (
                    "time/sts_max_dt_ratio=1.0",
                    "time/nlim=-1",
                    "mhd/nu_coll=1.0",
                ),
                comparison_group="finite_collision",
                comparison_role="sts",
            ),
            CaseSpec(
                "compare_collision_explicit",
                source,
                (
                    "mhd/cgl_heat_flux_integrator=explicit",
                    "time/sts_integrator=none",
                    "time/nlim=-1",
                    "mhd/nu_coll=1.0",
                ),
                comparison_group="finite_collision",
                comparison_role="explicit",
            ),
        ]
    if workflow == "full":
        return [
            CaseSpec(
                Path(source).stem,
                source,
                lf_active="cgl_lf_" in Path(source).stem,
                validation_output=True,
            )
            for source in FULL_INPUTS
        ]
    raise ValueError(f"unsupported workflow: {workflow}")


def executable_path(args: argparse.Namespace) -> tuple[Path, Path]:
    """Determine the build directory and AthenaK executable."""

    build_dir = resolve_from_root(args.build_dir)
    if args.athena_bin:
        executable = resolve_from_root(args.athena_bin)
    else:
        executable = build_dir / "src" / "athena"
    return build_dir, executable


def ensure_executable(args: argparse.Namespace, build_dir: Path,
                      executable: Path) -> None:
    """Build AthenaK in Release mode when no requested binary exists."""

    if executable.is_file() and os.access(executable, os.X_OK):
        return
    if args.no_build:
        raise RuntimeError(f"AthenaK executable is missing: {executable}")
    configure = [
        "cmake",
        "-S",
        str(ROOT_DIR),
        "-B",
        str(build_dir),
        "-DCMAKE_BUILD_TYPE=Release",
    ]
    build = ["cmake", "--build", str(build_dir), "-j", str(args.jobs)]
    subprocess.run(configure, cwd=ROOT_DIR, check=True)
    subprocess.run(build, cwd=ROOT_DIR, check=True)
    if not executable.is_file():
        raise RuntimeError(f"build did not produce AthenaK binary: {executable}")


def git_revision() -> str:
    """Read the source revision represented by a workflow bundle."""

    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT_DIR,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else "unknown"


def collect_outputs(spec: CaseSpec, paths: RunPaths,
                    case_dir: Path) -> dict[str, object]:
    """Copy analysis products into standard result-bundle directories."""

    outputs: dict[str, object] = {}
    for history in case_dir.glob("*.hst"):
        destination = paths.history / history.name
        shutil.copy2(history, destination)
        key = "user_history" if ".user." in history.name else "mhd_history"
        outputs[key] = display_path(destination, paths.root)
    tabs = sorted((case_dir / "tab").glob("*.tab"))
    if tabs:
        destination = paths.data / f"{spec.name}.mhd_w.tab"
        shutil.copy2(tabs[-1], destination)
        outputs["final_table"] = display_path(destination, paths.root)
    csv_files = sorted(paths.data.glob(f"{spec.name}.*.csv"))
    if csv_files:
        outputs["validation_data"] = [
            display_path(path, paths.root) for path in csv_files
        ]
    return outputs


def run_case(spec: CaseSpec, executable: Path,
             paths: RunPaths) -> dict[str, object]:
    """Run one case, retaining inputs, logs, and selected output products."""

    staged_input = stage_input(spec, paths)
    case_dir = paths.cases / spec.name
    case_dir.mkdir(parents=True, exist_ok=True)
    overrides = [f"job/basename={spec.name}", *spec.overrides]
    if spec.validation_output:
        overrides.extend(
            [
                "problem/validation_output=true",
                f"problem/validation_output_dir={paths.data}",
            ]
        )
    command = [str(executable), "-i", str(staged_input), *overrides]
    log_path = paths.logs / f"{spec.name}.log"
    print(f"  {spec.name}")
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            cwd=case_dir,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
        )
    outputs = collect_outputs(spec, paths, case_dir)
    case: dict[str, object] = {
        "name": spec.name,
        "input": spec.input_path,
        "execution_input": display_path(staged_input, paths.root),
        "command": command,
        "overrides": overrides,
        "log": display_path(log_path, paths.root),
        "status": "passed" if result.returncode == 0 else "failed",
        "returncode": result.returncode,
        "outputs": outputs,
        "lf_active": spec.lf_active,
        "amr": spec.amr,
    }
    if spec.comparison_group is not None:
        case["comparison_group"] = spec.comparison_group
        case["comparison_role"] = spec.comparison_role
    return case


def case_history_path(case: dict[str, object], root: Path,
                      key: str = "mhd_history") -> Path | None:
    """Resolve one recorded output path from a case manifest."""

    outputs = case.get("outputs", {})
    if not isinstance(outputs, dict) or key not in outputs:
        return None
    return root / str(outputs[key])


def evaluate_lf_case(case: dict[str, object], root: Path) -> dict[str, object]:
    """Evaluate LF history checks for a completed case."""

    history = case_history_path(case, root)
    if history is None:
        return {"clean": False, "error": "missing MHD history output"}
    try:
        return lf_diagnostics(history)
    except (OSError, ValueError) as error:
        return {"clean": False, "error": str(error)}


def evaluate_amr(case: dict[str, object], root: Path) -> dict[str, object]:
    """Evaluate supported AMR safety and conservation checks."""

    user_path = case_history_path(case, root, "user_history")
    mhd_path = case_history_path(case, root)
    if user_path is None or mhd_path is None:
        return {"passed": False, "error": "missing AMR history output"}
    try:
        user = parse_history(user_path)
        mhd = parse_history(mhd_path)
        ncell_initial = final_value({"value": [user["ncell"][0]]}, "value")
        ncell_final = final_value(user, "ncell")
        max_ndiv = max(user["max_ndiv"])
        bad_state = max(user["bad_state"])
        anisotropy_finite = all(math.isfinite(value)
                                for value in user["abs_anis"])
        energy_initial = mhd["tot-E"][0]
        energy_residual = abs(final_value(mhd, "tot-E") - energy_initial)
        energy_residual /= max(abs(energy_initial), 1.0e-30)
        passed = (
            ncell_final > ncell_initial
            and max_ndiv < 1.0e-12
            and bad_state == 0.0
            and anisotropy_finite
            and energy_residual < 5.0e-3
        )
        return {
            "passed": passed,
            "refined": ncell_final > ncell_initial,
            "ncell_initial": ncell_initial,
            "ncell_final": ncell_final,
            "max_normalized_divb": max_ndiv,
            "bad_state": bad_state,
            "anisotropy_finite": anisotropy_finite,
            "energy_residual": energy_residual,
        }
    except (OSError, ValueError, KeyError) as error:
        return {"passed": False, "error": str(error)}


def compare_tables(first: Path, second: Path) -> dict[str, object]:
    """Compare final primitive table data for two reference runs."""

    left = parse_table(first)
    right = parse_table(second)
    columns = sorted(set(left).intersection(right) - {"gid", "i", "x1v"})
    maximum = 0.0
    for column in columns:
        if len(left[column]) != len(right[column]):
            raise ValueError(f"table row counts differ for field {column}")
        maximum = max(
            maximum,
            max(
                abs(a - b)
                for a, b in zip(left[column], right[column])
            ),
        )
    return {
        "passed": maximum < 1.0e-3,
        "maximum_field_difference": maximum,
        "tolerance": 1.0e-3,
        "fields": columns,
    }


def evaluate_comparisons(cases: list[dict[str, object]],
                         root: Path) -> dict[str, object]:
    """Evaluate STS versus explicit comparisons represented in a manifest."""

    groups: dict[str, dict[str, dict[str, object]]] = {}
    for case in cases:
        group = case.get("comparison_group")
        role = case.get("comparison_role")
        if group and role:
            groups.setdefault(str(group), {})[str(role)] = case
    results: dict[str, object] = {}
    for name, pair in groups.items():
        try:
            paths = []
            for role in ("sts", "explicit"):
                outputs = pair[role]["outputs"]
                paths.append(root / str(outputs["final_table"]))
            results[name] = compare_tables(paths[0], paths[1])
        except (KeyError, OSError, ValueError) as error:
            results[name] = {"passed": False, "error": str(error)}
    return results


def evaluate_manifest(manifest: dict[str, object], root: Path) -> dict[str, object]:
    """Recompute analysis metrics from retained workflow products."""

    cases = manifest.get("cases", [])
    if not isinstance(cases, list):
        raise ValueError("manifest cases must be a list")
    lf_results: dict[str, object] = {}
    amr_results: dict[str, object] = {}
    for case in cases:
        if bool(case.get("lf_active")):
            lf_results[str(case["name"])] = evaluate_lf_case(case, root)
        if bool(case.get("amr")):
            amr_results[str(case["name"])] = evaluate_amr(case, root)
    comparisons = evaluate_comparisons(cases, root)
    processes_pass = all(case.get("status") == "passed" for case in cases)
    lf_pass = all(bool(value.get("clean")) for value in lf_results.values())
    amr_pass = all(bool(value.get("passed")) for value in amr_results.values())
    comparison_pass = all(
        bool(value.get("passed")) for value in comparisons.values()
    )
    return {
        "passed": processes_pass and lf_pass and amr_pass and comparison_pass,
        "lf": lf_results,
        "amr": amr_results,
        "comparisons": comparisons,
    }


def write_summary(manifest: dict[str, object], path: Path) -> None:
    """Write a concise human-readable result summary."""

    diagnostics = manifest["diagnostics"]
    lines = [
        f"# CGL-LF {manifest['workflow']} workflow summary",
        "",
        f"- Status: **{manifest['status']}**",
        f"- Created UTC: `{manifest['created_utc']}`",
        f"- Git revision: `{manifest['git_revision']}`",
        f"- Executable: `{manifest['executable']}`",
        "",
        "## Cases",
        "",
        "| Case | Status | LF safety |",
        "| --- | --- | --- |",
    ]
    lf_results = diagnostics.get("lf", {})
    for case in manifest["cases"]:
        lf = lf_results.get(case["name"])
        safety = "not applicable"
        if lf is not None:
            safety = "clean" if lf.get("clean") else "failed"
        lines.append(f"| `{case['name']}` | {case['status']} | {safety} |")
    comparisons = diagnostics.get("comparisons", {})
    if comparisons:
        lines.extend(["", "## Explicit Reference Comparisons", ""])
        for group, result in comparisons.items():
            lines.append(
                "- `{}:` max field difference = `{:.6e}` "
                "(limit `{:.6e}`), {}.".format(
                    group,
                    result.get("maximum_field_difference", float("nan")),
                    result.get("tolerance", float("nan")),
                    "passed" if result.get("passed") else "failed",
                )
            )
    amr_results = diagnostics.get("amr", {})
    if amr_results:
        lines.extend(["", "## AMR Diagnostics", ""])
        for name, result in amr_results.items():
            lines.append(
                "- `{}`: refined=`{}`, max normalized divB=`{:.6e}`, "
                "energy residual=`{:.6e}`, bad states=`{}`, {}.".format(
                    name,
                    result.get("refined", False),
                    result.get("max_normalized_divb", float("nan")),
                    result.get("energy_residual", float("nan")),
                    result.get("bad_state", "unknown"),
                    "passed" if result.get("passed") else "failed",
                )
            )
    if lf_results:
        lines.extend(["", "## LF Diagnostics", ""])
        for name, result in lf_results.items():
            if "counters" not in result:
                lines.append(f"- `{name}`: diagnostics unavailable.")
                continue
            counters = result["counters"]
            lines.append(
                "- `{}`: stages=`{:.0f}`, mirror occupancy=`{:.6e}`, "
                "firehose occupancy=`{:.6e}`, unsafe counters="
                "`{:.0f}/{:.0f}/{:.0f}/{:.0f}/{:.0f}`.".format(
                    name,
                    counters["lf_nstage"],
                    result["mirror_occupancy"],
                    result["firehose_occupancy"],
                    counters["lf_dfloor"],
                    counters["lf_pfloor"],
                    counters["lf_nonfin"],
                    counters["lf_nonpos"],
                    counters["lf_hardbd"],
                )
            )
    lines.extend(
        [
            "",
            "Generated output belongs to this result bundle and is not a "
            "source artifact.",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def write_manifest(manifest: dict[str, object], paths: RunPaths) -> None:
    """Write both machine-readable and concise human-readable results."""

    (paths.root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_summary(manifest, paths.root / "summary.md")


def plot_results(paths: RunPaths) -> None:
    """Invoke the plotting helper only for workflows with full data products."""

    command = [
        sys.executable,
        str(ROOT_DIR / "scripts" / "plot_cgl_lf_validation.py"),
        "--data-dir",
        str(paths.data),
        "--figure-dir",
        str(paths.figures),
    ]
    subprocess.run(command, cwd=ROOT_DIR, check=True)


def output_directory(args: argparse.Namespace) -> Path:
    """Resolve or choose a result-bundle directory."""

    if args.output_dir:
        return resolve_from_root(args.output_dir)
    if args.workflow in ("plot", "summarize"):
        raise ValueError("--output-dir is required for plot and summarize")
    build_dir = resolve_from_root(args.build_dir)
    return build_dir / "cgl_lf_runs" / f"{timestamp()}-{args.workflow}"


def execute_workflow(args: argparse.Namespace, paths: RunPaths) -> int:
    """Execute one AthenaK workflow and write a result bundle."""

    build_dir, executable = executable_path(args)
    ensure_executable(args, build_dir, executable)
    cases = workflow_cases(args.workflow)
    print(f"Running CGL-LF {args.workflow} workflow with {executable}")
    case_results = [run_case(spec, executable, paths) for spec in cases]
    manifest: dict[str, object] = {
        "schema_version": SCHEMA_VERSION,
        "workflow": args.workflow,
        "status": "pending",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "git_revision": git_revision(),
        "executable": str(executable),
        "build_dir": str(build_dir),
        "cases": case_results,
        "diagnostics": {},
    }
    manifest["diagnostics"] = evaluate_manifest(manifest, paths.root)
    if args.workflow == "full" and manifest["diagnostics"]["passed"]:
        plot_results(paths)
    manifest["status"] = (
        "passed" if manifest["diagnostics"]["passed"] else "failed"
    )
    write_manifest(manifest, paths)
    print(f"Wrote result bundle to {paths.root}")
    return 0 if manifest["status"] == "passed" else 1


def refresh_bundle(args: argparse.Namespace, paths: RunPaths) -> int:
    """Regenerate summaries or plots for a previously executed bundle."""

    manifest_path = paths.root / "manifest.json"
    with manifest_path.open(encoding="utf-8") as stream:
        manifest = json.load(stream)
    if args.workflow == "plot":
        plot_results(paths)
    manifest["diagnostics"] = evaluate_manifest(manifest, paths.root)
    manifest["status"] = (
        "passed" if manifest["diagnostics"]["passed"] else "failed"
    )
    manifest["last_action"] = args.workflow
    write_manifest(manifest, paths)
    print(f"Updated result bundle at {paths.root}")
    return 0 if manifest["status"] == "passed" else 1


def parser() -> argparse.ArgumentParser:
    """Build command-line argument parsing."""

    command = argparse.ArgumentParser(description=__doc__)
    command.add_argument(
        "workflow",
        choices=("quick", "compare", "amr", "full", "plot", "summarize"),
    )
    command.add_argument(
        "--build-dir",
        default=os.environ.get("BUILD_DIR", "build-cgl-implementation"),
    )
    command.add_argument("--athena-bin", default=os.environ.get("ATHENA_BIN"))
    command.add_argument("--output-dir", default=os.environ.get("OUTPUT_DIR"))
    command.add_argument(
        "--jobs",
        type=int,
        default=int(os.environ.get("JOBS", os.cpu_count() or 4)),
    )
    command.add_argument(
        "--no-build",
        action="store_true",
        help="Fail instead of building when the AthenaK executable is missing.",
    )
    return command


def main() -> int:
    """CLI entry point."""

    args = parser().parse_args()
    try:
        paths = RunPaths(output_directory(args))
        if args.workflow in ("plot", "summarize"):
            if not paths.root.is_dir():
                raise ValueError(f"result bundle does not exist: {paths.root}")
            return refresh_bundle(args, paths)
        paths.create()
        return execute_workflow(args, paths)
    except (OSError, RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"CGL-LF workflow failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
