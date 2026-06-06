"""Bounded outer-x1 shock-particle escape ledger and restart regression."""

from __future__ import annotations

import glob
import json
import logging
import math
import os
import re
import shlex
import shutil
import subprocess

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_parallel_shock_outer_x1_escape_restart.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_ESCAPE_FLOAT_FIELDS = [
    "ps_escape_last_audit_time",
    "ps_escaped_injected_cr_count_global",
    "ps_escaped_injected_cr_mass_global",
    "ps_escaped_injected_cr_momentum_x1_global",
    "ps_escaped_injected_cr_momentum_x2_global",
    "ps_escaped_injected_cr_momentum_x3_global",
    "ps_escaped_injected_cr_energy_global",
    "ps_escaped_initial_cr_count_global",
]
_STARTUP_FLOAT_FIELDS = [
    "ps_removed_cr_count_global",
    "ps_removed_cr_mass_global",
    "ps_removed_cr_momentum_x1_global",
    "ps_removed_cr_momentum_x2_global",
    "ps_removed_cr_momentum_x3_global",
    "ps_removed_cr_energy_global",
]
_MACRO_MASS = 1.0e-3
_CR_LIGHT_SPEED = 1.0e4
_PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE = 2
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_MPI_NPROC = int(os.environ.get("ATHENA_PIC_PARALLEL_SHOCK_ESCAPE_MPI_NPROC", "0"))
_TELEMETRY_PATTERN = re.compile(
    r"pic_parallel_shock escape_accounting_telemetry:"
    r" population_audit_calls=(?P<population>[0-9]+)"
    r" destruction_audit_calls=(?P<destruction>[0-9]+)"
    r" population_audit_policy=(?P<policy>[a-z_]+)"
    r" production_pilot_required=(?P<pilot>[01])"
)
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_ESCAPE_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _remove_outputs(basename):
    for path in glob.glob(os.path.join(_athena_exe_dir(), "rst", basename + ".*")):
        if os.path.isfile(path):
            os.remove(path)


def _latest_restart(basename):
    matches = sorted(
        glob.glob(os.path.join(_athena_exe_dir(), "rst", basename + ".*.rst"))
    )
    if not matches:
        raise RuntimeError("No restart files found for " + basename)
    return matches[-1]


def _parse_boolean(value):
    normalized = value.strip().lower()
    if normalized in {"1", "true"}:
        return True
    if normalized in {"0", "false"}:
        return False
    raise RuntimeError("Unable to parse restart boolean: " + value)


def _restart_problem_parameters(path):
    with open(path, "rb") as handle:
        text = handle.read(65536).decode("latin1", errors="ignore")
    if "<par_end>" not in text:
        raise RuntimeError("Restart parameter header is missing <par_end>: " + path)
    active_block = None
    parameters = {}
    for raw_line in text.split("<par_end>", 1)[0].splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("<") and line.endswith(">"):
            active_block = line[1:-1]
            continue
        if active_block != "problem" or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip()
        if key in parameters:
            raise RuntimeError("Duplicate restart problem parameter: " + key)
        parameters[key] = value.split("#", 1)[0].strip()
    return parameters


def _restart_metadata(path):
    parameters = _restart_problem_parameters(path)
    required = [
        "ps_escape_ledger_schema",
        "ps_escape_ledger_complete",
        "ps_escape_audit_calls",
        "ps_injected_cr_count_global",
        *_ESCAPE_FLOAT_FIELDS,
        *_STARTUP_FLOAT_FIELDS,
    ]
    missing = [name for name in required if name not in parameters]
    if missing:
        raise RuntimeError("Restart escape ledger is missing " + repr(missing))
    return {
        "ps_escape_ledger_schema": int(parameters["ps_escape_ledger_schema"]),
        "ps_escape_ledger_complete": _parse_boolean(
            parameters["ps_escape_ledger_complete"]
        ),
        "ps_escape_audit_calls": int(parameters["ps_escape_audit_calls"]),
        "ps_injected_cr_count_global": float(
            parameters["ps_injected_cr_count_global"]
        ),
        **{name: float(parameters[name]) for name in _ESCAPE_FLOAT_FIELDS},
        **{name: float(parameters[name]) for name in _STARTUP_FLOAT_FIELDS},
    }


def _fnv1a64(path):
    value = 14695981039346656037
    with open(path, "rb") as handle:
        while True:
            payload = handle.read(1024 * 1024)
            if not payload:
                break
            for byte in payload:
                value ^= byte
                value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def _write_completion_marker(path):
    with open(path + ".complete", "w", encoding="ascii") as handle:
        handle.write("ATHENAK_RESTART_COMPLETE_V1\n")
        handle.write("size=" + str(os.path.getsize(path)) + "\n")
        handle.write("fnv1a64=" + format(_fnv1a64(path), "016x") + "\n")


def _write_shared_restart_publication(path):
    _write_completion_marker(path)
    manifest_path = path + ".manifest"
    manifest = {
        "schema": "ATHENAK_RESTART_MANIFEST_V1",
        "members": [
            {
                "path": os.path.relpath(path, _athena_exe_dir()),
                "size": os.path.getsize(path),
                "fnv1a64": format(_fnv1a64(path), "016x"),
            }
        ],
    }
    with open(manifest_path, "w", encoding="ascii") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    _write_completion_marker(manifest_path)


def _duplicate_restart_parameter(source, basename, field, contradictory_value):
    _remove_outputs(basename)
    target = os.path.join(_athena_exe_dir(), "rst", basename + ".00000.rst")
    shutil.copyfile(source, target)
    with open(target, "rb") as handle:
        payload = handle.read()
    header_end = payload.find(b"<par_end>")
    if header_end < 0:
        raise RuntimeError("Restart parameter header is missing <par_end>")
    lines = payload[:header_end].splitlines(keepends=True)
    duplicate = (field + " = " + contradictory_value + "\n").encode("ascii")
    for index, line in enumerate(lines):
        if line.split(b"=", 1)[0].strip() == field.encode("ascii"):
            lines.insert(index, duplicate)
            break
    else:
        raise RuntimeError("Unable to duplicate restart parameter " + field)
    with open(target, "wb") as handle:
        handle.write(b"".join(lines) + payload[header_end:])
    _write_shared_restart_publication(target)
    return target


def _escape_telemetry(output):
    matches = list(_TELEMETRY_PATTERN.finditer(output))
    if not matches:
        raise RuntimeError("Escape-accounting telemetry is missing")
    values = matches[-1].groupdict()
    return {
        "population_audit_calls": int(values["population"]),
        "destruction_audit_calls": int(values["destruction"]),
        "population_audit_policy": values["policy"],
        "production_pilot_required": values["pilot"] == "1",
    }


def _execute(label, basename, nlim, restart_path=None, extra=(), nproc=1):
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", os.path.relpath(restart_path, _athena_exe_dir())]
    command += [
        "job/basename=" + basename,
        "time/nlim=" + str(nlim),
        "output1/dcycle=" + str(nlim),
        *extra,
    ]
    if nproc > 1:
        command = [*_MPIEXEC, "-n", str(nproc), *command]
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def _run_success(label, basename, nlim, restart_path=None, extra=(), nproc=1):
    code, output = _execute(
        label,
        basename,
        nlim,
        restart_path=restart_path,
        extra=extra,
        nproc=nproc,
    )
    if code != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _run_expect_reject(label, basename, nlim, restart_path, extra, expected):
    code, output = _execute(
        label, basename, nlim, restart_path=restart_path, extra=extra
    )
    if code == 0 or expected not in output:
        raise RuntimeError(
            "Escape-ledger corruption was not rejected: " + label + "\n" + output
        )


def _run_expect_recenter_reject():
    code, output = _execute(
        "paper recenter rejection",
        "pic_parallel_shock_escape_recenter_reject",
        1,
        extra=["problem/ps_enable_frame_tracking=true"],
    )
    expected = (
        "pic_parallel_shock paper-mode recentering with particle shifts "
        "is not qualified"
    )
    if code == 0 or expected not in output:
        raise RuntimeError("Paper-mode particle recentering was not rejected")


def _run_expect_initial_escape_reject():
    code, output = _execute(
        "initial-particle physical escape rejection",
        "pic_parallel_shock_initial_escape_reject",
        40,
        extra=[
            "particles/ppc=1",
            "particles/cr_vx0=100",
            "problem/ps_enable_injection=false",
        ],
    )
    expected = "pic_parallel_shock found an unaccounted or invalid particle destruction event"
    if code == 0 or expected not in output:
        raise RuntimeError("Initial-particle physical escape was not rejected\n" + output)


def _run_reflective_wall_success():
    basename = "pic_parallel_shock_initial_reflection"
    _remove_outputs(basename)
    _run_success(
        "initial-particle reflecting-wall retention",
        basename,
        1,
        extra=[
            "particles/ppc=1",
            "particles/cr_vx0=-100",
            "problem/ps_enable_injection=false",
        ],
    )
    metadata = _restart_metadata(_latest_restart(basename))
    return (
        metadata["ps_escape_audit_calls"] == _PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE
        and metadata["ps_escaped_initial_cr_count_global"] == 0.0
        and metadata["ps_escaped_injected_cr_count_global"] == 0.0
    )


def _kinetic_energy_lower_bound(mass, momentum_x1, momentum_x2, momentum_x3):
    if mass <= 0.0:
        return 0.0
    state_squared = (
        momentum_x1 * momentum_x1
        + momentum_x2 * momentum_x2
        + momentum_x3 * momentum_x3
    ) / (mass * mass)
    ratio = state_squared / (_CR_LIGHT_SPEED * _CR_LIGHT_SPEED)
    specific_energy = state_squared / (math.sqrt(1.0 + ratio) + 1.0)
    return mass * specific_energy


def _population_formula_is_valid(metadata, prefix):
    mass = metadata[prefix + "_mass_global"]
    energy = metadata[prefix + "_energy_global"]
    lower_bound = _kinetic_energy_lower_bound(
        mass,
        metadata[prefix + "_momentum_x1_global"],
        metadata[prefix + "_momentum_x2_global"],
        metadata[prefix + "_momentum_x3_global"],
    )
    return energy + 2.0e-12 * max(1.0, abs(energy)) >= lower_bound


def _validate_metadata(metadata):
    finite = all(
        math.isfinite(metadata[name])
        for name in ["ps_injected_cr_count_global", *_ESCAPE_FLOAT_FIELDS]
    )
    escaped = metadata["ps_escaped_injected_cr_count_global"]
    injected = metadata["ps_injected_cr_count_global"]
    removed = metadata["ps_removed_cr_count_global"]
    return {
        "schema_is_one": metadata["ps_escape_ledger_schema"] == 1,
        "ledger_complete": metadata["ps_escape_ledger_complete"],
        "finite": finite,
        "audit_calls_positive": metadata["ps_escape_audit_calls"] > 0,
        "injected_positive": injected > 0.0,
        "escaped_positive": escaped > 0.0,
        "escaped_not_above_injected": escaped <= injected,
        "removed_positive": removed > 0.0,
        "removed_plus_escaped_not_above_injected": removed + escaped <= injected,
        "mass_matches_count": math.isclose(
            metadata["ps_escaped_injected_cr_mass_global"],
            escaped * _MACRO_MASS,
            rel_tol=2.0e-14,
            abs_tol=2.0e-14,
        ),
        "energy_positive": metadata["ps_escaped_injected_cr_energy_global"] > 0.0,
        "escape_energy_respects_momentum_formula": _population_formula_is_valid(
            metadata, "ps_escaped_injected_cr"
        ),
        "removed_mass_matches_count": math.isclose(
            metadata["ps_removed_cr_mass_global"],
            removed * _MACRO_MASS,
            rel_tol=2.0e-14,
            abs_tol=2.0e-14,
        ),
        "removed_energy_respects_momentum_formula": _population_formula_is_valid(
            metadata, "ps_removed_cr"
        ),
        "initial_escape_is_zero": (
            metadata["ps_escaped_initial_cr_count_global"] == 0.0
        ),
    }


def _summary():
    segment = _RESULTS["segment"]
    continued = _RESULTS["continued"]
    return {
        "evidence_class": "bounded_local_host_regression",
        "not_qualification_evidence": True,
        "segment": _validate_metadata(segment),
        "continued": _validate_metadata(continued),
        "continuation_audit_calls_advance": (
            continued["ps_escape_audit_calls"] > segment["ps_escape_audit_calls"]
        ),
        "segment_has_exact_vl2_audits": (
            segment["ps_escape_audit_calls"]
            == _PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * 100
        ),
        "continuation_has_exact_vl2_audits": (
            continued["ps_escape_audit_calls"]
            == _PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * 120
        ),
        "continuation_time_advances": (
            continued["ps_escape_last_audit_time"]
            > segment["ps_escape_last_audit_time"]
        ),
        "continuation_injected_count_nondecreasing": (
            continued["ps_injected_cr_count_global"]
            >= segment["ps_injected_cr_count_global"]
        ),
        "continuation_escape_count_nondecreasing": (
            continued["ps_escaped_injected_cr_count_global"]
            >= segment["ps_escaped_injected_cr_count_global"]
        ),
        "corrupt_ledgers_rejected": _RESULTS["corrupt_ledgers_rejected"],
        "duplicate_ledgers_rejected_by_runtime_and_parser": _RESULTS[
            "duplicate_ledgers_rejected_by_runtime_and_parser"
        ],
        "initial_particle_escape_rejected": _RESULTS[
            "initial_particle_escape_rejected"
        ],
        "reflective_wall_does_not_enter_escape_ledger": _RESULTS[
            "reflective_wall_does_not_enter_escape_ledger"
        ],
        "paper_particle_recenter_rejected": (
            _RESULTS["paper_particle_recenter_rejected"]
        ),
        "segment_saw_escape_diagnostic": _RESULTS["segment_saw_escape_diagnostic"],
        "continuation_saw_escape_diagnostic": (
            _RESULTS["continuation_saw_escape_diagnostic"]
        ),
        "segment_population_audits_bounded": (
            1 <= _RESULTS["segment_telemetry"]["population_audit_calls"] <= 5
        ),
        "continuation_population_audits_bounded": (
            1 <= _RESULTS["continued_telemetry"]["population_audit_calls"] <= 5
        ),
        "telemetry_policy_is_checkpoint_restart_run_end": (
            _RESULTS["segment_telemetry"]["population_audit_policy"]
            == "checkpoint_restart_run_end"
            and _RESULTS["continued_telemetry"]["population_audit_policy"]
            == "checkpoint_restart_run_end"
        ),
        "telemetry_destruction_audits_match_ledgers": (
            _RESULTS["segment_telemetry"]["destruction_audit_calls"]
            == segment["ps_escape_audit_calls"]
            and _RESULTS["continued_telemetry"]["destruction_audit_calls"]
            == continued["ps_escape_audit_calls"]
        ),
        "production_performance_pilot_gate": {
            "status": "required_non_authorizing",
            "population_census_policy": "checkpoint_restart_run_end",
            "required_telemetry": "Q017 migration and checkpoint timings at production scale",
            "runtime_marks_gate_required": (
                _RESULTS["segment_telemetry"]["production_pilot_required"]
                and _RESULTS["continued_telemetry"]["production_pilot_required"]
            ),
        },
        "mpi2": _RESULTS["mpi2"],
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    segment_basename = "pic_parallel_shock_escape_segment"
    continued_basename = "pic_parallel_shock_escape_continued"
    for basename in [segment_basename, continued_basename]:
        _remove_outputs(basename)

    _run_expect_recenter_reject()
    _RESULTS["paper_particle_recenter_rejected"] = True
    _run_expect_initial_escape_reject()
    _RESULTS["initial_particle_escape_rejected"] = True
    _RESULTS["reflective_wall_does_not_enter_escape_ledger"] = (
        _run_reflective_wall_success()
    )

    segment_output = _run_success("segment", segment_basename, 100)
    segment_restart = _latest_restart(segment_basename)
    _RESULTS["segment"] = _restart_metadata(segment_restart)
    _RESULTS["segment_telemetry"] = _escape_telemetry(segment_output)
    _RESULTS["segment_saw_escape_diagnostic"] = (
        "pic_parallel_shock outer_x1_escape_sink:" in segment_output
    )

    numeric_error = "pic_parallel_shock restart CR ledger numeric metadata is invalid"
    chronology_error = "pic_parallel_shock restart paper-VL2 escape-audit chronology is invalid"
    corruptions = {
        "incomplete": (["problem/ps_escape_ledger_complete=false"], numeric_error),
        "zero_audit_calls": (["problem/ps_escape_audit_calls=0"], numeric_error),
        "inconsistent_mass": (
            ["problem/ps_escaped_injected_cr_mass_global=0"],
            numeric_error,
        ),
        "nonzero_initial_count": (
            ["problem/ps_escaped_initial_cr_count_global=1"],
            numeric_error,
        ),
        "wrong_positive_audit_calls": (
            ["problem/ps_escape_audit_calls=199"],
            chronology_error,
        ),
        "stale_audit_time": (
            ["problem/ps_escape_last_audit_time=0"],
            chronology_error,
        ),
        "future_audit_time": (
            ["problem/ps_escape_last_audit_time=1e99"],
            numeric_error,
        ),
        "nonfinite_momentum": (
            ["problem/ps_escaped_injected_cr_momentum_x1_global=nan"],
            numeric_error,
        ),
        "empty_escape_with_accumulated_state": (
            [
                "problem/ps_escaped_injected_cr_count_global=0",
                "problem/ps_escaped_injected_cr_mass_global=0",
            ],
            numeric_error,
        ),
    }
    for label, (extra, expected) in corruptions.items():
        _run_expect_reject(
            label,
            continued_basename + "_" + label,
            101,
            segment_restart,
            extra,
            expected,
        )
    _RESULTS["corrupt_ledgers_rejected"] = sorted(corruptions)

    duplicate_results = []
    duplicate_cases = {
        "duplicate_escape_ledger": ("ps_escape_ledger_complete", "false"),
        "duplicate_cr_ledger": ("ps_cr_ledger_complete", "false"),
    }
    for label, (field, value) in duplicate_cases.items():
        duplicate_restart = _duplicate_restart_parameter(
            segment_restart, continued_basename + "_" + label, field, value
        )
        try:
            _restart_metadata(duplicate_restart)
        except RuntimeError as exc:
            if "Duplicate restart problem parameter" not in str(exc):
                raise
        else:
            raise RuntimeError("Regression parser accepted " + label)
        _run_expect_reject(
            label,
            continued_basename + "_" + label + "_runtime",
            101,
            duplicate_restart,
            [],
            "pic_parallel_shock restart contains duplicate CR or escape ledger metadata",
        )
        duplicate_results.append(label)
    _RESULTS["duplicate_ledgers_rejected_by_runtime_and_parser"] = sorted(
        duplicate_results
    )

    continued_output = _run_success(
        "continuation", continued_basename, 120, restart_path=segment_restart
    )
    _RESULTS["continued"] = _restart_metadata(_latest_restart(continued_basename))
    _RESULTS["continued_telemetry"] = _escape_telemetry(continued_output)
    _RESULTS["continuation_saw_escape_diagnostic"] = (
        "pic_parallel_shock outer_x1_escape_sink:" in continued_output
    )

    if _MPI_NPROC >= 2:
        mpi_basename = "pic_parallel_shock_escape_mpi2"
        _remove_outputs(mpi_basename)
        mpi_output = _run_success(
            "mpi2 chronology",
            mpi_basename,
            20,
            extra=["mesh/nx1=16", "meshblock/nx1=8"],
            nproc=2,
        )
        mpi_metadata = _restart_metadata(_latest_restart(mpi_basename))
        mpi_telemetry = _escape_telemetry(mpi_output)
        _RESULTS["mpi2"] = {
            "requested": True,
            "passed": (
                mpi_metadata["ps_escape_audit_calls"]
                == _PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE * 20
                and mpi_telemetry["destruction_audit_calls"]
                == mpi_metadata["ps_escape_audit_calls"]
            ),
        }
    else:
        _RESULTS["mpi2"] = {"requested": False, "passed": False}


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("PIC parallel-shock outer-x1 escape metrics: %s", summary)
    return (
        all(summary["segment"].values())
        and all(summary["continued"].values())
        and summary["continuation_audit_calls_advance"]
        and summary["segment_has_exact_vl2_audits"]
        and summary["continuation_has_exact_vl2_audits"]
        and summary["continuation_time_advances"]
        and summary["continuation_injected_count_nondecreasing"]
        and summary["continuation_escape_count_nondecreasing"]
        and summary["corrupt_ledgers_rejected"]
        == [
            "empty_escape_with_accumulated_state",
            "future_audit_time",
            "incomplete",
            "inconsistent_mass",
            "nonfinite_momentum",
            "nonzero_initial_count",
            "stale_audit_time",
            "wrong_positive_audit_calls",
            "zero_audit_calls",
        ]
        and summary["duplicate_ledgers_rejected_by_runtime_and_parser"]
        == ["duplicate_cr_ledger", "duplicate_escape_ledger"]
        and summary["initial_particle_escape_rejected"]
        and summary["reflective_wall_does_not_enter_escape_ledger"]
        and summary["paper_particle_recenter_rejected"]
        and summary["segment_saw_escape_diagnostic"]
        and summary["continuation_saw_escape_diagnostic"]
        and summary["segment_population_audits_bounded"]
        and summary["continuation_population_audits_bounded"]
        and summary["telemetry_policy_is_checkpoint_restart_run_end"]
        and summary["telemetry_destruction_audits_match_ledgers"]
        and summary["production_performance_pilot_gate"]["runtime_marks_gate_required"]
        and (
            not summary["mpi2"]["requested"]
            or summary["mpi2"]["passed"]
        )
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
