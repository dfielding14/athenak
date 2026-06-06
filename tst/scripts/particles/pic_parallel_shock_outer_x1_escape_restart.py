"""Bounded outer-x1 shock-particle escape ledger and restart regression."""

from __future__ import annotations

import glob
import json
import logging
import math
import os
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
_MACRO_MASS = 1.0e-3
_PAPER_VL2_ESCAPE_AUDITS_PER_CYCLE = 2
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
        parameters[key.strip()] = value.split("#", 1)[0].strip()
    return parameters


def _restart_metadata(path):
    parameters = _restart_problem_parameters(path)
    required = [
        "ps_escape_ledger_schema",
        "ps_escape_ledger_complete",
        "ps_escape_audit_calls",
        "ps_injected_cr_count_global",
        *_ESCAPE_FLOAT_FIELDS,
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
    }


def _execute(label, basename, nlim, restart_path=None, extra=()):
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
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")


def _run_success(label, basename, nlim, restart_path=None):
    code, output = _execute(label, basename, nlim, restart_path=restart_path)
    if code != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _run_expect_reject(label, basename, nlim, restart_path, override):
    code, output = _execute(
        label, basename, nlim, restart_path=restart_path, extra=[override]
    )
    expected = "pic_parallel_shock restart CR ledger numeric metadata is invalid"
    if code == 0 or expected not in output:
        raise RuntimeError("Escape-ledger corruption was not rejected: " + label)


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


def _validate_metadata(metadata):
    finite = all(
        math.isfinite(metadata[name])
        for name in ["ps_injected_cr_count_global", *_ESCAPE_FLOAT_FIELDS]
    )
    escaped = metadata["ps_escaped_injected_cr_count_global"]
    injected = metadata["ps_injected_cr_count_global"]
    return {
        "schema_is_one": metadata["ps_escape_ledger_schema"] == 1,
        "ledger_complete": metadata["ps_escape_ledger_complete"],
        "finite": finite,
        "audit_calls_positive": metadata["ps_escape_audit_calls"] > 0,
        "injected_positive": injected > 0.0,
        "escaped_positive": escaped > 0.0,
        "escaped_not_above_injected": escaped <= injected,
        "mass_matches_count": math.isclose(
            metadata["ps_escaped_injected_cr_mass_global"],
            escaped * _MACRO_MASS,
            rel_tol=2.0e-14,
            abs_tol=2.0e-14,
        ),
        "energy_positive": metadata["ps_escaped_injected_cr_energy_global"] > 0.0,
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
        "paper_particle_recenter_rejected": (
            _RESULTS["paper_particle_recenter_rejected"]
        ),
        "segment_saw_escape_diagnostic": _RESULTS["segment_saw_escape_diagnostic"],
        "continuation_saw_escape_diagnostic": (
            _RESULTS["continuation_saw_escape_diagnostic"]
        ),
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

    segment_output = _run_success("segment", segment_basename, 100)
    segment_restart = _latest_restart(segment_basename)
    _RESULTS["segment"] = _restart_metadata(segment_restart)
    _RESULTS["segment_saw_escape_diagnostic"] = (
        "pic_parallel_shock outer_x1_escape_sink:" in segment_output
    )

    corruptions = {
        "incomplete": "problem/ps_escape_ledger_complete=false",
        "zero_audit_calls": "problem/ps_escape_audit_calls=0",
        "inconsistent_mass": "problem/ps_escaped_injected_cr_mass_global=0",
        "negative_initial_count": "problem/ps_escaped_initial_cr_count_global=-1",
    }
    for label, override in corruptions.items():
        _run_expect_reject(
            label,
            continued_basename + "_" + label,
            101,
            segment_restart,
            override,
        )
    _RESULTS["corrupt_ledgers_rejected"] = sorted(corruptions)

    continued_output = _run_success(
        "continuation", continued_basename, 120, restart_path=segment_restart
    )
    _RESULTS["continued"] = _restart_metadata(_latest_restart(continued_basename))
    _RESULTS["continuation_saw_escape_diagnostic"] = (
        "pic_parallel_shock outer_x1_escape_sink:" in continued_output
    )


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
        == ["incomplete", "inconsistent_mass", "negative_initial_count", "zero_audit_calls"]
        and summary["paper_particle_recenter_rejected"]
        and summary["segment_saw_escape_diagnostic"]
        and summary["continuation_saw_escape_diagnostic"]
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
