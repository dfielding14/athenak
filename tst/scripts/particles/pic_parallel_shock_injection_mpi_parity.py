"""Serial/MPI2 N=5 surface-smoothed subtraction parity at two macro-masses."""

from __future__ import annotations

import glob
import json
import logging
import os
import re
import shlex
import subprocess
import sys
import tempfile

import numpy as np
import scripts.utils.athena as athena

_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_PUBLICATION_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "publication")
)
sys.path.insert(0, _PUBLICATION_DIR)
sys.path.insert(0, os.path.join(_SOURCE_ROOT, "vis", "python"))
import bin_convert_new as bin_convert  # noqa: E402
from pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_parallel_shock_injection_mpi_parity.athinput"
_MPI_LAUNCHER_ENV = "ATHENA_PIC_PARALLEL_SHOCK_MPI_PARITY_MPI_LAUNCHER"
_RESTART_FLOAT_FIELDS = [
    "ps_mass_reservoir_global",
    "ps_injected_cr_count_global",
    "ps_injected_cr_mass_global",
    "ps_injected_cr_momentum_x1_global",
    "ps_injected_cr_momentum_x2_global",
    "ps_injected_cr_momentum_x3_global",
    "ps_injected_cr_energy_global",
]
_LEDGER_FIELDS = _RESTART_FLOAT_FIELDS[1:]
_REQUIRED_PARTICLE_SCALARS = {
    "gid",
    "ptag",
    "species",
    "cr_source",
    "macro_weight",
    "birth_time",
}
_MACRO_MASSES = {
    "production": 9.0e-4,
    "failed_pilot": 1.44e-2,
}
_METADATA_ABS_TOL = 2.0e-12
_TRANSACTION_ABS_TOL = 1.0e-5
_TRANSVERSE_ABS_TOL = 2.0e-12
_SUBTRACTION_STENCIL_CELLS = 5
_INJECTION_TIME_MIN = 4.0
_INJECTION_TIME_MAX = 4.2
_SOURCE_TRANSACTION_DIAG_RE = re.compile(
    r"pic_parallel_shock source_transaction_diag: cycle=([0-9]+) .*?"
    r"applied=\(([^,]+),([^,]+),([^,]+),([^,]+),([^)]+)\) .*?"
    r"expected=\(([^,]+),([^,]+),([^,]+),([^,]+),([^)]+)\)"
)
# Four blocks split both x1 and the ideal shock surface. Injection at
# t=[4,4.2] uses carrier cell 4, so each row's five downstream targets are
# cells 4 through 0 and cross the x1 MeshBlock boundary.
_CROSS_BLOCK_STENCIL_OVERRIDES = [
    "meshblock/nx1=4",
    "meshblock/nx2=4",
    "time/cfl_number=0.8",
    "time/nlim=128",
    "time/tlim=" + str(_INJECTION_TIME_MAX),
    "problem/ps_inject_t_start=" + str(_INJECTION_TIME_MIN),
    "problem/ps_inject_t_stop=" + str(_INJECTION_TIME_MAX),
    "problem/ps_enable_gas_subtraction=true",
    "problem/ps_feedback_diag_dcycle=1",
]
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        "ATHENA_PIC_PARALLEL_SHOCK_MPI_PARITY_EXE_DIR",
        os.path.join(os.getcwd(), "build", "src"),
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _stage_stencil_input():
    with open(_athena_input_path(), "r", encoding="utf-8") as source:
        payload = source.read()
    marker = "ps_enable_gas_subtraction     = false\n"
    if payload.count(marker) != 1:
        raise RuntimeError("Unable to locate unique gas-subtraction deck parameter")
    payload = payload.replace(
        marker,
        marker
        + "ps_subtract_stencil_cells     = "
        + str(_SUBTRACTION_STENCIL_CELLS)
        + "\n"
        + "ps_enable_surface_averaged_subtraction = true\n",
    )
    if "<output3>" in payload:
        raise RuntimeError("Parity input unexpectedly already defines output3")
    payload += """

<output3>
file_type   = bin
variable    = mhd_u
id          = mhd_u
dcycle      = 1
ghost_zones = false
"""
    descriptor, path = tempfile.mkstemp(
        prefix="pic_parallel_shock_injection_mpi_parity_",
        suffix=".athinput",
        dir=_athena_exe_dir(),
        text=True,
    )
    with os.fdopen(descriptor, "w", encoding="utf-8") as staged:
        staged.write(payload)
    return path


def _mpi_launcher(nproc):
    raw_launcher = os.environ.get(
        _MPI_LAUNCHER_ENV, os.environ.get("MPIEXEC", "mpiexec")
    )
    command = shlex.split(raw_launcher)
    if not command:
        raise RuntimeError(_MPI_LAUNCHER_ENV + " must not be empty")
    replaced = False
    for index, value in enumerate(command):
        if "{nproc}" in value:
            command[index] = value.replace("{nproc}", str(nproc))
            replaced = True
    if not replaced:
        command += ["-n", str(nproc)]
    return command


def _remove_outputs(basename):
    for dirname in ["bin", "pvtk", "rst"]:
        pattern = os.path.join(_athena_exe_dir(), dirname, basename + ".*")
        for path in glob.glob(pattern):
            if os.path.isfile(path):
                os.remove(path)


def _latest_file(dirname, pattern, label):
    paths = sorted(glob.glob(os.path.join(_athena_exe_dir(), dirname, pattern)))
    if not paths:
        raise RuntimeError("No " + label + " files found for pattern: " + pattern)
    return paths[-1]


def _restart_parameters(path):
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
        if active_block is None or "=" not in line:
            continue
        key, value = line.split("=", 1)
        parameters.setdefault(active_block, {})[key.strip()] = value.split(
            "#", 1
        )[0].strip()
    return parameters


def _restart_metadata(basename):
    path = _latest_file("rst", basename + ".*.rst", "restart")
    blocks = _restart_parameters(path)
    problem = blocks.get("problem", {})
    particles = blocks.get("particles", {})
    required = [
        "ps_next_tag",
        "ps_subtract_stencil_cells",
        "ps_enable_surface_averaged_subtraction",
        *_RESTART_FLOAT_FIELDS,
    ]
    missing = [name for name in required if name not in problem]
    if "deposit_qscale" not in particles:
        missing.append("particles/deposit_qscale")
    if missing:
        raise RuntimeError(
            "Restart metadata is missing injection-accounting keys "
            + repr(missing)
            + ": "
            + path
        )
    smoothing_control = problem[
        "ps_enable_surface_averaged_subtraction"
    ].lower()
    if smoothing_control not in {"true", "false"}:
        raise RuntimeError(
            "Restart smoothing control is not a strict boolean: " + path
        )
    return {
        "ps_next_tag": int(problem["ps_next_tag"]),
        "ps_subtract_stencil_cells": int(
            problem["ps_subtract_stencil_cells"]
        ),
        "ps_enable_surface_averaged_subtraction": smoothing_control == "true",
        "deposit_qscale": float(particles["deposit_qscale"]),
        **{name: float(problem[name]) for name in _RESTART_FLOAT_FIELDS},
    }


def _particle_snapshot(basename):
    path = _latest_file(
        "pvtk", basename + ".prtcl_all.*.part.vtk", "particle VTK"
    )
    data = read_particle_vtk(path)
    missing = sorted(_REQUIRED_PARTICLE_SCALARS - set(data.scalars))
    if missing:
        raise RuntimeError("Particle VTK is missing scalars " + repr(missing))
    if "vel" not in data.vectors:
        raise RuntimeError("Particle VTK is missing vector: vel")

    tags = data.scalars["ptag"]
    if np.unique(tags).size != tags.size:
        raise RuntimeError("Particle VTK ptag values must be unique")
    order = np.argsort(tags)
    return {
        "path": path,
        "points": np.asarray(data.points)[order],
        "scalars": {
            name: np.asarray(values)[order] for name, values in data.scalars.items()
        },
        "vectors": {
            name: np.asarray(values)[order] for name, values in data.vectors.items()
        },
    }


def _gas_snapshot(basename):
    path = _latest_file("bin", basename + ".mhd_u.*.bin", "MHD binary")
    data = bin_convert.read_binary_as_athdf(path, dtype=np.float64)
    names = ("dens", "mom1", "mom2", "mom3", "ener")
    missing = [name for name in names if name not in data]
    if missing:
        raise RuntimeError("MHD binary is missing fields " + repr(missing))
    return {
        "path": path,
        "fields": {name: np.asarray(data[name]) for name in names},
    }


def _source_transactions(output):
    transactions = []
    for match in _SOURCE_TRANSACTION_DIAG_RE.finditer(output):
        applied = np.array([float(value) for value in match.groups()[1:6]])
        expected = np.array([float(value) for value in match.groups()[6:11]])
        if expected[0] <= 0.0:
            continue
        transactions.append(
            {
                "cycle": int(match.group(1)),
                "applied": applied,
                "expected": expected,
            }
        )
    if not transactions:
        raise RuntimeError("No nonzero gas-subtraction transaction diagnostics found")
    return {
        "cycles": np.array([item["cycle"] for item in transactions], dtype=np.int64),
        "applied": np.stack([item["applied"] for item in transactions]),
        "expected": np.stack([item["expected"] for item in transactions]),
    }


def _run_case(name, nproc, input_path, macro_mass):
    basename = "pic_parallel_shock_injection_parity_" + name
    _remove_outputs(basename)
    command = [
        "./athena",
        "-i",
        input_path,
        "job/basename=" + basename,
        *_CROSS_BLOCK_STENCIL_OVERRIDES,
        "particles/deposit_qscale=" + repr(macro_mass),
    ]
    if nproc > 1:
        command = _mpi_launcher(nproc) + command
    logger.info("Executing %s: %s", name, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + name + "\n" + output)
    return {
        "basename": basename,
        "metadata": _restart_metadata(basename),
        "particles": _particle_snapshot(basename),
        "gas": _gas_snapshot(basename),
        "source_transactions": _source_transactions(output),
        "requested_macro_mass": macro_mass,
    }


def _array_map_parity(left, right):
    if set(left) != set(right):
        return False
    return all(np.array_equal(left[name], right[name]) for name in left)


def _metadata_absolute_errors(serial, mpi2):
    return {
        name: abs(serial[name] - mpi2[name]) for name in _RESTART_FLOAT_FIELDS
    }


def _case_invariants(case):
    metadata = case["metadata"]
    scalars = case["particles"]["scalars"]
    count = int(scalars["ptag"].size)
    next_tag = metadata["ps_next_tag"]
    ledger = np.array([metadata[name] for name in _LEDGER_FIELDS])
    transactions = case["source_transactions"]
    transaction_residual = transactions["applied"] - transactions["expected"]
    transaction_sum = np.sum(transactions["expected"], axis=0)
    birth_times = scalars["birth_time"]
    macro_mass = case["requested_macro_mass"]
    gas_fields = case["gas"]["fields"]
    gas_is_transversely_uniform = all(
        field.ndim == 3
        and field.shape[1] > 1
        and np.max(np.abs(field - field[:, :1, :])) <= _TRANSVERSE_ABS_TOL
        for field in gas_fields.values()
    )
    return {
        "particle_count": count,
        "next_tag": next_tag,
        "particle_count_matches_next_tag": count == next_tag,
        "tags_are_dense_from_zero": np.array_equal(
            scalars["ptag"], np.arange(next_tag, dtype=np.int64)
        ),
        "all_particles_are_shock_injected": bool(np.all(scalars["cr_source"] == 1)),
        "all_particles_are_species_zero": bool(np.all(scalars["species"] == 0)),
        "all_macro_weights_are_one": bool(np.all(scalars["macro_weight"] == 1.0)),
        "all_birth_times_are_in_cross_block_window": bool(
            np.all(birth_times >= _INJECTION_TIME_MIN)
            and np.all(birth_times <= _INJECTION_TIME_MAX)
        ),
        "subtraction_stencil_width_is_five": (
            metadata["ps_subtract_stencil_cells"] == _SUBTRACTION_STENCIL_CELLS
        ),
        "surface_averaged_subtraction_is_enabled": (
            metadata["ps_enable_surface_averaged_subtraction"] is True
        ),
        "macro_mass_matches_requested_case": (
            metadata["deposit_qscale"] == macro_mass
        ),
        "ledger_is_finite": bool(np.all(np.isfinite(ledger))),
        "ledger_count_matches_next_tag": bool(ledger[0] == next_tag),
        "ledger_mass_matches_macro_mass": bool(
            abs(ledger[1] - next_tag * macro_mass) <= _METADATA_ABS_TOL
        ),
        "reservoir_is_bounded": (
            0.0 <= metadata["ps_mass_reservoir_global"] < macro_mass
        ),
        "source_transactions_are_finite": bool(
            np.all(np.isfinite(transactions["applied"]))
            and np.all(np.isfinite(transactions["expected"]))
        ),
        "source_transactions_close": bool(
            np.max(np.abs(transaction_residual)) <= _TRANSACTION_ABS_TOL
        ),
        "transaction_sum_matches_injected_ledger": bool(
            np.max(np.abs(transaction_sum - ledger[1:]))
            <= _TRANSACTION_ABS_TOL * transactions["cycles"].size
        ),
        "gas_fields_are_finite": bool(
            all(np.all(np.isfinite(field)) for field in gas_fields.values())
        ),
        "gas_state_is_transversely_uniform": bool(
            gas_is_transversely_uniform
        ),
    }


def _pair_summary(macro_label):
    pair = _RESULTS[macro_label]
    serial = pair["serial"]
    mpi2 = pair["mpi2"]
    serial_metadata = serial["metadata"]
    mpi2_metadata = mpi2["metadata"]
    serial_particles = serial["particles"]
    mpi2_particles = mpi2["particles"]
    serial_transactions = serial["source_transactions"]
    mpi2_transactions = mpi2["source_transactions"]
    metadata_errors = _metadata_absolute_errors(serial_metadata, mpi2_metadata)
    transaction_shape_equal = (
        serial_transactions["applied"].shape == mpi2_transactions["applied"].shape
    )
    transaction_applied_error = float("inf")
    transaction_expected_error = float("inf")
    if transaction_shape_equal:
        transaction_applied_error = float(
            np.max(
                np.abs(
                    serial_transactions["applied"] - mpi2_transactions["applied"]
                )
            )
        )
        transaction_expected_error = float(
            np.max(
                np.abs(
                    serial_transactions["expected"] - mpi2_transactions["expected"]
                )
            )
        )
    return {
        "requested_macro_mass": _MACRO_MASSES[macro_label],
        "serial": _case_invariants(serial),
        "mpi2": _case_invariants(mpi2),
        "next_tag_equal": (
            serial_metadata["ps_next_tag"] == mpi2_metadata["ps_next_tag"]
        ),
        "subtraction_stencil_width_equal": (
            serial_metadata["ps_subtract_stencil_cells"]
            == mpi2_metadata["ps_subtract_stencil_cells"]
            == _SUBTRACTION_STENCIL_CELLS
        ),
        "surface_averaged_subtraction_equal": (
            serial_metadata["ps_enable_surface_averaged_subtraction"]
            is mpi2_metadata["ps_enable_surface_averaged_subtraction"]
            is True
        ),
        "macro_mass_equal": (
            serial_metadata["deposit_qscale"]
            == mpi2_metadata["deposit_qscale"]
            == _MACRO_MASSES[macro_label]
        ),
        "metadata_absolute_errors": metadata_errors,
        "metadata_within_tolerance": all(
            error <= _METADATA_ABS_TOL for error in metadata_errors.values()
        ),
        "particle_points_equal": np.array_equal(
            serial_particles["points"], mpi2_particles["points"]
        ),
        "particle_scalars_equal": _array_map_parity(
            serial_particles["scalars"], mpi2_particles["scalars"]
        ),
        "particle_vectors_equal": _array_map_parity(
            serial_particles["vectors"], mpi2_particles["vectors"]
        ),
        "gas_fields_equal": _array_map_parity(
            serial["gas"]["fields"], mpi2["gas"]["fields"]
        ),
        "source_transaction_cycles_equal": np.array_equal(
            serial_transactions["cycles"], mpi2_transactions["cycles"]
        ),
        "source_transaction_applied_max_absolute_error": (
            transaction_applied_error
        ),
        "source_transaction_expected_max_absolute_error": (
            transaction_expected_error
        ),
        "source_transactions_within_tolerance": (
            transaction_shape_equal
            and transaction_applied_error <= _TRANSACTION_ABS_TOL
            and transaction_expected_error <= _TRANSACTION_ABS_TOL
        ),
    }


def _summary():
    return {
        "evidence_class": "bounded_local_serial_vs_mpi2_host_regression",
        "not_qualification_evidence": True,
        "physical_setup": (
            "same_input_deck_four_meshblocks_surface_averaged_n5_downstream_"
            "subtraction_stencil_at_production_and_failed_pilot_macro_mass"
        ),
        "macro_mass_cases": {
            label: _pair_summary(label) for label in _MACRO_MASSES
        },
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    input_path = _stage_stencil_input()
    try:
        for label, macro_mass in _MACRO_MASSES.items():
            _RESULTS[label] = {
                "serial": _run_case(
                    label + "_serial", 1, input_path, macro_mass
                ),
                "mpi2": _run_case(
                    label + "_mpi2", 2, input_path, macro_mass
                ),
            }
    finally:
        os.remove(input_path)


def _pair_passes(summary):
    return (
        summary["serial"]["particle_count"] > 0
        and all(
            value
            for name, value in summary["serial"].items()
            if name not in {"particle_count", "next_tag"}
        )
        and all(
            value
            for name, value in summary["mpi2"].items()
            if name not in {"particle_count", "next_tag"}
        )
        and summary["next_tag_equal"]
        and summary["subtraction_stencil_width_equal"]
        and summary["surface_averaged_subtraction_equal"]
        and summary["macro_mass_equal"]
        and summary["metadata_within_tolerance"]
        and summary["particle_points_equal"]
        and summary["particle_scalars_equal"]
        and summary["particle_vectors_equal"]
        and summary["gas_fields_equal"]
        and summary["source_transaction_cycles_equal"]
        and summary["source_transactions_within_tolerance"]
    )


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("PIC parallel-shock injection MPI parity metrics: %s", summary)
    return all(
        _pair_passes(pair)
        for pair in summary["macro_mass_cases"].values()
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
