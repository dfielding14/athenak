"""Repeated MPI migration, physical destruction, restart, and ledger parity."""

from __future__ import annotations

import glob
import logging
import math
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena  # noqa: F401,E402


_SOURCE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_SOURCE_ROOT))
from tst.publication import q011_section54_restart as restart_layout  # noqa: E402
from tst.publication.pvtk_particles import read_particle_vtk  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = _SOURCE_ROOT / "inputs/tests/pic_migration_restart_ledger.athinput"
_MPIEXEC = shlex.split(os.environ.get("MPIEXEC", "mpiexec"))
_FINAL_CYCLE = 24
_SPLIT_CYCLE = 12
_INITIAL_COUNT = 192
_FINAL_COUNT = 128
_EVENT_CYCLES = {4, 10, 16, 22}
_LEDGER_NAMES = ("mass", "mom1", "mom2", "mom3", "energy")
_CYCLE_LEDGER_RE = re.compile(
    r"^pic_migration_ledger_cycle: cycle=(?P<cycle>[0-9]+)"
    r" time=(?P<time>\S+)"
    r" escaped_count=(?P<count>[0-9]+)"
    r" escape_mass=(?P<mass>\S+)"
    r" escape_mom1=(?P<mom1>\S+)"
    r" escape_mom2=(?P<mom2>\S+)"
    r" escape_mom3=(?P<mom3>\S+)"
    r" escape_energy=(?P<energy>\S+)"
    r" boundary_errors=(?P<errors>[0-9]+)$",
    re.MULTILINE,
)
_RESULTS = {}


def _athena_exe_dir():
    return Path(
        os.environ.get(
            "ATHENA_PIC_MIGRATION_RESTART_EXE_DIR",
            os.path.join(os.getcwd(), "build", "src"),
        )
    ).resolve()


def _athena_mpi_enabled():
    proc = subprocess.run(
        ["./athena", "-c"],
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration with -c")
    output = (proc.stdout or "") + (proc.stderr or "")
    return "MPI parallelism:            ON" in output


def _mpi_prefix(nproc):
    if not _MPIEXEC:
        raise RuntimeError("MPIEXEC must not be empty")
    return [*_MPIEXEC, "-n", str(nproc)]


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    patterns = [
        exe_dir / "pvtk" / (basename + ".*.part.vtk"),
        exe_dir / "rst" / (basename + ".*"),
    ]
    for pattern in patterns:
        for path in glob.glob(str(pattern)):
            if os.path.isfile(path):
                os.remove(path)


def _latest_restart(basename):
    paths = sorted(
        glob.glob(str(_athena_exe_dir() / "rst" / (basename + ".*.rst")))
    )
    if not paths:
        raise RuntimeError("No restart output found for " + basename)
    return Path(paths[-1])


def _execute(label, basename, nproc, nlim, restart_path=None, restart_dcycle=None):
    command = ["./athena"]
    if restart_path is None:
        command += ["-i", str(_INPUT_DECK)]
    else:
        command += ["-r", str(restart_path)]
    command += [
        "job/basename=" + basename,
        "time/nlim=" + str(nlim),
        "output2/dcycle=" + str(restart_dcycle or nlim),
    ]
    if nproc > 1:
        command = [*_mpi_prefix(nproc), *command]
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command,
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
        check=False,
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _cycle_ledgers(output):
    records = {}
    for match in _CYCLE_LEDGER_RE.finditer(output):
        cycle = int(match.group("cycle"))
        if cycle in records:
            raise RuntimeError("Duplicate migration cycle ledger: " + str(cycle))
        records[cycle] = {
            "time": float(match.group("time")),
            "count": int(match.group("count")),
            "errors": int(match.group("errors")),
            **{name: float(match.group(name)) for name in _LEDGER_NAMES},
        }
    return records


def _problem_parameters(payload, source):
    end = payload.find(b"<par_end>")
    if end < 0:
        raise RuntimeError(str(source) + ": restart parameter header is missing")
    text = payload[:end].decode("ascii")
    active = None
    values = {}
    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("<") and line.endswith(">"):
            active = line[1:-1]
            continue
        if active != "problem" or "=" not in line:
            continue
        key, value = (item.strip() for item in line.split("=", 1))
        if key in values:
            raise RuntimeError(str(source) + ": duplicate problem parameter " + key)
        values[key] = value
    return values


def _restart_snapshot(path):
    payload = path.read_bytes()
    probe = restart_layout.probe_current_restart_payload(
        payload, source=str(path)
    )
    if (
        probe.real_fields_per_particle != 26
        or probe.integer_fields_per_particle != 4
        or probe.state_kind != 1
        or probe.physical_mode != 4
    ):
        raise RuntimeError(str(path) + ": unexpected PIC restart schema contract")
    real = np.frombuffer(
        payload,
        dtype="<f8",
        count=probe.particle_count * probe.real_fields_per_particle,
        offset=probe.particle_real_offset,
    ).reshape(probe.particle_count, probe.real_fields_per_particle).copy()
    integer = np.frombuffer(
        payload,
        dtype="<i4",
        count=probe.particle_count * probe.integer_fields_per_particle,
        offset=probe.particle_integer_offset,
    ).reshape(probe.particle_count, probe.integer_fields_per_particle).copy()
    order = np.argsort(integer[:, 1], kind="stable")
    real = real[order]
    integer = integer[order]
    if np.unique(integer[:, 1]).size != probe.particle_count:
        raise RuntimeError(str(path) + ": restart particle tags are not unique")

    parameters = _problem_parameters(payload, path)
    required = {
        "mig1_ledger_schema",
        "mig1_ledger_complete",
        "mig1_committed_cycles",
        "mig1_committed_time",
        "mig1_escaped_count",
        "mig1_boundary_errors",
        *("mig1_escape_" + name for name in _LEDGER_NAMES),
        *("mig1_reflect_" + name for name in _LEDGER_NAMES),
    }
    missing = sorted(required - set(parameters))
    if missing:
        raise RuntimeError(str(path) + ": missing ledger fields " + repr(missing))
    ledger = {
        "schema": int(parameters["mig1_ledger_schema"]),
        "complete": parameters["mig1_ledger_complete"].lower() in {"1", "true"},
        "cycles": int(parameters["mig1_committed_cycles"]),
        "time": float(parameters["mig1_committed_time"]),
        "count": int(parameters["mig1_escaped_count"]),
        "errors": int(parameters["mig1_boundary_errors"]),
        "escape": np.array(
            [float(parameters["mig1_escape_" + name]) for name in _LEDGER_NAMES]
        ),
        "reflect": np.array(
            [float(parameters["mig1_reflect_" + name]) for name in _LEDGER_NAMES]
        ),
    }
    return {"real": real, "integer": integer, "ledger": ledger}


def _vtk_cycle(path):
    with open(path, "rb") as handle:
        header = handle.read(512).decode("ascii", errors="ignore")
    match = re.search(r"\bcycle=\s*([0-9]+)", header)
    if match is None:
        raise RuntimeError("Particle VTK cycle header is missing: " + str(path))
    return int(match.group(1))


def _trajectory(basename):
    paths = sorted(
        glob.glob(
            str(
                _athena_exe_dir()
                / "pvtk"
                / (basename + ".prtcl_all.*.part.vtk")
            )
        )
    )
    if not paths:
        raise RuntimeError("No particle trajectory outputs for " + basename)
    result = {}
    for path_string in paths:
        path = Path(path_string)
        cycle = _vtk_cycle(path)
        data = read_particle_vtk(str(path))
        required = {"ptag", "gid", "species", "cr_source"}
        if not required.issubset(data.scalars):
            raise RuntimeError("Particle VTK metadata is incomplete: " + str(path))
        tags = data.scalars["ptag"].astype(np.int64)
        order = np.argsort(tags, kind="stable")
        snapshot = np.stack(
            [
                tags[order],
                data.scalars["gid"][order].astype(np.int64),
                data.scalars["species"][order].astype(np.int64),
                data.scalars["cr_source"][order].astype(np.int64),
            ],
            axis=1,
        )
        if cycle in result and not np.array_equal(result[cycle], snapshot):
            raise RuntimeError("Conflicting trajectory snapshots at cycle " + str(cycle))
        result[cycle] = snapshot
    return result


def _combine_by_cycle(first, second, label):
    combined = dict(first)
    for cycle, value in second.items():
        if cycle in combined:
            if isinstance(value, np.ndarray):
                equal = np.array_equal(combined[cycle], value)
            else:
                equal = combined[cycle] == value
            if not equal:
                raise RuntimeError(label + " disagrees at cycle " + str(cycle))
        combined[cycle] = value
    return combined


def _run_full(case, nproc):
    basename = "pic_mig1_" + case + "_full"
    _remove_outputs(basename)
    output = _execute(case + " full", basename, nproc, _FINAL_CYCLE)
    restart = _latest_restart(basename)
    return {
        "final": _restart_snapshot(restart),
        "trajectory": _trajectory(basename),
        "cycle_ledgers": _cycle_ledgers(output),
    }


def _run_restarted(case, nproc):
    segment = "pic_mig1_" + case + "_segment"
    continued = "pic_mig1_" + case + "_continued"
    for basename in (segment, continued):
        _remove_outputs(basename)
    segment_output = _execute(
        case + " segment", segment, nproc, _SPLIT_CYCLE
    )
    segment_restart = _latest_restart(segment)
    segment_snapshot = _restart_snapshot(segment_restart)
    continued_output = _execute(
        case + " continuation",
        continued,
        nproc,
        _FINAL_CYCLE,
        restart_path=segment_restart,
        restart_dcycle=_FINAL_CYCLE,
    )
    final_restart = _latest_restart(continued)
    trajectory = _combine_by_cycle(
        _trajectory(segment), _trajectory(continued), case + " restart trajectory"
    )
    cycle_ledgers = _combine_by_cycle(
        _cycle_ledgers(segment_output),
        _cycle_ledgers(continued_output),
        case + " restart ledger",
    )
    return {
        "segment": segment_snapshot,
        "final": _restart_snapshot(final_restart),
        "trajectory": trajectory,
        "cycle_ledgers": cycle_ledgers,
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _RESULTS["serial_full"] = _run_full("serial", 1)
    _RESULTS["serial_restart"] = _run_restarted("serial", 1)

    mpi_enabled = _athena_mpi_enabled()
    launcher_available = bool(_MPIEXEC) and shutil.which(_MPIEXEC[0]) is not None
    _RESULTS["mpi_enabled"] = mpi_enabled
    if mpi_enabled and not launcher_available:
        raise RuntimeError(
            "MPI-enabled migration regression requires an available MPIEXEC launcher"
        )
    if mpi_enabled:
        _RESULTS["mpi2_full"] = _run_full("mpi2", 2)
        _RESULTS["mpi2_restart"] = _run_restarted("mpi2", 2)
    else:
        logger.info("Skipping MPI2 migration closure: Athena MPI support is disabled")


def _payload_equal(left, right):
    return np.array_equal(left["integer"], right["integer"]) and np.array_equal(
        left["real"].view(np.uint64), right["real"].view(np.uint64)
    )


def _ledger_close(left, right, tolerance=2.0e-13):
    scalar_equal = all(
        left[name] == right[name]
        for name in ("schema", "complete", "cycles", "count", "errors")
    )
    return (
        scalar_equal
        and math.isclose(left["time"], right["time"], rel_tol=0.0, abs_tol=tolerance)
        and np.allclose(left["escape"], right["escape"], rtol=0.0, atol=tolerance)
        and np.allclose(left["reflect"], right["reflect"], rtol=0.0, atol=tolerance)
    )


def _cycle_ledgers_close(left, right, tolerance=2.0e-13):
    if set(left) != set(right):
        return False
    for cycle in left:
        lrecord = left[cycle]
        rrecord = right[cycle]
        if lrecord["count"] != rrecord["count"]:
            return False
        if lrecord["errors"] != rrecord["errors"]:
            return False
        for name in ("time", *_LEDGER_NAMES):
            if not math.isclose(
                lrecord[name], rrecord[name], rel_tol=0.0, abs_tol=tolerance
            ):
                return False
    return True


def _canonical_trajectory_ok(full, restarted):
    expected_cycles = set(range(_FINAL_CYCLE + 1))
    return (
        set(full) == expected_cycles
        and set(restarted) == expected_cycles
        and all(np.array_equal(full[cycle], restarted[cycle]) for cycle in full)
    )


def _crossing_metrics(trajectory):
    initial = trajectory[0]
    final = trajectory[_FINAL_CYCLE]
    initial_by_tag = {int(row[0]): row for row in initial}
    final_tags = set(int(tag) for tag in final[:, 0])
    survivor_crossings = {}
    for tag in sorted(final_tags):
        owners = []
        for cycle in range(_FINAL_CYCLE + 1):
            rows = trajectory[cycle]
            matches = rows[rows[:, 0] == tag]
            if len(matches) != 1:
                return {"valid": False}
            owners.append(0 if int(matches[0, 1]) < 2 else 1)
        survivor_crossings[tag] = [
            cycle
            for cycle in range(1, _FINAL_CYCLE + 1)
            if owners[cycle] != owners[cycle - 1]
        ]

    destruction_by_cycle = {}
    simultaneous_each_rank = True
    for cycle in range(1, _FINAL_CYCLE + 1):
        previous = {int(tag) for tag in trajectory[cycle - 1][:, 0]}
        current = {int(tag) for tag in trajectory[cycle][:, 0]}
        removed = previous - current
        if not removed:
            continue
        destruction_by_cycle[cycle] = removed
        destroyed_owners = {
            0 if int(initial_by_tag[tag][1]) < 2 else 1 for tag in removed
        }
        crossing_origins = set()
        for tag, crossings in survivor_crossings.items():
            if cycle not in crossings:
                continue
            prior = trajectory[cycle - 1]
            row = prior[prior[:, 0] == tag][0]
            crossing_origins.add(0 if int(row[1]) < 2 else 1)
        simultaneous_each_rank &= destroyed_owners == {0, 1}
        simultaneous_each_rank &= crossing_origins == {0, 1}

    all_repeat = all(len(cycles) >= 2 for cycles in survivor_crossings.values())
    all_straddle = all(
        any(cycle <= _SPLIT_CYCLE for cycle in cycles)
        and any(cycle > _SPLIT_CYCLE for cycle in cycles)
        for cycles in survivor_crossings.values()
    )
    removed_tags = set().union(*destruction_by_cycle.values())
    removed_species = {int(initial_by_tag[tag][2]) for tag in removed_tags}
    return {
        "valid": True,
        "survivors": len(survivor_crossings),
        "minimum_crossings": min(map(len, survivor_crossings.values())),
        "all_repeat": all_repeat,
        "all_straddle_restart": all_straddle,
        "destruction_cycles": set(destruction_by_cycle),
        "destruction_counts": {
            cycle: len(tags) for cycle, tags in destruction_by_cycle.items()
        },
        "removed_species": removed_species,
        "simultaneous_each_rank": simultaneous_each_rank,
    }


def _final_inventory_ok(snapshot):
    integer = snapshot["integer"]
    return (
        integer.shape == (_FINAL_COUNT, 4)
        and np.unique(integer[:, 1]).size == _FINAL_COUNT
        and set(integer[:, 2]) == {0, 1}
        and np.count_nonzero(integer[:, 2] == 0) == _FINAL_COUNT // 2
        and np.count_nonzero(integer[:, 2] == 1) == _FINAL_COUNT // 2
        and np.all(integer[:, 3] == 0)
        and np.all(np.isfinite(snapshot["real"]))
    )


def _expected_final_ledger():
    qscale = 1.0 / 1024.0
    count = 64
    state_y = 1.6
    light_speed = 1.0e6
    kinetic = state_y * state_y / (
        math.sqrt(1.0 + state_y * state_y / (light_speed * light_speed)) + 1.0
    )
    mass = count * qscale
    return np.array([mass, 0.0, mass * state_y, 0.0, mass * kinetic])


def _ledger_contract_ok(ledger, expected_cycles, expected_count):
    return (
        ledger["schema"] == 1
        and ledger["complete"]
        and ledger["cycles"] == expected_cycles
        and ledger["count"] == expected_count
        and ledger["errors"] == 0
        and np.allclose(ledger["reflect"], 0.0, rtol=0.0, atol=0.0)
    )


def analyze():
    logger.debug("Analyzing test " + __name__)
    serial_full = _RESULTS["serial_full"]
    serial_restart = _RESULTS["serial_restart"]
    expected = _expected_final_ledger()
    serial_metrics = {
        "full_inventory": _final_inventory_ok(serial_full["final"]),
        "restart_inventory": _final_inventory_ok(serial_restart["final"]),
        "restart_payload_bitwise": _payload_equal(
            serial_full["final"], serial_restart["final"]
        ),
        "trajectory_exact": _canonical_trajectory_ok(
            serial_full["trajectory"], serial_restart["trajectory"]
        ),
        "ledger_exact": _ledger_close(
            serial_full["final"]["ledger"], serial_restart["final"]["ledger"]
        ),
        "cycle_ledgers_exact": _cycle_ledgers_close(
            serial_full["cycle_ledgers"], serial_restart["cycle_ledgers"]
        ),
        "final_ledger_contract": _ledger_contract_ok(
            serial_full["final"]["ledger"], _FINAL_CYCLE, 64
        ),
        "segment_ledger_contract": _ledger_contract_ok(
            serial_restart["segment"]["ledger"], _SPLIT_CYCLE, 32
        ),
        "final_ledger_expected": np.allclose(
            serial_full["final"]["ledger"]["escape"],
            expected,
            rtol=0.0,
            atol=2.0e-13,
        ),
    }
    if not _RESULTS["mpi_enabled"]:
        logger.info("MIG-1 serial restart metrics: %s", serial_metrics)
        return all(serial_metrics.values())

    mpi_full = _RESULTS["mpi2_full"]
    mpi_restart = _RESULTS["mpi2_restart"]
    crossing = _crossing_metrics(mpi_full["trajectory"])
    event_counts_ok = (
        crossing.get("destruction_cycles") == _EVENT_CYCLES
        and all(count == 16 for count in crossing.get("destruction_counts", {}).values())
    )
    cycle_events_ok = set(
        cycle
        for cycle, record in mpi_full["cycle_ledgers"].items()
        if record["count"] > 0
    ) == _EVENT_CYCLES
    metrics = {
        **{"serial_" + key: value for key, value in serial_metrics.items()},
        "mpi_full_inventory": _final_inventory_ok(mpi_full["final"]),
        "mpi_restart_inventory": _final_inventory_ok(mpi_restart["final"]),
        "all_final_payloads_bitwise": all(
            _payload_equal(serial_full["final"], candidate)
            for candidate in (
                serial_restart["final"],
                mpi_full["final"],
                mpi_restart["final"],
            )
        ),
        "mpi_restart_trajectory_exact": _canonical_trajectory_ok(
            mpi_full["trajectory"], mpi_restart["trajectory"]
        ),
        "all_final_ledgers_equal": all(
            _ledger_close(serial_full["final"]["ledger"], candidate["ledger"])
            for candidate in (
                serial_restart["final"],
                mpi_full["final"],
                mpi_restart["final"],
            )
        ),
        "all_cycle_ledgers_equal": all(
            _cycle_ledgers_close(serial_full["cycle_ledgers"], candidate)
            for candidate in (
                serial_restart["cycle_ledgers"],
                mpi_full["cycle_ledgers"],
                mpi_restart["cycle_ledgers"],
            )
        ),
        "mpi_final_ledger_expected": np.allclose(
            mpi_full["final"]["ledger"]["escape"],
            expected,
            rtol=0.0,
            atol=2.0e-13,
        ),
        "cycle_event_ledger_exact": cycle_events_ok,
        "destruction_events_exact": event_counts_ok,
        "removed_species_exact": crossing.get("removed_species") == {2},
        "all_survivors_repeat_cross": crossing.get("all_repeat", False),
        "all_survivors_cross_before_and_after_restart": crossing.get(
            "all_straddle_restart", False
        ),
        "simultaneous_send_destroy_each_rank": crossing.get(
            "simultaneous_each_rank", False
        ),
    }
    logger.info("MIG-1 repeated migration metrics: %s; crossing=%s", metrics, crossing)
    return all(metrics.values())
