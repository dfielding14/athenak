"""Disjoint adaptive restart byte-layout compatibility regression.

The executable directory and optional schema-absent legacy adaptive checkpoint
are supplied with ATHENA_DISJOINT_ADAPTIVE_RESTART_EXE_DIR and
ATHENA_DISJOINT_ADAPTIVE_RESTART_LEGACY_CHECKPOINT. Use
ATHENA_DISJOINT_ADAPTIVE_RESTART_REQUIRE_LEGACY=1 for a strict legacy run.
"""

from __future__ import annotations

import glob
import json
import logging
import os
import re
import shlex
import struct
import subprocess

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_disjoint_adaptive_restart_compatibility.athinput"
_SOURCE_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "..")
)
_EXE_DIR_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_EXE_DIR"
_LEGACY_CHECKPOINT_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_LEGACY_CHECKPOINT"
_REQUIRE_LEGACY_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_REQUIRE_LEGACY"
_NPROC_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_NPROC"
_LAUNCHER_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_MPI_LAUNCHER"
_OUTPUT_MODE_ENV = "ATHENA_DISJOINT_ADAPTIVE_RESTART_OUTPUT_MODE"
_MESH_METADATA_MAGIC = 0x4154484B4D455348
_RESTART_HEADER_LIMIT = 16 * 1024 * 1024
_RESULTS = {}


def _athena_exe_dir():
    return os.path.abspath(
        os.environ.get(_EXE_DIR_ENV, os.path.join(os.getcwd(), "build", "src"))
    )


def _athena_input_path():
    return os.path.join(_SOURCE_ROOT, "inputs", _INPUT_DECK)


def _positive_integer_env(name, default):
    raw = os.environ.get(name, str(default))
    try:
        value = int(raw)
    except ValueError as error:
        raise RuntimeError(name + " must be a positive integer") from error
    if value < 1:
        raise RuntimeError(name + " must be a positive integer")
    return value


def _boolean_env(name, default=False):
    raw = os.environ.get(name, str(default)).strip().lower()
    if raw in {"1", "true", "yes", "on"}:
        return True
    if raw in {"0", "false", "no", "off"}:
        return False
    raise RuntimeError(name + " must be a boolean")


def _rank_count():
    return _positive_integer_env(_NPROC_ENV, 1)


def _launcher_prefix():
    nproc = _rank_count()
    if nproc == 1:
        return []
    raw = os.environ.get(_LAUNCHER_ENV, os.environ.get("MPIEXEC", "mpiexec"))
    launcher = shlex.split(raw)
    if not launcher:
        raise RuntimeError(_LAUNCHER_ENV + " must not be empty")
    replaced = False
    for index, value in enumerate(launcher):
        if "{nproc}" in value:
            launcher[index] = value.replace("{nproc}", str(nproc))
            replaced = True
    if not replaced:
        launcher += ["-n", str(nproc)]
    return launcher


def _output_modes():
    raw = os.environ.get(_OUTPUT_MODE_ENV, "shared").strip().lower()
    if raw == "both":
        return ["shared", "per-rank"]
    if raw in {"shared", "per-rank"}:
        return [raw]
    raise RuntimeError(_OUTPUT_MODE_ENV + " must be shared, per-rank, or both")


def _athena_mpi_enabled():
    proc = subprocess.run(
        ["./athena", "-c"], cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise RuntimeError("Unable to query Athena configuration with -c")
    output = (proc.stdout or "") + (proc.stderr or "")
    return "MPI parallelism:            ON" in output


def _require_runtime_configuration():
    executable = os.path.join(_athena_exe_dir(), "athena")
    if not os.path.isfile(executable) or not os.access(executable, os.X_OK):
        raise RuntimeError(
            "Athena executable is missing or not executable; set "
            + _EXE_DIR_ENV
            + ": "
            + executable
        )
    if _rank_count() > 1 and not _athena_mpi_enabled():
        raise RuntimeError(_NPROC_ENV + " > 1 requires an MPI-enabled executable")


def _remove_outputs(basename):
    patterns = [
        os.path.join(_athena_exe_dir(), "rst", basename + ".*"),
        os.path.join(_athena_exe_dir(), "rst", "rank_*", basename + ".*"),
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            if os.path.isfile(path):
                os.remove(path)


def _run_athena(label, arguments, restart_path=None):
    command = _launcher_prefix() + ["./athena"]
    if restart_path is None:
        command += ["-i", _athena_input_path()]
    else:
        command += ["-r", restart_path]
    command += list(arguments)
    logger.info("Executing %s: %s", label, " ".join(command))
    proc = subprocess.run(
        command, cwd=_athena_exe_dir(), capture_output=True, text=True
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("Command failed for " + label + "\n" + output)
    return output


def _mode_arguments(mode):
    return ["output1/single_file_per_rank=" + str(mode == "per-rank").lower()]


def _seed_checkpoint(basename, mode):
    if mode == "shared":
        pattern = os.path.join(_athena_exe_dir(), "rst", basename + ".*.rst")
    else:
        pattern = os.path.join(
            _athena_exe_dir(), "rst", "rank_00000000", basename + ".*.rst"
        )
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError("No seed checkpoint found for pattern: " + pattern)
    return paths[-1]


def _relative_to_exe(path):
    return os.path.relpath(path, _athena_exe_dir())


def _resolved_checkpoint(path):
    if os.path.isabs(path):
        return path
    return os.path.join(_athena_exe_dir(), path)


def _contains_mesh_metadata(path):
    patterns = [
        struct.pack(byte_order + "Q", _MESH_METADATA_MAGIC)
        for byte_order in ["<", ">"]
    ]
    overlap = b""
    with open(_resolved_checkpoint(path), "rb") as handle:
        while True:
            payload = handle.read(1024 * 1024)
            if not payload:
                return False
            payload = overlap + payload
            if any(pattern in payload for pattern in patterns):
                return True
            overlap = payload[-7:]


def _restart_parameter(path, block, name):
    header = bytearray()
    with open(_resolved_checkpoint(path), "rb") as handle:
        while b"<par_end>" not in header:
            payload = handle.read(65536)
            if not payload:
                break
            header.extend(payload)
            if len(header) > _RESTART_HEADER_LIMIT:
                raise RuntimeError("Restart parameter header exceeds limit: " + path)
    text = header.decode("latin1", errors="ignore")
    if "<par_end>" not in text:
        raise RuntimeError("Restart parameter header is missing <par_end>: " + path)
    active_block = None
    for raw_line in text.split("<par_end>", 1)[0].splitlines():
        line = raw_line.strip()
        if line.startswith("<") and line.endswith(">"):
            active_block = line[1:-1]
            continue
        if active_block != block or "=" not in line:
            continue
        key, value = line.split("=", 1)
        if key.strip() == name:
            return value.split("#", 1)[0].strip()
    raise RuntimeError("Restart parameter is missing <" + block + ">/" + name)


def _created_meshblocks(output):
    matches = re.findall(r"(\d+) MeshBlocks created,\s*(\d+) deleted by AMR", output)
    if not matches:
        raise RuntimeError("Missing AMR transition telemetry\n" + output)
    return max(int(created) for created, _ in matches)


def _probe_reload(label, checkpoint, mode, refinement):
    basename = "pic_disjoint_adaptive_restart_" + label
    _remove_outputs(basename)
    arguments = [
        "job/basename=" + basename,
        "mesh_refinement/refinement=" + refinement,
        "time/nlim=0",
        "output1/dcycle=0",
        "output1/file_number=0",
    ] + _mode_arguments(mode)
    _run_athena(label, arguments, restart_path=checkpoint)


def _legacy_checkpoint():
    checkpoint = os.environ.get(_LEGACY_CHECKPOINT_ENV)
    if checkpoint is None or not checkpoint.strip():
        if _boolean_env(_REQUIRE_LEGACY_ENV):
            raise RuntimeError(
                _REQUIRE_LEGACY_ENV + " requires " + _LEGACY_CHECKPOINT_ENV
            )
        return None
    checkpoint = checkpoint.strip()
    resolved = _resolved_checkpoint(checkpoint)
    if not os.path.isfile(resolved):
        raise RuntimeError("Legacy checkpoint is missing: " + resolved)
    if _contains_mesh_metadata(checkpoint):
        raise RuntimeError(
            "Legacy checkpoint unexpectedly contains mesh metadata schema: " + resolved
        )
    if _restart_parameter(checkpoint, "mesh_refinement", "refinement") != "adaptive":
        raise RuntimeError("Legacy checkpoint is not adaptive: " + resolved)
    return checkpoint


def _run_mode(mode, legacy_checkpoint):
    seed_basename = "pic_disjoint_adaptive_restart_seed_" + mode.replace("-", "_")
    _remove_outputs(seed_basename)
    output = _run_athena(
        "seed_" + mode,
        [
            "job/basename=" + seed_basename,
            "time/nlim=1",
            "output1/file_number=0",
        ]
        + _mode_arguments(mode),
    )
    checkpoint = _seed_checkpoint(seed_basename, mode)
    if not _contains_mesh_metadata(checkpoint):
        raise RuntimeError("New checkpoint lacks mesh metadata schema: " + checkpoint)
    if _restart_parameter(checkpoint, "mesh_refinement", "refinement") != "adaptive":
        raise RuntimeError("New checkpoint is not adaptive: " + checkpoint)
    created = _created_meshblocks(output)
    if created <= 0:
        raise RuntimeError("Adaptive seed did not create refined MeshBlocks")

    probes = {}
    for refinement in ["adaptive", "none", "static"]:
        label = "new_" + mode.replace("-", "_") + "_" + refinement
        _probe_reload(label, _relative_to_exe(checkpoint), mode, refinement)
        probes[refinement] = True
    if legacy_checkpoint is not None:
        label = "legacy_" + mode.replace("-", "_")
        _probe_reload(label, legacy_checkpoint, mode, "adaptive")
        probes["legacy_schema_absent"] = True
    else:
        probes["legacy_schema_absent"] = None
    return {
        "new_checkpoint": checkpoint,
        "new_checkpoint_has_mesh_metadata_schema": True,
        "adaptive_seed_created_meshblocks": created,
        "reloads": probes,
    }


def _summary():
    return {
        "evidence_class": "bounded_restart_byte_layout_compatibility_regression",
        "configured_ranks": _RESULTS["configured_ranks"],
        "executable_dir": _RESULTS["executable_dir"],
        "legacy_checkpoint": _RESULTS["legacy_checkpoint"],
        "legacy_checkpoint_probe": (
            "passed" if _RESULTS["legacy_checkpoint"] is not None else "not configured"
        ),
        "output_modes": _RESULTS["output_modes"],
    }


def run(**kwargs):
    logger.debug("Running test " + __name__)
    _RESULTS.clear()
    _require_runtime_configuration()
    legacy_checkpoint = _legacy_checkpoint()
    _RESULTS["configured_ranks"] = _rank_count()
    _RESULTS["executable_dir"] = _athena_exe_dir()
    _RESULTS["legacy_checkpoint"] = legacy_checkpoint
    _RESULTS["output_modes"] = {
        mode: _run_mode(mode, legacy_checkpoint) for mode in _output_modes()
    }


def analyze():
    logger.debug("Analyzing test " + __name__)
    summary = _summary()
    logger.info("Disjoint adaptive restart compatibility metrics: %s", summary)
    for result in summary["output_modes"].values():
        if not result["new_checkpoint_has_mesh_metadata_schema"]:
            return False
        if result["adaptive_seed_created_meshblocks"] <= 0:
            return False
        if not all(
            result["reloads"][refinement]
            for refinement in ["adaptive", "none", "static"]
        ):
            return False
        if (
            summary["legacy_checkpoint"] is not None
            and not result["reloads"]["legacy_schema_absent"]
        ):
            return False
    return True


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
