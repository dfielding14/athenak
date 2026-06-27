"""Keep one particle while it crosses a 2D fine-to-coarse SMR corner."""

import glob
import logging
import os
import re
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger("athena" + __name__[7:])

_INPUT_DECK = "tests/pic_particle_2d_smr_corner_routing.athinput"
_BASENAME = "pic_particle_2d_smr_corner_routing"
_RESULT = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), "build", "src")


def _input_path():
    return "../../" + athena.athena_rel_path + "inputs/" + _INPUT_DECK


def _last_telemetry_value(output, field):
    matches = re.findall(
        rf"^q017\.telemetry\.{re.escape(field)}=([^\n]+)$", output, re.MULTILINE
    )
    if not matches:
        raise RuntimeError("Missing Q017 telemetry field: " + field)
    return float(matches[-1])


def _final_particle_count():
    pattern = os.path.join(
        _athena_exe_dir(), "pvtk", _BASENAME + ".prtcl_all.*.part.vtk"
    )
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError("Missing particle VTK output: " + pattern)
    with open(paths[-1], "rb") as handle:
        header = handle.read(1024)
    match = re.search(rb"\nPOINTS\s+([0-9]+)\s+float\n", header)
    if match is None:
        raise RuntimeError("Particle VTK POINTS header is missing")
    return int(match.group(1))


def run(**kwargs):
    logger.debug("Running test " + __name__)
    for path in glob.glob(
        os.path.join(_athena_exe_dir(), "pvtk", _BASENAME + ".*")
    ):
        os.remove(path)
    proc = subprocess.run(
        ["./athena", "-i", _input_path()],
        cwd=_athena_exe_dir(),
        capture_output=True,
        text=True,
    )
    output = (proc.stdout or "") + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError("2D SMR corner-routing run failed\n" + output)
    _RESULT.update(
        {
            "output": output,
            "particle_count": _final_particle_count(),
            "telemetry_particle_count": _last_telemetry_value(
                output, "particles.total"
            ),
            "coarse_level_count": _last_telemetry_value(
                output, "particle_memory.level.1.count"
            ),
            "fine_level_count": _last_telemetry_value(
                output, "particle_memory.level.2.count"
            ),
            "invalid_records": _last_telemetry_value(
                output, "particle_memory.invalid_records"
            ),
        }
    )


def analyze():
    logger.debug("Analyzing test " + __name__)
    output = _RESULT["output"]
    return (
        _RESULT["particle_count"] == 1
        and _RESULT["telemetry_particle_count"] == 1.0
        and _RESULT["coarse_level_count"] == 1.0
        and _RESULT["fine_level_count"] == 0.0
        and _RESULT["invalid_records"] == 0.0
        and "invalid particle destruction" not in output
        and "invalid_neighbor" not in output
        and "FATAL ERROR" not in output
    )
