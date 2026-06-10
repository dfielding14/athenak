#!/usr/bin/env python3
"""Run the CGL-LF Stage I matrix directly and continue clean partial runs.

This is intentionally a small execution harness, not a release controller.  It
pins the qualified executable and frozen inputs, gives every segment a unique
output directory, submits independent cases directly to Slurm, and performs a
minimal scientific sanity check before chaining a continuation.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex
import struct
import subprocess
import sys
import tempfile
from typing import Iterable


DEFAULT_ROOT = Path("/lustre/orion/ast207/proj-shared/dfielding/CGL")
FROZEN_SOURCE = Path("/autofs/nccs-svm1_home2/dfielding/athenak-cgl-e03-9e075422")
MATRIX_RELATIVE = Path("inputs/cgl_lf_paper/mks24_stage_i_manifest.json")
MATRIX_SHA256 = "bf31b88b985d1ad4ffe823108dd7c1132bdfa4d5e4a6abde51f66bb7778415c9"
ATHENA = (
    DEFAULT_ROOT
    / "build/frontier-hip-9e07542281e4-cpe25.09-cce20-rocm6.4.2/src/athena"
)
ATHENA_SHA256 = "68f243f9204df388b24365ae65a567f6f567dbe422a6d7a43b9fb4a499ef118c"
FAST_RUNS_RELATIVE = Path("runs/mks24-stage-i-fast/E03-forcing-policy")
LOGS_RELATIVE = Path("logs/slurm-fast")
TARGET_TIME = 10.0
RANKS_PER_NODE = 8
CPUS_PER_TASK = 7
REQUESTED_WALLTIME = "02:00:00"
ATHENA_WALLTIME = "01:50:00"
STRICT_FAILURE_COLUMNS = (
    "lf_dfloor",
    "lf_pfloor",
    "lf_nonfin",
    "lf_nonpos",
    "lf_hardbd",
)
CASE_NODES = {
    "R03": 1,
    **{f"R{number:02d}": 4 for number in range(4, 16)},
    "R16": 1,
    "R17": 8,
}
CASE_INPUT_SHA256 = {
    "R03": "997f449abb3c2e4d1de50509efffa127c6230e3ecb2632612200010b8fa6c0b0",
    "R04": "571eea2ccec5d069ccb1b49d132ba4c5b8bdac7ce15ee8d8a04c92327c453137",
    "R05": "6527a2072105d515287701c904c3af2cebb07e10ab0b45c3fb41f29d3467622c",
    "R06": "c6e038c198b23bf20a7cf1e2a90fa82e83e5ec544492ddc64e1cf2bb535a38ae",
    "R07": "72ed6a3b38342d00855f34a7c7eb67102b44cd22b6c12ff2448ddb6f141e5145",
    "R08": "3535aca5e3d47dc383e353cff262d79c56cf3c083a3ac6024d04ddaa5d353780",
    "R09": "2ec05b247d917259c8306911cfd31d6be4e1c2772d6cdea544db883edfc35f2c",
    "R10": "93d7ed9b0846019f59bce34feeb9b22b098afffdfc2dbaa59aaf5913655cd705",
    "R11": "ab0dba23f80ea2d6c173ecbb724dbcb332b13d3b1751b8776ad05313f9ccaf6c",
    "R12": "98ddea4b4f7fec18cc40abdbf5f7c8ba5b583a91f84f4dae00e1411f23d42e7c",
    "R13": "a190129ee34c46a4fe83a392a19120a58b9a9a4064742ac58511441370c59c35",
    "R14": "2b8d5837f8a7f3070ca2ef56f8b44e53d048839918185a8eef155cb084736578",
    "R15": "9a698a60bef4c558ccee4635d69c3acbf3fea478401bd633d09bbc0526d943d8",
    "R16": "c0ac4b54248e8f8dfb0f5fd34c0cfb4414b5330529cbf2836961c5277af3f2d1",
    "R17": "cc1092404b82129f807308a64f7585a6da31f45f41f1d2263acad0c8d30a7e04",
}
CASE_SEEDS = {
    "R03": {
        "time": 0.5,
        "restart": (
            DEFAULT_ROOT
            / "runs/mks24-stage-i/E03-forcing-policy/R03/"
            "s02_rankio_t0p312823_t0p5/output/rst/rank_00000000/"
            "E03_forcing_policy_paper_standard_active_alfvenic_beta100_"
            "s02_rankio_t0p312823_t0p5.00002.rst"
        ),
        "sha256": "89bf3613b4779f91ae47aed4a965677dd0059af46332bd5f77cb48b1faa2ed2d",
    },
    "R04": {
        "time": 0.25,
        "restart": (
            DEFAULT_ROOT
            / "runs/mks24-stage-i/E03-forcing-policy/R04/"
            "s01_rankio_t0_t0p25/output/rst/rank_00000000/"
            "E03_forcing_policy_paper_standard_active_random_beta10_"
            "s01_rankio_t0_t0p25.00001.rst"
        ),
        "sha256": "fb5c07d6247d5f6e6751dbb9a13b63221cee3d312466380397f52c1d206b57d2",
    },
    "R16": {
        "time": 1.5,
        "restart": (
            DEFAULT_ROOT
            / "runs/mks24-stage-i/E03-forcing-policy/R16/"
            "s00_rankio_t0_t1p5/output/rst/rank_00000000/"
            "E03_forcing_policy_paper_scale_separation_active_alfvenic_beta10_nperp96_"
            "s00_rankio_t0_t1p5.00002.rst"
        ),
        "sha256": "2b7f48e1a97b5e207a4081f5415c6642ca3041c84c7d6e2abd72d16385c3dcec",
    },
}
HISTORY_LABEL = re.compile(r"\[(\d+)\]=(\S+)")
SEGMENT_NAME = re.compile(r"fast_s(\d{3})_")


class FastRunError(RuntimeError):
    """A direct-run precondition or sanity check failed."""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_sha(path: Path, expected: str, label: str) -> None:
    if not path.is_file():
        raise FastRunError(f"missing {label}: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise FastRunError(f"{label} checksum mismatch: {path}: {actual}")


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(value, indent=2, sort_keys=True) + "\n"
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as stream:
        temporary = Path(stream.name)
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def load_json(path: Path) -> dict[str, object]:
    with path.open(encoding="utf-8") as stream:
        value = json.load(stream)
    if not isinstance(value, dict):
        raise FastRunError(f"expected JSON object: {path}")
    return value


def campaign_cases() -> dict[str, dict[str, object]]:
    matrix = FROZEN_SOURCE / MATRIX_RELATIVE
    require_sha(matrix, MATRIX_SHA256, "Stage I matrix")
    data = load_json(matrix)
    cases = data.get("cases")
    if not isinstance(cases, list):
        raise FastRunError("Stage I matrix has no cases")
    return {
        str(case["id"]): case
        for case in cases
        if isinstance(case, dict) and str(case.get("id")) in CASE_NODES
    }


def expand_cases(values: Iterable[str]) -> list[str]:
    result: list[str] = []
    for value in values:
        match = re.fullmatch(r"R(\d{2})-R(\d{2})", value)
        if match:
            start, stop = (int(item) for item in match.groups())
            if start > stop:
                raise FastRunError(f"descending case range is invalid: {value}")
            result.extend(f"R{number:02d}" for number in range(start, stop + 1))
        else:
            result.extend(item for item in value.split(",") if item)
    invalid = [case for case in result if case not in CASE_NODES]
    if invalid:
        raise FastRunError(f"unknown or unsupported cases: {', '.join(invalid)}")
    return list(dict.fromkeys(result))


def time_tag(value: float) -> str:
    text = f"{value:.12g}".replace(".", "p").replace("-", "m")
    return text


def case_root(root: Path, case_id: str) -> Path:
    return root / FAST_RUNS_RELATIVE / case_id


def segment_directories(root: Path, case_id: str) -> list[Path]:
    parent = case_root(root, case_id)
    if not parent.is_dir():
        return []
    return sorted(
        (path for path in parent.iterdir() if path.is_dir() and SEGMENT_NAME.match(path.name)),
        key=lambda path: int(SEGMENT_NAME.match(path.name).group(1)),  # type: ignore[union-attr]
    )


def segment_manifest(segment: Path) -> Path:
    return segment / "manifest/fast_run.json"


def parse_history(path: Path) -> dict[str, list[float]]:
    labels: list[str] | None = None
    rows: list[list[float]] = []
    with path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#"):
                found = HISTORY_LABEL.findall(line)
                if found:
                    labels = [name for _, name in sorted(found, key=lambda item: int(item[0]))]
                continue
            if line.strip():
                rows.append([float(item) for item in line.split()])
    if labels is None or not rows or any(len(row) != len(labels) for row in rows):
        raise FastRunError(f"malformed history: {path}")
    columns = {name: [row[index] for row in rows] for index, name in enumerate(labels)}
    if any(not math.isfinite(value) for values in columns.values() for value in values):
        raise FastRunError(f"non-finite history value: {path}")
    return columns


def restart_time(path: Path) -> float:
    with path.open("rb") as stream:
        payload = stream.read(128 * 1024)
    end = payload.find(b"<par_end>")
    if end < 0:
        raise FastRunError(f"restart lacks parameter terminator: {path}")
    text = payload[:end].decode("utf-8")
    block = ""
    values: list[float] = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
            continue
        if block == "time" and "=" in line:
            key, value = (item.strip() for item in line.split("=", 1))
            if key in {"time", "restart_time"}:
                values.append(float(value))
    if len(values) != 1 or not math.isfinite(values[0]):
        raise FastRunError(f"restart has ambiguous physical time: {path}")
    return values[0]


def binary_assignment(stream, path: Path, key: str) -> str:
    line = stream.readline(16 * 1024)
    if not line.endswith(b"\n"):
        raise FastRunError(f"binary product has truncated {key}: {path}")
    try:
        text = line.decode("ascii").strip()
    except UnicodeDecodeError as error:
        raise FastRunError(f"binary product has non-ASCII {key}: {path}") from error
    observed_key, separator, value = text.partition("=")
    if not separator or observed_key.strip() != key or not value.strip():
        raise FastRunError(f"binary product lacks {key}: {path}")
    return value.strip()


def binary_parameter_values(payload: bytes, path: Path) -> dict[str, str]:
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise FastRunError(f"binary product parameter header is not UTF-8: {path}") from error
    block = ""
    values: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line.startswith("<") and line.endswith(">"):
            block = line[1:-1].strip()
        elif block and "=" in line:
            key, value = (item.strip() for item in line.split("=", 1))
            values[f"{block}/{key}"] = value
    return values


def binary_profile(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        file_size = os.fstat(stream.fileno()).st_size
        if stream.readline(16 * 1024) != b"Athena binary output version=1.1\n":
            raise FastRunError(f"binary product has invalid magic/version: {path}")
        try:
            preheader = int(binary_assignment(stream, path, "size of preheader"))
            physical_time = float(binary_assignment(stream, path, "time"))
            cycle = int(binary_assignment(stream, path, "cycle"))
            location_size = int(binary_assignment(stream, path, "size of location"))
            variable_size = int(binary_assignment(stream, path, "size of variable"))
            variable_count = int(binary_assignment(stream, path, "number of variables"))
        except ValueError as error:
            raise FastRunError(f"binary product has invalid numeric header: {path}") from error
        if (
            preheader != 5
            or not math.isfinite(physical_time)
            or cycle < 0
            or location_size not in (4, 8)
            or variable_size not in (4, 8)
            or not 0 < variable_count <= 4096
        ):
            raise FastRunError(f"binary product has invalid header values: {path}")
        variable_line = stream.readline(64 * 1024)
        if not variable_line.endswith(b"\n"):
            raise FastRunError(f"binary product has truncated variable inventory: {path}")
        try:
            fields = variable_line.decode("ascii").strip().split()
        except UnicodeDecodeError as error:
            raise FastRunError(
                f"binary product has non-ASCII variable inventory: {path}"
            ) from error
        if (
            not fields
            or fields[0] != "variables:"
            or len(fields[1:]) != variable_count
            or len(set(fields[1:])) != variable_count
        ):
            raise FastRunError(f"binary product has invalid variable inventory: {path}")
        try:
            header_size = int(binary_assignment(stream, path, "header offset"))
        except ValueError as error:
            raise FastRunError(f"binary product has invalid header offset: {path}") from error
        if header_size <= 0 or stream.tell() + header_size >= file_size:
            raise FastRunError(f"binary product has invalid parameter header size: {path}")
        parameter_header = stream.read(header_size)
        if len(parameter_header) != header_size or b"<par_end>\n" not in parameter_header:
            raise FastRunError(f"binary product has invalid parameter header: {path}")
        parameters = binary_parameter_values(parameter_header, path)
        try:
            mesh_shape = tuple(int(parameters[f"mesh/nx{axis}"]) for axis in range(1, 4))
            meshblock_shape = tuple(
                int(parameters[f"meshblock/nx{axis}"]) for axis in range(1, 4)
            )
        except (KeyError, ValueError) as error:
            raise FastRunError(f"binary product lacks mesh decomposition: {path}") from error
        if (
            any(extent <= 0 for extent in mesh_shape + meshblock_shape)
            or any(
                mesh_extent % block_extent != 0
                for mesh_extent, block_extent in zip(mesh_shape, meshblock_shape)
            )
        ):
            raise FastRunError(f"binary product has invalid mesh decomposition: {path}")
        logical_shape = tuple(
            mesh_extent // block_extent
            for mesh_extent, block_extent in zip(mesh_shape, meshblock_shape)
        )
        expected_locations = frozenset(
            (lx1, lx2, lx3, 0)
            for lx3 in range(logical_shape[2])
            for lx2 in range(logical_shape[1])
            for lx1 in range(logical_shape[0])
        )
        locations: list[tuple[int, int, int, int]] = []
        while stream.tell() < file_size:
            if stream.tell() + 40 + 6 * location_size > file_size:
                raise FastRunError(f"binary product has truncated meshblock header: {path}")
            indices = struct.unpack("<10i", stream.read(40))
            shape = (
                indices[1] - indices[0] + 1,
                indices[3] - indices[2] + 1,
                indices[5] - indices[4] + 1,
            )
            if any(extent <= 0 for extent in shape):
                raise FastRunError(f"binary product has invalid meshblock shape: {path}")
            location = tuple(indices[6:10])
            if location not in expected_locations or location in locations:
                raise FastRunError(f"binary product has invalid logical location: {path}")
            locations.append(location)
            stream.seek(6 * location_size, os.SEEK_CUR)
            payload_size = math.prod(shape) * variable_count * variable_size
            if payload_size <= 0 or stream.tell() + payload_size > file_size:
                raise FastRunError(f"binary product has truncated meshblock payload: {path}")
            stream.seek(payload_size, os.SEEK_CUR)
        if stream.tell() != file_size or not locations:
            raise FastRunError(f"binary product lacks a complete meshblock payload: {path}")
    return {
        "physical_time": physical_time,
        "locations": tuple(locations),
        "expected_locations": expected_locations,
    }


def binary_time(path: Path) -> float:
    return float(binary_profile(path)["physical_time"])


def terminal_product_group(directory: Path, suffix: str, rank_count: int) -> dict[str, object]:
    rank_zero = directory / "rank_00000000"
    if not rank_zero.is_dir():
        raise FastRunError(f"missing rank-local product directory: {rank_zero}")
    candidates = sorted(rank_zero.glob(f"*{suffix}"))
    if not candidates:
        raise FastRunError(f"no {suffix} products: {directory}")
    if suffix == ".rst":
        terminal = max(candidates, key=restart_time)
        physical_time = restart_time(terminal)
    elif suffix == ".bin":
        terminal = max(candidates, key=binary_time)
        physical_time = binary_time(terminal)
    else:
        terminal = candidates[-1]
        physical_time = None
    siblings = sorted(directory.glob(f"rank_*/{terminal.name}"))
    expected_ranks = {f"rank_{rank:08d}" for rank in range(rank_count)}
    actual_ranks = {path.parent.name for path in siblings}
    if actual_ranks != expected_ranks or any(not path.is_file() or path.stat().st_size == 0 for path in siblings):
        raise FastRunError(f"incomplete terminal {suffix} rank group: {terminal.name}")
    if suffix in {".rst", ".bin"}:
        if suffix == ".rst":
            times = [restart_time(path) for path in siblings]
        else:
            profiles = [binary_profile(path) for path in siblings]
            times = [float(profile["physical_time"]) for profile in profiles]
            expected_locations = profiles[0]["expected_locations"]
            locations = [
                location
                for profile in profiles
                for location in profile["locations"]
            ]
            if (
                any(
                    profile["expected_locations"] != expected_locations
                    for profile in profiles
                )
                or len(locations) != len(set(locations))
                or frozenset(locations) != expected_locations
            ):
                raise FastRunError(
                    f"incomplete terminal {suffix} logical coverage: {terminal.name}"
                )
        if any(
            not math.isclose(
                time, float(physical_time), rel_tol=0.0, abs_tol=1.0e-12
            )
            for time in times
        ):
            raise FastRunError(f"inconsistent terminal {suffix} rank times: {terminal.name}")
    return {
        "rank_zero": str(terminal.resolve()),
        "rank_count": len(siblings),
        "physical_time": physical_time,
        "sha256": sha256_file(terminal),
    }


def analyze_segment(segment: Path, *, save: bool = True) -> dict[str, object]:
    manifest = load_json(segment_manifest(segment))
    output = segment / "output"
    rank_count = int(manifest["ranks"])
    mhd_paths = sorted(output.glob("*.mhd.hst"))
    user_paths = sorted(output.glob("*.user.hst"))
    if len(mhd_paths) != 1 or len(user_paths) != 1:
        raise FastRunError(f"expected one MHD and one user history: {output}")
    mhd = parse_history(mhd_paths[0])
    user = parse_history(user_paths[0])
    if "time" not in mhd or "time" not in user or "mass" not in mhd or "mass" not in user:
        raise FastRunError(f"histories lack required columns: {output}")
    if len(mhd["time"]) != len(user["time"]) or any(
        left != right for left, right in zip(mhd["time"], user["time"])
    ):
        raise FastRunError(f"history times are not synchronized: {output}")
    final_time = mhd["time"][-1]
    restart = terminal_product_group(output / "rst", ".rst", rank_count)
    snapshot = terminal_product_group(output / "bin", ".bin", rank_count)
    if not math.isclose(float(restart["physical_time"]), final_time, rel_tol=0.0, abs_tol=1.0e-12):
        raise FastRunError(f"terminal restart and history times differ: {output}")
    if not math.isclose(
        float(snapshot["physical_time"]), final_time, rel_tol=0.0, abs_tol=1.0e-12
    ):
        raise FastRunError(f"terminal snapshot and history times differ: {output}")
    strict_maxima = {
        name: max(abs(value) for value in mhd.get(name, [math.inf]))
        for name in STRICT_FAILURE_COLUMNS
    }
    initial_mass = mhd["mass"][0]
    mass_drift = max(abs(value - initial_mass) for value in mhd["mass"]) / max(
        abs(initial_mass), 1.0
    )
    mass_mismatch = max(
        abs(left - right) for left, right in zip(mhd["mass"], user["mass"])
    ) / max(abs(initial_mass), 1.0)
    passed = (
        all(value == 0.0 for value in strict_maxima.values())
        and mass_drift <= 1.0e-8
        and mass_mismatch <= 1.0e-8
        and float(restart["physical_time"]) > float(manifest["start_time"])
    )
    result: dict[str, object] = {
        "schema_version": 1,
        "analyzed_utc": utc_now(),
        "case_id": manifest["case_id"],
        "segment": segment.name,
        "start_time": manifest["start_time"],
        "target_time": manifest["target_time"],
        "final_time": final_time,
        "complete": final_time >= float(manifest["target_time"]) - 1.0e-12,
        "passed": passed,
        "strict_lf_failure_maxima": strict_maxima,
        "mass_relative_drift": mass_drift,
        "mhd_user_mass_relative_mismatch": mass_mismatch,
        "terminal_restart": restart,
        "terminal_snapshot": snapshot,
        "history_rows": len(mhd["time"]),
    }
    if save:
        write_json(segment / "manifest/fast_analysis.json", result)
    return result


def verify_restart_group(path: Path, rank_count: int, expected_sha: str, expected_time: float) -> None:
    require_sha(path, expected_sha, "restart rank zero")
    if not math.isclose(restart_time(path), expected_time, rel_tol=0.0, abs_tol=1.0e-12):
        raise FastRunError(f"restart seed time differs: {path}")
    siblings = sorted(path.parent.parent.glob(f"rank_*/{path.name}"))
    expected = {f"rank_{rank:08d}" for rank in range(rank_count)}
    if {item.parent.name for item in siblings} != expected:
        raise FastRunError(f"restart seed rank group is incomplete: {path}")


def batch_script_text(manifest: dict[str, object]) -> str:
    restart = str(manifest.get("restart") or "")
    restart_sha = str(manifest.get("restart_sha256") or "")
    input_path = str(manifest["input"])
    output = str(manifest["output_dir"])
    run_dir = str(manifest["run_dir"])
    log = str(manifest["slurm_log"])
    fast_script = str(manifest["fast_script"])
    fast_script_sha = str(manifest["fast_script_sha256"])
    run_args = f'-r {shlex.quote(restart)}' if restart else f'-i {shlex.quote(input_path)}'
    restart_check = (
        f'require_sha {restart_sha} "$RESTART" restart_rank_zero\n'
        if restart
        else ""
    )
    return f"""#!/bin/bash
#SBATCH -J cglf_{manifest['case_id']}_s{int(manifest['sequence']):03d}
#SBATCH -A AST207
#SBATCH -o {log}
#SBATCH -p batch
#SBATCH -t {REQUESTED_WALLTIME}
#SBATCH -N {manifest['nodes']}
#SBATCH --gpus-per-node=8
#SBATCH --threads-per-core=1

set -euo pipefail
ATHENA={shlex.quote(str(ATHENA))}
INPUT={shlex.quote(input_path)}
RESTART={shlex.quote(restart)}
OUT_DIR={shlex.quote(output)}
RUN_DIR={shlex.quote(run_dir)}
NNODES="${{SLURM_NNODES:?Missing SLURM_NNODES}}"
NRANKS="$((NNODES * {RANKS_PER_NODE}))"

require_sha() {{
  local expected="$1" path="$2" label="$3" actual
  test -f "$path" || {{ echo "missing $label: $path" >&2; exit 1; }}
  actual="$(sha256sum "$path" | awk '{{print $1}}')"
  test "$actual" = "$expected" || {{ echo "$label checksum mismatch: $path" >&2; exit 1; }}
}}

require_sha {ATHENA_SHA256} "$ATHENA" executable
require_sha {manifest['input_sha256']} "$INPUT" input
require_sha {fast_script_sha} {shlex.quote(fast_script)} fast_launcher
{restart_check}mkdir "$OUT_DIR"

module restore
module load PrgEnv-cray
module load craype-accel-amd-gfx90a
module load cpe/25.09 cray-mpich/9.0.1 rocm/6.4.2
module load cce/20.0.0
module unload darshan-runtime

export LD_LIBRARY_PATH=${{CRAY_LD_LIBRARY_PATH}}:${{LD_LIBRARY_PATH:-}}
export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_GPU_MANAGED_MEMORY_SUPPORT_ENABLED=0
export MPICH_OFI_NIC_POLICY=GPU
export MPICH_GPU_IPC_CACHE_MAX_SIZE=1000
export MPICH_MPIIO_HINTS="*:romio_cb_write=disable"
export MPICH_OFI_NUM_CQ_ENTRIES=131072
export FI_MR_CACHE_MONITOR=kdreg2
export FI_CXI_RX_MATCH_MODE=software
export HSA_XNACK=1
export OMP_NUM_THREADS=1

date -u +"started_utc=%Y-%m-%dT%H:%M:%SZ" > "$RUN_DIR/manifest/run_environment.txt"
module -t list 2>&1 >> "$RUN_DIR/manifest/run_environment.txt"

set +e
srun -N "$NNODES" -n "$NRANKS" --ntasks-per-node={RANKS_PER_NODE} \\
  -c {CPUS_PER_TASK} --threads-per-core=1 --cpu-bind=threads \\
  --gpus-per-task=1 --gpu-bind=closest \\
  "$ATHENA" {run_args} -d "$OUT_DIR" -t {ATHENA_WALLTIME} \\
  job/basename={manifest['run_basename']} time/tlim={manifest['target_time']}
run_rc=$?
set -e
echo "$run_rc" > "$RUN_DIR/manifest/run_exit_code"
date -u +"finished_utc=%Y-%m-%dT%H:%M:%SZ" >> "$RUN_DIR/manifest/run_environment.txt"
if (( run_rc != 0 )); then
  exit "$run_rc"
fi

require_sha {fast_script_sha} {shlex.quote(fast_script)} fast_launcher
/usr/bin/python3.11 -B {shlex.quote(fast_script)} --root {shlex.quote(str(manifest['root']))} \\
  advance --segment "$RUN_DIR" --submit
"""


def prepare_segment(
    root: Path,
    case: dict[str, object],
    sequence: int,
    start_time: float,
    restart: Path | None,
    restart_sha: str | None,
) -> Path:
    require_sha(ATHENA, ATHENA_SHA256, "qualified executable")
    case_id = str(case["id"])
    nodes = CASE_NODES[case_id]
    ranks = nodes * RANKS_PER_NODE
    input_path = FROZEN_SOURCE / str(case["input"])
    input_sha = CASE_INPUT_SHA256[case_id]
    require_sha(input_path, input_sha, f"{case_id} frozen input")
    if restart is not None:
        if restart_sha is None:
            restart_sha = sha256_file(restart)
        verify_restart_group(restart, ranks, restart_sha, start_time)
    parent = case_root(root, case_id)
    parent.mkdir(parents=True, exist_ok=True)
    name = f"fast_s{sequence:03d}_t{time_tag(start_time)}_to_t{time_tag(TARGET_TIME)}"
    run_dir = parent / name
    run_dir.mkdir()
    manifest_dir = run_dir / "manifest"
    manifest_dir.mkdir()
    logs = root / LOGS_RELATIVE
    logs.mkdir(parents=True, exist_ok=True)
    run_basename = f"E03_fast_{case['name']}_{name}"
    fast_script = Path(__file__).resolve()
    manifest: dict[str, object] = {
        "schema_version": 1,
        "prepared_utc": utc_now(),
        "root": str(root.resolve()),
        "case_id": case_id,
        "case_name": case["name"],
        "input": str(input_path.resolve()),
        "input_sha256": input_sha,
        "executable": str(ATHENA),
        "executable_sha256": ATHENA_SHA256,
        "matrix_sha256": MATRIX_SHA256,
        "fast_script": str(fast_script),
        "fast_script_sha256": sha256_file(fast_script),
        "nodes": nodes,
        "ranks_per_node": RANKS_PER_NODE,
        "ranks": ranks,
        "sequence": sequence,
        "start_time": start_time,
        "target_time": TARGET_TIME,
        "restart": str(restart.resolve()) if restart is not None else None,
        "restart_sha256": restart_sha,
        "run_basename": run_basename,
        "run_dir": str(run_dir.resolve()),
        "output_dir": str((run_dir / "output").resolve()),
        "slurm_log": str((logs / "%x.%j.log").resolve()),
        "job_id": None,
    }
    manifest_path = segment_manifest(run_dir)
    write_json(manifest_path, manifest)
    script = manifest_dir / "run.sbatch"
    script.write_text(batch_script_text(manifest), encoding="utf-8")
    script.chmod(0o755)
    return run_dir


def submit_segment(segment: Path) -> str:
    manifest_path = segment_manifest(segment)
    manifest = load_json(manifest_path)
    if manifest.get("job_id") is not None:
        raise FastRunError(f"segment already submitted: {segment}")
    completed = subprocess.run(
        ["/usr/bin/sbatch", "--parsable", str(segment / "manifest/run.sbatch")],
        check=True,
        text=True,
        capture_output=True,
    )
    response = completed.stdout.strip()
    if not re.fullmatch(r"[1-9][0-9]*(?:;[A-Za-z0-9_.-]+)?", response):
        raise FastRunError(f"unexpected sbatch response: {completed.stdout!r}")
    job_id = response.split(";", 1)[0]
    manifest["job_id"] = job_id
    manifest["submitted_utc"] = utc_now()
    write_json(manifest_path, manifest)
    print(f"submitted {manifest['case_id']} {segment.name}: {job_id}")
    return job_id


def job_state(job_id: str) -> str:
    queued = subprocess.run(
        ["/usr/bin/squeue", "-h", "-j", job_id, "-o", "%T"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    if queued:
        return queued.splitlines()[0]
    accounted = subprocess.run(
        ["/usr/bin/sacct", "-X", "-n", "-P", "-j", job_id, "-o", "State"],
        text=True,
        capture_output=True,
        check=False,
    ).stdout.strip()
    return accounted.split("|", 1)[0] if accounted else "UNKNOWN"


def seed_for_case(case_id: str) -> tuple[float, Path | None, str | None]:
    seed = CASE_SEEDS.get(case_id)
    if seed is None:
        return 0.0, None, None
    return float(seed["time"]), Path(seed["restart"]), str(seed["sha256"])


def next_from_segment(segment: Path, *, submit: bool) -> Path | None:
    manifest = load_json(segment_manifest(segment))
    analysis = analyze_segment(segment)
    if not analysis["passed"]:
        raise FastRunError(f"scientific sanity check failed: {segment}")
    if analysis["complete"]:
        print(f"complete {manifest['case_id']}: t={analysis['final_time']}")
        return None
    restart = analysis["terminal_restart"]
    if not isinstance(restart, dict):
        raise FastRunError(f"missing terminal restart: {segment}")
    cases = campaign_cases()
    next_segment = prepare_segment(
        Path(str(manifest["root"])),
        cases[str(manifest["case_id"])],
        int(manifest["sequence"]) + 1,
        float(analysis["final_time"]),
        Path(str(restart["rank_zero"])),
        str(restart["sha256"]),
    )
    if submit:
        submit_segment(next_segment)
    else:
        print(next_segment)
    return next_segment


def launch_case(root: Path, case_id: str, *, submit: bool) -> Path | None:
    cases = campaign_cases()
    segments = segment_directories(root, case_id)
    if segments:
        latest = segments[-1]
        manifest = load_json(segment_manifest(latest))
        job_id = manifest.get("job_id")
        if job_id is None:
            if submit:
                submit_segment(latest)
            else:
                print(latest)
            return latest
        if job_id is not None and job_state(str(job_id)) in {
            "PENDING",
            "RUNNING",
            "CONFIGURING",
            "COMPLETING",
        }:
            print(f"active {case_id}: {job_id}")
            return None
        return next_from_segment(latest, submit=submit)
    start_time, restart, restart_sha = seed_for_case(case_id)
    segment = prepare_segment(root, cases[case_id], 0, start_time, restart, restart_sha)
    if submit:
        submit_segment(segment)
    else:
        print(segment)
    return segment


def command_status(root: Path, cases: list[str]) -> None:
    for case_id in cases:
        segments = segment_directories(root, case_id)
        if not segments:
            seed_time, _, _ = seed_for_case(case_id)
            print(f"{case_id}\tnot-started\tseed_t={seed_time:g}")
            continue
        latest = segments[-1]
        manifest = load_json(segment_manifest(latest))
        job_id = manifest.get("job_id")
        state = job_state(str(job_id)) if job_id is not None else "PREPARED"
        analysis_path = latest / "manifest/fast_analysis.json"
        if analysis_path.is_file():
            analysis = load_json(analysis_path)
            progress = f"t={float(analysis['final_time']):.9g} passed={analysis['passed']}"
        else:
            progress = f"start_t={float(manifest['start_time']):.9g}"
        print(f"{case_id}\t{state}\tjob={job_id}\t{latest.name}\t{progress}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    commands = parser.add_subparsers(dest="command", required=True)

    launch = commands.add_parser("launch", help="prepare and optionally submit cases")
    launch.add_argument("cases", nargs="+")
    launch.add_argument("--submit", action="store_true")

    advance = commands.add_parser("advance", help="analyze a segment and continue it")
    advance.add_argument("--segment", required=True, type=Path)
    advance.add_argument("--submit", action="store_true")

    analyze = commands.add_parser("analyze", help="analyze one segment")
    analyze.add_argument("--segment", required=True, type=Path)

    status = commands.add_parser("status", help="show current case states")
    status.add_argument("cases", nargs="*", default=["R03-R17"])

    args = parser.parse_args()
    root = args.root.resolve()
    try:
        if args.command == "launch":
            for case_id in expand_cases(args.cases):
                launch_case(root, case_id, submit=args.submit)
        elif args.command == "advance":
            next_from_segment(args.segment.resolve(), submit=args.submit)
        elif args.command == "analyze":
            print(json.dumps(analyze_segment(args.segment.resolve()), indent=2, sort_keys=True))
        elif args.command == "status":
            command_status(root, expand_cases(args.cases))
        return 0
    except (FastRunError, OSError, subprocess.SubprocessError, ValueError, KeyError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
