"""Source-local paper-VL2 exact-conservation engineering self-consistency harness."""

from __future__ import annotations

from contextlib import contextmanager
import ctypes
from dataclasses import dataclass
import fcntl
import glob
import hashlib
import json
import logging
import os
from pathlib import Path
import re
import shlex
import stat
import struct
import subprocess
import sys
import tempfile
from typing import Iterator

import numpy as np

_SOURCE_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(_SOURCE_ROOT))
from tst.publication import q011_section54_exact_conservation_closure_successor_v1 as closure  # noqa: E402
from tst.publication import q011_section54_restart as restart_layout  # noqa: E402

logger = logging.getLogger("athena" + __name__[7:])

_LEDGER_DISABLED_PARITY_BASE_COMMIT = "1d72534619da277fdd6838cabdab0d465d0a8867"
_LEDGER_DISABLED_PARITY_BASE_TREE = "38c473644a42e01179a7b3140b2ad51e6345d667"
_LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE = "hip-mpi-release-paper-pic"
_LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION = (
    "930a04d1d39c873ea49abfcf500069011f6d5759240a8f5c5b3341a6d243b246"
)
_LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER_SHA256 = (
    "6053f190ed5bea093537ca5e6aef110212d54861f6726a6fa7294a1716eca2d5"
)
_LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER = (
    Path("/lustre/orion/ast207/proj-shared/dfielding/PIC/control_plane")
    / _LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION
    / "run_control_plane.py"
)
_LEDGER_DISABLED_PARITY_CONTROL_PLANE_PYTHON = (
    "/opt/cray/pe/python/3.11.7/bin/python3"
)
_LEDGER_DISABLED_PARITY_BASE_CANDIDATE_MANIFEST_ENV = (
    "ATHENA_PIC_EXACT_CONSERVATION_BASE_CANDIDATE_MANIFEST"
)
_LEDGER_DISABLED_PARITY_CANDIDATE_MANIFEST_ENV = (
    "ATHENA_PIC_EXACT_CONSERVATION_CANDIDATE_MANIFEST"
)
_LEDGER_DISABLED_PARITY_EXPECTED_CANDIDATE_COMMIT_ENV = (
    "ATHENA_PIC_EXACT_CONSERVATION_EXPECTED_CANDIDATE_COMMIT"
)
_DECK = (
    _SOURCE_ROOT
    / "inputs/tests/pic_parallel_shock_exact_conservation_closure_vl2_tsc.athinput"
)
_SOURCE_FILES = (
    "src/bvals/bvals_part.cpp",
    "src/particles/particles.cpp",
    "src/particles/particles.hpp",
    "src/particles/particles_pushers.cpp",
    "src/pgen/tests/pic_parallel_shock.cpp",
    "tst/publication/q011_section54_exact_conservation_closure_successor_v1.py",
    "tst/publication/q011_section54_restart.py",
    "tst/publication/q019_nonlinear_bell_particle_state.py",
)
_THREE_D_GEOMETRY_OVERRIDES = (
    "mesh/nx3=4",
    "mesh/x3max=4.0",
    "meshblock/nx3=4",
)
_LEDGER_DISABLED_NONPERIODIC_3D_OVERRIDES = _THREE_D_GEOMETRY_OVERRIDES + (
    "mesh/ix3_bc=outflow",
    "mesh/ox3_bc=outflow",
    "particles/pic_enable_2d3v=false",
    "particles/pic_boundary_conservation_ledger=false",
    "problem/ps_enable_conservation_ledger=false",
    "problem/user_hist=false",
)
_MAX_NORMALIZED_RESIDUAL = 2.0e-10
_STATE_ATOL = 2.0e-10
_RESULTS: dict[str, object] = {}
_RAW_RUN_RECORDS: dict[str, dict[str, object]] = {}
_INOTIFY_MUTATION_MASK = (
    0x00000002  # IN_MODIFY
    | 0x00000004  # IN_ATTRIB
    | 0x00000008  # IN_CLOSE_WRITE
    | 0x00000400  # IN_DELETE_SELF
    | 0x00000800  # IN_MOVE_SELF
)
_AT_EMPTY_PATH = 0x1000


@dataclass(frozen=True)
class _SealedPythonScript:
    source_path: Path
    descriptor: int
    sha256: str


def _read_regular_file_once(path: Path, *, label: str) -> bytes:
    resolved = path.resolve()
    descriptor = os.open(
        resolved,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise RuntimeError(f"{label} is not a regular file: {resolved}")
        chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            chunks.append(chunk)
        after = os.fstat(descriptor)
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
            raise RuntimeError(f"{label} changed while being captured")
        return b"".join(chunks)
    finally:
        os.close(descriptor)


@contextmanager
def _seal_python_script(path: Path, *, label: str) -> Iterator[_SealedPythonScript]:
    source_path = path.resolve()
    payload = _read_regular_file_once(source_path, label=label)
    if not hasattr(os, "memfd_create") or not hasattr(os, "MFD_ALLOW_SEALING"):
        raise RuntimeError("sealed executable staging requires Linux memfd sealing")
    descriptor = os.memfd_create(label, flags=os.MFD_ALLOW_SEALING)
    try:
        offset = 0
        while offset < len(payload):
            written = os.write(descriptor, payload[offset:])
            if written <= 0:
                raise RuntimeError(f"{label} sealed staging write made no progress")
            offset += written
        os.fchmod(descriptor, 0o500)
        seals = (
            fcntl.F_SEAL_WRITE
            | fcntl.F_SEAL_GROW
            | fcntl.F_SEAL_SHRINK
            | fcntl.F_SEAL_SEAL
        )
        fcntl.fcntl(descriptor, fcntl.F_ADD_SEALS, seals)
        if fcntl.fcntl(descriptor, fcntl.F_GET_SEALS) != seals:
            raise RuntimeError(f"{label} sealed staging is incomplete")
        if os.fstat(descriptor).st_size != len(payload):
            raise RuntimeError(f"{label} sealed staging size differs from captured bytes")
        os.lseek(descriptor, 0, os.SEEK_SET)
        staged_chunks = []
        while chunk := os.read(descriptor, 1024 * 1024):
            staged_chunks.append(chunk)
        if b"".join(staged_chunks) != payload:
            raise RuntimeError(f"{label} sealed staging bytes differ from captured bytes")
        yield _SealedPythonScript(
            source_path=source_path,
            descriptor=descriptor,
            sha256=_sha(payload),
        )
    finally:
        os.close(descriptor)


@dataclass(frozen=True)
class _BoundExecutable:
    source_path: Path
    descriptor: int
    watch_descriptor: int
    sha256: str
    byte_count: int
    identity: tuple[int, ...]


def _regular_identity(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _descriptor_payload(descriptor: int) -> bytes:
    size = os.fstat(descriptor).st_size
    chunks = []
    offset = 0
    while offset < size:
        chunk = os.pread(descriptor, min(1024 * 1024, size - offset), offset)
        if not chunk:
            raise RuntimeError("descriptor-bound executable read ended early")
        chunks.append(chunk)
        offset += len(chunk)
    return b"".join(chunks)


def _watch_descriptor(descriptor: int, *, label: str) -> int:
    libc = ctypes.CDLL(None, use_errno=True)
    watch_descriptor = libc.inotify_init1(os.O_NONBLOCK | os.O_CLOEXEC)
    if watch_descriptor < 0:
        error = ctypes.get_errno()
        raise RuntimeError(f"{label}: cannot initialize mutation watch: {os.strerror(error)}")
    result = libc.inotify_add_watch(
        watch_descriptor,
        f"/proc/self/fd/{descriptor}".encode("ascii"),
        _INOTIFY_MUTATION_MASK,
    )
    if result < 0:
        error = ctypes.get_errno()
        os.close(watch_descriptor)
        raise RuntimeError(f"{label}: cannot watch bound executable: {os.strerror(error)}")
    return watch_descriptor


def _require_no_mutation_events(descriptor: int, *, label: str) -> None:
    try:
        payload = os.read(descriptor, 1024 * 1024)
    except BlockingIOError:
        return
    if payload:
        raise RuntimeError(f"{label}: descriptor-bound executable changed while retained")


@contextmanager
def _bind_executable(path: Path, *, label: str) -> Iterator[_BoundExecutable]:
    source_path = path.resolve()
    descriptor = os.open(
        source_path,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
    )
    watch_descriptor: int | None = None
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode):
            raise RuntimeError(f"{label} is not a regular file: {source_path}")
        watch_descriptor = _watch_descriptor(descriptor, label=label)
        payload = _descriptor_payload(descriptor)
        after = os.fstat(descriptor)
        if _regular_identity(before) != _regular_identity(after):
            raise RuntimeError(f"{label} changed while being bound")
        _require_no_mutation_events(watch_descriptor, label=label)
        yield _BoundExecutable(
            source_path=source_path,
            descriptor=descriptor,
            watch_descriptor=watch_descriptor,
            sha256=_sha(payload),
            byte_count=len(payload),
            identity=_regular_identity(after),
        )
    finally:
        if watch_descriptor is not None:
            os.close(watch_descriptor)
        os.close(descriptor)


def _require_bound_executable_unchanged(
    executable: _BoundExecutable, *, label: str
) -> None:
    _require_no_mutation_events(executable.watch_descriptor, label=label)
    if _regular_identity(os.fstat(executable.descriptor)) != executable.identity:
        raise RuntimeError(f"{label}: descriptor-bound executable metadata changed")
    payload = _descriptor_payload(executable.descriptor)
    if len(payload) != executable.byte_count or _sha(payload) != executable.sha256:
        raise RuntimeError(f"{label}: descriptor-bound executable bytes changed")


def _bound_executables_have_identical_bytes(
    left: _BoundExecutable, right: _BoundExecutable
) -> bool:
    if left.byte_count != right.byte_count:
        return False
    offset = 0
    while offset < left.byte_count:
        size = min(1024 * 1024, left.byte_count - offset)
        left_chunk = os.pread(left.descriptor, size, offset)
        right_chunk = os.pread(right.descriptor, size, offset)
        if left_chunk != right_chunk or len(left_chunk) != size:
            return False
        offset += size
    return True


def _run_descriptor_bound_executable(
    *,
    executable: _BoundExecutable,
    arguments: list[str],
    cwd: Path,
) -> subprocess.CompletedProcess[str]:
    bootstrap = (
        "import ctypes,os,sys\n"
        "fd=int(sys.argv[1]); args=[item.encode() for item in sys.argv[2:]]\n"
        "env=[(key+'='+value).encode() for key,value in os.environ.items()]\n"
        "argv=(ctypes.c_char_p*(len(args)+1))(*args,None)\n"
        "envp=(ctypes.c_char_p*(len(env)+1))(*env,None)\n"
        "libc=ctypes.CDLL(None,use_errno=True)\n"
        "libc.execveat.argtypes=[ctypes.c_int,ctypes.c_char_p,"
        "ctypes.POINTER(ctypes.c_char_p),ctypes.POINTER(ctypes.c_char_p),ctypes.c_int]\n"
        "libc.execveat.restype=ctypes.c_int\n"
        f"libc.execveat(fd,b'',argv,envp,{_AT_EMPTY_PATH})\n"
        "error=ctypes.get_errno(); raise OSError(error,os.strerror(error))\n"
    )
    _require_bound_executable_unchanged(executable, label="pre-execution")
    completed = subprocess.run(
        [
            _LEDGER_DISABLED_PARITY_CONTROL_PLANE_PYTHON,
            "-I",
            "-c",
            bootstrap,
            str(executable.descriptor),
            *arguments,
        ],
        cwd=cwd,
        capture_output=True,
        text=True,
        pass_fds=(executable.descriptor,),
    )
    _require_bound_executable_unchanged(executable, label="post-execution")
    return completed


def _exe_dir() -> Path:
    return Path(
        os.environ.get(
            "ATHENA_PIC_EXACT_CONSERVATION_EXE_DIR",
            str(_SOURCE_ROOT / "build-q011-host/src"),
        )
    ).resolve()


def _mpi_exe_dir() -> Path:
    return Path(
        os.environ.get(
            "ATHENA_PIC_EXACT_CONSERVATION_MPI_EXE_DIR",
            str(_SOURCE_ROOT / "build-q011-mpi-host/src"),
        )
    ).resolve()


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _artifact(name: str, payload: bytes, lineage: dict[str, str]) -> dict[str, object]:
    return {
        "name": name,
        "payload": payload,
        "sha256": _sha(payload),
        "byte_count": len(payload),
        "lineage": dict(lineage),
    }


def _candidate_payload() -> bytes:
    manifest = {
        path: _sha((_SOURCE_ROOT / path).read_bytes())
        for path in _SOURCE_FILES
    }
    return (json.dumps(manifest, sort_keys=True, separators=(",", ":")) + "\n").encode(
        "ascii"
    )


def _candidate_commit() -> str:
    proc = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=_SOURCE_ROOT,
        capture_output=True,
        check=True,
        text=True,
    )
    return proc.stdout.strip()


def _mpi_launcher(nproc: int) -> list[str] | None:
    raw = os.environ.get("ATHENA_PIC_EXACT_CONSERVATION_MPI_LAUNCHER")
    if raw is None:
        return None
    command = shlex.split(raw)
    if not command:
        raise RuntimeError("ATHENA_PIC_EXACT_CONSERVATION_MPI_LAUNCHER is empty")
    replaced = False
    for index, value in enumerate(command):
        if "{nproc}" in value:
            command[index] = value.replace("{nproc}", str(nproc))
            replaced = True
    if not replaced:
        command += ["-n", str(nproc)]
    return command


def _remove_outputs(exe_dir: Path, basename: str) -> None:
    for path in glob.glob(str(exe_dir / f"{basename}.*")):
        if os.path.isfile(path):
            os.remove(path)
    for path in glob.glob(str(exe_dir / "rst" / f"{basename}.*")):
        if os.path.isfile(path):
            os.remove(path)


def _execute_athena(
    *,
    label: str,
    exe_dir: Path,
    basename: str,
    deck: Path | None = None,
    restart: Path | None = None,
    launcher: list[str] | None = None,
    overrides: tuple[str, ...] = (),
    bound_executable: _BoundExecutable | None = None,
) -> dict[str, object]:
    _remove_outputs(exe_dir, basename)
    if bound_executable is not None and launcher is not None:
        raise RuntimeError(
            label + ": descriptor-bound Athena execution does not accept an MPI launcher"
        )
    command = [
        "./athena" if bound_executable is None else str(bound_executable.source_path)
    ]
    if restart is None:
        if deck is None:
            raise RuntimeError(label + ": fresh execution requires a deck")
        command += ["-i", str(deck)]
    else:
        command += ["-r", os.path.relpath(restart, exe_dir)]
    command += ["job/basename=" + basename, *overrides]
    if launcher is not None:
        command = launcher + command
    logger.info("Executing %s: %s", label, " ".join(command))
    if bound_executable is None:
        proc = subprocess.run(command, cwd=exe_dir, capture_output=True, text=True)
    else:
        proc = _run_descriptor_bound_executable(
            executable=bound_executable,
            arguments=command,
            cwd=exe_dir,
        )
    stdout = proc.stdout or ""
    stderr = proc.stderr or ""
    record = {
        "argv": command,
        "cwd": str(exe_dir),
        "returncode": proc.returncode,
        "stdout": stdout,
        "stdout_sha256": _sha(stdout.encode("utf-8")),
        "stderr": stderr,
        "stderr_sha256": _sha(stderr.encode("utf-8")),
        "executable_sha256": (
            _sha((exe_dir / "athena").read_bytes())
            if bound_executable is None
            else bound_executable.sha256
        ),
        "deck_sha256": None if deck is None else _sha(deck.read_bytes()),
        "restart_sha256": None if restart is None else _sha(restart.read_bytes()),
    }
    _RAW_RUN_RECORDS[label] = record
    return record


def _run_athena(
    *,
    label: str,
    exe_dir: Path,
    basename: str,
    deck: Path | None = None,
    restart: Path | None = None,
    launcher: list[str] | None = None,
    overrides: tuple[str, ...] = (),
    bound_executable: _BoundExecutable | None = None,
) -> str:
    record = _execute_athena(
        label=label,
        exe_dir=exe_dir,
        basename=basename,
        deck=deck,
        restart=restart,
        launcher=launcher,
        overrides=overrides,
        bound_executable=bound_executable,
    )
    output = str(record["stdout"]) + str(record["stderr"])
    proc_returncode = int(record["returncode"])
    if proc_returncode != 0:
        raise RuntimeError(label + " failed\n" + output)
    return output


def _run_athena_expect_fail(
    *,
    label: str,
    exe_dir: Path,
    basename: str,
    reason: str,
    deck: Path | None = None,
    restart: Path | None = None,
    launcher: list[str] | None = None,
    overrides: tuple[str, ...] = (),
) -> dict[str, object]:
    record = _execute_athena(
        label=label,
        exe_dir=exe_dir,
        basename=basename,
        deck=deck,
        restart=restart,
        launcher=launcher,
        overrides=overrides,
    )
    output = str(record["stdout"]) + str(record["stderr"])
    if record["returncode"] == 0:
        raise RuntimeError(label + ": invalid exact-ledger input unexpectedly passed")
    if reason not in output:
        raise RuntimeError(label + ": missing fail-closed reason\n" + output)
    return {
        "returncode": record["returncode"],
        "reason": reason,
        "deck_sha256": record["deck_sha256"],
        "restart_sha256": record["restart_sha256"],
        "stdout_sha256": record["stdout_sha256"],
        "stderr_sha256": record["stderr_sha256"],
    }


def _history_payload(exe_dir: Path, basename: str, physics: str) -> bytes:
    path = exe_dir / f"{basename}.{physics}.hst"
    if not path.is_file():
        raise RuntimeError(f"{basename}: missing {physics} history")
    return path.read_bytes()


def _restart_paths(exe_dir: Path, basename: str) -> list[Path]:
    paths = [Path(path) for path in glob.glob(str(exe_dir / "rst" / f"{basename}.*.rst"))]
    if not paths:
        raise RuntimeError(basename + ": no restart checkpoints")
    return sorted(paths)


def _fnv1a64(path: Path) -> int:
    value = 14695981039346656037
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            for byte in chunk:
                value ^= byte
                value = (value * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return value


def _write_completion_marker(path: Path) -> None:
    path.with_name(path.name + ".complete").write_text(
        "ATHENAK_RESTART_COMPLETE_V1\n"
        f"size={path.stat().st_size}\n"
        f"fnv1a64={_fnv1a64(path):016x}\n",
        encoding="ascii",
    )


def _write_shared_restart_publication(path: Path, exe_dir: Path) -> None:
    _write_completion_marker(path)
    manifest_path = path.with_name(path.name + ".manifest")
    manifest = {
        "schema": "ATHENAK_RESTART_MANIFEST_V1",
        "members": [
            {
                "path": os.path.relpath(path, exe_dir),
                "size": path.stat().st_size,
                "fnv1a64": f"{_fnv1a64(path):016x}",
            }
        ],
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    _write_completion_marker(manifest_path)


def _restart_mesh_layout(payload: bytes, label: str) -> dict[str, object]:
    marker = b"<par_end>\n"
    header = payload.find(marker)
    if header < 0:
        raise RuntimeError(label + ": missing restart parameter terminator")
    header += len(marker)
    fixed_header_bytes = 2 * 4 + 9 * 8 + 2 * 19 * 4 + 2 * 8 + 2 * 4
    if header + fixed_header_bytes > len(payload):
        raise RuntimeError(label + ": truncated restart mesh header")
    nmb_total, root_level = struct.unpack_from("<ii", payload, header)
    original_nranks = struct.unpack_from("<i", payload, header + fixed_header_bytes - 4)[0]
    if nmb_total <= 0 or original_nranks <= 0:
        raise RuntimeError(label + ": invalid restart mesh counts")
    logical = header + fixed_header_bytes
    costs = logical + nmb_total * 4 * 4
    ranks = costs + nmb_total * 4
    metadata = ranks + nmb_total * 4 + original_nranks * 2 * 4
    if metadata > len(payload):
        raise RuntimeError(label + ": truncated restart mesh layout")
    levels = [
        struct.unpack_from("<i", payload, logical + gid * 4 * 4 + 3 * 4)[0]
        for gid in range(nmb_total)
    ]
    return {
        "nmb_total": nmb_total,
        "root_level": root_level,
        "logical_offset": logical,
        "cost_offset": costs,
        "metadata_offset": metadata,
        "levels": levels,
    }


def _publish_mutated_restart(
    *,
    source: Path,
    destination: Path,
    exe_dir: Path,
    mutate: object,
) -> Path:
    payload = bytearray(source.read_bytes())
    mutate(payload)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(payload)
    _write_shared_restart_publication(destination, exe_dir)
    return destination


def _run_restart_execution_negative_guards(
    *, exe_dir: Path, deck_payload: bytes, launcher: list[str] | None
) -> dict[str, object]:
    reason = "root_level=max_level, every reconstructed MeshBlock at root_level"
    results: dict[str, object] = {}
    with tempfile.TemporaryDirectory(prefix="q011-exact-guard-deck-", dir=exe_dir) as tmp:
        deck = Path(tmp) / _DECK.name
        deck.write_bytes(deck_payload)

        three_d_overrides = _THREE_D_GEOMETRY_OVERRIDES + (
            "time/nlim=1",
            "time/tlim=1.0e9",
            "output1/dt=1.0e9",
        )
        results["exact_ledger_rejects_non_2d3v"] = _run_athena_expect_fail(
            label="guard_exact_ledger_non_2d3v",
            exe_dir=exe_dir,
            basename="pic_parallel_shock_exact_guard_non_2d3v",
            reason="requires the 2D3V",
            deck=deck,
            launcher=launcher,
            overrides=three_d_overrides,
        )

        ledger_disabled_record = _execute_athena(
            label="guard_ledger_disabled_nonperiodic_3d",
            exe_dir=exe_dir,
            basename="pic_parallel_shock_exact_guard_ledger_disabled_3d",
            deck=deck,
            launcher=launcher,
            overrides=_LEDGER_DISABLED_NONPERIODIC_3D_OVERRIDES + (
                "time/nlim=1",
                "time/tlim=1.0e9",
                "output1/dt=1.0e9",
            ),
        )
        if ledger_disabled_record["returncode"] != 0:
            raise RuntimeError(
                "ledger-disabled nonperiodic 3D baseline unexpectedly failed\n"
                + str(ledger_disabled_record["stdout"])
                + str(ledger_disabled_record["stderr"])
            )
        results["ledger_disabled_nonperiodic_3d"] = {
            "returncode": ledger_disabled_record["returncode"],
            "deck_sha256": ledger_disabled_record["deck_sha256"],
            "stdout_sha256": ledger_disabled_record["stdout_sha256"],
            "stderr_sha256": ledger_disabled_record["stderr_sha256"],
        }

        uniform_basename = "pic_parallel_shock_exact_guard_uniform_seed"
        _run_athena(
            label="guard_uniform_seed",
            exe_dir=exe_dir,
            basename=uniform_basename,
            deck=deck,
            launcher=launcher,
            overrides=("time/nlim=1", "time/tlim=1.0e9", "output1/dt=1.0e-12"),
        )
        uniform_restart = _restart_paths(exe_dir, uniform_basename)[-1]
        nonunit_path = exe_dir / "rst" / "pic_parallel_shock_exact_guard_nonunit.00000.rst"

        def nonunit_cost(payload: bytearray) -> None:
            layout = _restart_mesh_layout(payload, "nonunit-cost restart")
            struct.pack_into("<f", payload, int(layout["cost_offset"]), 2.0)

        _publish_mutated_restart(
            source=uniform_restart,
            destination=nonunit_path,
            exe_dir=exe_dir,
            mutate=nonunit_cost,
        )
        results["nonunit_restored_cost"] = _run_athena_expect_fail(
            label="guard_nonunit_restored_cost",
            exe_dir=exe_dir,
            basename="pic_parallel_shock_exact_guard_nonunit_restart",
            reason=reason,
            restart=nonunit_path,
            launcher=launcher,
            overrides=("time/nlim=2",),
        )

        adaptive_basename = "pic_parallel_shock_exact_guard_adaptive_seed"
        _run_athena(
            label="guard_adaptive_seed",
            exe_dir=exe_dir,
            basename=adaptive_basename,
            deck=deck,
            launcher=launcher,
            overrides=(
                "time/nlim=2",
                "time/tlim=1.0e9",
                "output1/dt=1.0e-12",
                "mesh_refinement/refinement=adaptive",
                "mesh_refinement/num_levels=2",
                "problem/ps_enable_conservation_ledger=false",
                "problem/user_hist=false",
                "problem/ps_enable_curvature_amr=true",
                "problem/ps_refine_curv=-1.0",
                "problem/ps_derefine_curv=-2.0",
            ),
        )
        adaptive_restart = _restart_paths(exe_dir, adaptive_basename)[-1]
        adaptive_layout = _restart_mesh_layout(
            adaptive_restart.read_bytes(), "adaptive-seed restart"
        )
        if not any(
            level != adaptive_layout["root_level"] for level in adaptive_layout["levels"]
        ):
            raise RuntimeError("adaptive-seed restart did not retain refined MeshBlocks")
        results["restored_refined_or_adaptive_topology"] = _run_athena_expect_fail(
            label="guard_restored_refined_or_adaptive_topology",
            exe_dir=exe_dir,
            basename="pic_parallel_shock_exact_guard_adaptive_restart",
            reason=reason,
            restart=adaptive_restart,
            launcher=launcher,
            overrides=(
                "time/nlim=3",
                "mesh_refinement/refinement=none",
                "mesh_refinement/num_levels=1",
                "problem/ps_enable_conservation_ledger=true",
                "problem/ps_enable_curvature_amr=false",
            ),
        )
        results["adaptive_seed_mesh"] = {
            "nmb_total": adaptive_layout["nmb_total"],
            "root_level": adaptive_layout["root_level"],
            "levels": adaptive_layout["levels"],
        }

        adaptive_uniform_basename = "pic_parallel_shock_exact_guard_adaptive_uniform_seed"
        _run_athena(
            label="guard_adaptive_uniform_seed",
            exe_dir=exe_dir,
            basename=adaptive_uniform_basename,
            deck=deck,
            launcher=launcher,
            overrides=(
                "time/nlim=1",
                "time/tlim=1.0e9",
                "output1/dt=1.0e-12",
                "mesh_refinement/refinement=adaptive",
                "mesh_refinement/num_levels=2",
                "problem/ps_enable_conservation_ledger=false",
                "problem/user_hist=false",
                "problem/ps_enable_curvature_amr=false",
            ),
        )
        adaptive_uniform_restart = _restart_paths(exe_dir, adaptive_uniform_basename)[-1]
        adaptive_uniform_layout = _restart_mesh_layout(
            adaptive_uniform_restart.read_bytes(), "adaptive-uniform-seed restart"
        )
        if any(
            level != adaptive_uniform_layout["root_level"]
            for level in adaptive_uniform_layout["levels"]
        ):
            raise RuntimeError("adaptive-uniform-seed restart unexpectedly refined")
        results["restored_adaptive_metadata_without_refinement"] = (
            _run_athena_expect_fail(
                label="guard_restored_adaptive_metadata_without_refinement",
                exe_dir=exe_dir,
                basename="pic_parallel_shock_exact_guard_adaptive_uniform_restart",
                reason=reason,
                restart=adaptive_uniform_restart,
                launcher=launcher,
                overrides=(
                    "time/nlim=2",
                    "mesh_refinement/refinement=none",
                    "mesh_refinement/num_levels=1",
                    "problem/ps_enable_conservation_ledger=true",
                    "problem/ps_enable_curvature_amr=false",
                ),
            )
        )
        results["adaptive_uniform_seed_mesh"] = {
            "nmb_total": adaptive_uniform_layout["nmb_total"],
            "root_level": adaptive_uniform_layout["root_level"],
            "levels": adaptive_uniform_layout["levels"],
        }
    return results


def _committed_time(payload: bytes, label: str) -> float:
    parameters = restart_layout._problem_parameters(payload, label)
    try:
        value = float(parameters["ps_conservation_committed_time"])
    except (KeyError, ValueError) as error:
        raise RuntimeError(label + ": missing exact committed time") from error
    if not np.isfinite(value) or value < 0.0:
        raise RuntimeError(label + ": invalid exact committed time")
    return value


def _restart_payload_at(
    exe_dir: Path, basename: str, *, earliest_time: float
) -> tuple[Path, bytes]:
    candidates = []
    for path in _restart_paths(exe_dir, basename):
        payload = path.read_bytes()
        committed_time = _committed_time(payload, str(path))
        if committed_time >= earliest_time:
            candidates.append((committed_time, path, payload))
    if not candidates:
        raise RuntimeError(basename + ": no restart at requested time")
    _, path, payload = min(candidates, key=lambda value: value[0])
    return path, payload


def _reduce_case(
    *, exe_dir: Path, basename: str, deck_payload: bytes, executable: Path
) -> dict[str, object]:
    candidate = _candidate_payload()
    attempt = (
        json.dumps(
            {
                "attempt_id": basename,
                "basename": basename,
                "executable_sha256": _sha(executable.read_bytes()),
            },
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("ascii")
    lineage = {
        "attempt_id": basename,
        "attempt_sha256": _sha(attempt),
        "candidate_commit": _candidate_commit(),
        "candidate_sha256": _sha(candidate),
        "deck_sha256": _sha(deck_payload),
        "executable_sha256": _sha(executable.read_bytes()),
    }
    retained_restarts = []
    for path in _restart_paths(exe_dir, basename):
        payload = path.read_bytes()
        if _committed_time(payload, str(path)) >= 100.0:
            retained_restarts.append(_artifact(path.name, payload, lineage))
    return closure.reduce_exact_conservation_closure(
        lineage=lineage,
        candidate_artifact=_artifact("candidate-source-manifest.json", candidate, lineage),
        deck_artifact=_artifact(_DECK.name, deck_payload, lineage),
        attempt_artifact=_artifact("attempt.json", attempt, lineage),
        mhd_history_artifact=_artifact(
            basename + ".mhd.hst", _history_payload(exe_dir, basename, "mhd"), lineage
        ),
        user_history_artifact=_artifact(
            basename + ".user.hst", _history_payload(exe_dir, basename, "user"), lineage
        ),
        restart_artifacts=retained_restarts,
    )


def _terminal_summary(result: dict[str, object]) -> dict[str, object]:
    checkpoint = result["checkpoints"][-1]
    closure_result = checkpoint["conservation_closure"]
    terms = closure_result["cumulative_ledger_terms"]
    normalized = closure_result["absolute_normalized_residual"]["max_state_abs"]
    finite_normalized = [value for value in normalized.values() if value is not None]
    return {
        "time": checkpoint["time"],
        "cycle": checkpoint["cycle"],
        "absolute_residual": closure_result["absolute_residual"],
        "max_absolute_normalized_residual": max(finite_normalized),
        "mhd_boundary": terms["ps_cons_mhd_boundary"],
        "particle_reflection": terms["ps_cons_particle_reflect"],
        "particle_escape": terms["ps_cons_particle_escape"],
        "injected_cr": terms["injected_cr"],
        "gas_subtracted": terms["ps_cons_gas_subtracted"],
        "removed_cr": terms["removed_cr"],
        "post_history_mhd_delta":
            closure_result["mhd_state_change_after_history_before_restart"],
    }


def _assert_closure(case: str, summary: dict[str, object]) -> None:
    if summary["max_absolute_normalized_residual"] > _MAX_NORMALIZED_RESIDUAL:
        raise RuntimeError(case + ": combined closure residual exceeds tolerance")
    if not any(abs(value) > 0.0 for value in summary["mhd_boundary"].values()):
        raise RuntimeError(case + ": MHD boundary transport was not exercised")
    if summary["particle_reflection"]["momentum_x1"] == 0.0:
        raise RuntimeError(case + ": reflecting-wall impulse was not exercised")
    if not (
        summary["particle_escape"]["mass"] < 0.0
        and summary["particle_escape"]["energy"] < 0.0
    ):
        raise RuntimeError(case + ": physical CR escape was not exercised")
    for term in ("injected_cr", "gas_subtracted", "removed_cr"):
        if not (summary[term]["mass"] > 0.0 and summary[term]["energy"] > 0.0):
            raise RuntimeError(case + ": " + term + " was not exercised")


def _restart_state(payload: bytes, deck_payload: bytes, label: str) -> dict[str, object]:
    deck = closure._deck_physics(deck_payload)
    header, mhd = closure._restart_mhd_state(payload, deck, label)
    return {
        "header": header,
        "mhd": mhd,
        "cr": closure._particle_budget(payload, deck, label),
        "ledger": closure._conservation_ledger(payload, label),
    }


def _max_numeric_difference(left: object, right: object, label: str) -> float:
    if type(left) is dict and type(right) is dict:
        if set(left) != set(right):
            raise RuntimeError(label + ": restart object keys drifted")
        return max(
            (_max_numeric_difference(left[key], right[key], label + "/" + key)
             for key in left),
            default=0.0,
        )
    if type(left) is list and type(right) is list:
        if len(left) != len(right):
            raise RuntimeError(label + ": restart list length drifted")
        return max(
            (_max_numeric_difference(a, b, label) for a, b in zip(left, right)),
            default=0.0,
        )
    if type(left) is float and type(right) is float:
        return abs(left - right)
    if left != right:
        raise RuntimeError(label + ": restart discrete state drifted")
    return 0.0


def _compare_restart_continuation(
    *, full: bytes, restarted: bytes, deck_payload: bytes, label: str
) -> dict[str, object]:
    full_state = _restart_state(full, deck_payload, label + "/full")
    restarted_state = _restart_state(restarted, deck_payload, label + "/restarted")
    mhd_error = np.max(
        np.abs(np.asarray(full_state["mhd"]) - np.asarray(restarted_state["mhd"]))
    )
    cr_error = np.max(
        np.abs(np.asarray(full_state["cr"]) - np.asarray(restarted_state["cr"]))
    )
    ledger_error = _max_numeric_difference(
        full_state["ledger"], restarted_state["ledger"], label + "/ledger"
    )
    particle_state = _compare_particle_state(full, restarted, label + "/particles")
    if full_state["header"]["time"] != restarted_state["header"]["time"]:
        raise RuntimeError(label + ": restart continuation terminal time drifted")
    if full_state["header"]["cycle"] != restarted_state["header"]["cycle"]:
        raise RuntimeError(label + ": restart continuation terminal cycle drifted")
    if max(mhd_error, cr_error, ledger_error) > _STATE_ATOL:
        raise RuntimeError(label + ": restart continuation state drifted")
    return {
        "terminal_time": full_state["header"]["time"],
        "terminal_cycle": full_state["header"]["cycle"],
        "max_mhd_state_absolute_error": float(mhd_error),
        "max_cr_budget_absolute_error": float(cr_error),
        "max_ledger_absolute_error": ledger_error,
        "ledger_structure_exactly_equal": True,
        "particle_state": particle_state,
    }


def _run_exact_case(
    *,
    label: str,
    exe_dir: Path,
    deck_payload: bytes,
    launcher: list[str] | None,
) -> dict[str, object]:
    basename = "pic_parallel_shock_exact_" + label
    with tempfile.TemporaryDirectory(prefix="q011-exact-deck-", dir=exe_dir) as tmp:
        deck = Path(tmp) / _DECK.name
        deck.write_bytes(deck_payload)
        _run_athena(
            label=label + "_full",
            exe_dir=exe_dir,
            basename=basename,
            deck=deck,
            launcher=launcher,
        )
    reduced = _reduce_case(
        exe_dir=exe_dir,
        basename=basename,
        deck_payload=deck_payload,
        executable=exe_dir / "athena",
    )
    summary = _terminal_summary(reduced)
    _assert_closure(label, summary)

    restart_path, _ = _restart_payload_at(exe_dir, basename, earliest_time=100.0)
    restart_basename = basename + "_restart"
    _run_athena(
        label=label + "_restart",
        exe_dir=exe_dir,
        basename=restart_basename,
        restart=restart_path,
        launcher=launcher,
    )
    _, full_terminal = _restart_payload_at(exe_dir, basename, earliest_time=205.0)
    _, restarted_terminal = _restart_payload_at(
        exe_dir, restart_basename, earliest_time=205.0
    )
    return {
        "closure": summary,
        "restart_continuation": _compare_restart_continuation(
            full=full_terminal,
            restarted=restarted_terminal,
            deck_payload=deck_payload,
            label=label,
        ),
        "paper_vl2_boundary_stage_semantics":
            "nonzero_boundary_transport_closes_at_roundoff_with_cycle_local_stage2_overwrite",
    }


def _particle_payload(payload: bytes, label: str) -> tuple[bytes, bytes]:
    probe = restart_layout.probe_schema7_restart_payload(payload, source=label)
    real_bytes = probe.particle_count * probe.real_fields_per_particle * 8
    int_bytes = probe.particle_count * probe.integer_fields_per_particle * 4
    return (
        payload[probe.particle_real_offset:probe.particle_real_offset + real_bytes],
        payload[probe.particle_integer_offset:probe.particle_integer_offset + int_bytes],
    )


def _compare_particle_state(full: bytes, restarted: bytes, label: str) -> dict[str, object]:
    def arrays(payload: bytes, source: str) -> tuple[np.ndarray, np.ndarray]:
        probe = restart_layout.probe_schema7_restart_payload(payload, source=source)
        real = np.frombuffer(
            payload,
            dtype="<f8",
            count=probe.particle_count * probe.real_fields_per_particle,
            offset=probe.particle_real_offset,
        ).reshape(probe.particle_count, probe.real_fields_per_particle)
        integer = np.frombuffer(
            payload,
            dtype="<i4",
            count=probe.particle_count * probe.integer_fields_per_particle,
            offset=probe.particle_integer_offset,
        ).reshape(probe.particle_count, probe.integer_fields_per_particle)
        # Schema 7 integer columns are PGID, PTAG, species, source. PTAG is unique.
        if not np.all(np.isfinite(real)):
            raise RuntimeError(source + ": restart particle real state is nonfinite")
        if np.unique(integer[:, 1]).size != probe.particle_count:
            raise RuntimeError(source + ": restart particle tags are not unique")
        order = np.argsort(integer[:, 1], kind="stable")
        return real[order], integer[order]

    full_real, full_integer = arrays(full, label + "/full")
    restarted_real, restarted_integer = arrays(restarted, label + "/restarted")
    if full_real.shape != restarted_real.shape or full_integer.shape != restarted_integer.shape:
        raise RuntimeError(label + ": restart particle-state shape drifted")
    if not np.array_equal(full_integer, restarted_integer):
        raise RuntimeError(label + ": restart particle integer state drifted")
    real_error = float(np.max(np.abs(full_real - restarted_real), initial=0.0))
    if real_error > _STATE_ATOL:
        raise RuntimeError(label + ": restart particle real state drifted")
    return {
        "particle_count": int(full_real.shape[0]),
        "tag_sorted_integer_state_exactly_equal": True,
        "tag_sorted_real_state_max_absolute_error": real_error,
    }


def _run_sealed_python_script(
    *,
    sealed_script: _SealedPythonScript,
    original_path: Path,
    arguments: list[str],
    cwd: str,
    env: dict[str, str],
    timeout: int,
) -> subprocess.CompletedProcess[bytes]:
    if original_path.resolve() != sealed_script.source_path:
        raise RuntimeError("sealed Python script source path differs from execution filename")
    bootstrap = (
        "import os,sys\n"
        "fd=int(sys.argv[1]); filename=sys.argv[2]\n"
        "os.lseek(fd,0,os.SEEK_SET); chunks=[]\n"
        "while True:\n"
        " chunk=os.read(fd,1048576)\n"
        " if not chunk: break\n"
        " chunks.append(chunk)\n"
        "sys.argv=sys.argv[3:]\n"
        "scope={'__name__':'__main__','__file__':filename,'__package__':None,"
        "'__cached__':None,'__builtins__':__builtins__}\n"
        "exec(compile(b''.join(chunks),filename,'exec'),scope)\n"
    )
    command = [
        _LEDGER_DISABLED_PARITY_CONTROL_PLANE_PYTHON,
        "-I",
        "-c",
        bootstrap,
        str(sealed_script.descriptor),
        str(original_path),
        str(original_path),
        *arguments,
    ]
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        capture_output=True,
        check=False,
        timeout=timeout,
        pass_fds=(sealed_script.descriptor,),
    )


def _run_clean_candidate_revalidation(
    candidate_manifest: Path, manifest_sha256: str, expected_commit: str
) -> dict[str, object]:
    arguments = [
        "revalidate_clean_candidate.py",
        "--manifest",
        str(candidate_manifest),
        "--expected-manifest-sha256",
        manifest_sha256,
        "--expected-git-commit",
        expected_commit,
    ]
    with _seal_python_script(
        _LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER,
        label="installed-control-plane-runner",
    ) as sealed_runner:
        if (
            sealed_runner.sha256
            != _LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER_SHA256
        ):
            raise RuntimeError(
                "installed control-plane runner differs from pinned executable bytes"
            )
        completed = _run_sealed_python_script(
            sealed_script=sealed_runner,
            original_path=_LEDGER_DISABLED_PARITY_CONTROL_PLANE_RUNNER,
            arguments=arguments,
            cwd="/",
            env={"HOME": "/", "LANG": "C", "LC_ALL": "C", "PATH": "/usr/bin:/bin"},
            timeout=300,
        )
    if completed.returncode != 0 or completed.stderr:
        raise RuntimeError(
            "installed control-plane clean-candidate revalidation failed: "
            f"returncode={completed.returncode}"
        )
    try:
        report = json.loads(completed.stdout)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RuntimeError(
            "installed control-plane clean-candidate revalidation returned invalid JSON"
        ) from error
    if not isinstance(report, dict):
        raise RuntimeError(
            "installed control-plane clean-candidate revalidation returned no report object"
        )
    return report


def _authenticate_ledger_disabled_parity_base(
    *,
    candidate_executable: _BoundExecutable,
    base_executable: _BoundExecutable,
    base_manifest: Path,
) -> dict[str, object]:
    _require_bound_executable_unchanged(
        candidate_executable, label="parity candidate authentication"
    )
    _require_bound_executable_unchanged(
        base_executable, label="parity base authentication"
    )
    base_manifest = base_manifest.resolve()
    if candidate_executable.source_path == base_executable.source_path:
        raise RuntimeError(
            "ledger-disabled parity base and candidate executable paths are identical"
        )
    if (
        candidate_executable.sha256 == base_executable.sha256
        or _bound_executables_have_identical_bytes(
            candidate_executable, base_executable
        )
    ):
        raise RuntimeError(
            "ledger-disabled parity base and candidate executable bytes are identical"
        )
    if not base_manifest.is_file():
        raise RuntimeError(
            "ledger-disabled parity base clean-candidate manifest is not a regular file: "
            + str(base_manifest)
        )
    if base_executable.source_path != (base_manifest.parent / "athena").resolve():
        raise RuntimeError(
            "ledger-disabled parity baseline executable is not the clean-candidate "
            "manifest executable"
        )

    manifest_sha256 = _sha(
        _read_regular_file_once(
            base_manifest,
            label="base clean-candidate manifest",
        )
    )
    try:
        report = _run_clean_candidate_revalidation(
            base_manifest,
            manifest_sha256,
            _LEDGER_DISABLED_PARITY_BASE_COMMIT,
        )
    except Exception as error:
        raise RuntimeError(
            "ledger-disabled parity baseline clean-candidate build provenance "
            "authentication failed"
        ) from error
    _require_bound_executable_unchanged(
        candidate_executable, label="parity candidate authentication"
    )
    _require_bound_executable_unchanged(
        base_executable, label="parity base authentication"
    )
    if (
        report.get("status") != "passed"
        or report.get("source", {}).get("git_commit")
        != _LEDGER_DISABLED_PARITY_BASE_COMMIT
        or report.get("source", {}).get("git_tree")
        != _LEDGER_DISABLED_PARITY_BASE_TREE
        or report.get("build", {}).get("executable_sha256")
        != base_executable.sha256
        or report.get("build", {}).get("profile_id")
        != _LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE
        or report.get("current_control_plane_version")
        != _LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION
        or report.get("clean_candidate_manifest", {}).get("sha256")
        != manifest_sha256
        or Path(report.get("clean_candidate_manifest", {}).get("path", "")).resolve()
        != base_manifest
    ):
        raise RuntimeError(
            "ledger-disabled parity baseline authentication report does not match the "
            "pinned base source and production build profile or the bound executable"
        )
    return {
        "executable_source_path": str(base_executable.source_path),
        "source_commit": _LEDGER_DISABLED_PARITY_BASE_COMMIT,
        "source_tree": _LEDGER_DISABLED_PARITY_BASE_TREE,
        "clean_candidate_manifest_sha256": manifest_sha256,
        "build_profile_id": report["build"]["profile_id"],
        "build_receipt_control_plane_version": report["build"][
            "receipt_control_plane_version"
        ],
        "current_control_plane_version": _LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
        "executable_sha256": base_executable.sha256,
        "candidate_executable_sha256": candidate_executable.sha256,
        "execution_binding": "linux_execveat_descriptor_bound_exact_authenticated_bytes",
    }


def _authenticate_ledger_disabled_parity_candidate(
    *,
    candidate_executable: _BoundExecutable,
    candidate_manifest: Path,
    expected_commit: str,
) -> dict[str, object]:
    if re.fullmatch(r"[0-9a-f]{40}", expected_commit) is None:
        raise RuntimeError("ledger-disabled parity expected candidate commit is malformed")
    if expected_commit == _LEDGER_DISABLED_PARITY_BASE_COMMIT:
        raise RuntimeError(
            "ledger-disabled parity candidate commit is identical to the pinned base"
        )
    _require_bound_executable_unchanged(
        candidate_executable, label="parity candidate authentication"
    )
    candidate_manifest = candidate_manifest.resolve()
    if not candidate_manifest.is_file():
        raise RuntimeError(
            "ledger-disabled parity candidate clean-candidate manifest is not a "
            "regular file: " + str(candidate_manifest)
        )
    if candidate_executable.source_path != (candidate_manifest.parent / "athena").resolve():
        raise RuntimeError(
            "ledger-disabled parity candidate executable is not the clean-candidate "
            "manifest executable"
        )
    manifest_sha256 = _sha(
        _read_regular_file_once(
            candidate_manifest,
            label="candidate clean-candidate manifest",
        )
    )
    try:
        report = _run_clean_candidate_revalidation(
            candidate_manifest,
            manifest_sha256,
            expected_commit,
        )
    except Exception as error:
        raise RuntimeError(
            "ledger-disabled parity candidate clean-candidate build provenance "
            "authentication failed"
        ) from error
    _require_bound_executable_unchanged(
        candidate_executable, label="parity candidate authentication"
    )
    if (
        report.get("status") != "passed"
        or report.get("source", {}).get("git_commit") != expected_commit
        or report.get("build", {}).get("executable_sha256")
        != candidate_executable.sha256
        or report.get("build", {}).get("profile_id")
        != _LEDGER_DISABLED_PARITY_BASE_BUILD_PROFILE
        or report.get("current_control_plane_version")
        != _LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION
        or report.get("clean_candidate_manifest", {}).get("sha256")
        != manifest_sha256
        or Path(report.get("clean_candidate_manifest", {}).get("path", "")).resolve()
        != candidate_manifest
    ):
        raise RuntimeError(
            "ledger-disabled parity candidate authentication report does not match the "
            "explicit candidate commit, production build profile, or bound executable"
        )
    return {
        "executable_source_path": str(candidate_executable.source_path),
        "source_commit": expected_commit,
        "source_tree": report["source"]["git_tree"],
        "clean_candidate_manifest_sha256": manifest_sha256,
        "build_profile_id": report["build"]["profile_id"],
        "build_receipt_control_plane_version": report["build"][
            "receipt_control_plane_version"
        ],
        "current_control_plane_version": _LEDGER_DISABLED_PARITY_CONTROL_PLANE_VERSION,
        "executable_sha256": candidate_executable.sha256,
        "execution_binding": "linux_execveat_descriptor_bound_exact_authenticated_bytes",
    }


def _run_optional_ledger_disabled_parity(deck_payload: bytes) -> object:
    raw_base_manifest = os.environ.get(
        _LEDGER_DISABLED_PARITY_BASE_CANDIDATE_MANIFEST_ENV
    )
    raw_candidate_manifest = os.environ.get(
        _LEDGER_DISABLED_PARITY_CANDIDATE_MANIFEST_ENV
    )
    expected_candidate_commit = os.environ.get(
        _LEDGER_DISABLED_PARITY_EXPECTED_CANDIDATE_COMMIT_ENV
    )
    supplied = (raw_base_manifest, raw_candidate_manifest, expected_candidate_commit)
    if all(value is None for value in supplied):
        return "not_requested"
    if any(value is None for value in supplied):
        raise RuntimeError(
            "ledger-disabled parity requires base and candidate clean-candidate "
            "manifests plus an explicit expected candidate commit"
        )
    base_manifest = Path(str(raw_base_manifest)).resolve()
    candidate_manifest = Path(str(raw_candidate_manifest)).resolve()
    expected_candidate_commit = str(expected_candidate_commit)
    overrides = _LEDGER_DISABLED_NONPERIODIC_3D_OVERRIDES + (
        "time/tlim=5.0",
        "problem/ps_remove_birth_time_before=-1.0",
        "output1/dt=5.0",
        "output2/dt=5.0",
    )
    results = []
    with _bind_executable(
        candidate_manifest.parent / "athena",
        label="ledger-disabled-parity-candidate",
    ) as candidate_executable, _bind_executable(
        base_manifest.parent / "athena",
        label="ledger-disabled-parity-base",
    ) as base_executable:
        authentication = _authenticate_ledger_disabled_parity_base(
            candidate_executable=candidate_executable,
            base_executable=base_executable,
            base_manifest=base_manifest,
        )
        candidate_authentication = _authenticate_ledger_disabled_parity_candidate(
            candidate_executable=candidate_executable,
            candidate_manifest=candidate_manifest,
            expected_commit=expected_candidate_commit,
        )
        with tempfile.TemporaryDirectory(prefix="q011-disabled-parity-work-") as tmp:
            root = Path(tmp)
            for name, executable in (
                ("new", candidate_executable),
                ("base", base_executable),
            ):
                exe_dir = root / name
                exe_dir.mkdir()
                basename = "pic_parallel_shock_exact_disabled_" + name
                deck = exe_dir / _DECK.name
                deck.write_bytes(deck_payload)
                _run_athena(
                    label="ledger_disabled_" + name,
                    exe_dir=exe_dir,
                    basename=basename,
                    deck=deck,
                    overrides=overrides,
                    bound_executable=executable,
                )
                path = _restart_paths(exe_dir, basename)[-1]
                mhd_rows, _ = closure._parse_history(
                    _history_payload(exe_dir, basename, "mhd"),
                    closure._MHD_LABELS,
                    name,
                )
                results.append({
                    "executable_sha256": executable.sha256,
                    "mhd": mhd_rows[-1],
                    "particles": _particle_payload(path.read_bytes(), name),
                })
    mhd_equal = results[0]["mhd"] == results[1]["mhd"]
    particle_equal = results[0]["particles"] == results[1]["particles"]
    if not (mhd_equal and particle_equal):
        raise RuntimeError("ledger-disabled nonperiodic 3D new/base physics state differs")
    return {
        "case": "ledger_disabled_nonperiodic_3d",
        "new_executable_sha256": results[0]["executable_sha256"],
        "base_executable_sha256": results[1]["executable_sha256"],
        "base_authentication": {
            key: value for key, value in authentication.items()
        },
        "candidate_authentication": {
            key: value for key, value in candidate_authentication.items()
        },
        "mhd_history_state_exactly_equal": True,
        "particle_payload_exactly_equal": True,
    }


def _summary() -> dict[str, object]:
    return {
        "evidence_class": "bounded_source_local_engineering_self_consistency_only",
        "trusted_execution_provenance": "absent",
        "registered_science_evidence": False,
        "qualification_effect": "none",
        "authority": {
            "execution_authorized": False,
            "launch_authorized": False,
            "policy_mutation_authorized": False,
            "scientific_claim_authorized": False,
            "publication_authorized": False,
        },
        "raw_run_records": {
            label: {
                key: value
                for key, value in record.items()
                if key not in {"stdout", "stderr"}
            }
            for label, record in _RAW_RUN_RECORDS.items()
        },
        **_RESULTS,
    }


def run(**kwargs) -> None:
    del kwargs
    _RESULTS.clear()
    _RAW_RUN_RECORDS.clear()
    deck_payload = _DECK.read_bytes()
    _RESULTS["restart_execution_negative_guards"] = (
        _run_restart_execution_negative_guards(
            exe_dir=_exe_dir(), deck_payload=deck_payload, launcher=None
        )
    )
    _RESULTS["serial_inflow"] = _run_exact_case(
        label="serial_inflow",
        exe_dir=_exe_dir(),
        deck_payload=deck_payload,
        launcher=None,
    )
    outflow = deck_payload.replace(b"ox1_bc    = inflow", b"ox1_bc    = outflow")
    if outflow == deck_payload:
        raise RuntimeError("outflow deck mutation did not apply")
    _RESULTS["serial_outflow"] = _run_exact_case(
        label="serial_outflow",
        exe_dir=_exe_dir(),
        deck_payload=outflow,
        launcher=None,
    )
    launcher = _mpi_launcher(2)
    if launcher is None:
        _RESULTS["mpi2_inflow"] = "not_requested"
    else:
        _RESULTS["mpi2_inflow"] = _run_exact_case(
            label="mpi2_inflow",
            exe_dir=_mpi_exe_dir(),
            deck_payload=deck_payload,
            launcher=launcher,
        )
    _RESULTS["ledger_disabled_parity"] = _run_optional_ledger_disabled_parity(
        deck_payload
    )


def analyze() -> bool:
    summary = _summary()
    logger.info("PIC exact conservation closure runtime metrics: %s", summary)
    return (
        "serial_inflow" in _RESULTS
        and "serial_outflow" in _RESULTS
        and _RESULTS.get("mpi2_inflow") != "not_requested"
        and _RESULTS.get("ledger_disabled_parity") != "not_requested"
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    run()
    print(json.dumps(_summary(), indent=2, sort_keys=True))
    if not analyze():
        raise SystemExit(1)
