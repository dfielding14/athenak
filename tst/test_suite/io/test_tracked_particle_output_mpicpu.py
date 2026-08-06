"""MPI regressions for collective tracked-particle payload construction."""

import os
from pathlib import Path
import subprocess

import numpy as np
import pytest


def _require_part_random_mpi_build():
    config = Path.cwd().parent / "config.hpp"
    if not config.is_file():
        pytest.skip("run from a configured CMake build tree's src directory")
    text = config.read_text()
    if '#define PROBLEM_GENERATOR "part_random"' not in text:
        pytest.skip("tracked-particle tests require a PROBLEM=part_random build")
    if "#define MPI_PARALLEL_ENABLED 1" not in text:
        pytest.skip("MPI tracked-particle tests require an MPI-enabled build")


def _input_text(ntrack: int, *, nx1: int = 8, assign_tag: str = "index_order") -> str:
    return f"""
<job>
basename = tracked

<mesh>
nghost = 2
nx1 = {nx1}
x1min = -0.5
x1max = 0.5
ix1_bc = periodic
ox1_bc = periodic
nx2 = 4
x2min = -0.5
x2max = 0.5
ix2_bc = periodic
ox2_bc = periodic
nx3 = 4
x3min = -0.5
x3max = 0.5
ix3_bc = periodic
ox3_bc = periodic

<meshblock>
nx1 = 4
nx2 = 4
nx3 = 4

<mesh_refinement>
refinement = none

<time>
evolution = dynamic
integrator = rk2
cfl_number = 0.4
nlim = 0
tlim = 0.0
ndiag = 1

<particles>
particle_type = cosmic_ray
ppc = 0.0625
pusher = drift
assign_tag = {assign_tag}

<output1>
file_type = trk
variable = prtcl_all
nparticles = {ntrack}
dt = 1.0
"""


def _run(tmp_path: Path, name: str, text: str, *, tiny_chunks: bool = False):
    _require_part_random_mpi_build()
    input_file = tmp_path / f"{name}.athinput"
    input_file.write_text(text)
    run_dir = tmp_path / name
    run_dir.mkdir()
    env = os.environ.copy()
    if tiny_chunks:
        env["ATHENAK_TEST_MAX_MPI_BYTES"] = "7"
    proc = subprocess.run(
        ["mpirun", "-np", "2", "./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
        timeout=30,
        env=env,
    )
    return run_dir, proc


def _read_payloads(path: Path, ntrack: int):
    content = path.read_bytes()
    marker = f"ntracked_prtcls={ntrack}\n \n".encode()
    frame_prefix = b"\n# AthenaK tracked particle data"
    payloads = []
    search_begin = 0
    while True:
        marker_begin = content.find(marker, search_begin)
        if marker_begin < 0:
            break
        payload_begin = marker_begin + len(marker)
        payload_end = content.find(frame_prefix, payload_begin)
        if payload_end < 0:
            payload_end = len(content)
        payload = content[payload_begin:payload_end]
        assert len(payload) == 6 * ntrack * np.dtype("=f4").itemsize
        payloads.append(np.frombuffer(payload, dtype="=f4").reshape(ntrack, 6))
        search_begin = payload_end
    assert payloads
    return payloads


def test_collective_and_independent_writes_match_forced_tiny_chunks(tmp_path):
    text = _input_text(5)
    normal_dir, normal = _run(tmp_path, "normal", text)
    chunked_dir, chunked = _run(tmp_path, "chunked", text, tiny_chunks=True)
    assert normal.returncode == 0, normal.stdout + normal.stderr
    assert chunked.returncode == 0, chunked.stdout + chunked.stderr
    normal_payloads = _read_payloads(normal_dir / "trk" / "tracked.trk", 5)
    chunked_payloads = _read_payloads(chunked_dir / "trk" / "tracked.trk", 5)
    assert len(chunked_payloads) == len(normal_payloads)
    for normal_payload, chunked_payload in zip(normal_payloads, chunked_payloads):
        np.testing.assert_array_equal(chunked_payload, normal_payload)
        assert len({tuple(record) for record in normal_payload}) == 5


def test_zero_owning_rank_participates_without_deadlock(tmp_path):
    run_dir, proc = _run(tmp_path, "zero_owner", _input_text(4), tiny_chunks=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _read_payloads(run_dir / "trk" / "tracked.trk", 4)[0].shape == (4, 6)


def test_rank_order_gap_is_rejected_for_uneven_rank_populations(tmp_path):
    _, proc = _run(tmp_path, "gap", _input_text(12, nx1=12, assign_tag="rank_order"))
    assert proc.returncode != 0
    assert "must form the dense range [0, nparticles)" in (proc.stdout + proc.stderr)
