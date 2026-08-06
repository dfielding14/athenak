"""Serial regressions for tracked-particle output payload construction."""

from pathlib import Path
import subprocess

import numpy as np
import pytest


def _require_part_random_serial_build():
    config = Path.cwd().parent / "config.hpp"
    if not config.is_file():
        pytest.skip("run from a configured CMake build tree's src directory")
    text = config.read_text()
    if '#define PROBLEM_GENERATOR "part_random"' not in text:
        pytest.skip("tracked-particle tests require a PROBLEM=part_random build")
    if "#define MPI_PARALLEL_ENABLED 0" not in text:
        pytest.skip("serial tracked-particle tests require an MPI-disabled build")


def _input_text(
    ntrack: int,
    *,
    particles: bool = True,
    variable: bool = True,
    basename: str = "tracked",
) -> str:
    particle_block = """
<particles>
particle_type = cosmic_ray
ppc = 0.0625
pusher = drift
assign_tag = index_order
""" if particles else ""
    hydro_block = "" if particles else """
<hydro>
eos = ideal
reconstruct = plm
rsolver = llf
gamma = 1.4

<problem>
pgen_name = shock_tube
shock_dir = 1
xshock = 0.0
dl = 1.0
pl = 1.0
ul = 0.0
vl = 0.0
wl = 0.0
dr = 0.125
pr = 0.1
ur = 0.0
vr = 0.0
wr = 0.0
"""
    variable_line = "variable = prtcl_all\n" if variable else ""
    return f"""
<job>
basename = {basename}

<mesh>
nghost = 2
nx1 = 8
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
{particle_block}
{hydro_block}
<output1>
file_type = trk
{variable_line}nparticles = {ntrack}
dt = 1.0
"""


def _run_command(tmp_path: Path, text: str):
    input_file = tmp_path / "tracked.athinput"
    input_file.write_text(text)
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    proc = subprocess.run(
        ["./athena", "-i", str(input_file), "-d", str(run_dir)],
        capture_output=True,
        text=True,
        timeout=30,
    )
    return run_dir, proc


def _run(tmp_path: Path, text: str):
    _require_part_random_serial_build()
    return _run_command(tmp_path, text)


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


def test_subset_tracking_writes_distinct_dense_native_endian_records(tmp_path):
    run_dir, proc = _run(tmp_path, _input_text(3))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    payloads = _read_payloads(run_dir / "trk" / "tracked.trk", 3)
    for payload in payloads:
        assert np.isfinite(payload).all()
        assert len({tuple(record) for record in payload}) == 3


def test_tracked_output_does_not_require_an_irrelevant_variable(tmp_path):
    run_dir, proc = _run(tmp_path, _input_text(3, variable=False))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert _read_payloads(run_dir / "trk" / "tracked.trk", 3)[0].shape == (3, 6)


@pytest.mark.parametrize("ntrack", (-1, 9))
def test_invalid_tracked_particle_count_is_rejected(tmp_path, ntrack):
    _, proc = _run(tmp_path, _input_text(ntrack))
    assert proc.returncode != 0
    assert "requires nparticles between 0 and the total particle count (8)" in (
        proc.stdout + proc.stderr
    )


def test_tracked_output_without_particle_module_is_rejected(tmp_path):
    _, proc = _run_command(tmp_path, _input_text(0, particles=False))
    assert proc.returncode != 0
    assert (
        "requires <particles> block" in proc.stdout + proc.stderr
        or "requires a <particles> block" in proc.stdout + proc.stderr
    )
