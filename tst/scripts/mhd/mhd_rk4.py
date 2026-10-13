"""Unforced MHD register and temporal-order checks; use built_in_pgens."""

import logging
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

import numpy as np
import scripts.utils.athena as athena  # noqa: F401 (run_tests.py build entry)


ATHENA = None
logger = logging.getLogger("athena" + __name__[7:])
_RESULTS = []


def input_text(integrator, *, advection=False, cfl=0.3, isothermal=False):
    """Uniform dynamic MHD, or a sine advected by the linear upwind operator."""
    mesh = "<mesh>\nnghost = 2\n"
    blocks = "<meshblock>\n"
    for axis in (1, 2, 3):
        cells = (16 if axis == 1 else 1) if advection else 8
        width = (8 if axis == 1 else 1) if advection else 4
        mesh += (f"nx{axis} = {cells}\nx{axis}min = 0\nx{axis}max = 1\n"
                 f"ix{axis}_bc = periodic\nox{axis}_bc = periodic\n")
        blocks += f"nx{axis} = {width}\n"
    text = ("<job>\nbasename = mhd_rk\n" + mesh + blocks +
            f"<time>\nintegrator = {integrator}\ncfl_number = {cfl}\n"
            f"evolution = {'kinematic' if advection else 'dynamic'}\n"
            f"nlim = {-1 if advection else 4}\ntlim = 0.5\nndiag = 10000\n"
            "<mhd>\nreconstruct = dc\nnscalars = 1\n"
            "gamma = 1.666666666666667\niso_sound_speed = 1\n")
    if advection or isothermal:
        text += (f"eos = isothermal\nrsolver = {'advect' if advection else 'hlle'}\n"
                 "<problem>\npgen_name = advection\nflow_dir = 1\n"
                 f"iproblem = 1\nvelocity = 1\namplitude = {0.1 if advection else 0}\n")
    else:
        text += ("eos = ideal\nrsolver = hlle\n"
                 "<problem>\npgen_name = shock_tube\nshock_dir = 1\nxshock = 0.5\n")
        for side in ("l", "r"):
            for name, value in (("d", 1.3), ("u", 0.2), ("v", -0.1),
                                ("w", 0.3), ("p", 1.1), ("bx", 0.6),
                                ("by", -0.4), ("bz", 0.8), ("y", 0.37)):
                text += f"{name}{side} = {value}\n"
    # Tables retain double precision; binary output on older branches uses float32.
    text += ("<output1>\nfile_type = tab\nvariable = mhd_u_bcc\nid = state\n"
             "slice_x2 = 0.1\nslice_x3 = 0.1\ndata_format = %24.16e\n"
             f"dcycle = {100000 if advection else 1}\n"
             "<output2>\nfile_type = hst\ndcycle = 1\ndata_format = %24.16e\n")
    return text


def run_case(directory, text, nproc=1):
    directory.mkdir(parents=True)
    path = directory / "case.athinput"
    path.write_text(text)
    command = [str(ATHENA), "-i", str(path), "-d", str(directory)]
    if nproc > 1:
        command = shlex.split(os.environ.get("MPIEXEC", "mpiexec")) + [
            "-n", str(nproc)] + command
    result = subprocess.run(command,
                            capture_output=True, text=True, timeout=60, check=False)
    (directory / "run.log").write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    return sorted((directory / "tab").glob("*.tab"))


def read_table(path):
    names = path.read_text().splitlines()[1].split()[1:]
    values = np.loadtxt(path, ndmin=2)
    assert values.shape[1] == len(names)
    assert np.all(np.isfinite(values))
    return dict(zip(names, values.T))


def uniform_state(tmp_path, integrator, isothermal, nproc=1):
    """A 3D uniform state preserves energy, all B components, and a nonzero scalar."""
    tables = run_case(tmp_path / "uniform",
                      input_text(integrator, isothermal=isothermal), nproc)
    initial = read_table(tables[0])
    expected = {"dens": 1.3, "mom1": 0.26, "mom2": -0.13, "mom3": 0.39,
                "ener": 1.1 / (5 / 3 - 1) + 0.5 * 1.3 * 0.14 + 0.5 * 1.16,
                "r_00": 1.3 * 0.37, "bcc1": 0.6, "bcc2": -0.4, "bcc3": 0.8}
    if isothermal:
        expected = {"dens": 1, "mom1": 1, "mom2": 1, "mom3": 1,
                    "r_00": 1, "bcc1": 0, "bcc2": 1, "bcc3": 1}
    for name, value in expected.items():
        np.testing.assert_allclose(initial[name], value, rtol=0, atol=2e-14)
    assert len(tables) >= 5  # initial state and every completed cycle
    for path in tables[1:]:
        current = read_table(path)
        for name in expected:
            np.testing.assert_allclose(current[name], initial[name], rtol=0, atol=3e-13)
    history = np.loadtxt(next((tmp_path / "uniform").glob("*.hst")), ndmin=2)
    # Full-volume integrals complement the high-precision spatial slice.
    np.testing.assert_allclose(history[:, 2:], np.broadcast_to(history[0, 2:],
                               history[:, 2:].shape), rtol=0, atol=3e-13)


def temporal_order(tmp_path, integrator, order, nproc=1):
    """Compare B and scalar advection to the exact semidiscrete upwind solution."""
    errors = []
    for cfl in (0.8, 0.4, 0.2):
        directory = tmp_path / str(cfl)
        tables = run_case(directory,
                          input_text(integrator, advection=True, cfl=cfl), nproc)
        final = read_table(tables[-1])
        history = np.loadtxt(next(directory.glob("*.hst")), ndmin=2)
        time = history[-1, 0]
        assert abs(time - 0.5) < 1e-14
        dx = 1 / 16
        eigenvalue = (np.exp(-2j * np.pi * dx) - 1) / dx
        exact = 1 + 0.1 * np.exp(eigenvalue.real * time) * np.sin(
            2 * np.pi * final["x1v"] + eigenvalue.imag * time)
        # Both the cell-centered scalar and staggered transverse field evolve.
        errors.append([np.sqrt(np.mean((final[name] - exact)**2))
                       for name in ("r_00", "bcc2", "bcc3")])
    errors = np.asarray(errors)
    rates = np.log2(errors[:-1] / errors[1:])
    assert np.all(np.isfinite(rates))
    assert np.all(np.abs(rates - order) < 0.3), (errors, rates)
    logger.info("%s nproc=%d errors=%s rates=%s", integrator, nproc, errors, rates)


def run(**kwargs):
    global ATHENA
    ATHENA = Path.cwd() / "build" / "src" / "athena"
    _RESULTS.clear()
    directory = Path(tempfile.mkdtemp(prefix="mhd_rk4_", dir=ATHENA.parent))
    logger.info("RK register/order outputs: %s", directory)
    for integrator in ("rk1", "rk2", "rk3", "rk4"):
        for isothermal in (False, True):
            label = integrator + ("_isothermal" if isothermal else "_ideal")
            uniform_state(directory / label, integrator, isothermal)
            _RESULTS.append(label)
        temporal_order(directory / (integrator + "_advection"), integrator,
                       int(integrator[-1]))
        _RESULTS.append(integrator + "_advection")
    config = subprocess.check_output([str(ATHENA), "-c"], text=True)
    if "MPI parallelism:            ON" in config:
        uniform_state(directory / "rk4_mpi", "rk4", False, nproc=2)
        temporal_order(directory / "rk4_advection_mpi", "rk4", 4, nproc=2)
        _RESULTS.extend(("rk4_mpi", "rk4_advection_mpi"))


def analyze():
    logger.info("Passed %d unforced MHD register/order cases: %s",
                len(_RESULTS), ", ".join(_RESULTS))
    return len(_RESULTS) >= 12
