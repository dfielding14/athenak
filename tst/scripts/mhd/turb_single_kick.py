"""Exact impulse and per-cycle total-power checks; build tests/turb_forcing."""

import logging
import os
from pathlib import Path
import shlex
import subprocess

import numpy as np
import scripts.utils.athena as athena  # noqa: F401 (run_tests.py build entry)

logger = logging.getLogger('athena' + __name__[7:])
_RESULTS = []


def _input(fluids, integrator, direct, power):
    text = f"""<job>
basename = kick
<mesh>
nghost = 2
nx1 = 8
nx2 = 4
nx3 = 4
x1min = 0
x1max = 2
x2min = 0
x2max = 1
x3min = 0
x3max = 1
ix1_bc = periodic
ox1_bc = periodic
ix2_bc = periodic
ox2_bc = periodic
ix3_bc = periodic
ox3_bc = periodic
<meshblock>
nx1 = 4
nx2 = 4
nx3 = 4
<time>
evolution = dynamic
integrator = {integrator}
cfl_number = 0.2
nlim = {0 if direct else 4}
tlim = 1
ndiag = 1
<problem>
check_kick = {str(direct).lower()}
<turb_driving>
constant_edot = true
dedt = {power}
tcorr = 0.5
dt_turb_update = 0.01
nlow = 1
nhigh = 2
spect_form = 2
rseed = 2718
<output1>
file_type = hst
dcycle = 1
data_format = %.17e
"""
    for fluid, eos in fluids:
        text += f"""<{fluid}>
eos = {eos}
gamma = 1.6666666666666667
iso_sound_speed = 1
rsolver = llf
reconstruct = plm
nscalars = 1
"""
    if len(fluids) == 2:
        text += '<ion-neutral>\ndrag_coeff = 0\n'
    return text


def _case(exe_dir, label, fluids, integrator='rk2', direct=False,
          power=0.06, nproc=1):
    case = exe_dir / label
    case.mkdir(exist_ok=True)
    for history in case.glob('kick.*.hst'):
        history.unlink()
    path = case / 'athinput'
    path.write_text(_input(fluids, integrator, direct, power))
    command = [str(exe_dir / 'athena'), '-i', str(path)]
    if nproc > 1:
        command = shlex.split(os.environ.get('MPIEXEC', 'mpiexec')) + [
            '-n', str(nproc)] + command
    proc = subprocess.run(command, cwd=case, capture_output=True, text=True,
                          timeout=60)
    output = proc.stdout + proc.stderr
    (case / 'run.log').write_text(output)
    if proc.returncode:
        raise RuntimeError(label + ': ' + output)
    if direct:
        if 'max kick error=' not in output:
            raise RuntimeError('Direct impulse check did not execute: ' + label)
        _RESULTS.append((label, 0.0, 0.0))
        return

    histories = [np.atleast_2d(np.loadtxt(case / ('kick.' + fluid + '.hst')))
                 for fluid, _ in fluids]
    time = histories[0][:, 0]
    # Finalize can repeat the final history row; keep one row per time.
    _, indices = np.unique(time, return_index=True)
    histories = [h[indices] for h in histories]
    time = histories[0][:, 0]
    if len(time) != 5 or not all(np.isfinite(h).all() for h in histories):
        raise RuntimeError('Expected four finite completed cycles: ' + label)
    if not all(np.array_equal(h[:, 0], time) for h in histories):
        raise RuntimeError('Fluid history times differ: ' + label)
    energy = sum(h[:, 6] for h in histories)
    momentum = sum(h[:, 3:6] for h in histories)
    residual = np.diff(energy) - power*np.diff(time)
    _RESULTS.append((label, float(np.max(np.abs(residual))),
                     float(np.max(np.abs(momentum)))))


def run(**kwargs):
    _RESULTS.clear()
    exe_dir = Path.cwd() / 'build' / 'src'
    config = subprocess.check_output([str(exe_dir / 'athena'), '-c'], text=True)
    if 'Problem generator:          tests/turb_forcing' not in config:
        logger.info('Skipping: configure with -DPROBLEM=tests/turb_forcing')
        return
    for fluid in ('hydro', 'mhd'):
        for eos in ('ideal', 'isothermal'):
            _case(exe_dir, fluid + '_' + eos, [(fluid, eos)], direct=True)
        for integrator in ('rk1', 'rk2', 'rk3', 'rk4'):
            _case(exe_dir, fluid + '_' + integrator, [(fluid, 'ideal')],
                  integrator)
    for hydro_eos in ('ideal', 'isothermal'):
        for mhd_eos in ('ideal', 'isothermal'):
            _case(exe_dir, hydro_eos + '_' + mhd_eos,
                  [('hydro', hydro_eos), ('mhd', mhd_eos)], 'imex2', True)
    for integrator in ('imex2', 'imex+'):
        _case(exe_dir, integrator, [('hydro', 'ideal'), ('mhd', 'ideal')],
              integrator)
    _case(exe_dir, 'zero_power', [('mhd', 'ideal')], power=0.0)
    if 'MPI parallelism:            ON' in config:
        _case(exe_dir, 'mpi_power', [('mhd', 'ideal')], 'rk3', nproc=2)


def analyze():
    for label, work, momentum in _RESULTS:
        logger.info('%s max work residual=% .8e net momentum=% .8e',
                    label, work, momentum)
    return all(work < 2.0e-12 and momentum < 2.0e-12
               for _, work, momentum in _RESULTS)
