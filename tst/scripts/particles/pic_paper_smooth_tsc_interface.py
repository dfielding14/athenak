"""Deterministic static-SMR paper_smooth receiver-resolution TSC regression."""

import glob
import logging
import os
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena
from scripts.particles import pic_paper_smooth_tsc_oracle as oracle

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_paper_smooth_tsc_interface.athinput'
_INPUT_DECKS = {
    'j': 'tests/pic_paper_smooth_tsc_interface_3d.athinput',
}
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_PARTICLE_Y = -0.75
_DOMAIN_LENGTH_X = 8.0
_DOMAIN_LENGTH_Y = 4.0
_DOMAIN_LENGTH_Z = 4.0
_ATOL = 1.0e-6
_RESULTS = {}
_RUNTIME_PARTICLES = {
    label: (particle_x, _PARTICLE_Y, 0.0)
    for label, particle_x in oracle.PARTICLES.items()
}
_RUNTIME_PARTICLES['g'] = (-3.8, _PARTICLE_Y, 0.0)
_RUNTIME_PARTICLES['h'] = (1.0, -1.8, 0.0)
_RUNTIME_PARTICLES['i'] = (1.0, -1.8, 0.0)
_RUNTIME_PARTICLES['j'] = (3.8, _PARTICLE_Y, -1.8)
_RUNTIME_PARTICLES['k'] = (1.0, -1.8, 0.0)
_RUNTIME_PARTICLES['l'] = (1.0, -1.8, 0.0)
_RUNTIME_OVERRIDES = {
    'i': ['mesh/ix2_bc=reflect', 'mesh/ox2_bc=reflect'],
    'k': ['mesh/ix2_bc=outflow', 'mesh/ox2_bc=outflow'],
    'l': ['mesh/ix2_bc=inflow', 'mesh/ox2_bc=inflow'],
}
_PHYSICAL_BOUNDARY_CASES = {'i', 'k', 'l'}
_EXPECTED_TOTALS = dict(oracle.EXPECTED_RAW_TOTALS)
_EXPECTED_TOTALS['g'] = 1.14
_EXPECTED_TOTALS['h'] = 1.0
_EXPECTED_TOTALS['i'] = 1.0
_EXPECTED_TOTALS['j'] = 0.86
_EXPECTED_TOTALS['k'] = 1.0
_EXPECTED_TOTALS['l'] = 1.0
_NON_UNIT_TOTAL_PARTICLES = tuple(oracle.NON_UNIT_TOTAL_PARTICLES) + ('g',)


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(label):
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECKS.get(
        label, _INPUT_DECK)


def _athena_mpi_enabled():
    proc = subprocess.run(['./athena', '-c'], cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _latest_output_file(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           basename + '.prtcl_rho.*.bin')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _run_case(mode, label, nproc, require_split=False):
    basename = 'pic_paper_smooth_tsc_interface_' + mode + '_' + label
    args = [
        'job/basename=' + basename,
        'problem/particle_x=' + str(float(_RUNTIME_PARTICLES[label][0])),
        'problem/particle_y=' + str(float(_RUNTIME_PARTICLES[label][1])),
        'problem/particle_z=' + str(float(_RUNTIME_PARTICLES[label][2])),
    ]
    if require_split:
        args.append('problem/require_interface_mpi_split=true')
    args.extend(_RUNTIME_OVERRIDES.get(label, []))
    command = ['./athena', '-i', _athena_input_path(label)] + args
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    _remove_outputs(basename)
    logger.info('Executing %s_%s: %s', mode, label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + mode + '_' + label + '\n' + output)
    return _measure_case(basename, label)


def _measure_case(basename, label):
    data = bin_convert.read_binary(_latest_output_file(basename))
    particle_x, particle_y, particle_z = _RUNTIME_PARTICLES[label]
    collapsed = {}
    expected_collapsed = {}
    actual_cells = []
    expected_cells = []

    for geometry, rho in zip(data['mb_geometry'], data['mb_data']['prtcl_rho']):
        nx3, nx2, nx1 = rho.shape
        dx1 = (geometry[1] - geometry[0]) / nx1
        dx2 = (geometry[3] - geometry[2]) / nx2
        dx3 = (geometry[5] - geometry[4]) / nx3
        volume = dx1 * dx2 * dx3
        for k in range(nx3):
            zcenter = geometry[4] + (k + 0.5) * dx3
            wz = (1.0 if nx3 == 1 else max(
                oracle.raw_tsc_weight(particle_z + shift, zcenter, dx3)
                for shift in (-_DOMAIN_LENGTH_Z, 0.0, _DOMAIN_LENGTH_Z)))
            for j in range(nx2):
                ycenter = geometry[2] + (j + 0.5) * dx2
                yshifts = ((0.0,) if label in _PHYSICAL_BOUNDARY_CASES else
                           (-_DOMAIN_LENGTH_Y, 0.0, _DOMAIN_LENGTH_Y))
                wy = max(oracle.raw_tsc_weight(particle_y + shift, ycenter, dx2)
                         for shift in yshifts)
                if label in _PHYSICAL_BOUNDARY_CASES:
                    ny = int(round(_DOMAIN_LENGTH_Y / dx2))
                    ynorm = sum(
                        oracle.raw_tsc_weight(
                            particle_y, -0.5*_DOMAIN_LENGTH_Y + (n + 0.5)*dx2, dx2)
                        for n in range(ny))
                    wy /= ynorm
                for i in range(nx1):
                    xcenter = geometry[0] + (i + 0.5) * dx1
                    actual = float(rho[k, j, i] * volume)
                    wx = max(
                        oracle.raw_tsc_weight(particle_x + shift, xcenter, dx1)
                        for shift in (-_DOMAIN_LENGTH_X, 0.0, _DOMAIN_LENGTH_X))
                    expected = float(wx * wy * wz)
                    actual_cells.append(actual)
                    expected_cells.append(expected)
                    collapsed[xcenter] = collapsed.get(xcenter, 0.0) + actual
                    expected_collapsed[xcenter] = (
                        expected_collapsed.get(xcenter, 0.0) + expected)

    xcenters = sorted(collapsed)
    return {
        'actual_cells': np.asarray(actual_cells),
        'expected_cells': np.asarray(expected_cells),
        'actual_x': np.asarray([collapsed[x] for x in xcenters]),
        'expected_x': np.asarray([expected_collapsed[x] for x in xcenters]),
        'total': float(sum(collapsed.values())),
        'expected_total': float(_EXPECTED_TOTALS[label]),
    }


def _check_close(label, actual, expected):
    actual = np.asarray(actual)
    expected = np.asarray(expected)
    max_abs_err = float(np.max(np.abs(actual - expected)))
    logger.info('%s max_abs_err=% .8e atol=% .8e', label, max_abs_err, _ATOL)
    return bool(np.allclose(actual, expected, atol=_ATOL, rtol=0.0))


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    oracle.validate_oracle()
    _RESULTS.clear()
    for label in _RUNTIME_PARTICLES:
        _RESULTS['serial_' + label] = _run_case('serial', label, 1)

    if _athena_mpi_enabled():
        for label in _RUNTIME_PARTICLES:
            _RESULTS['mpi2_' + label] = _run_case(
                'mpi2', label, 2,
                require_split=(label in ('b', 'c', 'd', 'g')))
    else:
        logger.info('Skipping mpi2 cases: Athena build has MPI parallelism OFF')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True
    for case, result in _RESULTS.items():
        ok = _check_close(case + ':native_cells',
                          result['actual_cells'], result['expected_cells']) and ok
        ok = _check_close(case + ':collapsed_x',
                          result['actual_x'], result['expected_x']) and ok
        ok = _check_close(case + ':total',
                          result['total'], result['expected_total']) and ok

        label = case.rsplit('_', 1)[-1]
        if label in _NON_UNIT_TOTAL_PARTICLES:
            non_unit = abs(result['total'] - 1.0) > _ATOL
            logger.info('%s:non_unit_total measured=% .8e expected_non_unit=True',
                        case, result['total'])
            ok = non_unit and ok
    return ok
