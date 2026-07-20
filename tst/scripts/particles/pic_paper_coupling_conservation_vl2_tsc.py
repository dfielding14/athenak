import glob
import logging
import os
import re
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_paper_coupling_conservation_vl2_tsc.athinput'
_RESULTS = {}
_C = 3.0


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    for directory in ['bin', 'pvtk']:
        pattern = os.path.join(_athena_exe_dir(), directory, basename + '.*')
        for path in glob.glob(pattern):
            os.remove(path)
    history = os.path.join(_athena_exe_dir(), basename + '.mhd.hst')
    if os.path.isfile(history):
        os.remove(history)


def _run_case(basename, coeff, hall_mode='off'):
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'particles/couple_j_to_efield_coeff=' + str(coeff),
        'particles/pic_cr_hall_mode=' + hall_mode,
    ]
    if hall_mode == 'full':
        command.append('particles/pic_background_ion_q_over_mc=10.0')
    logger.info('Executing %s: %s', basename, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + basename + '\n' + output)
    induction = 'cr_hall_full' if hall_mode == 'full' else 'ideal_mhd_only'
    if 'induction=' + induction not in output:
        raise RuntimeError('Unexpected induction identity for ' + basename)


def _output_files(basename, directory, suffix):
    pattern = os.path.join(_athena_exe_dir(), directory,
                           basename + '.' + suffix)
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Expected initial and final outputs for ' + pattern)
    return paths


def _integrate(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _read_mhd_series(basename):
    result = {}
    for file_id, quantity in [
            ('mhd_u_m1', 'mom1'),
            ('mhd_u_m2', 'mom2'),
            ('mhd_u_m3', 'mom3'),
            ('mhd_u_e', 'ener')]:
        paths = _output_files(basename, 'bin', file_id + '.*.bin')
        initial = bin_convert.read_binary_as_athdf(paths[0])
        final = bin_convert.read_binary_as_athdf(paths[-1])
        result[quantity] = (_integrate(initial, quantity),
                            _integrate(final, quantity))

    paths = _output_files(basename, 'bin', 'mhd_bcc.*.bin')
    result['bcc_initial'] = bin_convert.read_binary_as_athdf(paths[0])
    result['bcc_final'] = bin_convert.read_binary_as_athdf(paths[-1])
    return result


def _find_marker(contents, pattern, offset, label):
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise RuntimeError('Could not find ' + label + ' in particle VTK output')
    return offset + match.end(), match


def _read_particle_state(path):
    with open(path, 'rb') as stream:
        contents = stream.read()

    offset, match = _find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS marker')
    npoint = int(match.group(1))
    offset += 4*3*npoint

    for name in ['gid', 'ptag', 'species']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' int\nLOOKUP_TABLE default\n')
        offset, _ = _find_marker(contents, pattern, offset, name + ' marker')
        offset += 4*npoint

    for name in ['deltaf_f0', 'deltaf_weight']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' float\nLOOKUP_TABLE default\n')
        offset, _ = _find_marker(contents, pattern, offset, name + ' marker')
        offset += 4*npoint

    offset, _ = _find_marker(contents, rb'\nVECTORS vel float\n',
                             offset, 'velocity marker')
    velocity = np.frombuffer(contents[offset:offset + 12*npoint],
                             dtype='>f4').astype(np.float64).reshape(npoint, 3)
    gamma = 1.0/np.sqrt(1.0 - np.sum(velocity*velocity, axis=1)/(_C*_C))
    momentum = gamma[:, None]*velocity
    return {
        'npoint': npoint,
        'momentum': np.sum(momentum, axis=0),
        'energy': float(np.sum((gamma - 1.0)*_C*_C)),
    }


def _measure_case(basename, hall_mode):
    vtk_paths = _output_files(basename, 'pvtk', 'prtcl_all.*.part.vtk')
    result = {
        'particles_initial': _read_particle_state(vtk_paths[0]),
        'particles_final': _read_particle_state(vtk_paths[-1]),
        'mhd': _read_mhd_series(basename),
    }
    if hall_mode == 'full':
        history_path = os.path.join(
            _athena_exe_dir(), basename + '.mhd.hst')
        result['hall_history'] = np.atleast_2d(np.loadtxt(history_path))[:, -2:]
    return result


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    cases = [('coeff0', 0.0, 'off'), ('coeff7', 7.0, 'off'),
             ('hall_full', 0.0, 'full')]
    for label, coeff, hall_mode in cases:
        basename = 'pic_paper_coupling_' + label
        _remove_outputs(basename)
        _run_case(basename, coeff, hall_mode)
        _RESULTS[label] = _measure_case(basename, hall_mode)


def _check_close(label, measured, expected, atol):
    error = float(np.max(np.abs(np.asarray(measured) - np.asarray(expected))))
    logger.info('%s max_abs_err=% .8e', label, error)
    return error <= atol


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True
    for label, result in _RESULTS.items():
        initial = result['particles_initial']
        final = result['particles_final']
        mhd = result['mhd']
        ok = initial['npoint'] == 64 and final['npoint'] == 64 and ok

        total_momentum_delta = final['momentum'] - initial['momentum']
        total_momentum_delta += np.array([
            mhd['mom1'][1] - mhd['mom1'][0],
            mhd['mom2'][1] - mhd['mom2'][0],
            mhd['mom3'][1] - mhd['mom3'][0],
        ])
        total_energy_delta = (final['energy'] - initial['energy'] +
                              mhd['ener'][1] - mhd['ener'][0])
        logger.info('%s total_momentum_delta=%s total_energy_delta=% .8e',
                    label, total_momentum_delta, total_energy_delta)
        ok = _check_close(label + ':momentum_conservation',
                          total_momentum_delta, np.zeros(3), 2.0e-6) and ok
        ok = abs(total_energy_delta) <= 3.0e-5 and ok
        ok = np.max(np.abs(final['momentum'] - initial['momentum'])) > 1.0e-3 and ok

    coeff0 = _RESULTS['coeff0']
    coeff7 = _RESULTS['coeff7']
    for quantity in ['bcc1', 'bcc2', 'bcc3']:
        ok = _check_close(
            'coefficient_invariant_induction:' + quantity,
            coeff0['mhd']['bcc_final'][quantity],
            coeff7['mhd']['bcc_final'][quantity],
            1.0e-12) and ok
    ok = _check_close('coefficient_invariant_particle_momentum',
                      coeff0['particles_final']['momentum'],
                      coeff7['particles_final']['momentum'], 1.0e-12) and ok
    hall = _RESULTS['hall_full']
    hall_history = hall['hall_history']
    ok = np.all(np.isfinite(hall_history)) and ok
    ok = hall_history[-1, 0] > 0.0 and hall_history[-1, 1] > 0.0 and ok
    hall_difference = np.max(np.abs(
        hall['particles_final']['momentum'] - coeff0['particles_final']['momentum']))
    logger.info('full-Hall particle-momentum difference=% .8e', hall_difference)
    ok = hall_difference > 1.0e-8 and ok
    return ok
