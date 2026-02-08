import glob
import logging
import os
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_mhd_passive_mode.athinput'
_COUPLED_DECK = 'tests/pic_mhd_current_coupling.athinput'
_NO_MHD_DECK = 'tests/pic_deposit_conservation.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_POSITIVE_RESULTS = {}
_NEGATIVE_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(input_deck):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_deck


def _athena_mpi_enabled():
    command = ['./athena', '-c']
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _deck_path_for_python(input_deck):
    return os.path.join('..', 'inputs', input_deck)


def _parse_override_value(raw_value):
    text = raw_value.strip()
    lower = text.lower()
    if lower == 'true':
        return True
    if lower == 'false':
        return False
    try:
        return int(text)
    except ValueError:
        pass
    try:
        return float(text)
    except ValueError:
        return text


def _apply_overrides(config, arguments):
    for arg in arguments:
        if '=' not in arg:
            continue
        lhs, rhs = arg.split('=', 1)
        if '/' not in lhs:
            continue
        block, param = lhs.split('/', 1)
        if block not in config:
            config[block] = {}
        config[block][param] = _parse_override_value(rhs)


def _expected_totals(arguments):
    config = bin_convert.athinput(_deck_path_for_python(_INPUT_DECK))
    _apply_overrides(config, arguments)

    mesh = config['mesh']
    particles = config['particles']
    species0 = config['species0']

    nx1 = int(mesh['nx1'])
    nx2 = int(mesh['nx2'])
    nx3 = int(mesh['nx3'])
    ppc = float(particles['ppc'])
    qscale = float(particles['deposit_qscale'])
    charge = float(species0['charge'])
    vx0 = float(particles['cr_vx0'])
    vy0 = float(particles['cr_vy0'])
    vz0 = float(particles['cr_vz0'])

    npart = int(ppc * nx1 * nx2 * nx3)
    return {
        'npart': float(npart),
        'Q': float(npart) * qscale * charge,
        'Jx': float(npart) * qscale * charge * vx0,
        'Jy': float(npart) * qscale * charge * vy0,
        'Jz': float(npart) * qscale * charge * vz0,
    }


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _latest_output_file(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _output_files(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) < 2:
        raise RuntimeError('Need at least two output snapshots for: ' + pattern)
    return matches


def _volume_weights(dataset):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    return dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]


def _integrate_quantity(dataset, quantity):
    dvol = _volume_weights(dataset)
    return float(np.sum(dataset[quantity] * dvol))


def _l2_quantity(dataset, quantity):
    dvol = _volume_weights(dataset)
    return float(np.sqrt(np.sum(dataset[quantity] * dataset[quantity] * dvol)))


def _l2_difference(first_dataset, last_dataset, quantity):
    dvol = _volume_weights(first_dataset)
    diff = last_dataset[quantity] - first_dataset[quantity]
    return float(np.sqrt(np.sum(diff * diff * dvol)))


def _measure_particle_totals(basename):
    rho_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_rho'))
    jx_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jx'))
    jy_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jy'))
    jz_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jz'))
    pdens_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_d'))

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'npart': float(np.sum(pdens_data['pdens'])),
    }


def _measure_bcc_drift(basename):
    files = _output_files(basename, 'mhd_bcc')
    first_data = bin_convert.read_binary_as_athdf(files[0])
    last_data = bin_convert.read_binary_as_athdf(files[-1])

    return {
        'bcc1_diff': _l2_difference(first_data, last_data, 'bcc1'),
        'bcc1_ref': _l2_quantity(first_data, 'bcc1'),
        'bcc2_diff': _l2_difference(first_data, last_data, 'bcc2'),
        'bcc2_ref': _l2_quantity(first_data, 'bcc2'),
        'bcc3_diff': _l2_difference(first_data, last_data, 'bcc3'),
        'bcc3_ref': _l2_quantity(first_data, 'bcc3'),
    }


def _measure_scalar_drift(basename, file_id, quantity, label):
    files = _output_files(basename, file_id)
    first_data = bin_convert.read_binary_as_athdf(files[0])
    last_data = bin_convert.read_binary_as_athdf(files[-1])

    return {
        f'{label}_diff': _l2_difference(first_data, last_data, quantity),
        f'{label}_ref': _l2_quantity(first_data, quantity),
    }


def _measure_case(basename):
    measured = {}
    measured.update(_measure_particle_totals(basename))
    measured.update(_measure_bcc_drift(basename))
    measured.update(_measure_scalar_drift(basename, 'mhd_u_m1', 'mom1', 'mom1'))
    measured.update(_measure_scalar_drift(basename, 'mhd_u_m2', 'mom2', 'mom2'))
    measured.update(_measure_scalar_drift(basename, 'mhd_u_m3', 'mom3', 'mom3'))
    measured.update(_measure_scalar_drift(basename, 'mhd_u_e', 'ener', 'ener'))
    return measured


def _run_command(label, nproc, arguments, expect_fail=False,
                 expected_message=None, input_deck=_INPUT_DECK):
    command = ['./athena', '-i', _athena_input_path(input_deck)] + list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')

    if expect_fail:
        if proc.returncode == 0:
            raise RuntimeError('Expected failure for ' + label + ', but run succeeded')
        if expected_message is not None and expected_message not in output:
            raise RuntimeError('Expected message not found for ' + label +
                               ': "' + expected_message + '"')
        _NEGATIVE_RESULTS[label] = {'matched_message': expected_message}
        logger.info('Expected failure observed for %s', label)
        return

    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + ': ' + ' '.join(command))


def _check_with_tolerance(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def _check_drift(label, diff, ref, abs_tol, rel_tol):
    threshold = max(abs_tol, rel_tol * max(ref, 1.0))
    logger.info('%s diff=% .8e ref=% .8e threshold=% .8e',
                label, diff, ref, threshold)
    return diff <= threshold


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    positive_cases = [
        {
            'name': 'serial',
            'basename': 'pic_mhd_passive_serial',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_passive_serial'],
        },
    ]
    if _athena_mpi_enabled():
        positive_cases.append({
            'name': 'mpi2',
            'basename': 'pic_mhd_passive_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_mhd_passive_mpi2'],
        })
    else:
        logger.info('Skipping mpi2 case: Athena build has MPI parallelism OFF')

    for case in positive_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _POSITIVE_RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(case['args']),
        }

    _run_command(
        'guard_passive_requires_mhd',
        1,
        ['particles/pic_background_mode=passive_mhd',
         'particles/pic_feedback_mode=test_particle',
         'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_background_mode=passive_mhd requires an active <mhd> block',
        input_deck=_NO_MHD_DECK,
    )
    _run_command(
        'guard_passive_requires_test_particle',
        1,
        ['particles/pic_background_mode=passive_mhd', 'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_background_mode=passive_mhd requires '
                         '<particles>/pic_feedback_mode=test_particle',
        input_deck=_COUPLED_DECK,
    )
    _run_command(
        'guard_test_particle_rejects_coupling_toggles',
        1,
        ['particles/couple_moments_to_mhd=true', 'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_feedback_mode=test_particle does not support '
                         'particle-to-MHD coupling toggles',
        input_deck=_INPUT_DECK,
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    drift_keys = ['bcc1', 'bcc2', 'bcc3', 'mom1', 'mom2', 'mom3', 'ener']

    for case_name in ['serial', 'mpi2']:
        if case_name not in _POSITIVE_RESULTS:
            continue
        result = _POSITIVE_RESULTS[case_name]
        measured = result['measured']
        expected = result['expected']

        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'], expected['Jx'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'], expected['Jy'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'], expected['Jz'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance(case_name + ':npart', measured['npart'],
                                   expected['npart'], 1.0e-6, 1.0e-10) and ok

        for key in drift_keys:
            diff_key = key + '_diff'
            ref_key = key + '_ref'
            ok = _check_drift(case_name + ':' + key,
                              measured[diff_key], measured[ref_key],
                              1.0e-8, 1.0e-10) and ok

    if 'mpi2' in _POSITIVE_RESULTS:
        serial = _POSITIVE_RESULTS['serial']['measured']
        mpi2 = _POSITIVE_RESULTS['mpi2']['measured']
        ok = _check_with_tolerance('serial_vs_mpi2:Q', mpi2['Q'], serial['Q'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jx', mpi2['Jx'], serial['Jx'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jy', mpi2['Jy'], serial['Jy'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jz', mpi2['Jz'], serial['Jz'],
                                   1.0e-6, 1.0e-10) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:npart',
                                   mpi2['npart'], serial['npart'],
                                   1.0e-6, 1.0e-10) and ok
        for key in drift_keys:
            diff_key = key + '_diff'
            ok = _check_with_tolerance('serial_vs_mpi2:' + key,
                                       mpi2[diff_key], serial[diff_key],
                                       1.0e-10, 1.0e-10) and ok

    ok = len(_NEGATIVE_RESULTS) == 3 and ok
    if len(_NEGATIVE_RESULTS) != 3:
        logger.warning('Expected 3 negative checks, got %d', len(_NEGATIVE_RESULTS))

    return ok
