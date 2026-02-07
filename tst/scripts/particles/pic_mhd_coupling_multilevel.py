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

_INPUT_DECK = 'tests/pic_mhd_coupling_multilevel.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}
_VX0 = 0.5
_VY0 = -0.25
_VZ0 = 0.125


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _athena_mpi_enabled():
    command = ['./athena', '-c']
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


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


def _integrate_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _l2_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sqrt(np.sum(dataset[quantity] * dataset[quantity] * dvol)))


def _measure_case(basename):
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
    bcc_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_bcc'))
    m1_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m1'))
    m2_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m2'))
    m3_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m3'))
    e_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_e'))

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'npart': float(np.sum(pdens_data['pdens'])),
        'bcc1_l2': _l2_quantity(bcc_data, 'bcc1'),
        'bcc2_l2': _l2_quantity(bcc_data, 'bcc2'),
        'bcc3_l2': _l2_quantity(bcc_data, 'bcc3'),
        'mom1': _integrate_quantity(m1_data, 'mom1'),
        'mom2': _integrate_quantity(m2_data, 'mom2'),
        'mom3': _integrate_quantity(m3_data, 'mom3'),
        'ener': _integrate_quantity(e_data, 'ener'),
    }


def _run_command(label, nproc, arguments):
    command = ['./athena', '-i', _athena_input_path()] + list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + ': ' + ' '.join(command))


def _check_with_tolerance(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    reps = [('cc', 'cell_centered'), ('edge', 'edge_staggered')]
    orders = [('mhd', 'mhd_src_terms'), ('efield', 'efield_src')]
    common = [
        'particles/cr_vx0=0.50',
        'particles/cr_vy0=-0.25',
        'particles/cr_vz0=0.125',
        'particles/couple_moments_momentum_coeff=1.0',
        'particles/couple_moments_energy_coeff=10.0',
    ]

    cases = []
    for rep_tag, rep_value in reps:
        base = 'pic_mhd_ml_serial_' + rep_tag + '_uncoupled'
        args = [
            'job/basename=' + base,
            'particles/couple_moments_to_mhd=false',
            'particles/couple_moments_momentum_to_mhd=false',
            'particles/couple_moments_energy_to_mhd=false',
            'particles/couple_j_to_efield_representation=' + rep_value,
        ] + common
        cases.append({'name': 'serial_' + rep_tag + '_uncoupled',
                      'basename': base, 'nproc': 1, 'args': args})

    for order_tag, order_value in orders:
        for rep_tag, rep_value in reps:
            base = 'pic_mhd_ml_serial_' + rep_tag + '_' + order_tag
            args = [
                'job/basename=' + base,
                'particles/couple_moments_to_mhd=true',
                'particles/couple_moments_momentum_to_mhd=true',
                'particles/couple_moments_energy_to_mhd=true',
                'particles/couple_j_to_efield_representation=' + rep_value,
                'particles/couple_fluid_feedback_order=' + order_value,
            ] + common
            cases.append({'name': 'serial_' + rep_tag + '_' + order_tag,
                          'basename': base, 'nproc': 1, 'args': args})

    if _athena_mpi_enabled():
        for order_tag, order_value in orders:
            for rep_tag, rep_value in reps:
                base = 'pic_mhd_ml_mpi2_' + rep_tag + '_' + order_tag
                args = [
                    'job/basename=' + base,
                    'particles/couple_moments_to_mhd=true',
                    'particles/couple_moments_momentum_to_mhd=true',
                    'particles/couple_moments_energy_to_mhd=true',
                    'particles/couple_j_to_efield_representation=' + rep_value,
                    'particles/couple_fluid_feedback_order=' + order_value,
                ] + common
                cases.append({'name': 'mpi2_' + rep_tag + '_' + order_tag,
                              'basename': base, 'nproc': 2, 'args': args})
    else:
        logger.info('MPI disabled: running serial-only multilevel checks')

    for case in cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _RESULTS[case['name']] = {'measured': _measure_case(case['basename'])}


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for name, result in _RESULTS.items():
        measured = result['measured']
        q = measured['Q']
        npart = measured['npart']
        if q <= 0.0:
            logger.error('%s has non-positive total deposited charge: %e', name, q)
            ok = False
        if npart <= 0.0:
            logger.error('%s has non-positive particle count: %e', name, npart)
            ok = False
        if q > 0.0:
            ok = _check_with_tolerance(name + ':vx_ratio',
                                       measured['Jx'] / q, _VX0,
                                       1.0e-8, 1.0e-8) and ok
            ok = _check_with_tolerance(name + ':vy_ratio',
                                       measured['Jy'] / q, _VY0,
                                       1.0e-8, 1.0e-8) and ok
            ok = _check_with_tolerance(name + ':vz_ratio',
                                       measured['Jz'] / q, _VZ0,
                                       1.0e-8, 1.0e-8) and ok

    for order_tag in ['mhd', 'efield']:
        cc_name = 'serial_cc_' + order_tag
        edge_name = 'serial_edge_' + order_tag
        if cc_name in _RESULTS and edge_name in _RESULTS:
            cc = _RESULTS[cc_name]['measured']
            edge = _RESULTS[edge_name]['measured']
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                ok = _check_with_tolerance(order_tag + ':cc_vs_edge:' + quantity,
                                           cc[quantity], edge[quantity],
                                           1.0e-6, 1.0e-8) and ok

    for rep_tag in ['cc', 'edge']:
        mhd_name = 'serial_' + rep_tag + '_mhd'
        ef_name = 'serial_' + rep_tag + '_efield'
        if mhd_name in _RESULTS and ef_name in _RESULTS:
            mhd = _RESULTS[mhd_name]['measured']
            efield = _RESULTS[ef_name]['measured']
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                ok = _check_with_tolerance(rep_tag + ':mhd_vs_efield:' + quantity,
                                           mhd[quantity], efield[quantity],
                                           1.0e-6, 1.0e-8) and ok

    for rep_tag in ['cc', 'edge']:
        unc_name = 'serial_' + rep_tag + '_uncoupled'
        for order_tag in ['mhd', 'efield']:
            coupled_name = 'serial_' + rep_tag + '_' + order_tag
            if unc_name not in _RESULTS or coupled_name not in _RESULTS:
                continue
            unc = _RESULTS[unc_name]['measured']
            coupled = _RESULTS[coupled_name]['measured']
            delta1 = abs(coupled['bcc1_l2'] - unc['bcc1_l2'])
            delta2 = abs(coupled['bcc2_l2'] - unc['bcc2_l2'])
            delta3 = abs(coupled['bcc3_l2'] - unc['bcc3_l2'])
            delta = max(delta1, delta2, delta3)
            logger.info(
                '%s coupled_field_delta bcc1=% .8e bcc2=% .8e bcc3=% .8e max=% .8e',
                coupled_name, delta1, delta2, delta3, delta)
            if delta <= 1.0e-8:
                logger.error(
                    '%s expected nonzero coupled-field response on multilevel mesh',
                    coupled_name)
                ok = False

    for order_tag in ['mhd', 'efield']:
        for rep_tag in ['cc', 'edge']:
            serial_name = 'serial_' + rep_tag + '_' + order_tag
            mpi_name = 'mpi2_' + rep_tag + '_' + order_tag
            if serial_name not in _RESULTS or mpi_name not in _RESULTS:
                continue
            serial = _RESULTS[serial_name]['measured']
            mpi = _RESULTS[mpi_name]['measured']
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                             'bcc1_l2', 'bcc2_l2', 'bcc3_l2',
                             'mom1', 'mom2', 'mom3', 'ener']:
                ok = _check_with_tolerance(
                    mpi_name + ':' + quantity + '_vs_serial',
                    mpi[quantity], serial[quantity], 1.0e-6, 1.0e-8) and ok

    return ok
