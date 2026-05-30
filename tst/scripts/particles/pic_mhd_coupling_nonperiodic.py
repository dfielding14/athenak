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

_INPUT_DECK = 'tests/pic_mhd_coupling_nonperiodic.athinput'
_DEFAULT_MODE_INPUT_DECK = (
    'tests/pic_mhd_coupling_nonperiodic_default_mode.athinput')
_UNSUPPORTED_DECK = 'tests/pic_mhd_coupling_guard_nonperiodic_unsupported.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}
_NEGATIVE_RESULTS = {}
_MOMENT_ABS_TOL = 1.0e-4
_MOMENT_REL_TOL = 1.0e-7


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

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'npart': float(np.sum(pdens_data['pdens'])),
        'bcc1_l2': _l2_quantity(bcc_data, 'bcc1'),
        'bcc2_l2': _l2_quantity(bcc_data, 'bcc2'),
        'bcc3_l2': _l2_quantity(bcc_data, 'bcc3'),
    }


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


def _tolerances_for_quantity(quantity):
    if quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
        return _MOMENT_ABS_TOL, _MOMENT_REL_TOL
    return 1.0e-6, 1.0e-8


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    common_args = [
        'particles/cr_vx0=0.50',
        'particles/cr_vy0=-0.25',
        'particles/cr_vz0=0.125',
    ]

    cases = [
        {
            'name': 'serial_uncoupled',
            'basename': 'pic_mhd_np_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_uncoupled',
                     'particles/couple_moments_to_mhd=false'] + common_args,
        },
        {
            'name': 'serial_coupled_cc',
            'basename': 'pic_mhd_np_coupled_cc',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_coupled_cc',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=cell_centered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge',
            'basename': 'pic_mhd_np_coupled_edge',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_coupled_edge',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_default_mode',
            'basename': 'pic_mhd_np_coupled_edge_default',
            'nproc': 1,
            'input_deck': _DEFAULT_MODE_INPUT_DECK,
            'args': ['job/basename=pic_mhd_np_coupled_edge_default',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct',
            'basename': 'pic_mhd_np_coupled_edge_direct',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_coupled_edge_direct',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_o2',
            'basename': 'pic_mhd_np_coupled_edge_direct_o2',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_coupled_edge_direct_o2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + common_args,
        },
        {
            'name': 'serial_zero_uncoupled',
            'basename': 'pic_mhd_np_zero_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_zero_uncoupled',
                     'particles/couple_moments_to_mhd=false',
                     'particles/cr_vx0=0.0',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'],
        },
        {
            'name': 'serial_zero_coupled',
            'basename': 'pic_mhd_np_zero_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_zero_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/cr_vx0=0.0',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'],
        },
        {
            'name': 'serial_coupled_inflow',
            'basename': 'pic_mhd_np_coupled_inflow',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_np_coupled_inflow',
                     'mesh/ix2_bc=inflow',
                     'mesh/ox2_bc=inflow',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=cell_centered'] +
                    common_args,
        },
    ]

    if _athena_mpi_enabled():
        cases.extend([
            {
                'name': 'mpi2_coupled_cc',
                'basename': 'pic_mhd_np_coupled_cc_mpi2',
                'nproc': 2,
                'args': ['job/basename=pic_mhd_np_coupled_cc_mpi2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=cell_centered'] +
                        common_args,
            },
            {
                'name': 'mpi2_coupled_edge',
                'basename': 'pic_mhd_np_coupled_edge_mpi2',
                'nproc': 2,
                'args': ['job/basename=pic_mhd_np_coupled_edge_mpi2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=edge_staggered',
                         'particles/couple_j_deposition_mode=cc_convert'] +
                        common_args,
            },
            {
                'name': 'mpi2_coupled_edge_default_mode',
                'basename': 'pic_mhd_np_coupled_edge_default_mpi2',
                'nproc': 2,
                'input_deck': _DEFAULT_MODE_INPUT_DECK,
                'args': ['job/basename=pic_mhd_np_coupled_edge_default_mpi2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=edge_staggered'] +
                        common_args,
            },
            {
                'name': 'mpi2_coupled_edge_direct',
                'basename': 'pic_mhd_np_coupled_edge_direct_mpi2',
                'nproc': 2,
                'args': ['job/basename=pic_mhd_np_coupled_edge_direct_mpi2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=edge_staggered',
                         'particles/couple_j_deposition_mode=direct_staggered'] +
                        common_args,
            },
            {
                'name': 'mpi2_coupled_edge_direct_o2',
                'basename': 'pic_mhd_np_coupled_edge_direct_o2_mpi2',
                'nproc': 2,
                'args': ['job/basename=pic_mhd_np_coupled_edge_direct_o2_mpi2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=edge_staggered',
                         'particles/couple_j_deposition_mode=direct_staggered',
                         'particles/deposit_order=2'] + common_args,
            },
        ])
    else:
        logger.info('MPI disabled: running serial-only non-periodic checks')

    for case in cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'],
                     input_deck=case.get('input_deck', _INPUT_DECK))
        _RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(case['args']),
        }

    _run_command(
        'guard_unsupported_boundary_diode',
        1,
        ['time/nlim=0'],
        expect_fail=True,
        expected_message='does not support mesh/ix1_bc=diode',
        input_deck=_UNSUPPORTED_DECK,
    )
    _run_command(
        'guard_unsupported_boundary_vacuum',
        1,
        ['mesh/ix3_bc=vacuum', 'mesh/ox3_bc=vacuum', 'time/nlim=0'],
        expect_fail=True,
        expected_message='does not support mesh/ix3_bc=vacuum',
    )
    _run_command(
        'guard_unsupported_boundary_user',
        1,
        ['mesh/ix1_bc=user', 'mesh/ox1_bc=user', 'time/nlim=0'],
        expect_fail=True,
        expected_message='does not support mesh/ix1_bc=user',
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name, result in _RESULTS.items():
        measured = result['measured']
        expected = result['expected']
        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'], expected['Jx'],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'], expected['Jy'],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'], expected['Jz'],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':npart', measured['npart'],
                                   expected['npart'], _MOMENT_ABS_TOL,
                                   _MOMENT_REL_TOL) and ok

    unc = _RESULTS['serial_uncoupled']['measured']
    cpl_cc = _RESULTS['serial_coupled_cc']['measured']
    cpl_edge = _RESULTS['serial_coupled_edge']['measured']
    cpl_edge_default = _RESULTS['serial_coupled_edge_default_mode']['measured']
    cpl_direct = _RESULTS['serial_coupled_edge_direct']['measured']
    cpl_direct_o2 = _RESULTS['serial_coupled_edge_direct_o2']['measured']

    delta_cc = max(abs(cpl_cc['bcc1_l2'] - unc['bcc1_l2']),
                   abs(cpl_cc['bcc2_l2'] - unc['bcc2_l2']),
                   abs(cpl_cc['bcc3_l2'] - unc['bcc3_l2']))
    delta_edge = max(abs(cpl_edge['bcc1_l2'] - unc['bcc1_l2']),
                     abs(cpl_edge['bcc2_l2'] - unc['bcc2_l2']),
                     abs(cpl_edge['bcc3_l2'] - unc['bcc3_l2']))
    delta_direct = max(abs(cpl_direct['bcc1_l2'] - unc['bcc1_l2']),
                       abs(cpl_direct['bcc2_l2'] - unc['bcc2_l2']),
                       abs(cpl_direct['bcc3_l2'] - unc['bcc3_l2']))
    delta_direct_o2 = max(abs(cpl_direct_o2['bcc1_l2'] - unc['bcc1_l2']),
                          abs(cpl_direct_o2['bcc2_l2'] - unc['bcc2_l2']),
                          abs(cpl_direct_o2['bcc3_l2'] - unc['bcc3_l2']))
    logger.info('nonperiodic_coupled_delta cc=% .8e edge=% .8e direct=% .8e '
                'direct_o2=% .8e',
                delta_cc, delta_edge, delta_direct, delta_direct_o2)
    if (delta_cc <= 1.0e-12 or delta_edge <= 1.0e-12 or
            delta_direct <= 1.0e-12 or delta_direct_o2 <= 1.0e-12):
        logger.warning('Expected nonzero coupled-field response in non-periodic run')
        ok = False

    for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
        ok = _check_with_tolerance('nonperiodic_default_mode_vs_direct:' + quantity,
                                   cpl_edge_default[quantity], cpl_direct[quantity],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance('nonperiodic_cc_vs_edge:' + quantity,
                                   cpl_cc[quantity], cpl_edge[quantity],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance('nonperiodic_cc_vs_direct:' + quantity,
                                   cpl_cc[quantity], cpl_direct[quantity],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
        ok = _check_with_tolerance('nonperiodic_cc_vs_direct_o2:' + quantity,
                                   cpl_cc[quantity], cpl_direct_o2[quantity],
                                   _MOMENT_ABS_TOL, _MOMENT_REL_TOL) and ok
    for quantity in ['bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
        ok = _check_with_tolerance('nonperiodic_default_mode_vs_direct:' + quantity,
                                   cpl_edge_default[quantity], cpl_direct[quantity],
                                   1.0e-6, 1.0e-8) and ok
    default_vs_convert_delta = max(
        abs(cpl_edge_default['bcc1_l2'] - cpl_edge['bcc1_l2']),
        abs(cpl_edge_default['bcc2_l2'] - cpl_edge['bcc2_l2']),
        abs(cpl_edge_default['bcc3_l2'] - cpl_edge['bcc3_l2']))
    logger.info('nonperiodic_default_edge_mode_delta_vs_cc max=% .8e',
                default_vs_convert_delta)
    if default_vs_convert_delta <= 1.0e-12:
        logger.warning('Default edge deposition mode appears identical to '
                       'cc_convert baseline')
        ok = False

    zero_unc = _RESULTS['serial_zero_uncoupled']['measured']
    zero_cpl = _RESULTS['serial_zero_coupled']['measured']
    ok = _check_with_tolerance('nonperiodic_zero_current:bcc1_l2',
                               zero_cpl['bcc1_l2'], zero_unc['bcc1_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('nonperiodic_zero_current:bcc2_l2',
                               zero_cpl['bcc2_l2'], zero_unc['bcc2_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('nonperiodic_zero_current:bcc3_l2',
                               zero_cpl['bcc3_l2'], zero_unc['bcc3_l2'],
                               1.0e-12, 1.0e-12) and ok

    if 'mpi2_coupled_cc' in _RESULTS:
        mpi_cc = _RESULTS['mpi2_coupled_cc']['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                         'bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
            abs_tol, rel_tol = _tolerances_for_quantity(quantity)
            ok = _check_with_tolerance('nonperiodic_mpi_cc:' + quantity,
                                       mpi_cc[quantity], cpl_cc[quantity],
                                       abs_tol, rel_tol) and ok

    if 'mpi2_coupled_edge' in _RESULTS:
        mpi_edge = _RESULTS['mpi2_coupled_edge']['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                         'bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
            abs_tol, rel_tol = _tolerances_for_quantity(quantity)
            ok = _check_with_tolerance('nonperiodic_mpi_edge:' + quantity,
                                       mpi_edge[quantity], cpl_edge[quantity],
                                       abs_tol, rel_tol) and ok

    if 'mpi2_coupled_edge_default_mode' in _RESULTS:
        mpi_edge_default = _RESULTS['mpi2_coupled_edge_default_mode']['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                         'bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
            abs_tol, rel_tol = _tolerances_for_quantity(quantity)
            ok = _check_with_tolerance('nonperiodic_mpi_default_mode:' + quantity,
                                       mpi_edge_default[quantity],
                                       cpl_edge_default[quantity],
                                       abs_tol, rel_tol) and ok

    if 'mpi2_coupled_edge_direct' in _RESULTS:
        mpi_direct = _RESULTS['mpi2_coupled_edge_direct']['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                         'bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
            abs_tol, rel_tol = _tolerances_for_quantity(quantity)
            ok = _check_with_tolerance('nonperiodic_mpi_direct:' + quantity,
                                       mpi_direct[quantity], cpl_direct[quantity],
                                       abs_tol, rel_tol) and ok

    if 'mpi2_coupled_edge_direct_o2' in _RESULTS:
        mpi_direct_o2 = _RESULTS['mpi2_coupled_edge_direct_o2']['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                         'bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
            abs_tol, rel_tol = _tolerances_for_quantity(quantity)
            ok = _check_with_tolerance('nonperiodic_mpi_direct_o2:' + quantity,
                                       mpi_direct_o2[quantity],
                                       cpl_direct_o2[quantity],
                                       abs_tol, rel_tol) and ok

    ok = (len(_NEGATIVE_RESULTS) == 3) and ok
    if len(_NEGATIVE_RESULTS) != 3:
        logger.warning('Expected 3 negative checks, got %d', len(_NEGATIVE_RESULTS))

    return ok
