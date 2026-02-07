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

_INPUT_DECK = 'tests/pic_mhd_current_coupling.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


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


def _deck_path_for_python():
    return os.path.join('..', 'inputs', _INPUT_DECK)


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
    config = bin_convert.athinput(_deck_path_for_python())
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

    common_args = [
        'particles/couple_moments_to_mhd=true',
        'particles/couple_moments_momentum_to_mhd=true',
        'particles/couple_moments_energy_to_mhd=true',
        'particles/couple_moments_momentum_coeff=1.0',
        'particles/couple_moments_energy_coeff=10.0',
        'particles/cr_vx0=0.50',
        'particles/cr_vy0=-0.25',
        'particles/cr_vz0=0.125',
    ]

    base_cases = [
        {
            'name': 'serial_mb444',
            'basename': 'pic_mhd_decomp_serial',
            'nproc': 1,
            'args': [
                'job/basename=pic_mhd_decomp_serial',
                'meshblock/nx1=4',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
            ],
        },
    ]

    if _athena_mpi_enabled():
        base_cases.extend([
            {
                'name': 'mpi2_mb444',
                'basename': 'pic_mhd_decomp_mpi2_same_mb',
                'nproc': 2,
                'args': [
                    'job/basename=pic_mhd_decomp_mpi2_same_mb',
                    'meshblock/nx1=4',
                    'meshblock/nx2=4',
                    'meshblock/nx3=4',
                ],
            },
            {
                'name': 'mpi4_mb444',
                'basename': 'pic_mhd_decomp_mpi4_same_mb',
                'nproc': 4,
                'args': [
                    'job/basename=pic_mhd_decomp_mpi4_same_mb',
                    'meshblock/nx1=4',
                    'meshblock/nx2=4',
                    'meshblock/nx3=4',
                ],
            },
        ])
    else:
        logger.info('MPI disabled: running serial-only decomposition baseline')

    cases = []
    representations = [('cc', 'cell_centered'), ('edge', 'edge_staggered')]
    feedback_orders = [('mhd', 'mhd_src_terms'), ('efield', 'efield_src')]
    for order_tag, order_value in feedback_orders:
        for rep_tag, rep_value in representations:
            for base_case in base_cases:
                combo_tag = rep_tag + '_' + order_tag
                rep_basename = base_case['basename'] + '_' + combo_tag
                cases.append({
                    'name': base_case['name'] + '_' + combo_tag,
                    'basename': rep_basename,
                    'nproc': base_case['nproc'],
                    'args': (['job/basename=' + rep_basename] +
                             list(base_case['args'][1:]) +
                             ['particles/couple_j_to_efield_representation=' + rep_value,
                              'particles/couple_fluid_feedback_order=' + order_value] +
                             common_args),
                })

    for case in cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(case['args']),
        }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name, result in _RESULTS.items():
        measured = result['measured']
        expected = result['expected']
        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'], expected['Jx'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'], expected['Jy'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'], expected['Jz'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':npart', measured['npart'],
                                   expected['npart'], 1.0e-6, 1.0e-8) and ok

    for order_tag in ['mhd', 'efield']:
        for rep_tag in ['cc', 'edge']:
            serial_case = 'serial_mb444_' + rep_tag + '_' + order_tag
            if serial_case not in _RESULTS:
                logger.warning('Missing serial baseline result for %s %s',
                               rep_tag, order_tag)
                return False

            baseline = _RESULTS[serial_case]['measured']
            for case_name in ['mpi2_mb444_' + rep_tag + '_' + order_tag,
                              'mpi4_mb444_' + rep_tag + '_' + order_tag]:
                if case_name not in _RESULTS:
                    continue
                measured = _RESULTS[case_name]['measured']
                ok = _check_with_tolerance(case_name + ':Q_vs_serial',
                                           measured['Q'], baseline['Q'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':Jx_vs_serial',
                                           measured['Jx'], baseline['Jx'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':Jy_vs_serial',
                                           measured['Jy'], baseline['Jy'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':Jz_vs_serial',
                                           measured['Jz'], baseline['Jz'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':npart_vs_serial',
                                           measured['npart'], baseline['npart'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':bcc1_l2_vs_serial',
                                           measured['bcc1_l2'], baseline['bcc1_l2'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':bcc2_l2_vs_serial',
                                           measured['bcc2_l2'], baseline['bcc2_l2'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':bcc3_l2_vs_serial',
                                           measured['bcc3_l2'], baseline['bcc3_l2'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':mom1_vs_serial',
                                           measured['mom1'], baseline['mom1'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':mom2_vs_serial',
                                           measured['mom2'], baseline['mom2'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':mom3_vs_serial',
                                           measured['mom3'], baseline['mom3'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance(case_name + ':ener_vs_serial',
                                           measured['ener'], baseline['ener'],
                                           1.0e-6, 1.0e-8) and ok

    for order_tag in ['mhd', 'efield']:
        serial_cc = 'serial_mb444_cc_' + order_tag
        serial_edge = 'serial_mb444_edge_' + order_tag
        if serial_cc in _RESULTS and serial_edge in _RESULTS:
            cc = _RESULTS[serial_cc]['measured']
            edge = _RESULTS[serial_edge]['measured']
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                ok = _check_with_tolerance(order_tag + ':cc_vs_edge:' + quantity,
                                           cc[quantity], edge[quantity],
                                           1.0e-6, 1.0e-8) and ok

    for rep_tag in ['cc', 'edge']:
        serial_mhd = 'serial_mb444_' + rep_tag + '_mhd'
        serial_efield = 'serial_mb444_' + rep_tag + '_efield'
        if serial_mhd in _RESULTS and serial_efield in _RESULTS:
            mhd = _RESULTS[serial_mhd]['measured']
            efield = _RESULTS[serial_efield]['measured']
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                ok = _check_with_tolerance(rep_tag + ':mhd_vs_efield:' + quantity,
                                           mhd[quantity], efield[quantity],
                                           1.0e-6, 1.0e-8) and ok

    return ok
