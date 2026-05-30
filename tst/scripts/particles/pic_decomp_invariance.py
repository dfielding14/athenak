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

_INPUT_DECK = 'tests/pic_deposit_conservation.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


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

    npart = int(ppc * nx1 * nx2 * nx3)
    return {
        'npart': float(npart),
        'Q': float(npart) * qscale * charge,
        'Jx': 0.0,
        'Jy': 0.0,
        'Jz': 0.0,
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

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'npart': float(np.sum(pdens_data['pdens'])),
        'rho_field': np.asarray(rho_data['prtcl_rho'], dtype=float),
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

    cases = [
        {
            'name': 'serial_mb444',
            'basename': 'pic_dep_decomp_serial',
            'nproc': 1,
            'args': [
                'job/basename=pic_dep_decomp_serial',
                'meshblock/nx1=4',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
            ],
        },
        {
            'name': 'mpi2_mb444',
            'basename': 'pic_dep_decomp_mpi2_same_mb',
            'nproc': 2,
            'args': [
                'job/basename=pic_dep_decomp_mpi2_same_mb',
                'meshblock/nx1=4',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
            ],
        },
        {
            'name': 'mpi2_mb842',
            'basename': 'pic_dep_decomp_mpi2_alt_mb',
            'nproc': 2,
            'args': [
                'job/basename=pic_dep_decomp_mpi2_alt_mb',
                'meshblock/nx1=8',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
            ],
        },
        {
            'name': 'serial_frac_random_mb444',
            'basename': 'pic_dep_decomp_frac_random_serial',
            'nproc': 1,
            'args': [
                'job/basename=pic_dep_decomp_frac_random_serial',
                'meshblock/nx1=4',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
                'particles/ppc=0.3',
                'particles/cr_distribution=random',
                'particles/pic_random_seed=17',
            ],
        },
        {
            'name': 'mpi2_frac_random_mb444',
            'basename': 'pic_dep_decomp_frac_random_mpi2',
            'nproc': 2,
            'args': [
                'job/basename=pic_dep_decomp_frac_random_mpi2',
                'meshblock/nx1=4',
                'meshblock/nx2=4',
                'meshblock/nx3=4',
                'particles/ppc=0.3',
                'particles/cr_distribution=random',
                'particles/pic_random_seed=17',
            ],
        },
    ]

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

    baseline = _RESULTS['serial_mb444']['measured']
    for case_name in ['mpi2_mb444', 'mpi2_mb842']:
        measured = _RESULTS[case_name]['measured']
        ok = _check_with_tolerance(case_name + ':Q_vs_serial',
                                   measured['Q'], baseline['Q'], 1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jx_vs_serial',
                                   measured['Jx'], baseline['Jx'], 1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jy_vs_serial',
                                   measured['Jy'], baseline['Jy'], 1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jz_vs_serial',
                                   measured['Jz'], baseline['Jz'], 1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':npart_vs_serial',
                                   measured['npart'], baseline['npart'],
                                   1.0e-6, 1.0e-8) and ok

    same_mb = _RESULTS['mpi2_mb444']['measured']
    rho_diff = np.max(np.abs(same_mb['rho_field'] - baseline['rho_field']))
    ok = _check_with_tolerance('mpi2_mb444:rho_field_vs_serial',
                               rho_diff, 0.0, 1.0e-12, 0.0) and ok

    frac_serial = _RESULTS['serial_frac_random_mb444']['measured']
    frac_mpi = _RESULTS['mpi2_frac_random_mb444']['measured']
    for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
        ok = _check_with_tolerance('mpi2_frac_random_mb444:' + quantity +
                                   '_vs_serial',
                                   frac_mpi[quantity], frac_serial[quantity],
                                   1.0e-10, 1.0e-10) and ok
    rho_diff = np.max(np.abs(frac_mpi['rho_field'] - frac_serial['rho_field']))
    ok = _check_with_tolerance('mpi2_frac_random_mb444:rho_field_vs_serial',
                               rho_diff, 0.0, 1.0e-12, 0.0) and ok

    return ok
