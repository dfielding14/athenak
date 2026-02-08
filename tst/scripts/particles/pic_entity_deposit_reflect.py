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

_INPUT_DECK = 'tests/pic_entity_deposit_reflect.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _deck_path_for_python():
    return os.path.join('..', 'inputs', _INPUT_DECK)


def _athena_mpi_enabled():
    proc = subprocess.run(['./athena', '-c'], cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


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


def _species_counts(npart, nspecies):
    base = npart // nspecies
    rem = npart % nspecies
    return [base + (1 if s < rem else 0) for s in range(nspecies)]


def _expected_totals(arguments):
    config = bin_convert.athinput(_deck_path_for_python())
    _apply_overrides(config, arguments)

    mesh = config['mesh']
    particles = config['particles']

    nx1 = int(mesh['nx1'])
    nx2 = int(mesh['nx2'])
    nx3 = int(mesh['nx3'])
    ppc = float(particles['ppc'])
    qscale = float(particles['deposit_qscale'])
    nspecies = int(particles['nspecies'])
    vx0 = float(particles['cr_vx0'])
    vy0 = float(particles['cr_vy0'])
    vz0 = float(particles['cr_vz0'])

    npart = int(ppc * nx1 * nx2 * nx3)
    counts = _species_counts(npart, nspecies)

    total_q = 0.0
    for s in range(nspecies):
        block = 'species' + str(s)
        charge = float(config[block]['charge'])
        total_q += counts[s] * charge
    total_q *= qscale

    return {
        'npart': float(npart),
        'Q': float(total_q),
        'Jx': float(total_q * vx0),
        'Jy': float(total_q * vy0),
        'Jz': float(total_q * vz0),
    }


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _latest_output_file(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if not matches:
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


def _read_quantity(basename, file_id, quantity):
    data = bin_convert.read_binary_as_athdf(_latest_output_file(basename, file_id))
    return np.array(data[quantity], copy=True)


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
        'jx_l2': _l2_quantity(jx_data, 'prtcl_jx'),
        'jy_l2': _l2_quantity(jy_data, 'prtcl_jy'),
        'jz_l2': _l2_quantity(jz_data, 'prtcl_jz'),
        'fields': {
            'rho': _read_quantity(basename, 'prtcl_rho', 'prtcl_rho'),
            'jx': _read_quantity(basename, 'prtcl_jx', 'prtcl_jx'),
            'jy': _read_quantity(basename, 'prtcl_jy', 'prtcl_jy'),
            'jz': _read_quantity(basename, 'prtcl_jz', 'prtcl_jz'),
        },
    }


def _run_athena(label, nproc, arguments):
    command = ['./athena', '-i', _athena_input_path()] + list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _check_with_tolerance(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    cases = [('np1', 1)]
    if _athena_mpi_enabled():
        cases += [('np2', 2), ('np4', 4)]
    else:
        logger.info('MPI disabled: running serial-only reflect-boundary check')

    for case_name, nproc in cases:
        basename = 'pic_entity_dep_reflect_' + case_name
        args = ['job/basename=' + basename]
        _remove_outputs(basename)
        _run_athena(case_name, nproc, args)
        _RESULTS[case_name] = {
            'measured': _measure_case(basename),
            'expected': _expected_totals(args),
        }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    base = _RESULTS['np1']

    for case_name, result in _RESULTS.items():
        measured = result['measured']
        expected = result['expected']
        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   1.0e-8, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'], expected['Jx'],
                                   1.0e-8, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'], expected['Jy'],
                                   1.0e-8, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'], expected['Jz'],
                                   1.0e-8, 1.0e-8) and ok
        logger.info('%s:npart measured=% .8e initial=% .8e retention=% .8e',
                    case_name, measured['npart'], expected['npart'],
                    measured['npart'] / max(expected['npart'], 1.0))
        ok = (measured['npart'] > 0.0) and ok

    base_measured = base['measured']
    for case_name, result in _RESULTS.items():
        if case_name == 'np1':
            continue
        measured = result['measured']
        for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart', 'jx_l2', 'jy_l2', 'jz_l2']:
            ok = _check_with_tolerance(
                case_name + ':' + quantity + '_vs_np1',
                measured[quantity], base_measured[quantity], 1.0e-8, 1.0e-8) and ok

        base_fields = base_measured['fields']
        fields = measured['fields']
        for field_name in ['rho', 'jx', 'jy', 'jz']:
            max_abs_diff = float(np.max(np.abs(fields[field_name] -
                                               base_fields[field_name])))
            logger.info('%s:%s max_abs_diff=% .8e',
                        case_name, field_name, max_abs_diff)
            ok = (max_abs_diff <= 1.0e-8) and ok

    return ok
