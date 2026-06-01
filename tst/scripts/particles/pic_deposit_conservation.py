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
_POSITIVE_RESULTS = {}
_NEGATIVE_RESULTS = {}


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
    nspecies = int(particles['nspecies'])
    cr_distribution = str(particles['cr_distribution'])
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
        'uniform_center_pdens': (
            cr_distribution == 'center' and nspecies == 1 and ppc == 1.0),
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
        'pdens_min': float(np.min(pdens_data['pdens'])),
        'pdens_max': float(np.max(pdens_data['pdens'])),
    }


def _run_command(label, nproc, arguments, expect_fail=False,
                 expected_message=None):
    command = ['./athena', '-i', _athena_input_path()] + list(arguments)
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


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    positive_cases = [
        {
            'name': 'serial',
            'basename': 'pic_dep_cons_serial',
            'nproc': 1,
            'args': ['job/basename=pic_dep_cons_serial'],
        },
        {
            'name': 'serial_nonunit_cell_volume',
            'basename': 'pic_dep_cons_nonunit_volume',
            'nproc': 1,
            'args': ['job/basename=pic_dep_cons_nonunit_volume',
                     'mesh/x1max=32.0',
                     'particles/cr_vx0=0.50',
                     'particles/cr_vy0=-0.25',
                     'particles/cr_vz0=0.125'],
        },
        {
            'name': 'serial_reflect_x1_crossing',
            'basename': 'pic_dep_cons_reflect_x1_crossing',
            'nproc': 1,
            'args': ['job/basename=pic_dep_cons_reflect_x1_crossing',
                     'mesh/ix1_bc=reflect',
                     'mesh/ox1_bc=reflect',
                     'particles/cr_vx0=8.0',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'],
            'check_current_totals': False,
            'expect_uniform_pdens': False,
        },
    ]
    if _athena_mpi_enabled():
        positive_cases.append({
            'name': 'mpi2',
            'basename': 'pic_dep_cons_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_dep_cons_mpi2'],
        })
    else:
        logger.info('Skipping mpi2 case: Athena build has MPI parallelism OFF')

    for case in positive_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _POSITIVE_RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(case['args']),
            'check_current_totals': case.get('check_current_totals', True),
            'expect_uniform_pdens': case.get('expect_uniform_pdens', True),
        }

    _run_command(
        'guard_deposit_moments_false',
        1,
        ['particles/deposit_moments=false', 'time/nlim=0'],
        expect_fail=True,
        expected_message='requires <particles>/deposit_moments=true',
    )
    _run_command(
        'guard_unsupported_deposit_order',
        1,
        ['particles/deposit_order=2', 'time/nlim=0'],
        expect_fail=True,
        expected_message='only deposit_order=1',
    )
    _run_command(
        'guard_unsupported_boundary_class',
        1,
        ['mesh/ix1_bc=diode', 'mesh/ox1_bc=diode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='does not support mesh/ix1_bc=diode',
    )
    _run_command(
        'guard_invalid_pic_background_mode',
        1,
        ['particles/pic_background_mode=bad_mode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_background_mode',
    )
    _run_command(
        'guard_invalid_pic_feedback_mode',
        1,
        ['particles/pic_feedback_mode=bad_mode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_feedback_mode',
    )
    _run_command(
        'guard_invalid_pic_interp_scheme',
        1,
        ['particles/pic_interp_scheme=linear', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_interp_scheme',
    )
    _run_command(
        'guard_invalid_pic_deltaf_mode',
        1,
        ['particles/pic_deltaf_mode=bad_mode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_deltaf_mode',
    )
    _run_command(
        'guard_invalid_pic_intermediate_arrays',
        1,
        ['particles/pic_intermediate_arrays=bad_mode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_intermediate_arrays',
    )
    _run_command(
        'guard_invalid_pic_expanding_box_mode',
        1,
        ['particles/pic_expanding_box_mode=bad_mode', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_expanding_box_mode',
    )
    _run_command(
        'guard_invalid_particle_type',
        1,
        ['particles/particle_type=bad_type', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/particle_type',
    )
    _run_command(
        'guard_invalid_ppc',
        1,
        ['particles/ppc=-1.0', 'time/nlim=0'],
        expect_fail=True,
        expected_message='<particles>/ppc must be >= 0',
    )
    _run_command(
        'guard_invalid_cr_distribution',
        1,
        ['particles/cr_distribution=bad_distribution', 'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/cr_distribution',
    )
    _run_command(
        'guard_invalid_pic_random_seed',
        1,
        ['particles/pic_random_seed=-1', 'time/nlim=0'],
        expect_fail=True,
        expected_message='<particles>/pic_random_seed must be >= 0',
    )
    _run_command(
        'guard_invalid_nspecies',
        1,
        ['particles/nspecies=0', 'time/nlim=0'],
        expect_fail=True,
        expected_message='<particles>/nspecies must be > 0',
    )
    _run_command(
        'guard_invalid_species_mass',
        1,
        ['species0/mass=0.0', 'time/nlim=0'],
        expect_fail=True,
        expected_message='<species0>/mass must be > 0',
    )
    _run_command(
        'guard_reserved_pic_cr_light_speed',
        1,
        ['particles/pic_cr_light_speed=2.0', 'time/nlim=0'],
        expect_fail=True,
        expected_message='<particles>/pic_cr_light_speed is reserved',
    )
    _run_command(
        'guard_pic_deltaf_requires_f0',
        1,
        ['particles/pic_deltaf_mode=on', 'time/nlim=0'],
        expect_fail=True,
        expected_message='requires <particles>/pic_deltaf_f0',
    )
    _run_command(
        'guard_invalid_pic_deltaf_f0',
        1,
        ['particles/pic_deltaf_mode=on',
         'particles/pic_deltaf_f0=not_a_supported_f0',
         'time/nlim=0'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/pic_deltaf_f0',
    )
    _run_command(
        'guard_pic_expansion_rate_requires_mode',
        1,
        ['particles/pic_expansion_rate_x1=1.0', 'time/nlim=0'],
        expect_fail=True,
        expected_message='require <particles>/pic_expanding_box_mode=on',
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name, result in _POSITIVE_RESULTS.items():
        measured = result['measured']
        expected = result['expected']
        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   1.0e-6, 1.0e-8) and ok
        if result['check_current_totals']:
            ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'],
                                       expected['Jx'], 1.0e-6, 1.0e-8) and ok
            ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'],
                                       expected['Jy'], 1.0e-6, 1.0e-8) and ok
            ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'],
                                       expected['Jz'], 1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':npart', measured['npart'],
                                   expected['npart'], 1.0e-6, 1.0e-8) and ok
        if expected['uniform_center_pdens'] and result['expect_uniform_pdens']:
            ok = _check_with_tolerance(case_name + ':pdens_min',
                                       measured['pdens_min'], 1.0,
                                       1.0e-6, 1.0e-8) and ok
            ok = _check_with_tolerance(case_name + ':pdens_max',
                                       measured['pdens_max'], 1.0,
                                       1.0e-6, 1.0e-8) and ok

    if 'mpi2' in _POSITIVE_RESULTS:
        serial = _POSITIVE_RESULTS['serial']['measured']
        mpi2 = _POSITIVE_RESULTS['mpi2']['measured']
        ok = _check_with_tolerance('serial_vs_mpi2:Q', mpi2['Q'], serial['Q'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jx', mpi2['Jx'], serial['Jx'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jy', mpi2['Jy'], serial['Jy'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jz', mpi2['Jz'], serial['Jz'],
                                   1.0e-6, 1.0e-8) and ok
    else:
        logger.info('No MPI run requested: serial-only non-MPI sanity path')

    ok = len(_NEGATIVE_RESULTS) == 19 and ok
    if len(_NEGATIVE_RESULTS) != 19:
        logger.warning('Expected 19 negative checks, got %d',
                       len(_NEGATIVE_RESULTS))

    return ok
