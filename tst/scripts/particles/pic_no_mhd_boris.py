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

_INPUT_DECK = 'tests/pic_no_mhd_boris.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_POSITIVE_RESULTS = {}
_NEGATIVE_RESULTS = {}
_TIMESTEP_RESULTS = {}
_PVTK_RESULTS = {}
_PDENS_RESULTS = {}
_AMR_RESULTS = {}
_MOMENT_TOTAL_ABS_TOL = 1.0e-4
_MOMENT_TOTAL_REL_TOL = 3.0e-7


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
    for pattern in [
            os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin'),
            os.path.join(_athena_exe_dir(), 'pvtk', basename + '.*.part.vtk')]:
        for fname in glob.glob(pattern):
            os.remove(fname)


def _latest_output_file(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches[-1]


def _latest_pvtk_file(basename, file_id='prtcl_all'):
    pattern = os.path.join(_athena_exe_dir(), 'pvtk',
                           f'{basename}.{file_id}.*.part.vtk')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No particle VTK files found for pattern: ' + pattern)
    return matches[-1]


def _read_big_endian_floats(contents, offset, count, label):
    payload = contents[offset:offset + 4 * count]
    if len(payload) != 4 * count:
        raise RuntimeError('Truncated ' + label + ' payload in particle VTK output')
    return np.frombuffer(payload, dtype='>f4').astype(np.float64)


def _read_big_endian_ints(contents, offset, count, label):
    payload = contents[offset:offset + 4 * count]
    if len(payload) != 4 * count:
        raise RuntimeError('Truncated ' + label + ' payload in particle VTK output')
    return np.frombuffer(payload, dtype='>i4').astype(np.int64)


def _find_marker(contents, pattern, offset, label):
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise RuntimeError('Could not find ' + label + ' marker in particle VTK output')
    return offset + match.start(), offset + match.end(), match


def _read_pvtk_snapshot(basename):
    path = _latest_pvtk_file(basename)
    with open(path, 'rb') as fp:
        contents = fp.read()

    _, offset, match = _find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS')
    npoint = int(match.group(1))
    points = _read_big_endian_floats(contents, offset, 3 * npoint, 'POINTS')
    points = points.reshape(npoint, 3)
    offset += 4 * 3 * npoint

    scalars = {}
    for name in ['gid', 'ptag', 'species']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' int\nLOOKUP_TABLE default\n')
        _, offset, _ = _find_marker(contents, pattern, offset, 'SCALARS ' + name)
        scalars[name] = _read_big_endian_ints(contents, offset, npoint, name)
        offset += 4 * npoint

    _, offset, _ = _find_marker(contents, rb'\nVECTORS vel float\n',
                                offset, 'VECTORS vel')
    velocity = _read_big_endian_floats(contents, offset, 3 * npoint, 'vel')
    velocity = velocity.reshape(npoint, 3)

    return {
        'path': path,
        'npoint': npoint,
        'points': points,
        'gid': scalars['gid'],
        'ptag': scalars['ptag'],
        'species': scalars['species'],
        'velocity': velocity,
    }


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
    }


def _measure_pdens_summary(basename):
    pdens_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_d'))
    pdens = pdens_data['pdens']
    return {
        'npart': float(np.sum(pdens)),
        'finite': bool(np.all(np.isfinite(pdens))),
        'max': float(np.max(pdens)),
    }


def _extract_cycle_dt(output, cycle):
    pattern = re.compile(r'cycle=\s*' + str(cycle) +
                         r'\s+time=\s*[0-9eE+\-.]+\s+dt=\s*([0-9eE+\-.]+)')
    for line in output.splitlines():
        match = pattern.search(line)
        if match is not None:
            return float(match.group(1))
    raise RuntimeError('Could not find cycle ' + str(cycle) + ' dt in Athena output')


def _parse_amr_created(output):
    match = re.search(r'(\d+) MeshBlocks created,\s+(\d+) deleted by AMR', output)
    if match is None:
        raise RuntimeError('Could not find AMR creation telemetry in Athena output')
    return int(match.group(1))


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
        return output

    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + ': ' +
                           ' '.join(command) + '\n' + output)
    return output


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
            'basename': 'pic_no_mhd_boris_serial',
            'nproc': 1,
            'args': ['job/basename=pic_no_mhd_boris_serial'],
        },
    ]
    if _athena_mpi_enabled():
        positive_cases.append({
            'name': 'mpi2',
            'basename': 'pic_no_mhd_boris_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_no_mhd_boris_mpi2'],
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

    dt_output = _run_command(
        'cell_cross_dt_limit',
        1,
        [
            'job/basename=pic_no_mhd_boris_dtlimit',
            'particles/cr_vx0=100.0',
            'particles/cr_vy0=0.0',
            'particles/cr_vz0=0.0',
            'particles/pic_max_cell_cross=1',
            'time/nlim=1',
            'output1/dcycle=0',
            'output2/dcycle=0',
            'output3/dcycle=0',
            'output4/dcycle=0',
            'output5/dcycle=0',
        ],
    )
    _TIMESTEP_RESULTS['cell_cross_dt_limit'] = {
        'measured_dt': _extract_cycle_dt(dt_output, 0),
        'expected_dt': 0.3 * 1.0 / 100.0,
    }

    theta_output = _run_command(
        'gyro_angle_dt_limit',
        1,
        [
            'job/basename=pic_no_mhd_boris_thetalimit',
            'particles/cr_vx0=0.0',
            'particles/cr_vy0=0.0',
            'particles/cr_vz0=0.0',
            'particles/pic_no_mhd_bz=50.0',
            'particles/pic_theta_max=0.1',
            'time/nlim=1',
            'output1/dcycle=0',
            'output2/dcycle=0',
            'output3/dcycle=0',
            'output4/dcycle=0',
            'output5/dcycle=0',
        ],
    )
    _TIMESTEP_RESULTS['gyro_angle_dt_limit'] = {
        'measured_dt': _extract_cycle_dt(theta_output, 0),
        'expected_dt': 0.3 * 0.1 / 50.0,
    }

    amr_output = _run_command(
        'no_mhd_amr_capacity',
        1,
        [
            'job/basename=pic_no_mhd_boris_amr_capacity',
            'mesh_refinement/refinement=adaptive',
            'mesh_refinement/num_levels=2',
            'mesh_refinement/ddens_max=1.0e-10',
            'mesh_refinement/ncycle_check=1',
            'mesh_refinement/refinement_interval=1',
            'mesh_refinement/max_nmb_per_rank=128',
            'time/nlim=3',
            'time/tlim=0.2',
            'output1/dcycle=0',
            'output2/dcycle=0',
            'output3/dcycle=0',
            'output4/dcycle=0',
            'output5/dcycle=0',
            'output6/dcycle=0',
        ],
    )
    _AMR_RESULTS['no_mhd_amr_capacity'] = {
        'n_created': float(_parse_amr_created(amr_output)),
    }

    base = 'pic_no_mhd_boris_periodic_upper_wrap'
    _remove_outputs(base)
    _run_command(
        'periodic_upper_exact_wrap',
        1,
        [
            'job/basename=' + base,
            'time/cfl_number=0.5',
            'time/nlim=1',
            'particles/cr_vx0=100.0',
            'particles/cr_vy0=0.0',
            'particles/cr_vz0=0.0',
            'particles/pic_max_cell_cross=1',
            'output1/dcycle=0',
            'output2/dcycle=0',
            'output3/dcycle=0',
            'output4/dcycle=0',
            'output6/dcycle=1',
        ],
    )
    _PVTK_RESULTS['periodic_upper_exact_wrap'] = _read_pvtk_snapshot(base)
    _PDENS_RESULTS['periodic_upper_exact_wrap'] = _measure_pdens_summary(base)

    base = 'pic_no_mhd_boris_2d3v_pvtk'
    _remove_outputs(base)
    _run_command(
        'pvtk_2d3v_velocity_z',
        1,
        [
            'job/basename=' + base,
            'mesh/nx3=1',
            'mesh/x3max=1.0',
            'meshblock/nx3=1',
            'particles/pic_enable_2d3v=true',
            'time/nlim=1',
            'output1/dcycle=0',
            'output2/dcycle=0',
            'output3/dcycle=0',
            'output4/dcycle=0',
            'output5/dcycle=0',
            'output6/dcycle=1',
        ],
    )
    _PVTK_RESULTS['2d3v_velocity_z'] = _read_pvtk_snapshot(base)

    if _athena_mpi_enabled():
        base = 'pic_no_mhd_boris_pvtk_low_ppc'
        _remove_outputs(base)
        _run_command(
            'pvtk_zero_particle_rank',
            2,
            [
                'job/basename=' + base,
                'particles/ppc=0.001',
                'time/nlim=1',
                'output1/dcycle=0',
                'output2/dcycle=0',
                'output3/dcycle=0',
                'output4/dcycle=0',
                'output5/dcycle=0',
                'output6/dcycle=1',
            ],
        )
        _PVTK_RESULTS['zero_particle_rank'] = _read_pvtk_snapshot(base)

    _run_command(
        'guard_boris_requires_field_carrier',
        1,
        ['particles/pic_background_mode=coupled', 'time/nlim=1'],
        expect_fail=True,
        expected_message='Boris pushers require MHD fields, or '
                         '<particles>/pic_background_mode=no_mhd',
    )
    _run_command(
        'guard_no_mhd_requires_test_particle',
        1,
        ['particles/pic_feedback_mode=coupled', 'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_background_mode=no_mhd requires '
                         '<particles>/pic_feedback_mode=test_particle',
    )
    _run_command(
        'guard_no_mhd_rejects_coupling_toggles',
        1,
        ['particles/couple_moments_to_mhd=true', 'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_feedback_mode=test_particle does not support '
                         'particle-to-MHD coupling toggles',
    )
    _run_command(
        'guard_pic_max_cell_cross_meshblock_bound',
        1,
        ['particles/pic_max_cell_cross=5', 'time/nlim=0'],
        expect_fail=True,
        expected_message='pic_max_cell_cross must not exceed the smallest '
                         'active MeshBlock dimension',
    )
    _run_command(
        'guard_2d_boris_requires_2d3v',
        1,
        ['mesh/nx3=1',
         'mesh/x3max=1.0',
         'meshblock/nx3=1',
         'time/nlim=0'],
        expect_fail=True,
        expected_message='Boris pushers in 2D require '
                         '<particles>/pic_enable_2d3v=true',
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name in ['serial', 'mpi2']:
        if case_name not in _POSITIVE_RESULTS:
            continue
        result = _POSITIVE_RESULTS[case_name]
        measured = result['measured']
        expected = result['expected']
        ok = _check_with_tolerance(case_name + ':Q', measured['Q'], expected['Q'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance(case_name + ':Jx', measured['Jx'], expected['Jx'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':Jy', measured['Jy'], expected['Jy'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':Jz', measured['Jz'], expected['Jz'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance(case_name + ':npart', measured['npart'],
                                   expected['npart'], 1.0e-6, 1.0e-8) and ok

    if 'mpi2' in _POSITIVE_RESULTS:
        serial = _POSITIVE_RESULTS['serial']['measured']
        mpi2 = _POSITIVE_RESULTS['mpi2']['measured']
        ok = _check_with_tolerance('serial_vs_mpi2:Q', mpi2['Q'], serial['Q'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jx', mpi2['Jx'], serial['Jx'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jy', mpi2['Jy'], serial['Jy'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:Jz', mpi2['Jz'], serial['Jz'],
                                   _MOMENT_TOTAL_ABS_TOL,
                                   _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_with_tolerance('serial_vs_mpi2:npart',
                                   mpi2['npart'], serial['npart'],
                                   1.0e-6, 1.0e-8) and ok

    for dt_label in ['cell_cross_dt_limit', 'gyro_angle_dt_limit']:
        dt_result = _TIMESTEP_RESULTS.get(dt_label)
        if dt_result is None:
            logger.warning('Missing %s result', dt_label)
            ok = False
        else:
            ok = _check_with_tolerance(
                dt_label + ':dt',
                dt_result['measured_dt'], dt_result['expected_dt'],
                1.0e-8, 1.0e-6) and ok

    amr_result = _AMR_RESULTS.get('no_mhd_amr_capacity')
    if amr_result is None:
        logger.warning('Missing no_mhd_amr_capacity result')
        ok = False
    else:
        ok = (amr_result['n_created'] >= 1.0) and ok

    if _PVTK_RESULTS:
        wrap = _PVTK_RESULTS['periodic_upper_exact_wrap']
        x1 = wrap['points'][:, 0]
        n_wrapped = int(np.count_nonzero(np.isclose(x1, 0.0, atol=1.0e-6)))
        logger.info('periodic_upper_exact_wrap path=%s npoint=%d '
                    'x1_min=% .8e x1_max=% .8e n_wrapped=%d',
                    wrap['path'], wrap['npoint'], float(np.min(x1)),
                    float(np.max(x1)), n_wrapped)
        ok = (wrap['npoint'] == 1024) and ok
        ok = bool(np.all(np.isfinite(wrap['points']))) and ok
        ok = bool(np.all(x1 >= -1.0e-6)) and ok
        ok = bool(np.all(x1 < 16.0 - 1.0e-6)) and ok
        ok = (n_wrapped == 64) and ok
        ok = bool(np.all(wrap['gid'] >= 0)) and ok
        ok = bool(np.all(wrap['species'] == 0)) and ok
        pdens = _PDENS_RESULTS.get('periodic_upper_exact_wrap')
        if pdens is None:
            logger.warning('Missing periodic_upper_exact_wrap particle-density result')
            ok = False
        else:
            logger.info('periodic_upper_exact_wrap pdens npart=% .8e max=% .8e',
                        pdens['npart'], pdens['max'])
            ok = pdens['finite'] and ok
            ok = _check_with_tolerance('periodic_upper_exact_wrap:pdens_npart',
                                       pdens['npart'], 1024.0,
                                       1.0e-6, 1.0e-8) and ok

    if '2d3v_velocity_z' in _PVTK_RESULTS:
        pvtk_2d3v = _PVTK_RESULTS['2d3v_velocity_z']
        expected_vel = np.array([0.4, -0.2, 0.1])
        logger.info('pvtk_2d3v_velocity_z path=%s npoint=%d '
                    'vz_min=% .8e vz_max=% .8e',
                    pvtk_2d3v['path'], pvtk_2d3v['npoint'],
                    float(np.min(pvtk_2d3v['velocity'][:, 2])),
                    float(np.max(pvtk_2d3v['velocity'][:, 2])))
        ok = (pvtk_2d3v['npoint'] == 128) and ok
        ok = bool(np.all(np.isfinite(pvtk_2d3v['velocity']))) and ok
        ok = bool(np.allclose(pvtk_2d3v['points'][:, 2], 0.0,
                              atol=1.0e-6)) and ok
        ok = bool(np.allclose(pvtk_2d3v['velocity'], expected_vel,
                              atol=1.0e-6)) and ok

    if 'zero_particle_rank' in _PVTK_RESULTS:
        pvtk = _PVTK_RESULTS['zero_particle_rank']
        logger.info('pvtk_zero_particle_rank path=%s npoint=%d',
                    pvtk['path'], pvtk['npoint'])
        ok = (pvtk['npoint'] == 1) and ok
        ok = bool(np.all(np.isfinite(pvtk['points']))) and ok
        ok = bool(np.all(np.isfinite(pvtk['velocity']))) and ok
        ok = bool(np.all((pvtk['points'] >= np.array([0.0, 0.0, 0.0])) &
                         (pvtk['points'] <= np.array([16.0, 8.0, 8.0])))) and ok
        ok = bool(np.array_equal(pvtk['ptag'], np.array([0], dtype=np.int64))) and ok
        ok = bool(np.array_equal(pvtk['species'], np.array([0], dtype=np.int64))) and ok
        ok = bool(np.all((pvtk['gid'] >= 0) & (pvtk['gid'] < 16))) and ok
        ok = bool(np.allclose(pvtk['velocity'], np.array([[0.4, -0.2, 0.1]]),
                              atol=1.0e-6)) and ok

    ok = len(_NEGATIVE_RESULTS) == 5 and ok
    if len(_NEGATIVE_RESULTS) != 5:
        logger.warning('Expected 5 negative checks, got %d', len(_NEGATIVE_RESULTS))

    return ok
