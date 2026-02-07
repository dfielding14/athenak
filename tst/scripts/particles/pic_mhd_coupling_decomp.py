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
_CONTINUITY_RESULTS = {}
_MESHBLOCK_CONFIGS = [
    ('mb444', 4, 4, 4),
    ('mb844', 8, 4, 4),
]
_REPRESENTATION_CASES = [
    ('cc', 'cell_centered', None),
    ('edge', 'edge_staggered', 'cc_convert'),
    ('edge_direct', 'edge_staggered', 'direct_staggered'),
]
_CC_CURRENT_IDS = ('prtcl_jx', 'prtcl_jy', 'prtcl_jz')
_EDGE_CURRENT_IDS = ('prtcl_jx_edge', 'prtcl_jy_edge', 'prtcl_jz_edge')


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


def _output_files(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No output files found for pattern: ' + pattern)
    return matches


def _volume_weights(dataset):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    return dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]


def _weighted_l2(field, dvol):
    return float(np.sqrt(np.sum(field * field * dvol)))


def _compute_divergence_from_cc_current(jx, jy, jz, x1f, x2f, x3f):
    div = np.zeros_like(jx)
    dx1 = np.diff(x1f)
    dx2 = np.diff(x2f)
    dx3 = np.diff(x3f)

    if dx1.size > 1:
        jx_face = 0.5 * (jx + np.roll(jx, -1, axis=2))
        div += (jx_face - np.roll(jx_face, 1, axis=2)) / dx1[None, None, :]
    if dx2.size > 1:
        jy_face = 0.5 * (jy + np.roll(jy, -1, axis=1))
        div += (jy_face - np.roll(jy_face, 1, axis=1)) / dx2[None, :, None]
    if dx3.size > 1:
        jz_face = 0.5 * (jz + np.roll(jz, -1, axis=0))
        div += (jz_face - np.roll(jz_face, 1, axis=0)) / dx3[:, None, None]

    return div


def _measure_continuity_residual(basename, current_ids):
    jx_id, jy_id, jz_id = current_ids
    rho_files = _output_files(basename, 'prtcl_rho')
    jx_files = _output_files(basename, jx_id)
    jy_files = _output_files(basename, jy_id)
    jz_files = _output_files(basename, jz_id)

    nfiles = min(len(rho_files), len(jx_files), len(jy_files), len(jz_files))
    if nfiles < 2:
        raise RuntimeError('Continuity residual requires at least two outputs: ' +
                           basename)

    curr_idx = nfiles - 1
    rho_curr = bin_convert.read_binary_as_athdf(rho_files[curr_idx])

    rho_prev = None
    prev_idx = curr_idx - 1
    while prev_idx >= 0:
        candidate = bin_convert.read_binary_as_athdf(rho_files[prev_idx])
        if (float(rho_curr['NumCycles']) > float(candidate['NumCycles']) or
                float(rho_curr['Time']) > float(candidate['Time']) + 1.0e-30):
            rho_prev = candidate
            break
        prev_idx -= 1
    if rho_prev is None:
        raise RuntimeError('Unable to find previous output with earlier time: ' +
                           basename)

    jx_curr = bin_convert.read_binary_as_athdf(jx_files[curr_idx])
    jy_curr = bin_convert.read_binary_as_athdf(jy_files[curr_idx])
    jz_curr = bin_convert.read_binary_as_athdf(jz_files[curr_idx])

    dt = float(rho_curr['Time'] - rho_prev['Time'])
    if (not np.isfinite(dt)) or dt <= 0.0:
        raise RuntimeError('Non-positive dt in continuity residual for ' + basename)

    drhodt = (rho_curr['prtcl_rho'] - rho_prev['prtcl_rho']) / dt
    divj = _compute_divergence_from_cc_current(jx_curr[jx_id],
                                               jy_curr[jy_id],
                                               jz_curr[jz_id],
                                               rho_curr['x1f'],
                                               rho_curr['x2f'],
                                               rho_curr['x3f'])
    residual = drhodt + divj
    dvol = _volume_weights(rho_curr)

    l2_drhodt = _weighted_l2(drhodt, dvol)
    l2_divj = _weighted_l2(divj, dvol)
    l2_residual = _weighted_l2(residual, dvol)
    linf_residual = float(np.max(np.abs(residual)))
    rel_residual = l2_residual / max(l2_drhodt + l2_divj, 1.0e-30)
    cycle_delta = float(rho_curr['NumCycles'] - rho_prev['NumCycles'])

    return {
        'dt': dt,
        'cycle_delta': cycle_delta,
        'l2_drhodt': l2_drhodt,
        'l2_divj': l2_divj,
        'l2_residual': l2_residual,
        'linf_residual': linf_residual,
        'rel_residual': rel_residual,
    }


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
    mpi_enabled = _athena_mpi_enabled()

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

    base_cases = []
    for mb_tag, mbx1, mbx2, mbx3 in _MESHBLOCK_CONFIGS:
        base_cases.append({
            'name': 'serial_' + mb_tag,
            'basename': 'pic_mhd_decomp_serial_' + mb_tag,
            'nproc': 1,
            'args': [
                'job/basename=pic_mhd_decomp_serial_' + mb_tag,
                'meshblock/nx1=' + str(mbx1),
                'meshblock/nx2=' + str(mbx2),
                'meshblock/nx3=' + str(mbx3),
            ],
        })

    if mpi_enabled:
        for mb_tag, mbx1, mbx2, mbx3 in _MESHBLOCK_CONFIGS:
            base_cases.extend([
                {
                    'name': 'mpi2_' + mb_tag,
                    'basename': 'pic_mhd_decomp_mpi2_' + mb_tag,
                    'nproc': 2,
                    'args': [
                        'job/basename=pic_mhd_decomp_mpi2_' + mb_tag,
                        'meshblock/nx1=' + str(mbx1),
                        'meshblock/nx2=' + str(mbx2),
                        'meshblock/nx3=' + str(mbx3),
                    ],
                },
                {
                    'name': 'mpi4_' + mb_tag,
                    'basename': 'pic_mhd_decomp_mpi4_' + mb_tag,
                    'nproc': 4,
                    'args': [
                        'job/basename=pic_mhd_decomp_mpi4_' + mb_tag,
                        'meshblock/nx1=' + str(mbx1),
                        'meshblock/nx2=' + str(mbx2),
                        'meshblock/nx3=' + str(mbx3),
                    ],
                },
            ])
    else:
        logger.info('MPI disabled: running serial-only decomposition baseline')

    cases = []
    feedback_orders = [('mhd', 'mhd_src_terms'), ('efield', 'efield_src')]
    for order_tag, order_value in feedback_orders:
        for rep_tag, rep_value, dep_mode in _REPRESENTATION_CASES:
            for base_case in base_cases:
                combo_tag = rep_tag + '_' + order_tag
                rep_basename = base_case['basename'] + '_' + combo_tag
                rep_args = ['particles/couple_j_to_efield_representation=' + rep_value,
                            'particles/couple_fluid_feedback_order=' + order_value]
                if dep_mode is not None:
                    rep_args.append('particles/couple_j_deposition_mode=' + dep_mode)
                cases.append({
                    'name': base_case['name'] + '_' + combo_tag,
                    'basename': rep_basename,
                    'nproc': base_case['nproc'],
                    'args': (['job/basename=' + rep_basename] +
                             list(base_case['args'][1:]) +
                             rep_args +
                             common_args),
                })

    for case in cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(case['args']),
        }

    continuity_cases = []
    cont_mb_tag, cont_mbx1, cont_mbx2, cont_mbx3 = _MESHBLOCK_CONFIGS[0]
    for rep_tag, rep_value, dep_mode in _REPRESENTATION_CASES:
        rep_args = ['particles/couple_j_to_efield_representation=' + rep_value]
        if dep_mode is not None:
            rep_args.append('particles/couple_j_deposition_mode=' + dep_mode)

        serial_basename = ('pic_mhd_decomp_cont_serial_' + cont_mb_tag +
                           '_' + rep_tag)
        continuity_cases.append({
            'name': 'serial_continuity_' + rep_tag,
            'basename': serial_basename,
            'nproc': 1,
            'current_ids': (_CC_CURRENT_IDS if rep_tag == 'cc'
                            else _EDGE_CURRENT_IDS),
            'args': ['job/basename=' + serial_basename,
                     'time/nlim=2',
                     'meshblock/nx1=' + str(cont_mbx1),
                     'meshblock/nx2=' + str(cont_mbx2),
                     'meshblock/nx3=' + str(cont_mbx3)] + rep_args + common_args,
        })

    if mpi_enabled:
        for nproc, mpi_tag in [(2, 'mpi2'), (4, 'mpi4')]:
            for rep_tag, rep_value, dep_mode in _REPRESENTATION_CASES:
                rep_args = ['particles/couple_j_to_efield_representation=' + rep_value]
                if dep_mode is not None:
                    rep_args.append('particles/couple_j_deposition_mode=' + dep_mode)

                mpi_basename = ('pic_mhd_decomp_cont_' + mpi_tag + '_' + cont_mb_tag +
                                '_' + rep_tag)
                continuity_cases.append({
                    'name': mpi_tag + '_continuity_' + rep_tag,
                    'basename': mpi_basename,
                    'nproc': nproc,
                    'current_ids': (_CC_CURRENT_IDS if rep_tag == 'cc'
                                    else _EDGE_CURRENT_IDS),
                    'args': (
                        ['job/basename=' + mpi_basename,
                         'time/nlim=2',
                         'meshblock/nx1=' + str(cont_mbx1),
                         'meshblock/nx2=' + str(cont_mbx2),
                         'meshblock/nx3=' + str(cont_mbx3)] +
                        rep_args + common_args
                    ),
                })

    for case in continuity_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _CONTINUITY_RESULTS[case['name']] = (
            _measure_continuity_residual(case['basename'], case['current_ids']))


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

    for mb_tag, _, _, _ in _MESHBLOCK_CONFIGS:
        for order_tag in ['mhd', 'efield']:
            for rep_tag, _, _ in _REPRESENTATION_CASES:
                serial_case = 'serial_' + mb_tag + '_' + rep_tag + '_' + order_tag
                if serial_case not in _RESULTS:
                    logger.warning('Missing serial baseline result for %s %s %s',
                                   mb_tag, rep_tag, order_tag)
                    return False

                baseline = _RESULTS[serial_case]['measured']
                for case_name in ['mpi2_' + mb_tag + '_' + rep_tag + '_' + order_tag,
                                  'mpi4_' + mb_tag + '_' + rep_tag + '_' + order_tag]:
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

    for mb_tag, _, _, _ in _MESHBLOCK_CONFIGS:
        for order_tag in ['mhd', 'efield']:
            serial_cc = 'serial_' + mb_tag + '_cc_' + order_tag
            serial_edge = 'serial_' + mb_tag + '_edge_' + order_tag
            serial_edge_direct = 'serial_' + mb_tag + '_edge_direct_' + order_tag
            if serial_cc in _RESULTS and serial_edge in _RESULTS:
                cc = _RESULTS[serial_cc]['measured']
                edge = _RESULTS[serial_edge]['measured']
                for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                    ok = _check_with_tolerance(mb_tag + ':' + order_tag +
                                               ':cc_vs_edge:' + quantity,
                                               cc[quantity], edge[quantity],
                                               1.0e-6, 1.0e-8) and ok
            if serial_edge in _RESULTS and serial_edge_direct in _RESULTS:
                edge = _RESULTS[serial_edge]['measured']
                edge_direct = _RESULTS[serial_edge_direct]['measured']
                for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                    ok = _check_with_tolerance(mb_tag + ':' + order_tag +
                                               ':edge_vs_edge_direct:' + quantity,
                                               edge[quantity], edge_direct[quantity],
                                               1.0e-6, 1.0e-8) and ok

    for mb_tag, _, _, _ in _MESHBLOCK_CONFIGS:
        for rep_tag, _, _ in _REPRESENTATION_CASES:
            serial_mhd = 'serial_' + mb_tag + '_' + rep_tag + '_mhd'
            serial_efield = 'serial_' + mb_tag + '_' + rep_tag + '_efield'
            if serial_mhd in _RESULTS and serial_efield in _RESULTS:
                mhd = _RESULTS[serial_mhd]['measured']
                efield = _RESULTS[serial_efield]['measured']
                for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
                    ok = _check_with_tolerance(mb_tag + ':' + rep_tag +
                                               ':mhd_vs_efield:' + quantity,
                                               mhd[quantity], efield[quantity],
                                               1.0e-6, 1.0e-8) and ok

    required_cont_cases = [
        'serial_continuity_cc',
        'serial_continuity_edge',
        'serial_continuity_edge_direct',
    ]
    for case_name in required_cont_cases:
        if case_name not in _CONTINUITY_RESULTS:
            logger.warning('Missing continuity-oracle case %s', case_name)
            return False

    cont_cc = _CONTINUITY_RESULTS['serial_continuity_cc']
    cont_edge = _CONTINUITY_RESULTS['serial_continuity_edge']
    cont_direct = _CONTINUITY_RESULTS['serial_continuity_edge_direct']
    continuity_sets = [('cc', cont_cc), ('edge', cont_edge),
                       ('edge_direct', cont_direct)]
    for label, cont in continuity_sets:
        logger.info('continuity_%s dt=% .8e cycle_delta=% .8e l2_drhodt=% .8e '
                    'l2_divj=% .8e l2_res=% .8e linf_res=% .8e rel=% .8e',
                    label, cont['dt'], cont['cycle_delta'], cont['l2_drhodt'],
                    cont['l2_divj'], cont['l2_residual'], cont['linf_residual'],
                    cont['rel_residual'])
        if (not np.isfinite(cont['l2_residual']) or
                not np.isfinite(cont['linf_residual']) or
                not np.isfinite(cont['rel_residual'])):
            logger.warning('Non-finite continuity metric for %s', label)
            ok = False
        if cont['dt'] <= 0.0 or cont['cycle_delta'] < 1.0:
            logger.warning('Invalid continuity time metadata for %s', label)
            ok = False

    if cont_edge['rel_residual'] > 0.0:
        direct_vs_edge_rel = cont_direct['rel_residual'] / cont_edge['rel_residual']
        logger.info('continuity_direct_vs_edge_rel_ratio=% .8e', direct_vs_edge_rel)
        if direct_vs_edge_rel > 1.25:
            logger.warning('Direct-edge continuity residual regressed vs edge mode')
            ok = False
    else:
        ok = _check_with_tolerance('continuity_edge_rel_zero',
                                   cont_direct['rel_residual'], 0.0,
                                   1.0e-12, 1.0e-12) and ok

    for mpi_tag in ['mpi2', 'mpi4']:
        if all((mpi_tag + '_continuity_' + rep_tag) in _CONTINUITY_RESULTS
               for rep_tag in ['cc', 'edge', 'edge_direct']):
            for rep_tag in ['cc', 'edge', 'edge_direct']:
                serial = _CONTINUITY_RESULTS['serial_continuity_' + rep_tag]
                mpi = _CONTINUITY_RESULTS[mpi_tag + '_continuity_' + rep_tag]
                ok = _check_with_tolerance('continuity_serial_vs_' + mpi_tag + ':' +
                                           rep_tag + ':l2_residual',
                                           mpi['l2_residual'], serial['l2_residual'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance('continuity_serial_vs_' + mpi_tag + ':' +
                                           rep_tag + ':linf_residual',
                                           mpi['linf_residual'],
                                           serial['linf_residual'],
                                           1.0e-6, 1.0e-8) and ok
                ok = _check_with_tolerance('continuity_serial_vs_' + mpi_tag + ':' +
                                           rep_tag + ':rel_residual',
                                           mpi['rel_residual'],
                                           serial['rel_residual'],
                                           1.0e-6, 1.0e-8) and ok

    return ok
