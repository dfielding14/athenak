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

_INPUT_DECK = 'tests/pic_mhd_restart_fidelity.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}
_REPRESENTATIONS = [
    ('cc', 'cell_centered', None),
    ('edge', 'edge_staggered', 'cc_convert'),
    ('edge_direct', 'edge_staggered', 'direct_staggered'),
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _athena_mpi_enabled():
    proc = subprocess.run(['./athena', '-c'], cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
            os.path.join(exe_dir, 'bin', basename + '.*.bin'),
            os.path.join(exe_dir, 'rst', basename + '.*.rst'),
            os.path.join(exe_dir, 'rst', 'rank_*', basename + '.*.rst')]:
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
    jx_edge_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jx_edge'))
    jy_edge_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jy_edge'))
    jz_edge_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jz_edge'))

    return {
        'Q': _integrate_quantity(rho_data, 'prtcl_rho'),
        'Jx': _integrate_quantity(jx_data, 'prtcl_jx'),
        'Jy': _integrate_quantity(jy_data, 'prtcl_jy'),
        'Jz': _integrate_quantity(jz_data, 'prtcl_jz'),
        'Jx_edge': _integrate_quantity(jx_edge_data, 'prtcl_jx_edge'),
        'Jy_edge': _integrate_quantity(jy_edge_data, 'prtcl_jy_edge'),
        'Jz_edge': _integrate_quantity(jz_edge_data, 'prtcl_jz_edge'),
        'jx_edge_l2': _l2_quantity(jx_edge_data, 'prtcl_jx_edge'),
        'jy_edge_l2': _l2_quantity(jy_edge_data, 'prtcl_jy_edge'),
        'jz_edge_l2': _l2_quantity(jz_edge_data, 'prtcl_jz_edge'),
        'npart': float(np.sum(pdens_data['pdens'])),
        'bcc1_l2': _l2_quantity(bcc_data, 'bcc1'),
        'bcc2_l2': _l2_quantity(bcc_data, 'bcc2'),
        'bcc3_l2': _l2_quantity(bcc_data, 'bcc3'),
        'mom1': _integrate_quantity(m1_data, 'mom1'),
        'mom2': _integrate_quantity(m2_data, 'mom2'),
        'mom3': _integrate_quantity(m3_data, 'mom3'),
        'ener': _integrate_quantity(e_data, 'ener'),
    }


def _run_athena(label, nproc, arguments, restart_file=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
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

    base_cases = [('serial', 1)]
    if _athena_mpi_enabled():
        base_cases.append(('mpi2', 2))
    else:
        logger.info('MPI disabled: running serial-only restart fidelity checks')

    for case_tag, nproc in base_cases:
        for rep_tag, rep_value, dep_mode in _REPRESENTATIONS:
            base = f'pic_mhd_rst_{case_tag}_{rep_tag}'
            base_full = base + '_full'
            base_seg = base + '_seg'
            base_rst = base + '_rst'

            for name in [base_full, base_seg, base_rst]:
                _remove_outputs(name)

            common_args = [
                'particles/couple_moments_to_mhd=true',
                'particles/couple_j_to_efield_representation=' + rep_value,
                'particles/cr_vx0=0.50',
                'particles/cr_vy0=-0.25',
                'particles/cr_vz0=0.125',
            ]
            if dep_mode is not None:
                common_args.append('particles/couple_j_deposition_mode=' + dep_mode)

            _run_athena(base + '_full_run', nproc,
                        ['job/basename=' + base_full, 'time/nlim=2'] + common_args)
            full_measured = _measure_case(base_full)

            _run_athena(base + '_seg_run', nproc,
                        ['job/basename=' + base_seg, 'time/nlim=1'] + common_args)

            rst_path = os.path.join('rst', base_seg + '.00000.rst')
            full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
            if not os.path.exists(full_rst_path):
                raise RuntimeError('Expected restart file not found: ' + full_rst_path)

            _run_athena(base + '_restart_run', nproc,
                        ['job/basename=' + base_rst,
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=' + rep_value] +
                        ([] if dep_mode is None else
                         ['particles/couple_j_deposition_mode=' + dep_mode]),
                        restart_file=rst_path)
            rst_measured = _measure_case(base_rst)

            _RESULTS[base] = {
                'full': full_measured,
                'restart': rst_measured,
            }

    for case_tag, nproc in base_cases:
        for rep_tag, rep_value, dep_mode in _REPRESENTATIONS:
            base = f'pic_mhd_rst_per_rank_{case_tag}_{rep_tag}'
            base_full = base + '_full'
            base_seg = base + '_seg'
            base_rst = base + '_rst'

            for name in [base_full, base_seg, base_rst]:
                _remove_outputs(name)

            common_args = [
                'particles/couple_moments_to_mhd=true',
                'particles/couple_j_to_efield_representation=' + rep_value,
                'particles/cr_vx0=0.50',
                'particles/cr_vy0=-0.25',
                'particles/cr_vz0=0.125',
                'output11/single_file_per_rank=true',
            ]
            if dep_mode is not None:
                common_args.append('particles/couple_j_deposition_mode=' + dep_mode)

            _run_athena(base + '_full_run', nproc,
                        ['job/basename=' + base_full, 'time/nlim=2'] + common_args)
            full_measured = _measure_case(base_full)

            _run_athena(base + '_seg_run', nproc,
                        ['job/basename=' + base_seg, 'time/nlim=1'] + common_args)

            rst_path = os.path.join('rst', 'rank_00000000', base_seg + '.00000.rst')
            full_rst_path = os.path.join(_athena_exe_dir(), rst_path)
            if not os.path.exists(full_rst_path):
                raise RuntimeError('Expected restart file not found: ' + full_rst_path)

            _run_athena(base + '_restart_run', nproc,
                        ['job/basename=' + base_rst,
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation=' + rep_value] +
                        ([] if dep_mode is None else
                         ['particles/couple_j_deposition_mode=' + dep_mode]),
                        restart_file=rst_path)
            rst_measured = _measure_case(base_rst)

            _RESULTS[base] = {
                'full': full_measured,
                'restart': rst_measured,
            }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    abs_tol = 1.0e-8
    rel_tol = 1.0e-8
    quantities = [
        'Q', 'Jx', 'Jy', 'Jz', 'npart', 'bcc1_l2', 'bcc2_l2', 'bcc3_l2',
        'mom1', 'mom2', 'mom3', 'ener',
        'Jx_edge', 'Jy_edge', 'Jz_edge',
        'jx_edge_l2', 'jy_edge_l2', 'jz_edge_l2'
    ]

    for case_name, case_data in _RESULTS.items():
        full = case_data['full']
        rst = case_data['restart']
        for quantity in quantities:
            ok = _check_with_tolerance(case_name + ':' + quantity,
                                       rst[quantity], full[quantity],
                                       abs_tol, rel_tol) and ok

    rep_groups = {}
    for case_name, case_data in _RESULTS.items():
        prefix = 'pic_mhd_rst_'
        if not case_name.startswith(prefix):
            continue
        suffix = case_name[len(prefix):]
        if suffix.startswith('serial_'):
            case_tag = 'serial'
            rep_tag = suffix[len('serial_'):]
        elif suffix.startswith('mpi2_'):
            case_tag = 'mpi2'
            rep_tag = suffix[len('mpi2_'):]
        else:
            continue
        rep_groups.setdefault(case_tag, {})[rep_tag] = case_data

    for case_tag, reps in rep_groups.items():
        if 'edge' not in reps or 'edge_direct' not in reps:
            continue
        edge = reps['edge']
        direct = reps['edge_direct']

        for state in ['full', 'restart']:
            edge_state = edge[state]
            direct_state = direct[state]
            for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart',
                             'Jx_edge', 'Jy_edge', 'Jz_edge']:
                ok = _check_with_tolerance(
                    case_tag + ':' + state + ':edge_vs_edge_direct:' + quantity,
                    direct_state[quantity], edge_state[quantity],
                    1.0e-6, 1.0e-8) and ok

        for quantity in ['bcc1_l2', 'bcc2_l2', 'bcc3_l2',
                         'mom1', 'mom2', 'mom3', 'ener',
                         'jx_edge_l2', 'jy_edge_l2', 'jz_edge_l2']:
            edge_full = edge['full'][quantity]
            edge_rst = edge['restart'][quantity]
            direct_full = direct['full'][quantity]
            direct_rst = direct['restart'][quantity]
            if max(abs(edge_full), abs(edge_rst)) <= 1.0e-12:
                logger.info('%s:%s edge baseline near zero; skipping direct/edge ratio',
                            case_tag, quantity)
                ok = _check_with_tolerance(
                    case_tag + ':direct_restart_vs_full:' + quantity,
                    direct_rst, direct_full, abs_tol, rel_tol) and ok
                continue
            full_ratio = direct_full / max(abs(edge_full), 1.0e-300)
            rst_ratio = direct_rst / max(abs(edge_rst), 1.0e-300)
            logger.info('%s:%s edge_direct_vs_edge ratio full=% .8e restart=% .8e',
                        case_tag, quantity, full_ratio, rst_ratio)
            if not np.isfinite(full_ratio) or not np.isfinite(rst_ratio):
                logger.error('%s:%s non-finite edge_direct/edge ratio',
                             case_tag, quantity)
                ok = False
                continue
            ok = _check_with_tolerance(
                case_tag + ':ratio_restart_vs_full:edge_direct_vs_edge:' + quantity,
                rst_ratio, full_ratio, 1.0e-8, 1.0e-8) and ok

    return ok
