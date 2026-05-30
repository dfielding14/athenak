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
_DEFAULT_MODE_INPUT_DECK = 'tests/pic_mhd_current_coupling_default_mode.athinput'
_BORIS_INPUT_DECK = 'tests/pic_boris_midpoint_eb.athinput'
_UNSUPPORTED_INPUT_DECK = 'tests/pic_deposit_conservation.athinput'
_RADIATION_GUARD_DECK = 'tests/pic_mhd_coupling_guard_radiation.athinput'
_NR_GUARD_DECK = 'tests/pic_mhd_coupling_guard_nr.athinput'
_HYDRO_GUARD_DECK = 'tests/pic_mhd_coupling_guard_hydro.athinput'
_ISOTHERMAL_GUARD_DECK = 'tests/pic_mhd_coupling_guard_isothermal.athinput'
_EDGE_REPR_REL_GUARD_DECK = 'tests/pic_mhd_coupling_guard_edge_relativistic.athinput'
_NONPERIODIC_UNSUPPORTED_GUARD_DECK = (
    'tests/pic_mhd_coupling_guard_nonperiodic_unsupported.athinput')
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_POSITIVE_RESULTS = {}
_NEGATIVE_RESULTS = {}
_PROJECTED_CURRENT_RESULTS = {}
_CC_CURRENT_IDS = ('prtcl_jx', 'prtcl_jy', 'prtcl_jz')
# Edge-current outputs are cell-centered projections of edge-centered E-field
# source arrays. They are useful output diagnostics, not raw staggered data.
_EDGE_CURRENT_IDS = ('prtcl_jx_edge', 'prtcl_jy_edge', 'prtcl_jz_edge')
_MOMENT_TOTAL_ABS_TOL = 2.5e-4
_MOMENT_TOTAL_REL_TOL = 3.0e-7
_EDGE_MOMENT_TOTAL_ABS_TOL = 1.0e-4
_COUNT_ABS_TOL = 1.0e-6
_COUNT_REL_TOL = 1.0e-8


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


def _resolved_config(arguments, input_deck=_INPUT_DECK):
    config = bin_convert.athinput(_deck_path_for_python(input_deck))
    _apply_overrides(config, arguments)
    return config


def _expected_totals(arguments, input_deck=_INPUT_DECK):
    config = _resolved_config(arguments, input_deck)

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


def _edge_projection_matches_moment_integral(arguments, input_deck=_INPUT_DECK):
    config = _resolved_config(arguments, input_deck)
    particles = config['particles']
    representation = particles['couple_j_to_efield_representation']
    return representation != 'edge_staggered'


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


def _latest_output_file_optional(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        return None
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
    """Diagnostic divergence of cell-centered or projected cell-centered J."""
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


def _measure_projected_current_residual(basename, current_ids):
    # This is a serial/MPI parity diagnostic for output projections. It is not a
    # discrete charge-conservation oracle for the raw edge-centered arrays.
    jx_id, jy_id, jz_id = current_ids
    rho_files = _output_files(basename, 'prtcl_rho')
    jx_files = _output_files(basename, jx_id)
    jy_files = _output_files(basename, jy_id)
    jz_files = _output_files(basename, jz_id)

    nfiles = min(len(rho_files), len(jx_files), len(jy_files), len(jz_files))
    if nfiles < 2:
        raise RuntimeError('Projected-current residual requires at least two '
                           'outputs: ' + basename)

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
        raise RuntimeError('Non-positive dt in projected-current residual for ' +
                           basename)

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


def _l2_difference_quantity(dataset_a, dataset_b, quantity):
    dx1 = np.diff(dataset_a['x1f'])
    dx2 = np.diff(dataset_a['x2f'])
    dx3 = np.diff(dataset_a['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    diff = dataset_a[quantity] - dataset_b[quantity]
    return float(np.sqrt(np.sum(diff * diff * dvol)))


def _l2_sum_difference_quantity(dataset_a, dataset_ref, dataset_b, quantity):
    dx1 = np.diff(dataset_a['x1f'])
    dx2 = np.diff(dataset_a['x2f'])
    dx3 = np.diff(dataset_a['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    sum_diff = (dataset_a[quantity] - dataset_ref[quantity] +
                dataset_b[quantity] - dataset_ref[quantity])
    return float(np.sqrt(np.sum(sum_diff * sum_diff * dvol)))


def _dot_difference_quantity(dataset_a, dataset_ref, dataset_b, quantity):
    dx1 = np.diff(dataset_a['x1f'])
    dx2 = np.diff(dataset_a['x2f'])
    dx3 = np.diff(dataset_a['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    diff_a = dataset_a[quantity] - dataset_ref[quantity]
    diff_b = dataset_b[quantity] - dataset_ref[quantity]
    return float(np.sum(diff_a * diff_b * dvol))


def _measure_dataset(basename, file_id):
    return bin_convert.read_binary_as_athdf(_latest_output_file(basename, file_id))


def _measure_bcc_dataset(basename):
    return _measure_dataset(basename, 'mhd_bcc')


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

    measured = {
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

    jx_edge_file = _latest_output_file_optional(basename, 'prtcl_jx_edge')
    jy_edge_file = _latest_output_file_optional(basename, 'prtcl_jy_edge')
    jz_edge_file = _latest_output_file_optional(basename, 'prtcl_jz_edge')
    if (jx_edge_file is not None and jy_edge_file is not None and
            jz_edge_file is not None):
        jx_edge_data = bin_convert.read_binary_as_athdf(jx_edge_file)
        jy_edge_data = bin_convert.read_binary_as_athdf(jy_edge_file)
        jz_edge_data = bin_convert.read_binary_as_athdf(jz_edge_file)
        measured.update({
            'Jx_edge': _integrate_quantity(jx_edge_data, 'prtcl_jx_edge'),
            'Jy_edge': _integrate_quantity(jy_edge_data, 'prtcl_jy_edge'),
            'Jz_edge': _integrate_quantity(jz_edge_data, 'prtcl_jz_edge'),
            'Jx_edge_l2': _l2_quantity(jx_edge_data, 'prtcl_jx_edge'),
            'Jy_edge_l2': _l2_quantity(jy_edge_data, 'prtcl_jy_edge'),
            'Jz_edge_l2': _l2_quantity(jz_edge_data, 'prtcl_jz_edge'),
        })

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


def _integrated_quantity_tolerances(quantity):
    if quantity == 'npart':
        return _COUNT_ABS_TOL, _COUNT_REL_TOL
    return _MOMENT_TOTAL_ABS_TOL, _MOMENT_TOTAL_REL_TOL


def _check_integrated_quantity(label, measured, expected, quantity):
    abs_tol, rel_tol = _integrated_quantity_tolerances(quantity)
    return _check_with_tolerance(label + ':' + quantity,
                                 measured[quantity], expected[quantity],
                                 abs_tol, rel_tol)


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    common_args = [
        'particles/cr_distribution=random',
        'particles/cr_vx0=0.50',
        'particles/cr_vy0=-0.25',
        'particles/cr_vz0=0.125',
    ]
    random_position_args = ['particles/cr_distribution=random']
    two_d_args = [
        'mesh/nx3=1',
        'mesh/x3max=1.0',
        'meshblock/nx3=1',
    ]
    boris_args = [
        'problem/vflow=0.3',
        'particles/pic_background_mode=coupled',
        'particles/pic_feedback_mode=coupled',
    ]

    positive_cases = [
        {
            'name': 'serial_default_optin_off',
            'basename': 'pic_mhd_default',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_default'] + common_args,
        },
        {
            'name': 'serial_explicit_optin_off',
            'basename': 'pic_mhd_explicit_off',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_explicit_off',
                     'particles/couple_moments_to_mhd=false'] + common_args,
        },
        {
            'name': 'serial_coupled',
            'basename': 'pic_mhd_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        },
        {
            'name': 'serial_coupled_edge_staggered',
            'basename': 'pic_mhd_coupled_edge',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_default_mode',
            'basename': 'pic_mhd_coupled_edge_default',
            'nproc': 1,
            'input_deck': _DEFAULT_MODE_INPUT_DECK,
            'args': ['job/basename=pic_mhd_coupled_edge_default',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_staggered',
            'basename': 'pic_mhd_coupled_edge_direct',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_staggered_o2',
            'basename': 'pic_mhd_coupled_edge_direct_o2',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_o2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + common_args,
        },
        {
            'name': 'serial_coupled_nonunit_volume',
            'basename': 'pic_mhd_coupled_nonunit_volume',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_nonunit_volume',
                     'mesh/x1max=32.0',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        },
        {
            'name': 'serial_coupled_edge_nonunit_volume',
            'basename': 'pic_mhd_coupled_edge_nonunit_volume',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_nonunit_volume',
                     'mesh/x1max=32.0',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_nonunit_volume',
            'basename': 'pic_mhd_coupled_edge_direct_nonunit_volume',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_nonunit_volume',
                     'mesh/x1max=32.0',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_o2_nonunit_volume',
            'basename': 'pic_mhd_coupled_edge_direct_o2_nonunit_volume',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_o2_nonunit_volume',
                     'mesh/x1max=32.0',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + common_args,
        },
        {
            'name': 'serial_coupled_edge_staggered_2d',
            'basename': 'pic_mhd_coupled_edge_2d',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_2d',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    two_d_args + common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_staggered_2d',
            'basename': 'pic_mhd_coupled_edge_direct_2d',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_2d',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    two_d_args + common_args,
        },
        {
            'name': 'serial_coupled_edge_direct_staggered_o2_2d',
            'basename': 'pic_mhd_coupled_edge_direct_o2_2d',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_o2_2d',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + two_d_args + common_args,
        },
        {
            'name': 'serial_feedback_momentum',
            'basename': 'pic_mhd_feedback_momentum',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_momentum',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=true',
                     'particles/couple_moments_momentum_coeff=1.0'] + common_args,
        },
        {
            'name': 'serial_feedback_energy',
            'basename': 'pic_mhd_feedback_energy',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_energy',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_energy_coeff=10.0'] + common_args,
        },
        {
            'name': 'serial_feedback_energy_rk1',
            'basename': 'pic_mhd_feedback_energy_rk1',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_energy_rk1',
                     'time/integrator=rk1',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_energy_coeff=10.0'] + common_args,
        },
        {
            'name': 'serial_feedback_energy_rk3',
            'basename': 'pic_mhd_feedback_energy_rk3',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_energy_rk3',
                     'time/integrator=rk3',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_energy_coeff=10.0'] + common_args,
        },
        {
            'name': 'serial_feedback_both',
            'basename': 'pic_mhd_feedback_both',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_both',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_momentum_coeff=1.0',
                     'particles/couple_moments_energy_coeff=10.0'] + common_args,
        },
        {
            'name': 'serial_feedback_both_efield_order',
            'basename': 'pic_mhd_feedback_both_efield',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_both_efield',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_momentum_coeff=1.0',
                     'particles/couple_moments_energy_coeff=10.0',
                     'particles/couple_fluid_feedback_order=efield_src'] +
                    common_args,
        },
        {
            'name': 'serial_feedback_both_efield_order_edge_staggered',
            'basename': 'pic_mhd_feedback_both_efield_edge',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_feedback_both_efield_edge',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert',
                     'particles/couple_moments_momentum_to_mhd=true',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_momentum_coeff=1.0',
                     'particles/couple_moments_energy_coeff=10.0',
                     'particles/couple_fluid_feedback_order=efield_src'] +
                    common_args,
        },
        {
            'name': 'serial_boris_coupled',
            'basename': 'pic_mhd_boris_coupled',
            'nproc': 1,
            'input_deck': _BORIS_INPUT_DECK,
            'args': ['job/basename=pic_mhd_boris_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=false',
                     'particles/couple_moments_energy_to_mhd=false'] +
                    boris_args,
        },
        {
            'name': 'serial_boris_feedback_energy',
            'basename': 'pic_mhd_boris_feedback_energy',
            'nproc': 1,
            'input_deck': _BORIS_INPUT_DECK,
            'args': ['job/basename=pic_mhd_boris_feedback_energy',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=false',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_energy_coeff=1.0'] +
                    boris_args,
        },
        {
            'name': 'serial_boris_feedback_energy_x2',
            'basename': 'pic_mhd_boris_feedback_energy_x2',
            'nproc': 1,
            'input_deck': _BORIS_INPUT_DECK,
            'args': ['job/basename=pic_mhd_boris_feedback_energy_x2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_moments_momentum_to_mhd=false',
                     'particles/couple_moments_energy_to_mhd=true',
                     'particles/couple_moments_energy_coeff=2.0'] +
                    boris_args,
        },
        {
            'name': 'serial_zero_uncoupled',
            'basename': 'pic_mhd_zero_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_zero_uncoupled',
                     'particles/couple_moments_to_mhd=false',
                     'particles/cr_vx0=0.0',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_zero_coupled',
            'basename': 'pic_mhd_zero_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_zero_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/cr_vx0=0.0',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_lin_lo_uncoupled',
            'basename': 'pic_mhd_lin_lo_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_lin_lo_uncoupled',
                     'particles/couple_moments_to_mhd=false',
                     'particles/deposit_qscale=1.0',
                     'particles/cr_vx0=0.01',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_lin_lo_coupled',
            'basename': 'pic_mhd_lin_lo_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_lin_lo_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/deposit_qscale=1.0',
                     'particles/cr_vx0=0.01',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_lin_lo_coupled_neg',
            'basename': 'pic_mhd_lin_lo_neg_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_lin_lo_neg_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/deposit_qscale=1.0',
                     'particles/couple_j_to_efield_coeff=-1.0',
                     'particles/cr_vx0=0.01',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_lin_hi_uncoupled',
            'basename': 'pic_mhd_lin_hi_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_lin_hi_uncoupled',
                     'particles/couple_moments_to_mhd=false',
                     'particles/deposit_qscale=2.0',
                     'particles/cr_vx0=0.01',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_lin_hi_coupled',
            'basename': 'pic_mhd_lin_hi_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_lin_hi_coupled',
                     'particles/couple_moments_to_mhd=true',
                     'particles/deposit_qscale=2.0',
                     'particles/cr_vx0=0.01',
                     'particles/cr_vy0=0.0',
                     'particles/cr_vz0=0.0'] + random_position_args,
        },
        {
            'name': 'serial_rk1_uncoupled',
            'basename': 'pic_mhd_rk1_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_rk1_uncoupled',
                     'time/integrator=rk1',
                     'particles/couple_moments_to_mhd=false'] + common_args,
        },
        {
            'name': 'serial_rk1_coupled',
            'basename': 'pic_mhd_rk1_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_rk1_coupled',
                     'time/integrator=rk1',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        },
        {
            'name': 'serial_rk3_uncoupled',
            'basename': 'pic_mhd_rk3_uncoupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_rk3_uncoupled',
                     'time/integrator=rk3',
                     'particles/couple_moments_to_mhd=false'] + common_args,
        },
        {
            'name': 'serial_rk3_coupled',
            'basename': 'pic_mhd_rk3_coupled',
            'nproc': 1,
            'args': ['job/basename=pic_mhd_rk3_coupled',
                     'time/integrator=rk3',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        },
    ]

    if _athena_mpi_enabled():
        positive_cases.append({
            'name': 'mpi2_coupled',
            'basename': 'pic_mhd_coupled_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_mhd_coupled_mpi2',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        })
        positive_cases.append({
            'name': 'mpi2_coupled_edge_staggered',
            'basename': 'pic_mhd_coupled_edge_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_mhd_coupled_edge_mpi2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    common_args,
        })
        positive_cases.append({
            'name': 'mpi2_coupled_edge_default_mode',
            'basename': 'pic_mhd_coupled_edge_default_mpi2',
            'nproc': 2,
            'input_deck': _DEFAULT_MODE_INPUT_DECK,
            'args': ['job/basename=pic_mhd_coupled_edge_default_mpi2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered'] +
                    common_args,
        })
        positive_cases.append({
            'name': 'mpi2_coupled_edge_direct_staggered',
            'basename': 'pic_mhd_coupled_edge_direct_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_mpi2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    common_args,
        })
        positive_cases.append({
            'name': 'mpi2_coupled_edge_direct_staggered_o2',
            'basename': 'pic_mhd_coupled_edge_direct_o2_mpi2',
            'nproc': 2,
            'args': ['job/basename=pic_mhd_coupled_edge_direct_o2_mpi2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + common_args,
        })
    else:
        logger.info('Skipping mpi2_coupled case: Athena build has MPI parallelism OFF')

    for case in positive_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'],
                     input_deck=case.get('input_deck', _INPUT_DECK))
        _POSITIVE_RESULTS[case['name']] = {
            'measured': _measure_case(case['basename']),
            'expected': _expected_totals(
                case['args'], case.get('input_deck', _INPUT_DECK)),
            'edge_projection_exact': _edge_projection_matches_moment_integral(
                case['args'], case.get('input_deck', _INPUT_DECK)),
        }

    projected_current_cases = [
        {
            'name': 'serial_projected_current_cc',
            'basename': 'pic_mhd_cont_cc',
            'nproc': 1,
            'current_ids': _CC_CURRENT_IDS,
            'args': ['job/basename=pic_mhd_cont_cc',
                     'time/nlim=2',
                     'particles/couple_moments_to_mhd=true'] + common_args,
        },
        {
            'name': 'serial_projected_current_edge',
            'basename': 'pic_mhd_cont_edge',
            'nproc': 1,
            'current_ids': _EDGE_CURRENT_IDS,
            'args': ['job/basename=pic_mhd_cont_edge',
                     'time/nlim=2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=cc_convert'] +
                    common_args,
        },
        {
            'name': 'serial_projected_current_edge_direct',
            'basename': 'pic_mhd_cont_edge_direct',
            'nproc': 1,
            'current_ids': _EDGE_CURRENT_IDS,
            'args': ['job/basename=pic_mhd_cont_edge_direct',
                     'time/nlim=2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered'] +
                    common_args,
        },
        {
            'name': 'serial_projected_current_edge_direct_o2',
            'basename': 'pic_mhd_cont_edge_direct_o2',
            'nproc': 1,
            'current_ids': _EDGE_CURRENT_IDS,
            'args': ['job/basename=pic_mhd_cont_edge_direct_o2',
                     'time/nlim=2',
                     'particles/couple_moments_to_mhd=true',
                     'particles/couple_j_to_efield_representation=edge_staggered',
                     'particles/couple_j_deposition_mode=direct_staggered',
                     'particles/deposit_order=2'] + common_args,
        },
    ]

    if _athena_mpi_enabled():
        projected_current_cases.extend([
            {
                'name': 'mpi2_projected_current_cc',
                'basename': 'pic_mhd_cont_cc_mpi2',
                'nproc': 2,
                'current_ids': _CC_CURRENT_IDS,
                'args': ['job/basename=pic_mhd_cont_cc_mpi2',
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true'] + common_args,
            },
            {
                'name': 'mpi2_projected_current_edge',
                'basename': 'pic_mhd_cont_edge_mpi2',
                'nproc': 2,
                'current_ids': _EDGE_CURRENT_IDS,
                'args': ['job/basename=pic_mhd_cont_edge_mpi2',
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation='
                         'edge_staggered',
                         'particles/couple_j_deposition_mode=cc_convert'] +
                        common_args,
            },
            {
                'name': 'mpi2_projected_current_edge_direct',
                'basename': 'pic_mhd_cont_edge_direct_mpi2',
                'nproc': 2,
                'current_ids': _EDGE_CURRENT_IDS,
                'args': ['job/basename=pic_mhd_cont_edge_direct_mpi2',
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation='
                         'edge_staggered',
                         'particles/couple_j_deposition_mode=direct_staggered'] +
                        common_args,
            },
            {
                'name': 'mpi2_projected_current_edge_direct_o2',
                'basename': 'pic_mhd_cont_edge_direct_o2_mpi2',
                'nproc': 2,
                'current_ids': _EDGE_CURRENT_IDS,
                'args': ['job/basename=pic_mhd_cont_edge_direct_o2_mpi2',
                         'time/nlim=2',
                         'particles/couple_moments_to_mhd=true',
                         'particles/couple_j_to_efield_representation='
                         'edge_staggered',
                         'particles/couple_j_deposition_mode=direct_staggered',
                         'particles/deposit_order=2'] + common_args,
            },
        ])

    for case in projected_current_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _PROJECTED_CURRENT_RESULTS[case['name']] = (
            _measure_projected_current_residual(case['basename'], case['current_ids']))

    _run_command(
        'guard_deposit_moments_false',
        1,
        ['particles/couple_moments_to_mhd=true',
         'particles/deposit_moments=false',
         'time/nlim=0'],
        expect_fail=True,
        expected_message='requires <particles>/deposit_moments=true',
    )

    _run_command(
        'guard_missing_mhd_composition',
        1,
        ['particles/couple_moments_to_mhd=true', 'time/nlim=0'],
        expect_fail=True,
        expected_message='requires an active <mhd> block in PR2',
        input_deck=_UNSUPPORTED_INPUT_DECK,
    )

    _run_command(
        'guard_radiation_mhd_composition',
        1,
        ['time/nlim=0', 'particles/couple_j_to_efield_representation=edge_staggered'],
        expect_fail=True,
        expected_message='does not support radiation+MHD compositions in PR2',
        input_deck=_RADIATION_GUARD_DECK,
    )

    _run_command(
        'guard_nr_mhd_composition',
        1,
        ['time/nlim=0', 'particles/couple_j_to_efield_representation=edge_staggered'],
        expect_fail=True,
        expected_message='does not support numerical relativity task paths in PR2',
        input_deck=_NR_GUARD_DECK,
    )

    _run_command(
        'guard_hydro_mhd_composition',
        1,
        ['time/nlim=0', 'particles/couple_j_to_efield_representation=edge_staggered'],
        expect_fail=True,
        expected_message='does not support ion-neutral/two-fluid compositions in PR2',
        input_deck=_HYDRO_GUARD_DECK,
    )

    _run_command(
        'guard_feedback_requires_coupled_mode',
        1,
        ['particles/couple_moments_to_mhd=false',
         'particles/couple_moments_momentum_to_mhd=true',
         'time/nlim=0'],
        expect_fail=True,
        expected_message='requires <particles>/couple_moments_to_mhd=true',
    )

    _run_command(
        'guard_energy_feedback_requires_ideal_eos',
        1,
        ['time/nlim=0'],
        expect_fail=True,
        expected_message='requires <mhd>/eos=ideal',
        input_deck=_ISOTHERMAL_GUARD_DECK,
    )

    _run_command(
        'guard_edge_representation_relativistic',
        1,
        ['time/nlim=0', 'particles/couple_j_to_efield_representation=edge_staggered'],
        expect_fail=True,
        expected_message='requires non-relativistic Cartesian MHD in PR2',
        input_deck=_EDGE_REPR_REL_GUARD_DECK,
    )

    _run_command(
        'guard_invalid_feedback_order',
        1,
        ['time/nlim=0', 'particles/couple_fluid_feedback_order=bad_order'],
        expect_fail=True,
        expected_message='Unsupported value for <particles>/'
                         'couple_fluid_feedback_order: bad_order',
    )

    _run_command(
        'guard_direct_staggered_requires_edge_representation',
        1,
        ['time/nlim=0',
         'particles/couple_moments_to_mhd=true',
         'particles/couple_j_to_efield_representation=cell_centered',
         'particles/couple_j_deposition_mode=direct_staggered'],
        expect_fail=True,
        expected_message=('requires '
                          '<particles>/couple_j_to_efield_representation='
                          'edge_staggered'),
    )

    _run_command(
        'guard_direct_staggered_large_step',
        1,
        ['job/basename=pic_mhd_direct_large_step',
         'time/cfl_number=1.0',
         'time/nlim=1',
         'particles/couple_moments_to_mhd=true',
         'particles/couple_j_to_efield_representation=edge_staggered',
         'particles/couple_j_deposition_mode=direct_staggered',
         'particles/cr_distribution=random',
         'particles/cr_vx0=4.0',
         'particles/cr_vy0=0.0',
         'particles/cr_vz0=0.0'],
        expect_fail=True,
        expected_message=('requires particle shape support to move by at most '
                          'one cell per step'),
    )

    _run_command(
        'guard_nonperiodic_unsupported_boundary',
        1,
        ['time/nlim=0'],
        expect_fail=True,
        expected_message='does not support mesh/ix1_bc=diode',
        input_deck=_NONPERIODIC_UNSUPPORTED_GUARD_DECK,
    )


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case_name, result in _POSITIVE_RESULTS.items():
        measured = result['measured']
        expected = result['expected']
        is_boris_case = case_name.startswith('serial_boris_')
        ok = _check_integrated_quantity(case_name, measured, expected, 'Q') and ok
        if not is_boris_case:
            ok = _check_integrated_quantity(case_name, measured, expected, 'Jx') and ok
            ok = _check_integrated_quantity(case_name, measured, expected, 'Jy') and ok
            ok = _check_integrated_quantity(case_name, measured, expected, 'Jz') and ok
            if result['edge_projection_exact']:
                edge_expectations = [('Jx_edge', 'Jx'), ('Jy_edge', 'Jy'),
                                     ('Jz_edge', 'Jz')]
                for edge_quantity, expected_quantity in edge_expectations:
                    if edge_quantity in measured:
                        expected_edge = expected[expected_quantity]
                        edge_abs_tol = (_EDGE_MOMENT_TOTAL_ABS_TOL
                                        if abs(expected_edge) > 1.0 else 1.0e-8)
                        ok = _check_with_tolerance(
                            case_name + ':' + edge_quantity,
                            measured[edge_quantity], expected_edge,
                            edge_abs_tol, _MOMENT_TOTAL_REL_TOL) and ok
        ok = _check_integrated_quantity(case_name, measured, expected, 'npart') and ok

    default_case = _POSITIVE_RESULTS['serial_default_optin_off']['measured']
    explicit_off = _POSITIVE_RESULTS['serial_explicit_optin_off']['measured']
    ok = _check_with_tolerance('default_vs_explicit_off:bcc1_l2',
                               explicit_off['bcc1_l2'], default_case['bcc1_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('default_vs_explicit_off:bcc2_l2',
                               explicit_off['bcc2_l2'], default_case['bcc2_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('default_vs_explicit_off:bcc3_l2',
                               explicit_off['bcc3_l2'], default_case['bcc3_l2'],
                               1.0e-12, 1.0e-12) and ok

    coupled = _POSITIVE_RESULTS['serial_coupled']['measured']
    delta_b1 = abs(coupled['bcc1_l2'] - default_case['bcc1_l2'])
    delta_b2 = abs(coupled['bcc2_l2'] - default_case['bcc2_l2'])
    delta_b3 = abs(coupled['bcc3_l2'] - default_case['bcc3_l2'])
    max_delta = max(delta_b1, delta_b2, delta_b3)
    logger.info('coupled_field_delta bcc1=% .8e bcc2=% .8e bcc3=% .8e max=% .8e',
                delta_b1, delta_b2, delta_b3, max_delta)
    if max_delta <= 1.0e-12:
        logger.warning('Coupled and uncoupled B-field norms are too similar')
        ok = False

    coupled_edge = _POSITIVE_RESULTS['serial_coupled_edge_staggered']['measured']
    coupled_edge_default = _POSITIVE_RESULTS[
        'serial_coupled_edge_default_mode']['measured']
    delta_edge_b1 = abs(coupled_edge['bcc1_l2'] - default_case['bcc1_l2'])
    delta_edge_b2 = abs(coupled_edge['bcc2_l2'] - default_case['bcc2_l2'])
    delta_edge_b3 = abs(coupled_edge['bcc3_l2'] - default_case['bcc3_l2'])
    max_delta_edge = max(delta_edge_b1, delta_edge_b2, delta_edge_b3)
    logger.info('edge_repr_coupled_field_delta bcc1=% .8e bcc2=% .8e '
                'bcc3=% .8e max=% .8e',
                delta_edge_b1, delta_edge_b2, delta_edge_b3, max_delta_edge)
    if max_delta_edge <= 1.0e-12:
        logger.warning('Edge-staggered coupled and uncoupled B-field norms are too '
                       'similar')
        ok = False

    coupled_edge_direct = _POSITIVE_RESULTS[
        'serial_coupled_edge_direct_staggered']['measured']
    coupled_edge_direct_o2 = _POSITIVE_RESULTS[
        'serial_coupled_edge_direct_staggered_o2']['measured']
    delta_direct_b1 = abs(coupled_edge_direct['bcc1_l2'] - default_case['bcc1_l2'])
    delta_direct_b2 = abs(coupled_edge_direct['bcc2_l2'] - default_case['bcc2_l2'])
    delta_direct_b3 = abs(coupled_edge_direct['bcc3_l2'] - default_case['bcc3_l2'])
    max_delta_direct = max(delta_direct_b1, delta_direct_b2, delta_direct_b3)
    logger.info('direct_edge_coupled_field_delta bcc1=% .8e bcc2=% .8e '
                'bcc3=% .8e max=% .8e',
                delta_direct_b1, delta_direct_b2, delta_direct_b3,
                max_delta_direct)
    if max_delta_direct <= 1.0e-12:
        logger.warning('Direct-edge coupled and uncoupled B-field norms are too '
                       'similar')
        ok = False

    for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
        ok = _check_integrated_quantity('default_edge_mode_vs_direct_staggered',
                                        coupled_edge_default, coupled_edge_direct,
                                        quantity) and ok
        ok = _check_integrated_quantity('cell_centered_vs_edge_staggered',
                                        coupled_edge, coupled, quantity) and ok
        ok = _check_integrated_quantity('edge_staggered_vs_direct_staggered',
                                        coupled_edge_direct, coupled_edge,
                                        quantity) and ok
        ok = _check_integrated_quantity('edge_staggered_vs_direct_staggered_o2',
                                        coupled_edge_direct_o2, coupled_edge,
                                        quantity) and ok
    for quantity in ['bcc1_l2', 'bcc2_l2', 'bcc3_l2']:
        ok = _check_with_tolerance('default_edge_mode_vs_direct_staggered:' +
                                   quantity,
                                   coupled_edge_default[quantity],
                                   coupled_edge_direct[quantity],
                                   1.0e-6, 1.0e-8) and ok
    default_vs_convert_delta = max(
        abs(coupled_edge_default['bcc1_l2'] - coupled_edge['bcc1_l2']),
        abs(coupled_edge_default['bcc2_l2'] - coupled_edge['bcc2_l2']),
        abs(coupled_edge_default['bcc3_l2'] - coupled_edge['bcc3_l2']))
    logger.info('default_edge_mode_delta_vs_cc_convert max=% .8e',
                default_vs_convert_delta)
    if default_vs_convert_delta <= 1.0e-12:
        logger.warning('Default edge deposition mode appears identical to '
                       'cc_convert baseline')
        ok = False

    required_projected_current_cases = [
        'serial_projected_current_cc',
        'serial_projected_current_edge',
        'serial_projected_current_edge_direct',
        'serial_projected_current_edge_direct_o2',
    ]
    for case_name in required_projected_current_cases:
        if case_name not in _PROJECTED_CURRENT_RESULTS:
            logger.warning('Missing projected-current diagnostic case %s', case_name)
            return False

    projected_cc = _PROJECTED_CURRENT_RESULTS['serial_projected_current_cc']
    projected_edge = _PROJECTED_CURRENT_RESULTS['serial_projected_current_edge']
    projected_direct = _PROJECTED_CURRENT_RESULTS[
        'serial_projected_current_edge_direct']
    projected_direct_o2 = _PROJECTED_CURRENT_RESULTS[
        'serial_projected_current_edge_direct_o2']
    projected_current_sets = [('cc', projected_cc), ('edge', projected_edge),
                              ('edge_direct', projected_direct),
                              ('edge_direct_o2', projected_direct_o2)]
    for label, projected in projected_current_sets:
        logger.info('projected_current_%s dt=% .8e cycle_delta=% .8e l2_drhodt=% .8e '
                    'l2_divj=% .8e l2_res=% .8e linf_res=% .8e rel=% .8e',
                    label, projected['dt'], projected['cycle_delta'],
                    projected['l2_drhodt'], projected['l2_divj'],
                    projected['l2_residual'], projected['linf_residual'],
                    projected['rel_residual'])
        if (not np.isfinite(projected['l2_residual']) or
                not np.isfinite(projected['linf_residual']) or
                not np.isfinite(projected['rel_residual'])):
            logger.warning('Non-finite projected-current metric for %s', label)
            ok = False
        if projected['dt'] <= 0.0 or projected['cycle_delta'] < 1.0:
            logger.warning('Invalid projected-current time metadata for %s', label)
            ok = False

    coupled_m1 = _measure_dataset('pic_mhd_coupled', 'mhd_u_m1')
    coupled_m2 = _measure_dataset('pic_mhd_coupled', 'mhd_u_m2')
    coupled_m3 = _measure_dataset('pic_mhd_coupled', 'mhd_u_m3')
    coupled_e = _measure_dataset('pic_mhd_coupled', 'mhd_u_e')
    coupled_edge_m1 = _measure_dataset('pic_mhd_coupled_edge', 'mhd_u_m1')
    coupled_edge_m2 = _measure_dataset('pic_mhd_coupled_edge', 'mhd_u_m2')
    coupled_edge_m3 = _measure_dataset('pic_mhd_coupled_edge', 'mhd_u_m3')
    coupled_edge_e = _measure_dataset('pic_mhd_coupled_edge', 'mhd_u_e')
    fb_mom_m1 = _measure_dataset('pic_mhd_feedback_momentum', 'mhd_u_m1')
    fb_mom_m2 = _measure_dataset('pic_mhd_feedback_momentum', 'mhd_u_m2')
    fb_mom_m3 = _measure_dataset('pic_mhd_feedback_momentum', 'mhd_u_m3')
    fb_eng_e = _measure_dataset('pic_mhd_feedback_energy', 'mhd_u_e')
    rk1_coupled_e = _measure_dataset('pic_mhd_rk1_coupled', 'mhd_u_e')
    rk3_coupled_e = _measure_dataset('pic_mhd_rk3_coupled', 'mhd_u_e')
    fb_eng_rk1_e = _measure_dataset('pic_mhd_feedback_energy_rk1', 'mhd_u_e')
    fb_eng_rk3_e = _measure_dataset('pic_mhd_feedback_energy_rk3', 'mhd_u_e')
    fb_both_m1 = _measure_dataset('pic_mhd_feedback_both', 'mhd_u_m1')
    fb_both_m2 = _measure_dataset('pic_mhd_feedback_both', 'mhd_u_m2')
    fb_both_m3 = _measure_dataset('pic_mhd_feedback_both', 'mhd_u_m3')
    fb_both_e = _measure_dataset('pic_mhd_feedback_both', 'mhd_u_e')
    fb_both_ef_m1 = _measure_dataset('pic_mhd_feedback_both_efield', 'mhd_u_m1')
    fb_both_ef_m2 = _measure_dataset('pic_mhd_feedback_both_efield', 'mhd_u_m2')
    fb_both_ef_m3 = _measure_dataset('pic_mhd_feedback_both_efield', 'mhd_u_m3')
    fb_both_ef_e = _measure_dataset('pic_mhd_feedback_both_efield', 'mhd_u_e')
    fb_both_ef_edge_m1 = _measure_dataset('pic_mhd_feedback_both_efield_edge',
                                          'mhd_u_m1')
    fb_both_ef_edge_m2 = _measure_dataset('pic_mhd_feedback_both_efield_edge',
                                          'mhd_u_m2')
    fb_both_ef_edge_m3 = _measure_dataset('pic_mhd_feedback_both_efield_edge',
                                          'mhd_u_m3')
    fb_both_ef_edge_e = _measure_dataset('pic_mhd_feedback_both_efield_edge',
                                         'mhd_u_e')
    boris_coupled_e = _measure_dataset('pic_mhd_boris_coupled', 'mhd_u_e')
    boris_fb_eng_e = _measure_dataset('pic_mhd_boris_feedback_energy', 'mhd_u_e')
    boris_fb_eng_x2_e = _measure_dataset('pic_mhd_boris_feedback_energy_x2',
                                         'mhd_u_e')

    delta_mom_m1 = _l2_difference_quantity(fb_mom_m1, coupled_m1, 'mom1')
    delta_mom_m2 = _l2_difference_quantity(fb_mom_m2, coupled_m2, 'mom2')
    delta_mom_m3 = _l2_difference_quantity(fb_mom_m3, coupled_m3, 'mom3')
    delta_mom = float(np.sqrt(delta_mom_m1 * delta_mom_m1 +
                              delta_mom_m2 * delta_mom_m2 +
                              delta_mom_m3 * delta_mom_m3))
    delta_eng = _l2_difference_quantity(fb_eng_e, coupled_e, 'ener')
    delta_eng_rk1 = _l2_difference_quantity(fb_eng_rk1_e, rk1_coupled_e, 'ener')
    delta_eng_rk3 = _l2_difference_quantity(fb_eng_rk3_e, rk3_coupled_e, 'ener')

    delta_both_m1 = _l2_difference_quantity(fb_both_m1, coupled_m1, 'mom1')
    delta_both_m2 = _l2_difference_quantity(fb_both_m2, coupled_m2, 'mom2')
    delta_both_m3 = _l2_difference_quantity(fb_both_m3, coupled_m3, 'mom3')
    delta_both_mom = float(np.sqrt(delta_both_m1 * delta_both_m1 +
                                   delta_both_m2 * delta_both_m2 +
                                   delta_both_m3 * delta_both_m3))
    delta_both_eng = _l2_difference_quantity(fb_both_e, coupled_e, 'ener')
    delta_both_ef_m1 = _l2_difference_quantity(fb_both_ef_m1, coupled_m1, 'mom1')
    delta_both_ef_m2 = _l2_difference_quantity(fb_both_ef_m2, coupled_m2, 'mom2')
    delta_both_ef_m3 = _l2_difference_quantity(fb_both_ef_m3, coupled_m3, 'mom3')
    delta_both_ef_mom = float(np.sqrt(delta_both_ef_m1 * delta_both_ef_m1 +
                                      delta_both_ef_m2 * delta_both_ef_m2 +
                                      delta_both_ef_m3 * delta_both_ef_m3))
    delta_both_ef_eng = _l2_difference_quantity(fb_both_ef_e, coupled_e, 'ener')
    delta_both_ef_edge_m1 = _l2_difference_quantity(fb_both_ef_edge_m1,
                                                    coupled_edge_m1, 'mom1')
    delta_both_ef_edge_m2 = _l2_difference_quantity(fb_both_ef_edge_m2,
                                                    coupled_edge_m2, 'mom2')
    delta_both_ef_edge_m3 = _l2_difference_quantity(fb_both_ef_edge_m3,
                                                    coupled_edge_m3, 'mom3')
    delta_both_ef_edge_mom = float(np.sqrt(
        delta_both_ef_edge_m1 * delta_both_ef_edge_m1 +
        delta_both_ef_edge_m2 * delta_both_ef_edge_m2 +
        delta_both_ef_edge_m3 * delta_both_ef_edge_m3))
    delta_both_ef_edge_eng = _l2_difference_quantity(fb_both_ef_edge_e,
                                                     coupled_edge_e, 'ener')
    logger.info('feedback_delta momentum_only=% .8e energy_only=% .8e '
                'energy_rk1=% .8e energy_rk3=% .8e '
                'both_mom=% .8e both_eng=% .8e '
                'both_ef_mom=% .8e both_ef_eng=% .8e '
                'both_ef_edge_mom=% .8e both_ef_edge_eng=% .8e',
                delta_mom, delta_eng, delta_eng_rk1, delta_eng_rk3,
                delta_both_mom, delta_both_eng,
                delta_both_ef_mom, delta_both_ef_eng,
                delta_both_ef_edge_mom, delta_both_ef_edge_eng)
    if delta_mom <= 1.0e-12:
        logger.warning('Momentum-feedback branch did not change fluid momentum')
        ok = False
    if delta_eng <= 1.0e-12:
        logger.warning('Energy-feedback branch did not change fluid energy')
        ok = False
    if delta_eng_rk1 <= 1.0e-12 or delta_eng_rk3 <= 1.0e-12:
        logger.warning('Integrator feedback scaling branches did not change fluid energy')
        ok = False
    if delta_eng > 1.0e-12 and delta_eng_rk1 > 1.0e-12 and delta_eng_rk3 > 1.0e-12:
        for label, delta in [('rk1', delta_eng_rk1), ('rk3', delta_eng_rk3)]:
            ratio = delta / delta_eng
            logger.info('feedback_integrator_ratio_%s measured=% .8e expected=% .8e',
                        label, ratio, 1.0)
            if ratio < 0.8 or ratio > 1.25:
                logger.warning('Feedback source strength changes with %s integrator',
                               label)
                ok = False
    if delta_both_mom <= 1.0e-12 or delta_both_eng <= 1.0e-12:
        logger.warning('Combined feedback branch did not affect both targets')
        ok = False
    if delta_both_ef_mom <= 1.0e-12 or delta_both_ef_eng <= 1.0e-12:
        logger.warning('EField-ordered combined feedback branch did not affect both '
                       'targets')
        ok = False
    if delta_both_ef_edge_mom <= 1.0e-12 or delta_both_ef_edge_eng <= 1.0e-12:
        logger.warning('EField-ordered edge-staggered feedback branch did not affect '
                       'both targets')
        ok = False

    boris_delta_eng = _l2_difference_quantity(boris_fb_eng_e, boris_coupled_e,
                                              'ener')
    boris_delta_eng_x2 = _l2_difference_quantity(boris_fb_eng_x2_e,
                                                 boris_coupled_e, 'ener')
    logger.info('boris_feedback_delta energy_coeff1=% .8e energy_coeff2=% .8e',
                boris_delta_eng, boris_delta_eng_x2)
    if boris_delta_eng <= 1.0e-12:
        logger.warning('Boris energy-feedback branch did not change fluid energy')
        ok = False
    if boris_delta_eng_x2 <= 1.0e-12:
        logger.warning('Boris energy-feedback coeff=2 branch did not change fluid '
                       'energy')
        ok = False
    if boris_delta_eng > 1.0e-12 and boris_delta_eng_x2 > 1.0e-12:
        boris_coeff_ratio = boris_delta_eng_x2 / boris_delta_eng
        logger.info('boris_feedback_coeff_ratio measured=% .8e expected=% .8e',
                    boris_coeff_ratio, 2.0)
        if boris_coeff_ratio < 1.5 or boris_coeff_ratio > 2.5:
            logger.warning('Boris energy-feedback coefficient sensitivity is outside '
                           'expected range')
            ok = False

    zero_uncoupled = _POSITIVE_RESULTS['serial_zero_uncoupled']['measured']
    zero_coupled = _POSITIVE_RESULTS['serial_zero_coupled']['measured']
    ok = _check_with_tolerance('zero_current:bcc1_l2',
                               zero_coupled['bcc1_l2'], zero_uncoupled['bcc1_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('zero_current:bcc2_l2',
                               zero_coupled['bcc2_l2'], zero_uncoupled['bcc2_l2'],
                               1.0e-12, 1.0e-12) and ok
    ok = _check_with_tolerance('zero_current:bcc3_l2',
                               zero_coupled['bcc3_l2'], zero_uncoupled['bcc3_l2'],
                               1.0e-12, 1.0e-12) and ok

    lin_lo_unc_bcc = _measure_bcc_dataset('pic_mhd_lin_lo_uncoupled')
    lin_lo_cpl_bcc = _measure_bcc_dataset('pic_mhd_lin_lo_coupled')
    lin_lo_neg_cpl_bcc = _measure_bcc_dataset('pic_mhd_lin_lo_neg_coupled')
    lin_hi_unc_bcc = _measure_bcc_dataset('pic_mhd_lin_hi_uncoupled')
    lin_hi_cpl_bcc = _measure_bcc_dataset('pic_mhd_lin_hi_coupled')

    delta_lo_b1 = _l2_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc, 'bcc1')
    delta_lo_b2 = _l2_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc, 'bcc2')
    delta_lo_b3 = _l2_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc, 'bcc3')
    delta_lo = float(np.sqrt(delta_lo_b1 * delta_lo_b1 +
                             delta_lo_b2 * delta_lo_b2 +
                             delta_lo_b3 * delta_lo_b3))

    delta_hi_b1 = _l2_difference_quantity(lin_hi_cpl_bcc, lin_hi_unc_bcc, 'bcc1')
    delta_hi_b2 = _l2_difference_quantity(lin_hi_cpl_bcc, lin_hi_unc_bcc, 'bcc2')
    delta_hi_b3 = _l2_difference_quantity(lin_hi_cpl_bcc, lin_hi_unc_bcc, 'bcc3')
    delta_hi = float(np.sqrt(delta_hi_b1 * delta_hi_b1 +
                             delta_hi_b2 * delta_hi_b2 +
                             delta_hi_b3 * delta_hi_b3))

    logger.info('linearity_delta low=% .8e high=% .8e', delta_lo, delta_hi)
    if delta_lo <= 1.0e-12 or delta_hi <= 1.0e-12:
        logger.warning('Linearity check deltas are too small')
        ok = False
    else:
        ratio = delta_hi / delta_lo
        expected_ratio = 2.0
        ratio_abs_err = abs(ratio - expected_ratio)
        logger.info('linearity_ratio measured=% .8e expected=% .8e abs_err=% .8e',
                    ratio, expected_ratio, ratio_abs_err)
        if ratio_abs_err > 2.0e-1:
            logger.warning('Linearity ratio is outside tolerance')
            ok = False

    for quantity in ['bcc1', 'bcc2', 'bcc3']:
        delta_pos = _l2_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc, quantity)
        delta_neg = _l2_difference_quantity(lin_lo_neg_cpl_bcc,
                                            lin_lo_unc_bcc, quantity)
        anti_delta = _l2_sum_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc,
                                                 lin_lo_neg_cpl_bcc, quantity)
        dot = _dot_difference_quantity(lin_lo_cpl_bcc, lin_lo_unc_bcc,
                                       lin_lo_neg_cpl_bcc, quantity)
        if delta_pos <= 1.0e-12 or delta_neg <= 1.0e-12:
            logger.warning('Sign oracle deltas are too small for %s', quantity)
            ok = False
            continue
        cosine = dot / (delta_pos * delta_neg)
        anti_rel = anti_delta / max(delta_pos, delta_neg)
        logger.info('sign_oracle_%s cos=% .8e anti_rel=% .8e',
                    quantity, cosine, anti_rel)
        if cosine >= -9.5e-1:
            logger.warning('Sign oracle cosine check failed for %s', quantity)
            ok = False
        if anti_rel >= 2.5e-1:
            logger.warning('Sign oracle anti-symmetry check failed for %s', quantity)
            ok = False

    rk1_unc = _POSITIVE_RESULTS['serial_rk1_uncoupled']['measured']
    rk1_cpl = _POSITIVE_RESULTS['serial_rk1_coupled']['measured']
    rk3_unc = _POSITIVE_RESULTS['serial_rk3_uncoupled']['measured']
    rk3_cpl = _POSITIVE_RESULTS['serial_rk3_coupled']['measured']

    rk1_delta = max(abs(rk1_cpl['bcc1_l2'] - rk1_unc['bcc1_l2']),
                    abs(rk1_cpl['bcc2_l2'] - rk1_unc['bcc2_l2']),
                    abs(rk1_cpl['bcc3_l2'] - rk1_unc['bcc3_l2']))
    rk3_delta = max(abs(rk3_cpl['bcc1_l2'] - rk3_unc['bcc1_l2']),
                    abs(rk3_cpl['bcc2_l2'] - rk3_unc['bcc2_l2']),
                    abs(rk3_cpl['bcc3_l2'] - rk3_unc['bcc3_l2']))
    logger.info('integrator_delta rk1=% .8e rk3=% .8e', rk1_delta, rk3_delta)
    if rk1_delta <= 1.0e-12 or rk3_delta <= 1.0e-12:
        logger.warning('Coupled field delta missing in rk1/rk3 cases')
        ok = False

    for quantity in ['Q', 'Jx', 'Jy', 'Jz', 'npart']:
        ok = _check_integrated_quantity('rk1_vs_rk2_coupled', rk1_cpl,
                                        coupled, quantity) and ok
        ok = _check_integrated_quantity('rk3_vs_rk2_coupled', rk3_cpl,
                                        coupled, quantity) and ok

    if 'mpi2_coupled' in _POSITIVE_RESULTS:
        mpi2 = _POSITIVE_RESULTS['mpi2_coupled']['measured']
        ok = _check_integrated_quantity('serial_vs_mpi2_coupled', mpi2,
                                        coupled, 'Q') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_coupled', mpi2,
                                        coupled, 'Jx') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_coupled', mpi2,
                                        coupled, 'Jy') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_coupled', mpi2,
                                        coupled, 'Jz') and ok
        ok = _check_with_tolerance('serial_vs_mpi2_coupled:bcc1_l2',
                                   mpi2['bcc1_l2'], coupled['bcc1_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_coupled:bcc2_l2',
                                   mpi2['bcc2_l2'], coupled['bcc2_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_coupled:bcc3_l2',
                                   mpi2['bcc3_l2'], coupled['bcc3_l2'],
                                   1.0e-6, 1.0e-8) and ok

    if 'mpi2_coupled_edge_staggered' in _POSITIVE_RESULTS:
        mpi2_edge = _POSITIVE_RESULTS['mpi2_coupled_edge_staggered']['measured']
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_staggered', mpi2_edge,
                                        coupled_edge, 'Q') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_staggered', mpi2_edge,
                                        coupled_edge, 'Jx') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_staggered', mpi2_edge,
                                        coupled_edge, 'Jy') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_staggered', mpi2_edge,
                                        coupled_edge, 'Jz') and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_staggered:bcc1_l2',
                                   mpi2_edge['bcc1_l2'], coupled_edge['bcc1_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_staggered:bcc2_l2',
                                   mpi2_edge['bcc2_l2'], coupled_edge['bcc2_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_staggered:bcc3_l2',
                                   mpi2_edge['bcc3_l2'], coupled_edge['bcc3_l2'],
                                   1.0e-6, 1.0e-8) and ok

    if 'mpi2_coupled_edge_default_mode' in _POSITIVE_RESULTS:
        mpi2_edge_default = _POSITIVE_RESULTS[
            'mpi2_coupled_edge_default_mode']['measured']
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_default_mode',
                                        mpi2_edge_default, coupled_edge_default,
                                        'Q') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_default_mode',
                                        mpi2_edge_default, coupled_edge_default,
                                        'Jx') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_default_mode',
                                        mpi2_edge_default, coupled_edge_default,
                                        'Jy') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_edge_default_mode',
                                        mpi2_edge_default, coupled_edge_default,
                                        'Jz') and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_default_mode:bcc1_l2',
                                   mpi2_edge_default['bcc1_l2'],
                                   coupled_edge_default['bcc1_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_default_mode:bcc2_l2',
                                   mpi2_edge_default['bcc2_l2'],
                                   coupled_edge_default['bcc2_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_edge_default_mode:bcc3_l2',
                                   mpi2_edge_default['bcc3_l2'],
                                   coupled_edge_default['bcc3_l2'],
                                   1.0e-6, 1.0e-8) and ok

    if 'mpi2_coupled_edge_direct_staggered' in _POSITIVE_RESULTS:
        mpi2_edge_direct = _POSITIVE_RESULTS[
            'mpi2_coupled_edge_direct_staggered']['measured']
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered',
                                        mpi2_edge_direct, coupled_edge_direct,
                                        'Q') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered',
                                        mpi2_edge_direct, coupled_edge_direct,
                                        'Jx') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered',
                                        mpi2_edge_direct, coupled_edge_direct,
                                        'Jy') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered',
                                        mpi2_edge_direct, coupled_edge_direct,
                                        'Jz') and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered:bcc1_l2',
                                   mpi2_edge_direct['bcc1_l2'],
                                   coupled_edge_direct['bcc1_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered:bcc2_l2',
                                   mpi2_edge_direct['bcc2_l2'],
                                   coupled_edge_direct['bcc2_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered:bcc3_l2',
                                   mpi2_edge_direct['bcc3_l2'],
                                   coupled_edge_direct['bcc3_l2'],
                                   1.0e-6, 1.0e-8) and ok
    if 'mpi2_coupled_edge_direct_staggered_o2' in _POSITIVE_RESULTS:
        mpi2_edge_direct_o2 = _POSITIVE_RESULTS[
            'mpi2_coupled_edge_direct_staggered_o2']['measured']
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered_o2',
                                        mpi2_edge_direct_o2,
                                        coupled_edge_direct_o2,
                                        'Q') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered_o2',
                                        mpi2_edge_direct_o2,
                                        coupled_edge_direct_o2,
                                        'Jx') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered_o2',
                                        mpi2_edge_direct_o2,
                                        coupled_edge_direct_o2,
                                        'Jy') and ok
        ok = _check_integrated_quantity('serial_vs_mpi2_direct_staggered_o2',
                                        mpi2_edge_direct_o2,
                                        coupled_edge_direct_o2,
                                        'Jz') and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered_o2:bcc1_l2',
                                   mpi2_edge_direct_o2['bcc1_l2'],
                                   coupled_edge_direct_o2['bcc1_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered_o2:bcc2_l2',
                                   mpi2_edge_direct_o2['bcc2_l2'],
                                   coupled_edge_direct_o2['bcc2_l2'],
                                   1.0e-6, 1.0e-8) and ok
        ok = _check_with_tolerance('serial_vs_mpi2_direct_staggered_o2:bcc3_l2',
                                   mpi2_edge_direct_o2['bcc3_l2'],
                                   coupled_edge_direct_o2['bcc3_l2'],
                                   1.0e-6, 1.0e-8) and ok

    if ('mpi2_projected_current_cc' in _PROJECTED_CURRENT_RESULTS and
            'mpi2_projected_current_edge' in _PROJECTED_CURRENT_RESULTS and
            'mpi2_projected_current_edge_direct' in _PROJECTED_CURRENT_RESULTS and
            'mpi2_projected_current_edge_direct_o2' in _PROJECTED_CURRENT_RESULTS):
        for rep_tag in ['cc', 'edge', 'edge_direct', 'edge_direct_o2']:
            serial = _PROJECTED_CURRENT_RESULTS['serial_projected_current_' + rep_tag]
            mpi2 = _PROJECTED_CURRENT_RESULTS['mpi2_projected_current_' + rep_tag]
            ok = _check_with_tolerance('projected_current_serial_vs_mpi2:' + rep_tag +
                                       ':l2_residual',
                                       mpi2['l2_residual'], serial['l2_residual'],
                                       1.0e-6, 1.0e-8) and ok
            ok = _check_with_tolerance('projected_current_serial_vs_mpi2:' + rep_tag +
                                       ':linf_residual',
                                       mpi2['linf_residual'], serial['linf_residual'],
                                       1.0e-6, 1.0e-8) and ok
            ok = _check_with_tolerance('projected_current_serial_vs_mpi2:' + rep_tag +
                                       ':rel_residual',
                                       mpi2['rel_residual'], serial['rel_residual'],
                                       1.0e-6, 1.0e-8) and ok

    ok = len(_NEGATIVE_RESULTS) == 12 and ok
    if len(_NEGATIVE_RESULTS) != 12:
        logger.warning('Expected 12 negative checks, got %d', len(_NEGATIVE_RESULTS))

    return ok
