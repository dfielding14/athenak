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

_INPUT_DECK = 'tests/pic_boris_midpoint_eb.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_B0 = np.array([1.0, np.sqrt(2.0), 0.5])
_V0 = np.array([0.2, 0.1, 0.0])
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


def _boris_midpoint_step(v0, dt, vflow):
    u = np.array([vflow, 0.0, 0.0])
    c_e = -np.cross(u, _B0)
    vminus = v0 + 0.5 * dt * c_e
    t = 0.5 * dt * _B0
    s = 2.0 * t / (1.0 + np.dot(t, t))
    vprime = vminus + np.cross(vminus, t)
    vplus = vminus + np.cross(vprime, s)
    vnew = vplus + 0.5 * dt * c_e
    return vnew


def _run_command(label, nproc, arguments):
    command = ['./athena', '-i', _athena_input_path()] + list(arguments)
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + ': ' + ' '.join(command))


def _measure_case(basename):
    rho_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_rho'))
    jx_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jx'))
    jy_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jy'))
    jz_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_jz'))
    d_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_d'))
    dpxdt_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_dpxdt'))
    dpydt_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_dpydt'))
    dpzdt_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_dpzdt'))
    dedt_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_dedt'))
    ebdot_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'prtcl_ebdot'))
    m1_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m1'))
    m2_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m2'))
    m3_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_m3'))
    e_data = bin_convert.read_binary_as_athdf(
        _latest_output_file(basename, 'mhd_u_e'))

    qtot = _integrate_quantity(rho_data, 'prtcl_rho')
    jtot = np.array([
        _integrate_quantity(jx_data, 'prtcl_jx'),
        _integrate_quantity(jy_data, 'prtcl_jy'),
        _integrate_quantity(jz_data, 'prtcl_jz'),
    ])
    npart = _integrate_quantity(d_data, 'pdens')
    if abs(qtot) < 1.0e-14:
        raise RuntimeError(
            'Integrated particle charge is too small for velocity recovery'
        )

    vavg = jtot / qtot
    pke = 0.5 * npart * np.dot(vavg, vavg)

    return {
        'time': float(rho_data['Time']),
        'qtot': qtot,
        'jtot': jtot,
        'dpdt_tot': np.array([
            _integrate_quantity(dpxdt_data, 'prtcl_dpxdt'),
            _integrate_quantity(dpydt_data, 'prtcl_dpydt'),
            _integrate_quantity(dpzdt_data, 'prtcl_dpzdt'),
        ]),
        'dedt_tot': _integrate_quantity(dedt_data, 'prtcl_dedt'),
        'ebdot_tot': _integrate_quantity(ebdot_data, 'prtcl_ebdot'),
        'vavg': vavg,
        'npart': npart,
        'pke': pke,
        'gas_mom': np.array([
            _integrate_quantity(m1_data, 'mom1'),
            _integrate_quantity(m2_data, 'mom2'),
            _integrate_quantity(m3_data, 'mom3'),
        ]),
        'gas_e': _integrate_quantity(e_data, 'ener'),
    }


def _check_close_vec(label, measured, expected, abs_tol):
    err = np.linalg.norm(measured - expected)
    logger.info('%s |meas-exp|=% .8e', label, err)
    return err <= abs_tol


def _check_close_scalar(label, measured, expected, abs_tol):
    err = abs(measured - expected)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e',
                label, measured, expected, err)
    return err <= abs_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    serial_cases = [
        {
            'name': 'passive_zero',
            'basename': 'pic_boris_mid_passive_zero',
            'nproc': 1,
            'args': [
                'job/basename=pic_boris_mid_passive_zero',
                'problem/vflow=0.0',
                'particles/pic_background_mode=passive_mhd',
                'particles/pic_feedback_mode=test_particle',
                'particles/couple_moments_to_mhd=false',
                'particles/couple_moments_momentum_to_mhd=false',
                'particles/couple_moments_energy_to_mhd=false',
            ],
        },
        {
            'name': 'passive_flow',
            'basename': 'pic_boris_mid_passive_flow',
            'nproc': 1,
            'args': [
                'job/basename=pic_boris_mid_passive_flow',
                'problem/vflow=0.3',
                'particles/pic_background_mode=passive_mhd',
                'particles/pic_feedback_mode=test_particle',
                'particles/couple_moments_to_mhd=false',
                'particles/couple_moments_momentum_to_mhd=false',
                'particles/couple_moments_energy_to_mhd=false',
            ],
        },
        {
            'name': 'coupled_uncoupled',
            'basename': 'pic_boris_mid_coupled_uncoupled',
            'nproc': 1,
            'args': [
                'job/basename=pic_boris_mid_coupled_uncoupled',
                'problem/vflow=0.3',
                'particles/pic_background_mode=coupled',
                'particles/pic_feedback_mode=test_particle',
                'particles/couple_moments_to_mhd=false',
                'particles/couple_moments_momentum_to_mhd=false',
                'particles/couple_moments_energy_to_mhd=false',
            ],
        },
        {
            'name': 'coupled_feedback',
            'basename': 'pic_boris_mid_coupled_feedback',
            'nproc': 1,
            'args': [
                'job/basename=pic_boris_mid_coupled_feedback',
                'problem/vflow=0.3',
                'particles/pic_background_mode=coupled',
                'particles/pic_feedback_mode=coupled',
                'particles/couple_moments_to_mhd=true',
                'particles/couple_moments_momentum_to_mhd=true',
                'particles/couple_moments_energy_to_mhd=true',
                'particles/couple_moments_momentum_coeff=1.0',
                'particles/couple_moments_energy_coeff=1.0',
            ],
        },
    ]

    for case in serial_cases:
        _remove_outputs(case['basename'])
        _run_command(case['name'], case['nproc'], case['args'])
        _RESULTS[case['name']] = _measure_case(case['basename'])

    if _athena_mpi_enabled():
        mpi_cases = [
            {
                'name': 'passive_flow_mpi2',
                'basename': 'pic_boris_mid_passive_flow_mpi2',
                'nproc': 2,
                'args': [
                    'job/basename=pic_boris_mid_passive_flow_mpi2',
                    'problem/vflow=0.3',
                    'particles/pic_background_mode=passive_mhd',
                    'particles/pic_feedback_mode=test_particle',
                    'particles/couple_moments_to_mhd=false',
                    'particles/couple_moments_momentum_to_mhd=false',
                    'particles/couple_moments_energy_to_mhd=false',
                ],
            },
            {
                'name': 'coupled_feedback_mpi2',
                'basename': 'pic_boris_mid_coupled_feedback_mpi2',
                'nproc': 2,
                'args': [
                    'job/basename=pic_boris_mid_coupled_feedback_mpi2',
                    'problem/vflow=0.3',
                    'particles/pic_background_mode=coupled',
                    'particles/pic_feedback_mode=coupled',
                    'particles/couple_moments_to_mhd=true',
                    'particles/couple_moments_momentum_to_mhd=true',
                    'particles/couple_moments_energy_to_mhd=true',
                    'particles/couple_moments_momentum_coeff=1.0',
                    'particles/couple_moments_energy_coeff=1.0',
                ],
            },
        ]
        for case in mpi_cases:
            _remove_outputs(case['basename'])
            _run_command(case['name'], case['nproc'], case['args'])
            _RESULTS[case['name']] = _measure_case(case['basename'])


def analyze():
    logger.debug('Analyzing test ' + __name__)

    ok = True

    dt0 = _RESULTS['passive_zero']['time']
    dt1 = _RESULTS['passive_flow']['time']
    vexp_zero = _boris_midpoint_step(_V0, dt0, 0.0)
    vexp_flow = _boris_midpoint_step(_V0, dt1, 0.3)

    ok = _check_close_vec('passive_zero:vavg',
                          _RESULTS['passive_zero']['vavg'], vexp_zero,
                          1.0e-6) and ok
    ok = _check_close_vec('passive_flow:vavg',
                          _RESULTS['passive_flow']['vavg'], vexp_flow,
                          1.0e-6) and ok

    speed_delta = abs(np.linalg.norm(_RESULTS['passive_flow']['vavg']) -
                      np.linalg.norm(_RESULTS['passive_zero']['vavg']))
    logger.info('passive_speed_delta=% .8e', speed_delta)
    ok = (speed_delta > 1.0e-3) and ok

    c_e = -np.cross(np.array([0.3, 0.0, 0.0]), _B0)
    ok = _check_close_scalar('frozen_in_cE_dot_B', float(np.dot(c_e, _B0)),
                             0.0, 1.0e-14) and ok
    ok = _check_close_scalar('frozen_in_ebdot_output',
                             _RESULTS['passive_flow']['ebdot_tot'], 0.0,
                             1.0e-12) and ok

    unc = _RESULTS['coupled_uncoupled']
    cpl = _RESULTS['coupled_feedback']
    # A cycle-fixed source applied through the RK stages integrates to one full
    # dt over the step; stage weights are an implementation detail of the RK
    # update, not a half-strength physical feedback split.
    feedback_integral = 1.0

    d_particle_mom = unc['jtot'] - unc['qtot'] * _V0
    d_gas_mom = cpl['gas_mom'] - unc['gas_mom']
    ok = _check_close_vec('coupled_momentum_exchange',
                          d_gas_mom + feedback_integral * d_particle_mom,
                          np.zeros(3), 1.0e-5) and ok

    dt_cpl = cpl['time']
    ok = _check_close_vec('coupled_momentum_exchange_diag',
                          d_gas_mom + feedback_integral * dt_cpl * cpl['dpdt_tot'],
                          np.zeros(3), 1.0e-5) and ok

    ke0 = 0.5 * unc['npart'] * np.dot(_V0, _V0)
    d_particle_ke = unc['pke'] - ke0
    d_gas_e = cpl['gas_e'] - unc['gas_e']
    energy_tol = 5.0e-5
    ok = _check_close_scalar('coupled_energy_exchange',
                             d_gas_e + feedback_integral * d_particle_ke, 0.0,
                             energy_tol) and ok
    ok = _check_close_scalar('coupled_energy_exchange_diag',
                             d_gas_e + feedback_integral * dt_cpl * cpl['dedt_tot'],
                             0.0, energy_tol) and ok

    if 'passive_flow_mpi2' in _RESULTS:
        ok = _check_close_vec('serial_vs_mpi2:passive_flow_vavg',
                              _RESULTS['passive_flow_mpi2']['vavg'],
                              _RESULTS['passive_flow']['vavg'],
                              5.0e-12) and ok
    if 'coupled_feedback_mpi2' in _RESULTS:
        ok = _check_close_vec('serial_vs_mpi2:coupled_feedback_jtot',
                              _RESULTS['coupled_feedback_mpi2']['jtot'],
                              _RESULTS['coupled_feedback']['jtot'],
                              5.0e-11) and ok
        ok = _check_close_vec('serial_vs_mpi2:coupled_feedback_gas_mom',
                              _RESULTS['coupled_feedback_mpi2']['gas_mom'],
                              _RESULTS['coupled_feedback']['gas_mom'],
                              5.0e-11) and ok

    return ok
