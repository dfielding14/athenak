import glob
import logging
import os
import re
import subprocess
import sys

import numpy as np

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_GYRO_INPUT = 'tests/pic_relativistic_gyro_paper.athinput'
_CPAW_INPUT = 'tests/pic_mhd_expanding_box_cpaw_local.athinput'
_GYRO_RATE = 0.01
_GYRO_TEND = 6.0
_GYRO_C = 3.0
_GYRO_P0 = 1.0
_GYRO_THETA_VALUES = (0.2, 0.1, 0.05)
_CPAW_RATES = {
    'static': 0.0,
    'expanding': 0.01,
    'compressing': -0.01,
}
_CPAW_GAMMA = 1.66666666667
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get('ATHENA_Q008_EXE_DIR',
                          os.path.join(os.getcwd(), 'build', 'src'))


def _athena_input_path(input_deck):
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..',
                                        '..', 'inputs', input_deck))


def _run_athena(label, input_deck, arguments):
    command = ['./athena', '-i', _athena_input_path(input_deck)] + list(arguments)
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
            os.path.join(exe_dir, 'pvtk', basename + '.*.part.vtk'),
            os.path.join(exe_dir, 'bin', basename + '.*.bin'),
            os.path.join(exe_dir, basename + '-errs.dat')]:
        for path in glob.glob(pattern):
            os.remove(path)


def _read_big_endian_floats(contents, offset, count, label):
    payload = contents[offset:offset + 4*count]
    if len(payload) != 4*count:
        raise RuntimeError('Truncated ' + label + ' payload in particle VTK output')
    return np.frombuffer(payload, dtype='>f4').astype(np.float64)


def _find_marker(contents, pattern, offset, label):
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise RuntimeError('Could not find ' + label + ' marker in particle VTK output')
    return offset + match.start(), offset + match.end(), match


def _read_pvtk(path):
    with open(path, 'rb') as fp:
        contents = fp.read()

    header = re.search(
        rb'# AthenaK particle data at time=\s*([^ ]+)\s+nranks=.*cycle=([0-9]+)',
        contents)
    if header is None:
        raise RuntimeError('Could not read particle VTK cycle header in ' + path)

    _, offset, match = _find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS')
    npoint = int(match.group(1))
    offset += 4*3*npoint
    _, offset, _ = _find_marker(contents, rb'\nVECTORS vel float\n',
                                offset, 'VECTORS vel')
    velocity = _read_big_endian_floats(contents, offset, 3*npoint, 'vel')
    return {
        'time': float(header.group(1)),
        'cycle': int(header.group(2)),
        'npoint': npoint,
        'velocity': velocity.reshape(npoint, 3),
    }


def _particle_momentum(velocity):
    velocity2 = np.sum(velocity*velocity, axis=1)
    gamma = 1.0/np.sqrt(1.0 - velocity2/(_GYRO_C*_GYRO_C))
    return velocity*gamma[:, None], gamma


def _gyro_phase(time):
    scale = 1.0 + _GYRO_RATE*time
    p0_over_c = _GYRO_P0/_GYRO_C
    return (np.arcsinh(p0_over_c)
            - np.arcsinh(p0_over_c/scale)) / (_GYRO_RATE*p0_over_c)


def _measure_gyro_case(basename):
    pattern = os.path.join(_athena_exe_dir(), 'pvtk',
                           basename + '.prtcl_all.*.part.vtk')
    snapshots = [_read_pvtk(path) for path in sorted(glob.glob(pattern))]
    if len(snapshots) < 3:
        raise RuntimeError('Need at least three gyro snapshots for ' + basename)

    measured_phase = []
    expected_phase = []
    pnorm_errors = []
    gamma_errors = []
    momentum_spreads = []
    for snapshot in snapshots:
        momentum, gamma = _particle_momentum(snapshot['velocity'])
        mean_momentum = np.mean(momentum, axis=0)
        scale = 1.0 + _GYRO_RATE*snapshot['time']
        expected_pnorm = _GYRO_P0/scale
        expected_gamma = np.sqrt(1.0 + expected_pnorm**2/(_GYRO_C*_GYRO_C))
        measured_phase.append(np.arctan2(-mean_momentum[1], mean_momentum[0]))
        expected_phase.append(_gyro_phase(snapshot['time']))
        pnorm_errors.append(
            np.max(np.abs(np.linalg.norm(momentum, axis=1) - expected_pnorm))
            / expected_pnorm)
        gamma_errors.append(np.max(np.abs(gamma - expected_gamma)))
        momentum_spreads.append(np.max(np.abs(momentum - mean_momentum[None, :])))

    measured_phase = np.unwrap(np.asarray(measured_phase))
    expected_phase = np.asarray(expected_phase)
    phase_errors = np.abs(measured_phase - expected_phase)
    final = snapshots[-1]
    return {
        'time': final['time'],
        'cycle': final['cycle'],
        'dt': final['time']/final['cycle'],
        'snapshot_count': len(snapshots),
        'npoint': final['npoint'],
        'history_pnorm_rel_error': float(np.max(pnorm_errors)),
        'history_gamma_abs_error': float(np.max(gamma_errors)),
        'history_phase_abs_error': float(np.max(phase_errors)),
        'endpoint_phase_abs_error': float(phase_errors[-1]),
        'momentum_spread': float(np.max(momentum_spreads)),
    }


def _output_files(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Need initial and final CPAW snapshots for ' + pattern)
    return paths


def _read_endpoints(basename, file_id):
    paths = _output_files(basename, file_id)
    return (bin_convert.read_binary_as_athdf(paths[0]),
            bin_convert.read_binary_as_athdf(paths[-1]))


def _relative_error(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-14)
    return float(np.max(np.abs(actual - expected))) / scale


def _wrapped_abs(angle):
    return float(abs(np.arctan2(np.sin(angle), np.cos(angle))))


def _mode_coefficient(dataset, y_name, z_name):
    x1 = dataset['x1v']
    wavenumber = 2.0*np.pi/(float(dataset['x1f'][-1]) - float(dataset['x1f'][0]))
    phase = np.exp(1j*wavenumber*x1)[None, None, :]
    return complex(np.mean((dataset[y_name] + 1j*dataset[z_name])*phase))


def _measure_cpaw_case(basename):
    wfirst, wlast = _read_endpoints(basename, 'mhd_w')
    bfirst, blast = _read_endpoints(basename, 'mhd_bcc')
    return {
        'time0': float(wfirst['Time']),
        'time1': float(wlast['Time']),
        'wfirst': wfirst,
        'wlast': wlast,
        'bfirst': bfirst,
        'blast': blast,
        'v_mode_initial': _mode_coefficient(wfirst, 'vely', 'velz'),
        'v_mode_final': _mode_coefficient(wlast, 'vely', 'velz'),
        'b_mode_initial': _mode_coefficient(bfirst, 'bcc2', 'bcc3'),
        'b_mode_final': _mode_coefficient(blast, 'bcc2', 'bcc3'),
    }


def _summarize_cpaw(reference, measured, rate):
    time_error = abs(measured['time1'] - reference['time1'])
    scale = 1.0 + rate*measured['time1']
    ratio = 1.0/scale
    expected = {
        'dens': reference['wlast']['dens']*ratio,
        'velx': reference['wlast']['velx']*ratio,
        'vely': reference['wlast']['vely'],
        'velz': reference['wlast']['velz'],
        'eint': reference['wlast']['eint']*ratio**_CPAW_GAMMA,
        'bcc1': reference['blast']['bcc1'],
        'bcc2': reference['blast']['bcc2']*ratio,
        'bcc3': reference['blast']['bcc3']*ratio,
    }
    actual = {
        'dens': measured['wlast']['dens'],
        'velx': measured['wlast']['velx'],
        'vely': measured['wlast']['vely'],
        'velz': measured['wlast']['velz'],
        'eint': measured['wlast']['eint'],
        'bcc1': measured['blast']['bcc1'],
        'bcc2': measured['blast']['bcc2'],
        'bcc3': measured['blast']['bcc3'],
    }
    return {
        'time': measured['time1'],
        'time_abs_error_vs_static': time_error,
        'source_map_errors': {
            name: _relative_error(actual[name], expected[name]) for name in expected
        },
        'v_mode_amplitude_rel_error': abs(
            abs(measured['v_mode_final']/reference['v_mode_final']) - 1.0),
        'b_mode_amplitude_rel_error': abs(
            abs(measured['b_mode_final']/reference['b_mode_final']) - ratio),
        'v_mode_phase_abs_error': _wrapped_abs(
            np.angle(measured['v_mode_final']/reference['v_mode_final'])),
        'b_mode_phase_abs_error': _wrapped_abs(
            np.angle(measured['b_mode_final']/reference['b_mode_final'])),
    }


def _summarize_static_cpaw(measured):
    dt = measured['time1'] - measured['time0']
    wavenumber = 2.0*np.pi
    phase_delta = np.angle(measured['v_mode_final']/measured['v_mode_initial'])
    translation_phase_error = min(
        _wrapped_abs(phase_delta - wavenumber*dt),
        _wrapped_abs(phase_delta + wavenumber*dt),
    )
    initial_lock = measured['v_mode_initial']/measured['b_mode_initial']
    final_lock = measured['v_mode_final']/measured['b_mode_final']
    return {
        'time': measured['time1'],
        'v_mode_amplitude_rel_drift': abs(
            abs(measured['v_mode_final']/measured['v_mode_initial']) - 1.0),
        'b_mode_amplitude_rel_drift': abs(
            abs(measured['b_mode_final']/measured['b_mode_initial']) - 1.0),
        'translation_phase_abs_error': translation_phase_error,
        'vb_phase_lock_abs_error': _wrapped_abs(np.angle(final_lock/initial_lock)),
    }


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    gyro_cases = []
    for theta_max in _GYRO_THETA_VALUES:
        tag = str(theta_max).replace('.', 'p')
        basename = 'pic_box_appendix_local_gyro_theta' + tag
        _remove_outputs(basename)
        _run_athena(
            'gyro_theta_' + tag, _GYRO_INPUT,
            ['job/basename=' + basename,
             'particles/pic_expanding_box_mode=on',
             'particles/pic_expansion_law=linear',
             'particles/pic_expansion_rate_x1=' + str(_GYRO_RATE),
             'particles/pic_expansion_rate_x2=' + str(_GYRO_RATE),
             'particles/pic_expansion_rate_x3=0.0',
             'particles/pic_theta_max=' + str(theta_max),
             'time/tlim=' + str(_GYRO_TEND/4.0),
             'time/nlim=1000',
             'output2/dcycle=0'])
        gyro_cases.append(_measure_gyro_case(basename))
    _RESULTS['gyro'] = gyro_cases

    cpaw_cases = {}
    for label, rate in _CPAW_RATES.items():
        basename = 'pic_box_appendix_local_cpaw_' + label
        _remove_outputs(basename)
        _run_athena(
            'cpaw_' + label, _CPAW_INPUT,
            ['job/basename=' + basename,
             'particles/pic_expansion_rate_x1=' + str(rate)])
        cpaw_cases[label] = _measure_cpaw_case(basename)
    reference = cpaw_cases['static']
    _RESULTS['cpaw'] = {
        'static': _summarize_static_cpaw(reference),
        'expanding': _summarize_cpaw(
            reference, cpaw_cases['expanding'], _CPAW_RATES['expanding']),
        'compressing': _summarize_cpaw(
            reference, cpaw_cases['compressing'], _CPAW_RATES['compressing']),
    }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True
    gyro_cases = _RESULTS['gyro']
    for theta_max, measured in zip(_GYRO_THETA_VALUES, gyro_cases):
        logger.info('gyro theta=% .3e metrics=%s', theta_max, measured)
        ok = measured['npoint'] == 64 and ok
        ok = abs(measured['time'] - _GYRO_TEND) <= 1.0e-12 and ok
        ok = measured['history_pnorm_rel_error'] <= 3.0e-6 and ok
        ok = measured['history_gamma_abs_error'] <= 1.0e-6 and ok
        ok = measured['momentum_spread'] <= 2.0e-7 and ok
    phase_errors = np.asarray(
        [case['endpoint_phase_abs_error'] for case in gyro_cases])
    timesteps = np.asarray([case['dt'] for case in gyro_cases])
    phase_orders = np.log(phase_errors[:-1]/phase_errors[1:]) / np.log(
        timesteps[:-1]/timesteps[1:])
    logger.info('gyro endpoint_phase_errors=%s phase_orders=%s',
                phase_errors, phase_orders)
    ok = bool(np.all(np.diff(timesteps) < 0.0)) and ok
    ok = bool(np.all(np.diff(phase_errors) < 0.0)) and ok
    ok = float(np.min(phase_orders)) >= 1.7 and ok
    ok = gyro_cases[-1]['history_phase_abs_error'] <= 2.0e-3 and ok

    static = _RESULTS['cpaw']['static']
    logger.info('cpaw static metrics=%s', static)
    ok = static['time'] > 0.0 and ok
    ok = static['v_mode_amplitude_rel_drift'] <= 2.0e-4 and ok
    ok = static['b_mode_amplitude_rel_drift'] <= 2.0e-4 and ok
    ok = static['translation_phase_abs_error'] <= 1.0e-3 and ok
    ok = static['vb_phase_lock_abs_error'] <= 1.0e-5 and ok
    for label in ['expanding', 'compressing']:
        measured = _RESULTS['cpaw'][label]
        logger.info('cpaw %s metrics=%s', label, measured)
        ok = measured['time_abs_error_vs_static'] <= 1.0e-14 and ok
        ok = max(measured['source_map_errors'].values()) <= 1.0e-6 and ok
        ok = measured['v_mode_amplitude_rel_error'] <= 1.0e-6 and ok
        ok = measured['b_mode_amplitude_rel_error'] <= 1.0e-6 and ok
        ok = measured['v_mode_phase_abs_error'] <= 1.0e-6 and ok
        ok = measured['b_mode_phase_abs_error'] <= 1.0e-6 and ok
    return ok
