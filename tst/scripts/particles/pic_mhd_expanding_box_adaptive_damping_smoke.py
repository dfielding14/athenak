import glob
import logging
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = '../../../inputs/tests/pic_mhd_expanding_box_adaptive_damping_smoke.athinput'
_BASENAME = 'pic_mhd_expanding_box_adaptive_damping'
_BACKGROUND_RHO = 1.0
_NU_IN = 0.5
_EXPANSION_RATES = (0.01, 0.02, 0.03)
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _remove_outputs(basename):
    hst = Path(_athena_exe_dir()) / (basename + '.mhd.hst')
    if hst.exists():
        hst.unlink()
    for path in glob.glob(os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')):
        os.remove(path)


def _last_output(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.' + file_id + '.*.bin')
    paths = sorted(glob.glob(pattern))
    if len(paths) < 2:
        raise RuntimeError('Need initial and final outputs for: ' + pattern)
    return bin_convert.read_binary_as_athdf(paths[-1])


def _run_case(label, background_rho, damping_mode, collision_rate):
    basename = _BASENAME + '_' + label
    _remove_outputs(basename)
    command = [
        './athena', '-i', _INPUT_DECK,
        'job/basename=' + basename,
        'particles/pic_deltaf_background_rho=' + str(background_rho),
        'particles/pic_wave_damping_mode=' + damping_mode,
        'particles/pic_ion_neutral_collision_rate=' + str(collision_rate),
    ]
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)
    for token in [
            'background=coupled',
            'feedback=coupled',
            'deltaf_adapt=global_bikappa_moments_experimental',
            'expanding_box=on',
            'wave_damping=' + damping_mode,
            'PIC adaptive delta-f fit:']:
        if token not in output:
            raise RuntimeError(label + ' missing runtime token: ' + token)
    rows = np.loadtxt(Path(_athena_exe_dir()) / (basename + '.mhd.hst'))
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    return {
        'history': rows,
        'u': _last_output(basename, 'mhd_u'),
        'w': _last_output(basename, 'mhd_w'),
        'bcc': _last_output(basename, 'mhd_bcc'),
    }


def _max_abs(actual, expected):
    return float(np.max(np.abs(np.asarray(actual) - np.asarray(expected))))


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _RESULTS['control_off'] = _run_case('control_off', 0.0, 'off', 0.0)
    _RESULTS['source_off'] = _run_case('source_off', _BACKGROUND_RHO, 'off', 0.0)
    _RESULTS['control_damped'] = _run_case(
        'control_damped', 0.0, 'ion_neutral_friction', _NU_IN)
    _RESULTS['source_damped'] = _run_case(
        'source_damped', _BACKGROUND_RHO, 'ion_neutral_friction', _NU_IN)


def _expected_damped_state(off):
    dt = float(off['u']['Time'])
    factor = float(np.exp(-_NU_IN*dt))
    u = off['u']
    return {
        'factor': factor,
        'dens': u['dens'],
        'mom1': u['mom1'],
        'mom2': factor*u['mom2'],
        'mom3': factor*u['mom3'],
        'ener': (u['ener'] - 0.5*(1.0 - factor*factor)*
                 (u['mom2']*u['mom2'] + u['mom3']*u['mom3'])/u['dens']),
    }


def analyze():
    control_off = _RESULTS['control_off']
    source_off = _RESULTS['source_off']
    control_damped = _RESULTS['control_damped']
    source_damped = _RESULTS['source_damped']
    dt = float(control_off['u']['Time'])
    volume_scale = float(np.exp(sum(_EXPANSION_RATES)*dt))
    physical_rho = _BACKGROUND_RHO/volume_scale

    w = control_off['w']
    b = control_off['bcc']
    ex = -(w['vely']*b['bcc3'] - w['velz']*b['bcc2'])
    ey = -(w['velz']*b['bcc1'] - w['velx']*b['bcc3'])
    ez = -(w['velx']*b['bcc2'] - w['vely']*b['bcc1'])
    expected_source_delta = {
        'mom1': -dt*physical_rho*ex,
        'mom2': -dt*physical_rho*ey,
        'mom3': -dt*physical_rho*ez,
        'ener': np.zeros_like(control_off['u']['ener']),
    }
    measured_source_delta = {
        name: source_off['u'][name] - control_off['u'][name]
        for name in expected_source_delta
    }
    source_errors = {
        name: _max_abs(measured_source_delta[name], expected_source_delta[name])
        for name in expected_source_delta
    }

    damping_errors = {}
    for prefix, off, damped in [
            ('control', control_off, control_damped),
            ('source', source_off, source_damped)]:
        expected = _expected_damped_state(off)
        for name in ['dens', 'mom1', 'mom2', 'mom3', 'ener']:
            damping_errors[prefix + ':' + name] = _max_abs(
                damped['u'][name], expected[name])

    factor = _expected_damped_state(source_off)['factor']
    expected_ordered_delta = {
        'mom1': measured_source_delta['mom1'],
        'mom2': factor*measured_source_delta['mom2'],
        'mom3': factor*measured_source_delta['mom3'],
    }
    ordered_errors = {
        name: _max_abs(
            source_damped['u'][name] - control_damped['u'][name],
            expected_ordered_delta[name])
        for name in expected_ordered_delta
    }
    histories_finite = all(
        case['history'].shape[0] >= 2
        and case['history'][-1, 0] > case['history'][0, 0]
        and np.all(np.isfinite(case['history']))
        for case in _RESULTS.values())
    source_norm = max(float(np.max(np.abs(value)))
                      for value in expected_source_delta.values())
    measured = {
        'dt': dt,
        'volume_scale': volume_scale,
        'physical_background_rho': physical_rho,
        'damping_factor': factor,
        'source_norm': source_norm,
        'source_errors': source_errors,
        'damping_errors': damping_errors,
        'ordered_errors': ordered_errors,
        'histories_finite': histories_finite,
    }
    logger.info('expanding adaptive damping manufactured-source metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (
        dt > 0.0
        and source_norm > 1.0e-5
        and max(source_errors.values()) <= 2.0e-6
        and max(damping_errors.values()) <= 2.0e-6
        and max(ordered_errors.values()) <= 2.0e-6
        and histories_finite
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_mhd_expanding_box_adaptive_damping_smoke: FAIL')
    print('pic_mhd_expanding_box_adaptive_damping_smoke: PASS')
