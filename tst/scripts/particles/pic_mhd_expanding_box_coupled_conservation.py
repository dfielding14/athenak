import logging
import os
import subprocess

import numpy as np

from scripts.particles import pic_paper_coupling_conservation as conservation

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_paper_coupling_conservation.athinput'
_RATES = np.asarray([0.05, 0.10, 0.15])
_C = 3.0
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../../inputs/' + _INPUT_DECK


def _run_case(label, feedback_coeff):
    basename = 'pic_mhd_box_coupled_' + label
    conservation._remove_outputs(basename)
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'time/nlim=1',
        'particles/pic_expanding_box_mode=on',
        'particles/pic_expansion_rate_x1=' + str(_RATES[0]),
        'particles/pic_expansion_rate_x2=' + str(_RATES[1]),
        'particles/pic_expansion_rate_x3=' + str(_RATES[2]),
        'particles/couple_moments_momentum_coeff=' + str(feedback_coeff),
        'particles/couple_moments_energy_coeff=' + str(feedback_coeff),
    ]
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)
    measured = conservation._measure_case(basename)
    measured['time1'] = float(measured['mhd']['bcc_final']['Time'])
    return measured


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _RESULTS['control'] = _run_case('control', 0.0)
    _RESULTS['feedback'] = _run_case('feedback', 1.0)


def _physical_mhd_delta(feedback, control, quantity, volume_scale):
    return volume_scale*(feedback['mhd'][quantity][1] -
                         control['mhd'][quantity][1])


def analyze():
    logger.debug('Analyzing test ' + __name__)
    control = _RESULTS['control']
    feedback = _RESULTS['feedback']
    time1 = feedback['time1']
    if abs(time1 - control['time1']) > 1.0e-14 or time1 <= 0.0:
        logger.warning('Control and feedback endpoint times disagree: %s', _RESULTS)
        return False

    scale = 1.0 + _RATES*time1
    ratio = 1.0/scale
    volume_scale = float(np.prod(scale))
    initial = feedback['particles_initial']
    final = feedback['particles_final']
    baseline_particle_momentum = initial['momentum']*ratio
    particle_delta = final['momentum'] - baseline_particle_momentum
    gas_delta = np.asarray([
        _physical_mhd_delta(feedback, control, 'mom1', volume_scale),
        _physical_mhd_delta(feedback, control, 'mom2', volume_scale),
        _physical_mhd_delta(feedback, control, 'mom3', volume_scale),
    ])
    momentum_residual = particle_delta + gas_delta

    baseline_p = initial['momentum']*ratio/initial['npoint']
    baseline_gamma = np.sqrt(1.0 + np.sum(baseline_p*baseline_p)/(_C*_C))
    baseline_energy = initial['npoint']*(baseline_gamma - 1.0)*_C*_C
    particle_energy_delta = final['energy'] - baseline_energy
    gas_energy_delta = _physical_mhd_delta(feedback, control, 'ener', volume_scale)
    energy_residual = particle_energy_delta + gas_energy_delta

    logger.info('expanding coupled particle_delta=%s gas_delta=%s residual=%s',
                particle_delta, gas_delta, momentum_residual)
    logger.info('expanding coupled particle_energy_delta=% .8e '
                'gas_energy_delta=% .8e residual=% .8e',
                particle_energy_delta, gas_energy_delta, energy_residual)
    return (
        np.max(np.abs(particle_delta)) > 1.0e-4
        and np.max(np.abs(momentum_residual)) <= 3.0e-6
        and abs(particle_energy_delta) > 1.0e-5
        and abs(energy_residual) <= 4.0e-5
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_mhd_expanding_box_coupled_conservation: FAIL')
    print('pic_mhd_expanding_box_coupled_conservation: PASS')
