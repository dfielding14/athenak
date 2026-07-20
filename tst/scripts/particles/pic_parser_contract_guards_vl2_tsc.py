import logging
import os
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_parser_contract_guards.athinput'
_RESULTS = {}

_VL2 = [
    'time/integrator=vl2',
    'particles/deposit_order=2',
]

_TEST_PARTICLE = [
    'particles/pic_background_mode=passive_mhd',
    'particles/pic_feedback_mode=test_particle',
    'particles/deposit_moments=false',
    'particles/couple_moments_to_mhd=false',
    'particles/couple_moments_momentum_to_mhd=false',
    'particles/couple_moments_energy_to_mhd=false',
]

_ADAPTIVE_DAMPED_BOX = [
    'particles/pic_deltaf_mode=physical',
    'particles/pic_deltaf_f0=kappa_aniso',
    'particles/pic_deltaf_kappa=2.0',
    'particles/pic_deltaf_adapt_mode=global_bikappa_moments_experimental',
    'particles/pic_deltaf_adapt_interval=1.0',
    'particles/pic_expanding_box_mode=on',
    'particles/pic_expansion_rate_x1=0.01',
    'particles/pic_expansion_rate_x2=0.02',
    'particles/pic_expansion_rate_x3=0.03',
    'particles/pic_wave_damping_mode=ion_neutral_friction',
    'particles/pic_ion_neutral_collision_rate=0.5',
]

_POSITIVE_CASES = [
    ('coupled_rk1', [], [
        'integrator=rk1',
        'state=momentum_p_over_m',
        'background=coupled',
        'feedback=coupled',
    ]),
    ('coupled_vl2', _VL2, [
        'integrator=vl2',
        'state=momentum_p_over_m',
        'deposition=tsc',
        'background=coupled',
        'feedback=coupled',
    ]),
    ('test_particle', _TEST_PARTICLE, [
        'background=passive_mhd',
        'feedback=test_particle',
    ]),
    ('adaptive_damped_box', _ADAPTIVE_DAMPED_BOX, [
        'deltaf_adapt=global_bikappa_moments_experimental',
        'expanding_box=on',
        'wave_damping=ion_neutral_friction',
    ]),
]

_REJECTION_CASES = [
    ('unsupported_initial_state',
     ['particles/pic_cr_initial_state=bad_state'],
     'Unsupported value for <particles>/pic_cr_initial_state: bad_state'),
    ('unsupported_hall_mode',
     ['particles/pic_cr_hall_mode=bad_mode'],
     'Unsupported value for <particles>/pic_cr_hall_mode: bad_mode'),
    ('vl2_requires_boris_tsc',
     _VL2 + [
         'particles/pusher=drift',
         'particles/pic_cr_initial_state=velocity',
         'particles/pic_cr_light_speed=1.0',
     ],
     '<time>/integrator=vl2 with particles requires '
     '<particles>/particle_type=cosmic_ray and <particles>/pusher=boris_tsc'),
    ('vl2_requires_conservative_feedback',
     _VL2 + ['particles/couple_moments_momentum_to_mhd=false'],
     '<time>/integrator=vl2 requires coupled MHD background, coupled feedback, '
     'moment deposition, and conservative momentum feedback'),
    ('vl2_requires_tsc_deposition',
     ['time/integrator=vl2'],
     '<time>/integrator=vl2 requires <particles>/deposit_order=2'),
    ('vl2_rejects_direct_current_ct',
     _VL2 + [
         'particles/couple_j_to_efield_representation=edge_staggered',
         'particles/couple_j_deposition_mode=direct_staggered',
     ],
     '<time>/integrator=vl2 rejects direct-current CT induction options'),
    ('vl2_requires_mhd_source_feedback',
     _VL2 + ['particles/couple_fluid_feedback_order=efield_src'],
     '<time>/integrator=vl2 requires '
     '<particles>/couple_fluid_feedback_order=mhd_src_terms'),
    ('vl2_requires_unit_feedback',
     _VL2 + ['particles/couple_moments_momentum_coeff=0.5'],
     '<time>/integrator=vl2 requires unit conservative momentum and enabled '
     'energy feedback coefficients'),
    ('vl2_rejects_expanding_box',
     _VL2 + ['particles/pic_expanding_box_mode=on'],
     '<time>/integrator=vl2 with <particles>/pic_expanding_box_mode=on is not '
     'yet supported'),
    ('full_hall_requires_vl2',
     ['particles/pic_cr_hall_mode=full'],
     '<particles>/pic_cr_hall_mode=full requires <time>/integrator=vl2'),
    ('wave_damping_requires_positive_rate',
     ['particles/pic_wave_damping_mode=ion_neutral_friction'],
     '<particles>/pic_wave_damping_mode=ion_neutral_friction requires '
     '<particles>/pic_ion_neutral_collision_rate > 0'),
    ('adaptive_deltaf_requires_bounded_configuration',
     ['particles/pic_deltaf_adapt_mode=global_bikappa_moments_experimental'],
     '<particles>/pic_deltaf_adapt_mode=global_bikappa_moments_experimental '
     'requires physical kappa_aniso delta-f'),
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _execute(label, arguments):
    command = ['./athena', '-i', _athena_input_path(), 'time/nlim=0'] + arguments
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    return proc.returncode, (proc.stdout or '') + (proc.stderr or '')


def _run_success(label, arguments, expected_tokens):
    code, output = _execute(label, arguments)
    if code != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)
    for token in expected_tokens:
        if token not in output:
            raise RuntimeError(label + ' missing runtime token: ' + token)


def _run_expect_fail(label, arguments, expected_message):
    code, output = _execute(label, arguments)
    if code == 0:
        raise RuntimeError('Guard unexpectedly passed: ' + label)
    if expected_message not in output:
        raise RuntimeError(label + ' missing guard reason\n'
                           'Expected substring: ' + expected_message + '\n'
                           'Output:\n' + output)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for label, arguments, expected_tokens in _POSITIVE_CASES:
        _run_success(label, arguments, expected_tokens)
    for label, arguments, expected_message in _REJECTION_CASES:
        _run_expect_fail(label, arguments, expected_message)
    _RESULTS['positive_cases'] = len(_POSITIVE_CASES)
    _RESULTS['rejection_cases'] = len(_REJECTION_CASES)


def analyze():
    logger.info('PIC parser contract guards: %s', _RESULTS)
    return (_RESULTS.get('positive_cases') == len(_POSITIVE_CASES)
            and _RESULTS.get('rejection_cases') == len(_REJECTION_CASES))
