import logging
import os
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_parser_contract_guards.athinput'
_Q006_RUNTIME_LOCAL_INPUT_DECK = (
    'tests/pic_q006_paper_multispecies_oscillation_uniform_runtime_local.athinput'
)
_Q006_RUNTIME_LOCAL_CASES = {
    'paper_mhd_pic_isothermal_fullf_momentum_only',
    'paper_mhd_pic_isothermal_rejects_energy_feedback',
    'paper_mhd_pic_isothermal_fullf_rejects_non_q006_generator',
}
_RESULTS = {}

_PAPER_TEST_PARTICLE = [
    'particles/pic_physical_mode=paper_test_particle',
    'particles/pic_background_mode=no_mhd',
    'particles/pic_feedback_mode=test_particle',
    'particles/couple_moments_to_mhd=false',
    'particles/couple_moments_momentum_to_mhd=false',
    'particles/couple_moments_energy_to_mhd=false',
]

_ADAPTIVE_DELTAF = [
    'particles/pic_background_mode=passive_mhd',
    'particles/pic_feedback_mode=test_particle',
    'particles/deposit_moments=false',
    'particles/couple_moments_to_mhd=false',
    'particles/couple_moments_momentum_to_mhd=false',
    'particles/couple_moments_energy_to_mhd=false',
    'particles/pic_deltaf_mode=physical',
    'particles/pic_deltaf_f0=kappa_aniso',
    'particles/pic_deltaf_kappa=2.0',
    'particles/pic_deltaf_adapt_mode=global_bikappa_moments_experimental',
    'particles/pic_deltaf_adapt_interval=1.0',
    'particles/pic_expanding_box_mode=on',
]

_EXPANDING_COUPLED_PHYSICAL_DAMPING = [
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
    ('extended_mhd_pic', [], ['physical_mode=extended_mhd_pic']),
    ('paper_mhd_pic',
     ['particles/pic_physical_mode=paper_mhd_pic'],
     ['physical_mode=paper_mhd_pic']),
    ('paper_mhd_pic_isothermal_fullf_momentum_only',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'mhd/eos=isothermal',
      'particles/couple_moments_energy_to_mhd=false'],
     ['physical_mode=paper_mhd_pic']),
    ('paper_test_particle',
     _PAPER_TEST_PARTICLE,
     ['physical_mode=paper_test_particle']),
    ('extended_adaptive_deltaf',
     _ADAPTIVE_DELTAF,
     ['physical_mode=extended_mhd_pic',
      'deltaf_adapt=global_bikappa_moments_experimental']),
    ('extended_expanding_coupled_physical_damping',
     _EXPANDING_COUPLED_PHYSICAL_DAMPING,
     ['physical_mode=extended_mhd_pic',
      'background=coupled',
      'feedback=coupled',
      'deltaf_adapt=global_bikappa_moments_experimental',
      'expanding_box=on',
      'wave_damping=ion_neutral_friction']),
]

_REJECTION_CASES = [
    ('unsupported_physical_mode',
     ['particles/pic_physical_mode=bad_mode'],
     'Unsupported value for <particles>/pic_physical_mode: bad_mode'),
    ('unsupported_initial_state',
     ['particles/pic_cr_initial_state=bad_state'],
     'Unsupported value for <particles>/pic_cr_initial_state: bad_state'),
    ('unsupported_hall_mode',
     ['particles/pic_cr_hall_mode=bad_mode'],
     'Unsupported value for <particles>/pic_cr_hall_mode: bad_mode'),
    ('unsupported_wave_damping_mode',
     ['particles/pic_wave_damping_mode=bad_mode'],
     'Unsupported value for <particles>/pic_wave_damping_mode: bad_mode'),
    ('unsupported_deltaf_adapt_mode',
     ['particles/pic_deltaf_adapt_mode=bad_mode'],
     'Unsupported value for <particles>/pic_deltaf_adapt_mode: bad_mode'),
    ('unsupported_expansion_law',
     ['particles/pic_expansion_law=bad_law'],
     'Unsupported value for <particles>/pic_expansion_law: bad_law'),
    ('unsupported_current_representation',
     ['particles/couple_j_to_efield_representation=bad_representation'],
     'Unsupported value for <particles>/couple_j_to_efield_representation: '
     'bad_representation'),
    ('unsupported_current_deposition',
     ['particles/couple_j_deposition_mode=bad_deposition'],
     'Unsupported value for <particles>/couple_j_deposition_mode: '
     'bad_deposition'),
    ('unsupported_feedback_order',
     ['particles/couple_fluid_feedback_order=bad_order'],
     'Unsupported value for <particles>/couple_fluid_feedback_order: '
     'bad_order'),
    ('engineering_requires_velocity',
     ['particles/pic_physical_mode=engineering',
      'particles/pic_cr_light_speed=1.0',
      'particles/pic_cr_initial_state=momentum'],
     '<particles>/pic_physical_mode=engineering requires '
     '<particles>/pic_cr_initial_state=velocity'),
    ('physical_mode_requires_boris_tsc',
     _PAPER_TEST_PARTICLE + ['particles/pusher=drift'],
     '<particles>/pic_physical_mode=paper_test_particle requires '
     '<particles>/pusher=boris_tsc'),
    ('paper_test_particle_rejects_feedback',
     ['particles/pic_physical_mode=paper_test_particle'],
     '<particles>/pic_physical_mode=paper_test_particle requires '
     'test-particle feedback with all particle-to-MHD coupling toggles disabled'),
    ('paper_mhd_pic_requires_conservative_feedback',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'particles/couple_moments_momentum_to_mhd=false'],
     '<particles>/pic_physical_mode=paper_mhd_pic requires coupled MHD '
     'background, coupled feedback, moment deposition, and conservative '
     'momentum feedback'),
    ('paper_mhd_pic_isothermal_rejects_energy_feedback',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'mhd/eos=isothermal',
      'particles/couple_moments_energy_to_mhd=true',
      'particles/couple_moments_energy_coeff=1.0'],
     '<particles>/couple_moments_energy_to_mhd=true requires <mhd>/eos=ideal'),
    ('paper_mhd_pic_isothermal_fullf_rejects_non_q006_generator',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'mhd/eos=isothermal',
      'particles/couple_moments_energy_to_mhd=false',
      'problem/pgen_name=linear_wave'],
     '<particles>/pic_physical_mode=paper_mhd_pic requires coupled MHD '
     'background, coupled feedback, moment deposition, and conservative '
     'momentum feedback'),
    ('paper_mhd_pic_ideal_fullf_rejects_momentum_only',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'particles/couple_moments_energy_to_mhd=false',
      'particles/couple_moments_energy_coeff=0.0'],
     '<particles>/pic_physical_mode=paper_mhd_pic requires coupled MHD '
     'background, coupled feedback, moment deposition, and conservative '
     'momentum feedback'),
    ('paper_mhd_pic_rejects_direct_ct',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'particles/couple_j_to_efield_representation=edge_staggered',
      'particles/couple_j_deposition_mode=direct_staggered'],
     '<particles>/pic_physical_mode=paper_mhd_pic rejects direct-current CT '
     'induction options'),
    ('paper_mhd_pic_rejects_relativistic_mhd_background',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'coord/special_rel=true'],
     'Fluid momentum/energy feedback is limited to non-relativistic MHD in PR2'),
    ('paper_mode_rejects_hall_extension',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'particles/pic_cr_hall_mode=current_to_ct_experimental'],
     '<particles>/pic_cr_hall_mode=current_to_ct_experimental requires '
     '<particles>/pic_physical_mode=extended_mhd_pic'),
    ('paper_mode_rejects_wave_damping_extension',
     ['particles/pic_physical_mode=paper_mhd_pic',
      'particles/pic_wave_damping_mode=ion_neutral_friction',
      'particles/pic_ion_neutral_collision_rate=0.5'],
     '<particles>/pic_wave_damping_mode=ion_neutral_friction requires '
     '<particles>/pic_physical_mode=extended_mhd_pic'),
    ('paper_mode_rejects_adaptive_deltaf_extension',
     ['particles/pic_physical_mode=paper_mhd_pic'] + _ADAPTIVE_DELTAF,
     '<particles>/pic_deltaf_adapt_mode='
     'global_bikappa_moments_experimental requires '
     '<particles>/pic_physical_mode=extended_mhd_pic'),
    ('hall_extension_requires_coupled_moments',
     ['particles/pic_cr_hall_mode=current_to_ct_experimental',
      'particles/couple_moments_to_mhd=false',
      'particles/couple_moments_momentum_to_mhd=false',
      'particles/couple_moments_energy_to_mhd=false'],
     '<particles>/pic_cr_hall_mode=current_to_ct_experimental requires '
     '<particles>/couple_moments_to_mhd=true'),
    ('wave_damping_requires_positive_rate',
     ['particles/pic_wave_damping_mode=ion_neutral_friction'],
     '<particles>/pic_wave_damping_mode=ion_neutral_friction requires '
     '<particles>/pic_ion_neutral_collision_rate > 0'),
    ('wave_damping_requires_coupled_mhd',
     ['particles/pic_wave_damping_mode=ion_neutral_friction',
      'particles/pic_ion_neutral_collision_rate=0.5',
      'particles/pic_background_mode=no_mhd',
      'particles/pic_feedback_mode=test_particle',
      'particles/couple_moments_to_mhd=false',
      'particles/couple_moments_momentum_to_mhd=false',
      'particles/couple_moments_energy_to_mhd=false'],
     '<particles>/pic_wave_damping_mode=ion_neutral_friction requires '
     'an active coupled <mhd> background'),
    ('adaptive_deltaf_requires_bounded_configuration',
     ['particles/pic_deltaf_adapt_mode=global_bikappa_moments_experimental'],
     '<particles>/pic_deltaf_adapt_mode=global_bikappa_moments_experimental '
     'requires physical kappa_aniso delta-f'),
    ('expanding_box_rejects_direct_current_deposition',
     ['particles/pic_expanding_box_mode=on',
      'particles/couple_j_to_efield_representation=edge_staggered',
      'particles/couple_j_deposition_mode=direct_staggered'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     'staggered or direct CR-current deposition'),
    ('expanding_box_rejects_hall_current_induction',
     ['particles/pic_expanding_box_mode=on',
      'particles/pic_cr_hall_mode=current_to_ct_experimental'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     'particle feedback outside the admitted cell-centered source splits'),
    ('expanding_box_rejects_efield_feedback_order',
     ['particles/pic_expanding_box_mode=on',
      'particles/couple_fluid_feedback_order=efield_src'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     'particle feedback outside the admitted cell-centered source splits'),
    ('expanding_box_rejects_nonadaptive_deltaf_feedback',
     ['particles/pic_expanding_box_mode=on',
      'particles/pic_deltaf_mode=physical',
      'particles/pic_deltaf_f0=kappa_aniso'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     'particle feedback outside the admitted cell-centered source splits'),
    ('expanding_box_rejects_no_mhd_guard_bypass',
     ['particles/pic_expanding_box_mode=on',
      'particles/pic_background_mode=no_mhd',
      'particles/pic_feedback_mode=test_particle',
      'particles/deposit_moments=false',
      'particles/couple_moments_to_mhd=false',
      'particles/couple_moments_momentum_to_mhd=false',
      'particles/couple_moments_energy_to_mhd=false'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     '<particles>/pic_background_mode=no_mhd with an active <mhd> block'),
    ('expanding_box_rejects_unmapped_deltaf_background_current',
     ['particles/pic_expanding_box_mode=on',
      'particles/pic_deltaf_mode=physical',
      'particles/pic_deltaf_f0=kappa_aniso',
      'particles/pic_deltaf_background_jx=0.1'],
     '<particles>/pic_expanding_box_mode=on with active MHD does not support '
     'nonzero delta-f background current'),
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(label):
    input_deck = (_Q006_RUNTIME_LOCAL_INPUT_DECK
                  if label in _Q006_RUNTIME_LOCAL_CASES else _INPUT_DECK)
    return '../../' + athena.athena_rel_path + 'inputs/' + input_deck


def _execute(label, arguments):
    command = ['./athena', '-i', _athena_input_path(label), 'time/nlim=0'] + arguments
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
