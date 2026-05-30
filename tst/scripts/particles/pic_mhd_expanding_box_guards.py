import logging
import os
from pathlib import Path
import subprocess
import tempfile

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_mhd_expanding_box_uniform.athinput'
_RESULTS = {}

_REJECTION_CASES = [
    (
        'initial_turb',
        '''
<initial_turb>
''',
        'does not support coupled fluid, radiation, relativity, or turbulence blocks',
    ),
    (
        'user_history_callback',
        '''
<problem>
user_hist = true
''',
        'does not support user-defined history callbacks',
    ),
    (
        'drift_pusher',
        '''
<particles>
pusher = drift
''',
        'does not support particle types or pushers other than cosmic-ray Boris '
        'pushers',
    ),
    (
        'contracting_particles',
        '''
<particles>
ppc                   = 1.0
pic_expansion_rate_x1 = -0.05
''',
        'does not support contracting scale factors with active particles',
    ),
]


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return Path('../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK)


def _expect_rejection(label, appended_input, expected_message):
    source = (Path(_athena_exe_dir()) / _athena_input_path()).resolve()
    with tempfile.TemporaryDirectory() as tmpdir:
        deck = Path(tmpdir) / ('pic_mhd_expanding_box_' + label + '.athinput')
        deck.write_text(source.read_text(encoding='ascii') + appended_input,
                        encoding='ascii')
        command = ['./athena', '-i', str(deck), 'time/nlim=0']
        logger.info('Executing: %s', ' '.join(command))
        proc = subprocess.run(command, cwd=_athena_exe_dir(),
                              capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode == 0:
        raise RuntimeError('Guard unexpectedly passed: ' + label)
    if expected_message not in output:
        raise RuntimeError(label + ' missing guard reason\n'
                           'Expected substring: ' + expected_message + '\n'
                           'Output:\n' + output)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    for label, appended_input, expected_message in _REJECTION_CASES:
        _expect_rejection(label, appended_input, expected_message)
    _RESULTS['rejection_cases'] = len(_REJECTION_CASES)


def analyze():
    logger.info('expanding-MHD parser guards: %s', _RESULTS)
    return _RESULTS.get('rejection_cases') == len(_REJECTION_CASES)
