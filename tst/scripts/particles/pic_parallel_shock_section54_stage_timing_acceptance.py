import logging
import os
import subprocess

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = (
    'tests/pic_parallel_shock_section54_stage_timing_acceptance.athinput'
)
_EXPECTED_MARKERS = (
    'physical_mode=paper_mhd_pic',
    'C=10000',
    'cycle=1 time=5.721748e-03',
    'q017.telemetry.particles.total=2.00000000000000000e+00',
)
_RESULTS = {}


def _athena_exe_dir():
    return os.environ.get(
        'ATHENA_PIC_PARALLEL_SHOCK_STAGE_TIMING_EXE_DIR',
        os.path.join(os.getcwd(), 'build', 'src'),
    )


def _athena_input_path():
    return os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', '..', '..', 'inputs', _INPUT_DECK,
    ))


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    command = ['./athena', '-i', _athena_input_path()]
    logger.info('Executing Section 5.4 stage-timing acceptance: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    _RESULTS.update({
        'returncode': proc.returncode,
        'output': output,
    })


def analyze():
    logger.debug('Analyzing test ' + __name__)
    output = _RESULTS.get('output', '')
    matched_markers = all(marker in output for marker in _EXPECTED_MARKERS)
    saw_floor_rejection = 'gas subtraction would violate a fluid floor' in output
    logger.info(
        'Section 5.4 stage-timing acceptance returncode=%s '
        'matched_markers=%s saw_floor_rejection=%s',
        _RESULTS.get('returncode'), matched_markers, saw_floor_rejection)
    return (
        _RESULTS.get('returncode') == 0
        and matched_markers
        and not saw_floor_rejection
    )
