import glob
import logging
import os
import subprocess

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_parallel_shock_coarse_uniform.athinput'
_RESULTS = {}
_MISMATCH_REASON = (
    'pic_parallel_shock restart continuation-control fingerprint mismatch')


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    for subdir in ['bin', 'pvtk', 'rst']:
        pattern = os.path.join(_athena_exe_dir(), subdir, basename + '.*')
        for fname in glob.glob(pattern):
            os.remove(fname)


def _latest_restart_file(basename):
    pattern = os.path.join(_athena_exe_dir(), 'rst', basename + '.*.rst')
    matches = sorted(glob.glob(pattern))
    if len(matches) == 0:
        raise RuntimeError('No restart files found for pattern: ' + pattern)
    return os.path.relpath(matches[-1], _athena_exe_dir())


def _args(basename, nlim, rst_dcycle):
    return [
        'job/basename=' + basename,
        'mesh/nx1=16',
        'mesh/nx2=8',
        'mesh/nx3=1',
        'meshblock/nx1=8',
        'meshblock/nx2=8',
        'meshblock/nx3=1',
        'particles/ppc=0.0',
        'particles/deposit_moments=true',
        'particles/deposit_qscale=1.0e-3',
        'problem/ps_eta=0.1',
        'problem/ps_enable_injection=true',
        'problem/ps_enable_gas_subtraction=false',
        'problem/ps_enable_frame_tracking=true',
        'problem/ps_frame_mode=velocity',
        'problem/ps_frame_t_start=0.0',
        'problem/ps_frame_t_ramp=0.0',
        'problem/ps_frame_vfrac=0.0',
        'time/nlim=' + str(nlim),
        'time/tlim=1.0',
        'output1/dcycle=0',
        'output2/dcycle=0',
        'output3/dcycle=0',
        'output4/dcycle=' + str(rst_dcycle),
    ]


def _execute(label, arguments, restart_file=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    return proc.returncode, (proc.stdout or '') + (proc.stderr or '')


def _run_success(label, arguments, restart_file=None):
    code, output = _execute(label, arguments, restart_file=restart_file)
    if code != 0:
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _run_expect_mismatch(label, arguments, restart_file):
    code, output = _execute(label, arguments, restart_file=restart_file)
    if code == 0:
        raise RuntimeError('Changed restart control unexpectedly passed: ' +
                           label)
    if _MISMATCH_REASON not in output:
        raise RuntimeError(label + ' missing restart mismatch reason\n' + output)


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    basenames = [
        'pic_parallel_shock_restart_controls_seed',
        'pic_parallel_shock_restart_controls_same',
        'pic_parallel_shock_restart_controls_eta',
        'pic_parallel_shock_restart_controls_frame',
        'pic_parallel_shock_restart_controls_stencil',
        'pic_parallel_shock_restart_controls_surface_average',
    ]
    for basename in basenames:
        _remove_outputs(basename)

    _run_success('seed', _args(basenames[0], 1, 1))
    restart_file = _latest_restart_file(basenames[0])
    _run_success('unchanged_restart', _args(basenames[1], 2, 0),
                 restart_file=restart_file)
    _run_expect_mismatch(
        'changed_injection_control',
        _args(basenames[2], 2, 0) + ['problem/ps_eta=0.2'],
        restart_file)
    _run_expect_mismatch(
        'changed_frame_control',
        _args(basenames[3], 2, 0) + ['problem/ps_frame_vfrac=0.25'],
        restart_file)
    _run_expect_mismatch(
        'changed_subtraction_stencil_control',
        _args(basenames[4], 2, 0) +
        ['problem/ps_subtract_stencil_cells=3'],
        restart_file)
    _run_expect_mismatch(
        'changed_surface_averaged_subtraction_control',
        _args(basenames[5], 2, 0) +
        ['problem/ps_enable_surface_averaged_subtraction=true'],
        restart_file)

    _RESULTS['unchanged_restart'] = True
    _RESULTS['rejected_changed_injection_control'] = True
    _RESULTS['rejected_changed_frame_control'] = True
    _RESULTS['rejected_changed_subtraction_stencil_control'] = True
    _RESULTS['rejected_changed_surface_averaged_subtraction_control'] = True


def analyze():
    logger.info('PIC parallel-shock restart control guards: %s', _RESULTS)
    return all(_RESULTS.values()) and len(_RESULTS) == 5
