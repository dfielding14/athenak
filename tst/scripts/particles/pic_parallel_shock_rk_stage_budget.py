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

_INPUT_DECK = 'tests/pic_parallel_shock_rk_stage_budget.athinput'
_INTEGRATORS = ('rk1', 'rk2', 'rk3')
_MACRO_MASS = 1.0e-3
_MASS_FLUX = 0.10 * 1.0 * (3.0 + 1.0) * 8.0
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    for subdir in ['bin', 'rst']:
        pattern = os.path.join(_athena_exe_dir(), subdir, basename + '.*')
        for fname in glob.glob(pattern):
            os.remove(fname)
    history = os.path.join(_athena_exe_dir(), basename + '.mhd.hst')
    if os.path.exists(history):
        os.remove(history)


def _latest_file(subdir, pattern):
    matches = sorted(glob.glob(os.path.join(_athena_exe_dir(), subdir, pattern)))
    if not matches:
        raise RuntimeError('No output files found for: ' + pattern)
    return matches[-1]


def _restart_parameter(path, block, name, converter):
    with open(path, 'rb') as handle:
        text = handle.read(65536).decode('latin1', errors='ignore')
    if '<par_end>' not in text:
        raise RuntimeError('Restart parameter header is missing <par_end>: ' + path)
    active_block = None
    for raw_line in text.split('<par_end>', 1)[0].splitlines():
        line = raw_line.strip()
        if not line or line.startswith('#'):
            continue
        if line.startswith('<') and line.endswith('>'):
            active_block = line[1:-1]
            continue
        if active_block != block or '=' not in line:
            continue
        key, value = line.split('=', 1)
        if key.strip() == name:
            return converter(value.split('#', 1)[0].strip())
    raise RuntimeError(f'Missing restart parameter <{block}>/{name} in ' + path)


def _run_case(integrator, subtraction):
    suffix = 'on' if subtraction else 'off'
    basename = f'pic_parallel_shock_rk_stage_budget_{integrator}_{suffix}'
    _remove_outputs(basename)
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'time/integrator=' + integrator,
        'problem/ps_enable_gas_subtraction=' + str(subtraction).lower(),
    ]
    logger.info('Executing %s: %s', basename, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if proc.returncode != 0:
        raise RuntimeError('Command failed for ' + basename + '\n' + output)

    restart = _latest_file('rst', basename + '.*.rst')
    snapshot = _latest_file('bin', basename + '.mhd_u.*.bin')
    density = np.asarray(bin_convert.read_binary_as_athdf(snapshot)['dens'])
    history = np.loadtxt(os.path.join(_athena_exe_dir(), basename + '.mhd.hst'))
    if history.ndim == 1:
        history = history.reshape(1, -1)
    consumed_dt = float(history[-1, 0] - history[0, 0])
    if consumed_dt <= 0.0:
        raise RuntimeError('History did not record a positive consumed timestep')
    return {
        'next_tag': _restart_parameter(restart, 'problem', 'ps_next_tag', int),
        'reservoir': _restart_parameter(
            restart, 'problem', 'ps_mass_reservoir_global', float),
        'gas_mass': float(np.sum(density)),
        'dt': consumed_dt,
    }


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for integrator in _INTEGRATORS:
        off = _run_case(integrator, False)
        on = _run_case(integrator, True)
        _RESULTS[integrator] = {
            'off': off,
            'on': on,
            'gas_mass_removed': off['gas_mass'] - on['gas_mass'],
        }


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True
    reference = None
    for integrator in _INTEGRATORS:
        result = _RESULTS[integrator]
        off = result['off']
        on = result['on']
        created_mass = off['next_tag'] * _MACRO_MASS
        budget = created_mass + off['reservoir']
        expected_budget = _MASS_FLUX * off['dt']
        gas_mass_removed = result['gas_mass_removed']
        logger.info(
            '%s dt=% .12e next_tag=%d reservoir=% .12e budget=% .12e '
            'expected_budget=% .12e gas_mass_removed=% .12e',
            integrator, off['dt'], off['next_tag'], off['reservoir'], budget,
            expected_budget, gas_mass_removed)
        ok = abs(budget - expected_budget) <= 2.0e-12 and ok
        ok = abs(gas_mass_removed - created_mass) <= 1.0e-5 and ok
        ok = on['next_tag'] == off['next_tag'] and ok
        ok = abs(on['reservoir'] - off['reservoir']) <= 2.0e-12 and ok
        if reference is None:
            reference = (off['next_tag'], off['reservoir'], gas_mass_removed)
        else:
            ok = off['next_tag'] == reference[0] and ok
            ok = abs(off['reservoir'] - reference[1]) <= 2.0e-12 and ok
            ok = abs(gas_mass_removed - reference[2]) <= 2.0e-12 and ok
    return ok
