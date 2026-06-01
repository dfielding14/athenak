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

_INPUT_DECK = 'tests/pic_parallel_shock_rk_stage_budget.athinput'
_INTEGRATORS = ('rk1', 'rk2', 'rk3')
_MACRO_MASS = 1.0e-3
_MASS_FLUX = 0.10 * 1.0 * (3.0 + 1.0) * 8.0
_RESULTS = {}
_FLOOR_REJECTION = {}
_INTEGRATOR_REJECTION = {}
_PARAMETER_REJECTIONS = {}
_FEEDBACK_DIAG_RE = re.compile(
    r'pic_parallel_shock feedback_diag: .*?'
    r'j_rms=\(([^,]+),([^,]+),([^)]+)\).*?'
    r'dpdt_rms=\(([^,]+),([^,]+),([^)]+)\).*?'
    r'dedt_rms=([^ ]+)'
)
_SOURCE_TRANSACTION_DIAG_RE = re.compile(
    r'pic_parallel_shock source_transaction_diag: .*?'
    r'applied=\(([^,]+),([^,]+),([^,]+),([^,]+),([^)]+)\) .*?'
    r'expected=\(([^,]+),([^,]+),([^,]+),([^,]+),([^)]+)\)'
)


def _athena_exe_dir():
    return os.environ.get(
        'ATHENA_PIC_PARALLEL_SHOCK_RK_EXE_DIR',
        os.path.join(os.getcwd(), 'build', 'src'),
    )


def _athena_input_path():
    return os.path.abspath(os.path.join(
        os.path.dirname(__file__), '..', '..', '..', 'inputs', _INPUT_DECK,
    ))


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
    fields = bin_convert.read_binary_as_athdf(snapshot)
    density = np.asarray(fields['dens'])
    momentum = np.array([
        np.sum(np.asarray(fields[name])) for name in ('mom1', 'mom2', 'mom3')
    ])
    energy = float(np.sum(np.asarray(fields['ener'])))
    history = np.loadtxt(os.path.join(_athena_exe_dir(), basename + '.mhd.hst'))
    if history.ndim == 1:
        history = history.reshape(1, -1)
    consumed_dt = float(history[-1, 0] - history[0, 0])
    if consumed_dt <= 0.0:
        raise RuntimeError('History did not record a positive consumed timestep')
    match = _FEEDBACK_DIAG_RE.search(output)
    if match is None:
        raise RuntimeError('Missing birth-cycle feedback diagnostic for ' + basename)
    source_match = _SOURCE_TRANSACTION_DIAG_RE.search(output)
    if subtraction and source_match is None:
        raise RuntimeError('Missing source-transaction diagnostic for ' + basename)
    source_transaction = None
    if source_match is not None:
        values = np.array([float(value) for value in source_match.groups()])
        source_transaction = {
            'applied': values[:5],
            'expected': values[5:],
        }
    return {
        'next_tag': _restart_parameter(restart, 'problem', 'ps_next_tag', int),
        'reservoir': _restart_parameter(
            restart, 'problem', 'ps_mass_reservoir_global', float),
        'injected_ledger': np.array([
            _restart_parameter(restart, 'problem', name, float)
            for name in (
                'ps_injected_cr_count_global',
                'ps_injected_cr_mass_global',
                'ps_injected_cr_momentum_x1_global',
                'ps_injected_cr_momentum_x2_global',
                'ps_injected_cr_momentum_x3_global',
                'ps_injected_cr_energy_global',
            )
        ]),
        'gas_mass': float(np.sum(density)),
        'gas_momentum': momentum,
        'gas_energy': energy,
        'dt': consumed_dt,
        'feedback_diag': np.array([float(value) for value in match.groups()]),
        'source_transaction': source_transaction,
    }


def _run_expected_floor_rejection():
    basename = 'pic_parallel_shock_rk_stage_budget_floor_reject'
    _remove_outputs(basename)
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'time/integrator=rk1',
        'problem/ps_p0=0.10',
        'problem/ps_eta=0.30',
        'problem/ps_enable_gas_subtraction=true',
    ]
    logger.info('Executing expected rejection: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    _FLOOR_REJECTION.update({
        'returncode': proc.returncode,
        'saw_floor_rejection': (
            'gas subtraction would violate a fluid floor' in output
            and 'process is stopping before clipping or checkpoint publication' in output
        ),
    })


def _run_expected_integrator_rejection():
    basename = 'pic_parallel_shock_rk_stage_budget_rk4_reject'
    _remove_outputs(basename)
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'time/integrator=rk4',
    ]
    logger.info('Executing expected rejection: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    _INTEGRATOR_REJECTION.update({
        'returncode': proc.returncode,
        'saw_integrator_rejection': (
            'injection is qualified only with time/integrator=rk1, rk2, or rk3'
            in output
        ),
    })


def _run_expected_parameter_rejection(label, overrides, reason):
    basename = 'pic_parallel_shock_rk_stage_budget_' + label
    _remove_outputs(basename)
    command = [
        './athena', '-i', _athena_input_path(),
        'job/basename=' + basename,
        'time/integrator=rk1',
    ] + overrides
    logger.info('Executing expected rejection: %s', ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    _PARAMETER_REJECTIONS[label] = {
        'returncode': proc.returncode,
        'saw_rejection': reason in output,
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
            'gas_momentum_removed': off['gas_momentum'] - on['gas_momentum'],
            'gas_energy_removed': off['gas_energy'] - on['gas_energy'],
        }
    _run_expected_floor_rejection()
    _run_expected_integrator_rejection()
    _run_expected_parameter_rejection(
        'negative_floor_reject',
        ['problem/ps_rho_floor_frac=-1.0'],
        'floor fractions must be finite and non-negative',
    )
    _run_expected_parameter_rejection(
        'nan_floor_reject',
        ['problem/ps_p_floor_frac=nan'],
        'floor fractions must be finite and non-negative',
    )
    _run_expected_parameter_rejection(
        'passive_subtraction_reject',
        [
            'particles/pic_background_mode=passive_mhd',
            'particles/pic_feedback_mode=test_particle',
            'particles/couple_moments_to_mhd=false',
            'particles/couple_moments_momentum_to_mhd=false',
            'particles/couple_moments_energy_to_mhd=false',
        ],
        'gas subtraction requires <particles>/pic_background_mode=coupled',
    )
    _run_expected_parameter_rejection(
        'source_transaction_relative_bound_reject',
        [
            'problem/ps_test_source_transaction_terms_override=true',
            'problem/ps_test_source_transaction_terms=3.0e15',
        ],
        'applied gas-subtraction transaction does not match the injected-particle ledger',
    )
    _run_expected_parameter_rejection(
        'source_transaction_nonfinite_terms_reject',
        [
            'problem/ps_test_source_transaction_terms_override=true',
            'problem/ps_test_source_transaction_terms=nan',
        ],
        'applied gas-subtraction transaction does not match the injected-particle ledger',
    )


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
        gas_momentum_removed = result['gas_momentum_removed']
        gas_energy_removed = result['gas_energy_removed']
        logger.info(
            '%s dt=% .12e next_tag=%d reservoir=% .12e budget=% .12e '
            'expected_budget=% .12e gas_mass_removed=% .12e '
            'gas_momentum_removed=%s gas_energy_removed=% .12e '
            'injected_ledger=%s',
            integrator, off['dt'], off['next_tag'], off['reservoir'], budget,
            expected_budget, gas_mass_removed, gas_momentum_removed,
            gas_energy_removed, off['injected_ledger'])
        # The history endpoint prints dt with fewer digits than restart metadata.
        ok = abs(budget - expected_budget) <= 2.0e-7 and ok
        ok = abs(gas_mass_removed - created_mass) <= 1.0e-5 and ok
        ok = on['next_tag'] == off['next_tag'] and ok
        ok = abs(on['reservoir'] - off['reservoir']) <= 2.0e-12 and ok
        ok = np.array_equal(on['injected_ledger'], off['injected_ledger']) and ok
        ok = off['injected_ledger'][0] == off['next_tag'] and ok
        ok = abs(off['injected_ledger'][1] - created_mass) <= 2.0e-12 and ok
        ok = np.all(np.isfinite(off['injected_ledger'])) and ok
        ok = off['injected_ledger'][5] > 0.0 and ok
        ok = np.linalg.norm(off['feedback_diag'][:3]) > 0.0 and ok
        ok = np.all(np.isfinite(gas_momentum_removed)) and ok
        ok = np.isfinite(gas_energy_removed) and gas_energy_removed > 0.0 and ok
        source_transaction = on['source_transaction']
        ok = source_transaction is not None and ok
        ok = np.max(np.abs(
            source_transaction['applied'] -
            source_transaction['expected'])) <= 1.0e-5 and ok
        ok = np.max(np.abs(
            source_transaction['expected'] - off['injected_ledger'][1:])) <= 1.0e-5 and ok
        if reference is None:
            reference = (
                off['next_tag'], off['reservoir'], off['injected_ledger'],
            )
        else:
            ok = off['next_tag'] == reference[0] and ok
            ok = abs(off['reservoir'] - reference[1]) <= 2.0e-12 and ok
            ok = np.array_equal(off['injected_ledger'], reference[2]) and ok
    ok = _FLOOR_REJECTION.get('returncode', 0) != 0 and ok
    ok = _FLOOR_REJECTION.get('saw_floor_rejection', False) and ok
    ok = _INTEGRATOR_REJECTION.get('returncode', 0) != 0 and ok
    ok = _INTEGRATOR_REJECTION.get('saw_integrator_rejection', False) and ok
    for rejection in _PARAMETER_REJECTIONS.values():
        ok = rejection.get('returncode', 0) != 0 and ok
        ok = rejection.get('saw_rejection', False) and ok
    return ok
