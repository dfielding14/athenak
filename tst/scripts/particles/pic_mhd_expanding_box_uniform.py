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

_INPUT_DECK = 'tests/pic_mhd_expanding_box_uniform.athinput'
_LAWS = ('linear', 'reciprocal_linear', 'exponential')
_RATES = (0.05, 0.10, 0.15)
_GAMMA = 1.66666666667
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)
    history = os.path.join(_athena_exe_dir(), basename + '.mhd.hst')
    if os.path.exists(history):
        os.remove(history)


def _output_files(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin',
                           f'{basename}.{file_id}.*.bin')
    matches = sorted(glob.glob(pattern))
    if len(matches) < 2:
        raise RuntimeError('Need at least two output snapshots for: ' + pattern)
    return matches


def _read_endpoints(basename, file_id):
    files = _output_files(basename, file_id)
    return (bin_convert.read_binary_as_athdf(files[0]),
            bin_convert.read_binary_as_athdf(files[-1]))


def _scale_factor(law, rate, time):
    if law == 'linear':
        return 1.0 + rate * time
    if law == 'reciprocal_linear':
        return 1.0 / (1.0 + rate * time)
    if law == 'exponential':
        return np.exp(rate * time)
    raise ValueError('Unsupported law: ' + law)


def _relative_error(actual, expected):
    scale = max(float(np.max(np.abs(expected))), 1.0e-14)
    return float(np.max(np.abs(actual - expected))) / scale


def _history_rows(basename):
    path = os.path.join(_athena_exe_dir(), basename + '.mhd.hst')
    rows = np.loadtxt(path)
    if rows.ndim == 1:
        rows = rows[np.newaxis, :]
    if rows.shape[0] < 2:
        raise RuntimeError('Need at least two history rows for: ' + path)
    return rows


def _measure_case(basename, law):
    ufirst, ulast = _read_endpoints(basename, 'mhd_u')
    wfirst, wlast = _read_endpoints(basename, 'mhd_w')
    bfirst, blast = _read_endpoints(basename, 'mhd_bcc')

    time0 = float(ufirst['Time'])
    time1 = float(ulast['Time'])
    a0 = np.asarray([_scale_factor(law, rate, time0) for rate in _RATES])
    a1 = np.asarray([_scale_factor(law, rate, time1) for rate in _RATES])
    ratio = a0 / a1
    rvol = float(np.prod(ratio))

    expected = {
        'dens': ufirst['dens'] * rvol,
        'mom1': ufirst['mom1'] * rvol * ratio[0],
        'mom2': ufirst['mom2'] * rvol * ratio[1],
        'mom3': ufirst['mom3'] * rvol * ratio[2],
        'velx': wfirst['velx'] * ratio[0],
        'vely': wfirst['vely'] * ratio[1],
        'velz': wfirst['velz'] * ratio[2],
        'eint': wfirst['eint'] * rvol**_GAMMA,
        'bcc1': bfirst['bcc1'] * rvol / ratio[0],
        'bcc2': bfirst['bcc2'] * rvol / ratio[1],
        'bcc3': bfirst['bcc3'] * rvol / ratio[2],
    }
    expected['ener'] = (
        expected['eint']
        + 0.5 * (expected['mom1']**2 + expected['mom2']**2
                 + expected['mom3']**2) / expected['dens']
        + 0.5 * (expected['bcc1']**2 + expected['bcc2']**2
                 + expected['bcc3']**2)
    )

    actual = {
        'dens': ulast['dens'],
        'mom1': ulast['mom1'],
        'mom2': ulast['mom2'],
        'mom3': ulast['mom3'],
        'ener': ulast['ener'],
        'velx': wlast['velx'],
        'vely': wlast['vely'],
        'velz': wlast['velz'],
        'eint': wlast['eint'],
        'bcc1': blast['bcc1'],
        'bcc2': blast['bcc2'],
        'bcc3': blast['bcc3'],
    }
    history = _history_rows(basename)
    physical_cell_volume = float(np.prod(a1))
    history_expected = np.asarray([
        physical_cell_volume * np.sum(actual['dens']),
        physical_cell_volume * np.sum(actual['mom1']),
        physical_cell_volume * np.sum(actual['mom2']),
        physical_cell_volume * np.sum(actual['mom3']),
        physical_cell_volume * np.sum(actual['ener']),
        physical_cell_volume * np.sum(0.5 * actual['mom1']**2 / actual['dens']),
        physical_cell_volume * np.sum(0.5 * actual['mom2']**2 / actual['dens']),
        physical_cell_volume * np.sum(0.5 * actual['mom3']**2 / actual['dens']),
        physical_cell_volume * np.sum(0.5 * actual['bcc1']**2),
        physical_cell_volume * np.sum(0.5 * actual['bcc2']**2),
        physical_cell_volume * np.sum(0.5 * actual['bcc3']**2),
    ])
    return {
        'time0': time0,
        'time1': time1,
        'errors': {name: _relative_error(actual[name], expected[name])
                   for name in expected},
        'history_errors': {
            'physical_integrals': _relative_error(history[-1, 2:], history_expected),
            'mass_conservation': _relative_error(history[-1, 2], history[0, 2]),
        },
    }


def _run_case(law, basename):
    args = [
        'job/basename=' + basename,
        'particles/pic_expansion_law=' + law,
    ]
    command = ['./athena', '-i', _athena_input_path()] + args
    logger.info('Executing %s: %s', law, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + law + '\n' + output)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for law in _LAWS:
        basename = 'pic_mhd_box_uniform_' + law
        _remove_outputs(basename)
        _run_case(law, basename)
        _RESULTS[law] = _measure_case(basename, law)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True
    tolerance = 5.0e-7
    for law in _LAWS:
        measured = _RESULTS[law]
        if measured['time1'] <= measured['time0']:
            logger.warning('%s did not advance time: %s', law, measured)
            ok = False
        for name, error in measured['errors'].items():
            logger.info('%s:%s relative_error=% .8e tolerance=% .8e',
                        law, name, error, tolerance)
            ok = (error <= tolerance) and ok
        for name, error in measured['history_errors'].items():
            logger.info('%s:history:%s relative_error=% .8e tolerance=% .8e',
                        law, name, error, tolerance)
            ok = (error <= tolerance) and ok
    return ok
