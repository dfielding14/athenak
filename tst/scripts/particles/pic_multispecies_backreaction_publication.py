"""Legacy-named multispecies engineering proxy; not reproduction evidence."""

import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUTS = {
    'uniform': 'tests/pic_multispecies_osc_uniform_publication.athinput',
    'smr': 'tests/pic_multispecies_osc_smr_publication.athinput',
    'amr_publication': 'tests/pic_multispecies_osc_amr_publication.athinput',
}
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}
_AMPLITUDE_MIN = {
    'uniform': 10.0,
    'smr': 10.0,
    'amr_publication': 10.0,
}
_TURN_MIN = {
    'uniform': 4.0,
    'smr': 4.0,
}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path(input_rel):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_rel


def _athena_mpi_enabled():
    proc = subprocess.run(['./athena', '-c'], cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration with -c')
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _remove_outputs(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.*.bin')
    for fname in glob.glob(pattern):
        os.remove(fname)


def _integrate_quantity(dataset, quantity):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    dvol = dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]
    return float(np.sum(dataset[quantity] * dvol))


def _dominant_frequency_nonuniform(time, signal, fmin, fmax, nsample=4000):
    t = np.asarray(time, dtype=float)
    y = np.asarray(signal, dtype=float)
    if t.size != y.size or t.size < 8:
        raise RuntimeError('insufficient samples for frequency fit')

    y = y - np.mean(y)
    freqs = np.linspace(fmin, fmax, nsample)
    best_f = float(freqs[0])
    best_amp = -1.0

    for f in freqs:
        omega_t = 2.0 * np.pi * f * t
        design = np.column_stack((np.sin(omega_t), np.cos(omega_t)))
        coeff, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
        amp = float(np.sqrt(coeff[0] * coeff[0] + coeff[1] * coeff[1]))
        if amp > best_amp:
            best_amp = amp
            best_f = float(f)

    return best_f, best_amp


def _zero_crossings(signal):
    values = np.asarray(signal, dtype=float)
    shifted = values - np.mean(values)
    signs = np.sign(shifted)
    signs = signs[signs != 0.0]
    if signs.size < 2:
        return 0.0
    return float(np.count_nonzero(signs[1:] != signs[:-1]))


def _load_series(basename):
    pattern_m2 = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_u_m2.*.bin')
    files = sorted(glob.glob(pattern_m2))
    if len(files) < 20:
        raise RuntimeError('Not enough outputs for oscillation fit: ' + basename)

    times = []
    mom2 = []
    ener = []

    for m2_file in files:
        stem = os.path.basename(m2_file)
        cycle = stem.split('.')[-2]
        e_file = os.path.join(_athena_exe_dir(), 'bin',
                              basename + '.mhd_u_e.' + cycle + '.bin')

        m2_data = bin_convert.read_binary_as_athdf(m2_file)
        e_data = bin_convert.read_binary_as_athdf(e_file)

        times.append(float(m2_data['Time']))
        mom2.append(_integrate_quantity(m2_data, 'mom2'))
        ener.append(_integrate_quantity(e_data, 'ener'))

    t = np.asarray(times, dtype=float)
    m = np.asarray(mom2, dtype=float)
    e = np.asarray(ener, dtype=float)

    finite = np.all(np.isfinite(m)) and np.all(np.isfinite(e))
    if finite:
        freq, _ = _dominant_frequency_nonuniform(t, m, 0.08, 0.12)
        amp = float(np.max(m) - np.min(m))
        turns = _zero_crossings(m)
        edrift = float(np.max(np.abs(e - e[0])) / max(abs(e[0]), 1.0))
    else:
        freq = -1.0
        amp = -1.0
        turns = -1.0
        edrift = float('inf')

    return {
        'freq': float(freq),
        'amp': float(amp),
        'turns': float(turns),
        'edrift': float(edrift),
    }


def _run_case(label, nproc, input_rel, basename):
    args = ['job/basename=' + basename]
    command = ['./athena', '-i', _athena_input_path(input_rel)] + args
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _check_range(label, measured, lower, upper):
    logger.info('%s measured=% .8e range=[% .8e,% .8e]',
                label, measured, lower, upper)
    return (measured >= lower) and (measured <= upper)


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    for tag, input_rel in _INPUTS.items():
        base = 'pic_mso_pub_' + tag + '_np1'
        _remove_outputs(base)
        _run_case(tag + '_serial', 1, input_rel, base)
        _RESULTS[tag + '_np1'] = _load_series(base)

    if _athena_mpi_enabled():
        for tag, input_rel in _INPUTS.items():
            base = 'pic_mso_pub_' + tag + '_np2'
            _remove_outputs(base)
            _run_case(tag + '_mpi2', 2, input_rel, base)
            _RESULTS[tag + '_np2'] = _load_series(base)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for tag in _INPUTS:
        res = _RESULTS[tag + '_np1']
        ok = _check_lower(tag + ':amp_np1',
                          res['amp'], _AMPLITUDE_MIN[tag]) and ok
        if tag in _TURN_MIN:
            ok = _check_lower(tag + ':turns_np1',
                              res['turns'], _TURN_MIN[tag]) and ok
        ok = _check_upper(tag + ':edrift_np1', res['edrift'], 5.0e-2) and ok

    if 'uniform_np2' in _RESULTS:
        for tag in _INPUTS:
            np1 = _RESULTS[tag + '_np1']
            np2 = _RESULTS[tag + '_np2']
            ok = _check_lower(tag + ':amp_np2',
                              np2['amp'], _AMPLITUDE_MIN[tag]) and ok
            if tag in _TURN_MIN:
                ok = _check_lower(tag + ':turns_np2',
                                  np2['turns'], _TURN_MIN[tag]) and ok
            ok = _check_upper(tag + ':edrift_np2', np2['edrift'], 5.0e-2) and ok
            ok = _check_close(tag + ':amp_np1_vs_np2',
                              np2['amp'], np1['amp'], 1.0, 1.0e-2) and ok
            ok = _check_close(tag + ':turns_np1_vs_np2',
                              np2['turns'], np1['turns'], 0.0, 0.0) and ok

    return ok
