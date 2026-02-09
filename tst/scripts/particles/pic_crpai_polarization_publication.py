import glob
import logging
import os
import subprocess

import numpy as np
import scripts.utils.athena as athena
from scripts.particles.pic_analysis_utils import fit_exponential_growth_windowed
from scripts.particles.pic_analysis_utils import polarization_split

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_CASES = {
    'prolate': 'tests/pic_crpai_prolate_publication.athinput',
    'oblate': 'tests/pic_crpai_oblate_publication.athinput',
}
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


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


def _mode1_complex(dataset, field):
    values = np.asarray(dataset[field], dtype=float)
    x_mode = np.mean(values, axis=(0, 1))
    fluc = x_mode - np.mean(x_mode)

    x1f = np.asarray(dataset['x1f'], dtype=float)
    x1v = np.asarray(dataset['x1v'], dtype=float)
    length = float(x1f[-1] - x1f[0])
    phase = np.exp(-2.0j * np.pi * (x1v - x1f[0]) / length)
    return np.sum(fluc * phase) / x_mode.size


def _fit_branch(time, amp):
    floor = max(1.0e-30, 1.0e-4 * float(np.max(amp)))
    fit = fit_exponential_growth_windowed(
        time, amp, min_points=6, floor=floor, min_growth_factor=1.2
    )
    return {
        'gamma': float(fit['gamma']),
        'r2': float(fit['r2']),
        'tmin': float(fit['tmin']),
        'tmax': float(fit['tmax']),
        'growth_factor': float(fit['growth_factor']),
    }


def _load_growth_metrics(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.mhd_bcc.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 10:
        raise RuntimeError('Not enough mhd_bcc outputs for CRPAI publication fit')

    times = []
    by_mode = []
    bz_mode = []
    for fname in files:
        data = bin_convert.read_binary_as_athdf(fname)
        times.append(float(data['Time']))
        by_mode.append(_mode1_complex(data, 'bcc2'))
        bz_mode.append(_mode1_complex(data, 'bcc3'))

    t = np.asarray(times, dtype=float)
    by = np.asarray(by_mode, dtype=complex)
    bz = np.asarray(bz_mode, dtype=complex)
    right, left = polarization_split(by, bz)

    amp_r = np.abs(right)
    amp_l = np.abs(left)

    fit_r = _fit_branch(t, amp_r)
    fit_l = _fit_branch(t, amp_l)

    dom = 'right' if amp_r[-1] >= amp_l[-1] else 'left'
    asym = float((amp_r[-1] - amp_l[-1]) / max(amp_r[-1] + amp_l[-1], 1.0e-30))

    return {
        'dom_branch': dom,
        'asymmetry': asym,
        'gamma_r': float(fit_r['gamma']),
        'gamma_l': float(fit_l['gamma']),
        'r2_r': float(fit_r['r2']),
        'r2_l': float(fit_l['r2']),
        'growth_r': float(fit_r['growth_factor']),
        'growth_l': float(fit_l['growth_factor']),
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


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_log_ratio_close(label, measured, expected, tol_log):
    if measured <= 0.0 or expected <= 0.0:
        logger.error('%s requires positive values (measured=%e expected=%e)',
                     label, measured, expected)
        return False
    log_ratio = abs(np.log(measured / expected))
    logger.info('%s measured=% .8e expected=% .8e |log-ratio|=% .8e',
                label, measured, expected, log_ratio)
    return log_ratio <= tol_log


def run(**kwargs):
    logger.debug('Running test ' + __name__)

    for label, input_rel in _CASES.items():
        base = 'pic_crpai_pub_' + label + '_np1'
        _remove_outputs(base)
        _run_case(label + '_serial', 1, input_rel, base)
        _RESULTS[label + '_np1'] = _load_growth_metrics(base)

    if _athena_mpi_enabled():
        for label, input_rel in _CASES.items():
            base = 'pic_crpai_pub_' + label + '_np2'
            _remove_outputs(base)
            _run_case(label + '_mpi2', 2, input_rel, base)
            _RESULTS[label + '_np2'] = _load_growth_metrics(base)

        for label, input_rel in _CASES.items():
            base = 'pic_crpai_pub_' + label + '_np4'
            _remove_outputs(base)
            _run_case(label + '_mpi4', 4, input_rel, base)
            _RESULTS[label + '_np4'] = _load_growth_metrics(base)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for label in _CASES:
        case = _RESULTS[label + '_np1']
        ok = _check_lower(label + ':max_gamma_np1',
                          max(case['gamma_r'], case['gamma_l']), 5.0e-2) and ok
        ok = _check_lower(label + ':max_r2_np1',
                          max(case['r2_r'], case['r2_l']), 0.5) and ok
        ok = _check_lower(label + ':branch_split_np1',
                          abs(case['asymmetry']), 1.0e-2) and ok

    pro = _RESULTS['prolate_np1']
    obl = _RESULTS['oblate_np1']
    ok = _check_lower('np1:prolate_oblate_asym_separation',
                      abs(pro['asymmetry'] - obl['asymmetry']),
                      2.0e-2) and ok

    if 'prolate_np2' in _RESULTS:
        for label in _CASES:
            np1 = _RESULTS[label + '_np1']
            np2 = _RESULTS[label + '_np2']
            ok = _check_lower(label + ':max_gamma_np2',
                              max(np2['gamma_r'], np2['gamma_l']), 5.0e-2) and ok
            ok = _check_lower(label + ':max_r2_np2',
                              max(np2['r2_r'], np2['r2_l']), 0.4) and ok
            ok = _check_lower(label + ':branch_split_np2',
                              abs(np2['asymmetry']), 2.0e-2) and ok
            ok = _check_log_ratio_close(label + ':max_gamma_np1_vs_np2',
                                        max(np2['gamma_r'], np2['gamma_l']),
                                        max(np1['gamma_r'], np1['gamma_l']), 1.6) and ok

    if 'prolate_np4' in _RESULTS:
        for label in _CASES:
            np1 = _RESULTS[label + '_np1']
            np4 = _RESULTS[label + '_np4']
            ok = _check_lower(label + ':max_gamma_np4',
                              max(np4['gamma_r'], np4['gamma_l']), 5.0e-2) and ok
            ok = _check_lower(label + ':max_r2_np4',
                              max(np4['r2_r'], np4['r2_l']), 0.3) and ok
            ok = _check_lower(label + ':branch_split_np4',
                              abs(np4['asymmetry']), 2.0e-2) and ok
            ok = _check_log_ratio_close(label + ':max_gamma_np1_vs_np4',
                                        max(np4['gamma_r'], np4['gamma_l']),
                                        max(np1['gamma_r'], np1['gamma_l']), 1.8) and ok

    return ok
