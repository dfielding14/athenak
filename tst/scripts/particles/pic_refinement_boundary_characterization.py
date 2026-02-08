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

_CASES = {
    'smr': 'tests/pic_refinement_boundary_smr_proxy.athinput',
    'amr_proxy': 'tests/pic_refinement_boundary_amr_proxy.athinput',
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


def _volume_weights(dataset):
    dx1 = np.diff(dataset['x1f'])
    dx2 = np.diff(dataset['x2f'])
    dx3 = np.diff(dataset['x3f'])
    return dx3[:, None, None] * dx2[None, :, None] * dx1[None, None, :]


def _integrate(dataset, quantity):
    dvol = _volume_weights(dataset)
    return float(np.sum(dataset[quantity] * dvol))


def _rho_metrics(dataset):
    rho = np.asarray(dataset['prtcl_rho'], dtype=float)
    dvol = _volume_weights(dataset)

    total_w = float(np.sum(dvol))
    mean_rho = float(np.sum(rho * dvol) / max(total_w, 1.0e-30))
    var_rho = float(np.sum(((rho - mean_rho) ** 2) * dvol) / max(total_w, 1.0e-30))
    cv = float(np.sqrt(max(var_rho, 0.0)) / max(abs(mean_rho), 1.0e-12))

    wx = np.sum(dvol, axis=(0, 1))
    profile = np.sum(rho * dvol, axis=(0, 1)) / np.maximum(wx, 1.0e-30)
    max_adj = float(np.max(np.abs(np.diff(profile))))
    jump = float(max_adj / max(abs(np.mean(profile)), 1.0e-12))

    return cv, jump


def _load_metrics(basename):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.prtcl_rho.*.bin')
    files = sorted(glob.glob(pattern))
    if len(files) < 6:
        raise RuntimeError('Not enough outputs for refinement-boundary metrics')

    data0 = bin_convert.read_binary_as_athdf(files[0])
    dataf = bin_convert.read_binary_as_athdf(files[-1])

    q0 = _integrate(data0, 'prtcl_rho')
    qf = _integrate(dataf, 'prtcl_rho')
    q_drift = abs(qf - q0) / max(abs(q0), 1.0)

    cv, jump = _rho_metrics(dataf)

    return {
        'q_drift': float(q_drift),
        'cv': float(cv),
        'jump': float(jump),
    }


def _run_case(label, nproc, input_rel, basename, extra_args=None):
    args = ['job/basename=' + basename]
    if extra_args is not None:
        args.extend(extra_args)
    command = ['./athena', '-i', _athena_input_path(input_rel)] + args
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _check_close(label, measured, expected, abs_tol, rel_tol):
    abs_err = abs(measured - expected)
    rel_err = abs_err / max(abs(expected), 1.0)
    logger.info('%s measured=% .8e expected=% .8e abs_err=% .8e rel_err=% .8e',
                label, measured, expected, abs_err, rel_err)
    return abs_err <= abs_tol or rel_err <= rel_tol


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    alt_decomp_args = [
        'meshblock/nx1=8',
        'meshblock/nx2=4',
        'meshblock/nx3=4',
    ]

    for case, input_rel in _CASES.items():
        base = 'pic_ref_boundary_' + case + '_np1'
        _remove_outputs(base)
        _run_case(case + '_serial', 1, input_rel, base)
        _RESULTS[case + '_np1'] = _load_metrics(base)

        base_alt = 'pic_ref_boundary_' + case + '_np1_alt'
        _remove_outputs(base_alt)
        _run_case(case + '_serial_alt', 1, input_rel, base_alt, alt_decomp_args)
        _RESULTS[case + '_np1_alt'] = _load_metrics(base_alt)

    if _athena_mpi_enabled():
        for case, input_rel in _CASES.items():
            base = 'pic_ref_boundary_' + case + '_np2'
            _remove_outputs(base)
            _run_case(case + '_mpi2', 2, input_rel, base)
            _RESULTS[case + '_np2'] = _load_metrics(base)

            base = 'pic_ref_boundary_' + case + '_np4'
            _remove_outputs(base)
            _run_case(case + '_mpi4', 4, input_rel, base)
            _RESULTS[case + '_np4'] = _load_metrics(base)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    for case in _CASES:
        metrics = _RESULTS[case + '_np1']
        ok = _check_upper(case + ':q_drift_np1', metrics['q_drift'], 3.0e4) and ok
        ok = _check_upper(case + ':cv_np1', metrics['cv'], 1.2) and ok
        ok = _check_upper(case + ':jump_np1', metrics['jump'], 6.0) and ok

    smr = _RESULTS['smr_np1']
    amr = _RESULTS['amr_proxy_np1']
    ok = _check_close('smr_vs_amr:q_drift_np1', amr['q_drift'], smr['q_drift'],
                      1.0e4, 0.0) and ok
    ok = _check_close('smr_vs_amr:cv_np1', amr['cv'], smr['cv'], 5.0e-1, 0.0) and ok
    ok = _check_close('smr_vs_amr:jump_np1', amr['jump'], smr['jump'], 4.0, 0.0) and ok

    if 'smr_np2' in _RESULTS:
        for case in _CASES:
            np1 = _RESULTS[case + '_np1']
            np2 = _RESULTS[case + '_np2']
            ok = _check_upper(case + ':q_drift_np2', np2['q_drift'], 3.0e4) and ok
            ok = _check_upper(case + ':cv_np2', np2['cv'], 1.2) and ok
            ok = _check_upper(case + ':jump_np2', np2['jump'], 6.0) and ok
            ok = _check_close(case + ':q_drift_np1_vs_np2', np2['q_drift'],
                              np1['q_drift'], 1.0e3, 3.0e-2) and ok
            ok = _check_close(case + ':cv_np1_vs_np2', np2['cv'], np1['cv'],
                              1.0e-1, 0.0) and ok
            ok = _check_close(case + ':jump_np1_vs_np2', np2['jump'], np1['jump'],
                              2.0e-1, 0.0) and ok

    for case in _CASES:
        alt = _RESULTS[case + '_np1_alt']
        ok = _check_upper(case + ':q_drift_np1_alt', alt['q_drift'], 3.0e4) and ok
        ok = _check_upper(case + ':cv_np1_alt', alt['cv'], 1.2) and ok
        ok = _check_upper(case + ':jump_np1_alt', alt['jump'], 6.0) and ok

    if 'smr_np4' in _RESULTS:
        for case in _CASES:
            np1 = _RESULTS[case + '_np1']
            np4 = _RESULTS[case + '_np4']
            ok = _check_upper(case + ':q_drift_np4', np4['q_drift'], 3.0e4) and ok
            ok = _check_upper(case + ':cv_np4', np4['cv'], 1.2) and ok
            ok = _check_upper(case + ':jump_np4', np4['jump'], 6.0) and ok
            ok = _check_close(case + ':q_drift_np1_vs_np4', np4['q_drift'],
                              np1['q_drift'], 1.5e3, 5.0e-2) and ok
            ok = _check_close(case + ':cv_np1_vs_np4', np4['cv'], np1['cv'],
                              1.5e-1, 0.0) and ok
            ok = _check_close(case + ':jump_np1_vs_np4', np4['jump'], np1['jump'],
                              3.0e-1, 0.0) and ok

    if 'smr_np2' in _RESULTS and 'smr_np4' in _RESULTS:
        for case in _CASES:
            np2 = _RESULTS[case + '_np2']
            np4 = _RESULTS[case + '_np4']
            ok = _check_close(case + ':q_drift_np2_vs_np4', np4['q_drift'],
                              np2['q_drift'], 1.5e3, 5.0e-2) and ok
            ok = _check_close(case + ':cv_np2_vs_np4', np4['cv'], np2['cv'],
                              1.5e-1, 0.0) and ok
            ok = _check_close(case + ':jump_np2_vs_np4', np4['jump'], np2['jump'],
                              3.0e-1, 0.0) and ok

    return ok
