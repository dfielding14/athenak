import glob
import logging
import os
import re
import subprocess

import numpy as np
import scripts.utils.athena as athena

import sys
sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_amr_shock_lb_smoke.athinput'
_MPIEXEC = os.environ.get('MPIEXEC', 'mpiexec')
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


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


def _output_files_by_cycle(basename, file_id):
    pattern = os.path.join(_athena_exe_dir(), 'bin', basename + '.' + file_id + '.*.bin')
    files = sorted(glob.glob(pattern))
    out = {}
    for fname in files:
        cycle = os.path.basename(fname).split('.')[-2]
        out[cycle] = fname
    return out


def _load_smoke_metrics(basename):
    bfiles = _output_files_by_cycle(basename, 'mhd_bcc')
    jxfiles = _output_files_by_cycle(basename, 'prtcl_jx')
    jyfiles = _output_files_by_cycle(basename, 'prtcl_jy')
    jzfiles = _output_files_by_cycle(basename, 'prtcl_jz')

    shared_cycles = sorted(set(bfiles) & set(jxfiles) & set(jyfiles) & set(jzfiles))
    if len(shared_cycles) < 8:
        raise RuntimeError('Not enough matched outputs for AMR shock smoke metrics')

    b1std_series = []
    jrms_series = []
    tail_series = []

    for cyc in shared_cycles:
        bdat = bin_convert.read_binary_as_athdf(bfiles[cyc])
        jxdat = bin_convert.read_binary_as_athdf(jxfiles[cyc])
        jydat = bin_convert.read_binary_as_athdf(jyfiles[cyc])
        jzdat = bin_convert.read_binary_as_athdf(jzfiles[cyc])

        b1 = np.asarray(bdat['bcc1'], dtype=float)
        b1std_series.append(float(np.std(b1)))

        jx = np.asarray(jxdat['prtcl_jx'], dtype=float)
        jy = np.asarray(jydat['prtcl_jy'], dtype=float)
        jz = np.asarray(jzdat['prtcl_jz'], dtype=float)
        jmag = np.sqrt(jx * jx + jy * jy + jz * jz)
        jrms = float(np.sqrt(np.mean(jmag * jmag)))
        tail = float(np.quantile(jmag, 0.99) / max(jrms, 1.0e-30))
        jrms_series.append(jrms)
        tail_series.append(tail)

    b1std = np.asarray(b1std_series, dtype=float)
    jrms = np.asarray(jrms_series, dtype=float)
    tail = np.asarray(tail_series, dtype=float)

    return {
        'n_samples': float(len(shared_cycles)),
        'b1_std_final_amp': float(b1std[-1] / max(b1std[0], 1.0e-30)),
        'b1_std_peak_amp': float(np.max(b1std) / max(b1std[0], 1.0e-30)),
        'jrms_final': float(jrms[-1]),
        'jrms_peak': float(np.max(jrms)),
        'tail_final': float(tail[-1]),
        'tail_peak': float(np.max(tail)),
        'tail_delta': float(tail[-1] - tail[0]),
    }


def _parse_telemetry(output):
    out = {
        'n_created': float('nan'),
        'n_deleted': float('nan'),
        'n_communicated': float('nan'),
        'lb_efficiency': float('nan'),
    }

    amr_match = re.search(r'(\d+) MeshBlocks created,\s+(\d+) deleted by AMR', output)
    if amr_match is not None:
        out['n_created'] = float(amr_match.group(1))
        out['n_deleted'] = float(amr_match.group(2))

    lb_match = re.search(
        r'(\d+) communicated for load balancing,\s*load balancing efficiency =\s*'
        r'([0-9eE+.-]+)',
        output)
    if lb_match is not None:
        out['n_communicated'] = float(lb_match.group(1))
        out['lb_efficiency'] = float(lb_match.group(2))

    return out


def _run_case(label, basename, nproc, extra_args=None):
    args = ['job/basename=' + basename]
    if extra_args is not None:
        args.extend(extra_args)
    command = ['./athena', '-i', _athena_input_path()] + args
    if nproc > 1:
        command = [_MPIEXEC, '-n', str(nproc)] + command

    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Command failed for ' + label + '\n' + output)

    output = (proc.stdout or '') + (proc.stderr or '')
    metrics = _load_smoke_metrics(basename)
    metrics.update(_parse_telemetry(output))
    _RESULTS[label] = metrics


def _check_lower(label, measured, lower):
    logger.info('%s measured=% .8e lower=% .8e margin=% .8e',
                label, measured, lower, measured - lower)
    return measured >= lower


def _check_upper(label, measured, upper):
    logger.info('%s measured=% .8e upper=% .8e margin=% .8e',
                label, measured, upper, measured - upper)
    return measured <= upper


def _check_log_ratio(label, measured, expected, upper):
    ratio = abs(np.log(max(measured, 1.0e-30) / max(expected, 1.0e-30)))
    logger.info('%s measured=% .8e expected=% .8e |log-ratio|=% .8e upper=% .8e',
                label, measured, expected, ratio, upper)
    return ratio <= upper


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    common_args = [
        'time/nlim=10',
        'time/tlim=0.08',
        'mesh_refinement/ddens_max=0.05',
    ]
    alt_decomp_args = common_args + [
        'meshblock/nx1=4',
        'meshblock/nx2=8',
        'meshblock/nx3=4',
    ]

    _remove_outputs('pic_amr_shock_lb_serial')
    _run_case('serial', 'pic_amr_shock_lb_serial', 1, common_args)

    _remove_outputs('pic_amr_shock_lb_serial_alt')
    _run_case('serial_alt', 'pic_amr_shock_lb_serial_alt', 1, alt_decomp_args)

    if _athena_mpi_enabled():
        _remove_outputs('pic_amr_shock_lb_mpi2')
        _run_case('mpi2', 'pic_amr_shock_lb_mpi2', 2, common_args)

        _remove_outputs('pic_amr_shock_lb_mpi4')
        _run_case('mpi4', 'pic_amr_shock_lb_mpi4', 4, common_args)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    ok = True

    serial = _RESULTS['serial']
    ok = _check_lower('serial:n_samples', serial['n_samples'], 8.0) and ok
    ok = _check_lower('serial:b1_std_final_amp', serial['b1_std_final_amp'], 1.002) and ok
    ok = _check_lower('serial:b1_std_peak_amp', serial['b1_std_peak_amp'], 1.002) and ok
    ok = _check_lower('serial:jrms_peak', serial['jrms_peak'], 1.0e-6) and ok
    ok = _check_lower('serial:tail_peak', serial['tail_peak'], 1.5) and ok
    ok = _check_lower('serial:tail_delta', serial['tail_delta'], 0.0) and ok
    ok = _check_lower('serial:n_created', serial['n_created'], 1.0) and ok

    serial_alt = _RESULTS['serial_alt']
    ok = _check_log_ratio('serial_vs_serial_alt:b1_std_peak_amp',
                          serial_alt['b1_std_peak_amp'], serial['b1_std_peak_amp'],
                          5.0e-2) and ok
    ok = _check_log_ratio('serial_vs_serial_alt:tail_peak',
                          serial_alt['tail_peak'], serial['tail_peak'], 3.0e-1) and ok

    if 'mpi2' in _RESULTS:
        mpi2 = _RESULTS['mpi2']
        ok = _check_lower('mpi2:n_samples', mpi2['n_samples'], 8.0) and ok
        ok = _check_lower('mpi2:b1_std_final_amp', mpi2['b1_std_final_amp'], 1.002) and ok
        ok = _check_lower('mpi2:b1_std_peak_amp', mpi2['b1_std_peak_amp'], 1.002) and ok
        ok = _check_lower('mpi2:jrms_peak', mpi2['jrms_peak'], 1.0e-6) and ok
        ok = _check_lower('mpi2:tail_peak', mpi2['tail_peak'], 1.5) and ok
        ok = _check_lower('mpi2:tail_delta', mpi2['tail_delta'], 0.0) and ok
        ok = _check_lower('mpi2:n_created', mpi2['n_created'], 1.0) and ok
        ok = _check_log_ratio('serial_vs_mpi2:b1_std_peak_amp',
                              mpi2['b1_std_peak_amp'], serial['b1_std_peak_amp'],
                              1.0e-3) and ok
        ok = _check_log_ratio('serial_vs_mpi2:tail_peak',
                              mpi2['tail_peak'], serial['tail_peak'], 0.2) and ok

    if 'mpi4' in _RESULTS:
        mpi4 = _RESULTS['mpi4']
        ok = _check_lower('mpi4:n_samples', mpi4['n_samples'], 8.0) and ok
        ok = _check_lower('mpi4:b1_std_final_amp', mpi4['b1_std_final_amp'], 1.002) and ok
        ok = _check_lower('mpi4:b1_std_peak_amp', mpi4['b1_std_peak_amp'], 1.002) and ok
        ok = _check_lower('mpi4:jrms_peak', mpi4['jrms_peak'], 1.0e-6) and ok
        ok = _check_lower('mpi4:tail_peak', mpi4['tail_peak'], 1.5) and ok
        ok = _check_lower('mpi4:tail_delta', mpi4['tail_delta'], 0.0) and ok
        ok = _check_lower('mpi4:n_created', mpi4['n_created'], 1.0) and ok
        ok = _check_lower('mpi4:n_communicated', mpi4['n_communicated'], 1.0) and ok
        ok = _check_lower('mpi4:lb_efficiency_lower', mpi4['lb_efficiency'], 0.1) and ok
        ok = _check_upper('mpi4:lb_efficiency_upper', mpi4['lb_efficiency'], 1.0) and ok

        ok = _check_log_ratio('serial_vs_mpi4:b1_std_peak_amp',
                              mpi4['b1_std_peak_amp'], serial['b1_std_peak_amp'],
                              1.0e-3) and ok
        ok = _check_log_ratio('serial_vs_mpi4:tail_peak',
                              mpi4['tail_peak'], serial['tail_peak'], 0.2) and ok

    if 'mpi2' in _RESULTS and 'mpi4' in _RESULTS:
        ok = _check_log_ratio('mpi2_vs_mpi4:b1_std_peak_amp',
                              _RESULTS['mpi4']['b1_std_peak_amp'],
                              _RESULTS['mpi2']['b1_std_peak_amp'], 1.0e-3) and ok
        ok = _check_log_ratio('mpi2_vs_mpi4:tail_peak',
                              _RESULTS['mpi4']['tail_peak'],
                              _RESULTS['mpi2']['tail_peak'], 0.2) and ok

    return ok
