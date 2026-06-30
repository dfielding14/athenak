"""Regression tests for one-time Fourier initial perturbations."""

import glob
import logging
import os
import shlex
import shutil
import subprocess
import sys

import numpy as np
import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import bin_convert_new as bin_convert  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_3D = 'tests/initial_perturbations.athinput'
_INPUT_2D = 'tests/initial_perturbations_2d.athinput'
_MPIEXEC = shlex.split(os.environ.get('MPIEXEC', 'mpiexec'))
_RESULTS = {}
_SKIPPED = False


def _exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _input_path(input_deck):
    return '../../' + athena.athena_rel_path + 'inputs/' + input_deck


def _configured_problem():
    proc = subprocess.run(
        ['./athena', '-c'], cwd=_exe_dir(), capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError('Unable to query Athena configuration')
    output = (proc.stdout or '') + (proc.stderr or '')
    for line in output.splitlines():
        if 'Problem generator:' in line:
            return line.split(':', 1)[1].strip()
    return ''


def _mpi_enabled():
    proc = subprocess.run(
        ['./athena', '-c'], cwd=_exe_dir(), capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    return 'MPI parallelism:            ON' in output


def _mpi_launcher_available():
    return bool(_MPIEXEC and shutil.which(_MPIEXEC[0]))


def _remove_outputs(basename):
    for pattern in (
            os.path.join(_exe_dir(), 'bin', basename + '.*.bin'),
            os.path.join(_exe_dir(), basename + '.*.hst'),
            os.path.join(_exe_dir(), basename + '-errs.dat')):
        for path in glob.glob(pattern):
            os.remove(path)


def _run(label, input_deck, basename, nproc, extra_args=None):
    args = ['job/basename=' + basename]
    if extra_args:
        args += list(extra_args)
    command = ['./athena', '-i', _input_path(input_deck)] + args
    if nproc > 1:
        command = _MPIEXEC + ['-n', str(nproc)] + command
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(
        command, cwd=_exe_dir(), capture_output=True, text=True)
    if proc.returncode != 0:
        output = (proc.stdout or '') + (proc.stderr or '')
        raise RuntimeError('Initial perturbation run failed:\n' + output)


def _latest(basename, file_id):
    pattern = os.path.join(
        _exe_dir(), 'bin', basename + '.' + file_id + '.*.bin')
    paths = sorted(glob.glob(pattern))
    if not paths:
        raise RuntimeError('No output found for ' + pattern)
    return paths[-1]


def _measure(basename, nlow, nhigh):
    fields = bin_convert.read_binary_as_athdf(
        _latest(basename, 'mhd_w_bcc'))
    divb_data = bin_convert.read_binary_as_athdf(
        _latest(basename, 'mhd_divb'))

    density = np.asarray(fields['dens'], dtype=float)
    density_delta = density/np.mean(density) - 1.0
    velocity = np.stack([
        np.asarray(fields['velx'], dtype=float),
        np.asarray(fields['vely'], dtype=float),
        np.asarray(fields['velz'], dtype=float),
    ])
    magnetic = np.stack([
        np.asarray(fields['bcc1'], dtype=float),
        np.asarray(fields['bcc2'], dtype=float),
        np.asarray(fields['bcc3'], dtype=float),
    ])
    magnetic_fluctuation = magnetic - np.mean(
        magnetic, axis=tuple(range(1, magnetic.ndim)), keepdims=True)

    power = np.abs(np.fft.fftn(density_delta))**2
    kz = np.fft.fftfreq(density_delta.shape[0])*density_delta.shape[0]
    ky = np.fft.fftfreq(density_delta.shape[1])*density_delta.shape[1]
    kx = np.fft.fftfreq(density_delta.shape[2])*density_delta.shape[2]
    ksq = (kz[:, None, None]**2 + ky[None, :, None]**2
            + kx[None, None, :]**2)
    support = (ksq >= nlow*nlow) & (ksq <= nhigh*nhigh)
    total_power = float(np.sum(power))

    return {
        'density_mean': float(np.mean(density_delta)),
        'density_rms': float(np.sqrt(np.mean(density_delta**2))),
        'zero_fraction': float(power[0, 0, 0]/total_power),
        'leakage_fraction': float(np.sum(power[~support])/total_power),
        'velocity_rms': float(np.sqrt(np.mean(np.sum(velocity**2, axis=0)))),
        'magnetic_rms': float(np.sqrt(
            np.mean(np.sum(magnetic_fluctuation**2, axis=0)))),
        'magnetic_mean': np.mean(
            magnetic, axis=tuple(range(1, magnetic.ndim))),
        'max_divb': float(np.max(np.abs(divb_data['divb']))),
        'density': density_delta,
        'velocity': velocity,
        'magnetic_fluctuation': magnetic_fluctuation,
    }


def run(**kwargs):
    global _SKIPPED
    logger.debug('Running test ' + __name__)
    if _configured_problem() != 'built_in_pgens':
        logger.info('Skipping: this regression requires the built-in pgen build')
        _SKIPPED = True
        return

    cases = [
        ('3d_np1', _INPUT_3D, 'InitialPerturbations_np1', 1, 1, 4, None),
        ('2d_np1', _INPUT_2D, 'InitialPerturbations2D_np1', 1, 1, 8, None),
    ]
    if _mpi_enabled() and _mpi_launcher_available():
        cases += [
            ('3d_np2', _INPUT_3D, 'InitialPerturbations_np2',
             2, 1, 4, None),
            ('3d_np4_mb8', _INPUT_3D, 'InitialPerturbations_np4_mb8',
             4, 1, 4,
             ['meshblock/nx1=8', 'meshblock/nx2=8',
              'meshblock/nx3=8']),
        ]
    elif _mpi_enabled():
        logger.info('Skipping MPI cases: launcher %s not found',
                    _MPIEXEC[0] if _MPIEXEC else '<empty>')

    for name, deck, basename, nproc, nlow, nhigh, extra in cases:
        _remove_outputs(basename)
        _run(name, deck, basename, nproc, extra)
        _RESULTS[name] = _measure(basename, nlow, nhigh)


def _check_metrics(name, result):
    ok = True
    checks = [
        ('density mean', abs(result['density_mean']), 1.0e-12),
        ('density rms error', abs(result['density_rms'] - 1.0e-2),
         5.0e-6),
        ('zero-mode fraction', result['zero_fraction'], 1.0e-20),
        # Binary outputs are written in single precision, so roundoff in the
        # stored density produces O(1e-11) off-band FFT power.
        ('spectral leakage', result['leakage_fraction'], 1.0e-10),
        ('velocity rms error', abs(result['velocity_rms'] - 1.0e-3),
         5.0e-6),
        ('magnetic rms error', abs(result['magnetic_rms'] - 1.0e-3),
         5.0e-6),
        ('max divB', result['max_divb'], 1.0e-10),
    ]
    for metric, value, tolerance in checks:
        logger.info('%s %s=% .8e tolerance=% .8e',
                    name, metric, value, tolerance)
        ok = (value <= tolerance) and ok
    return ok


def analyze():
    if _SKIPPED:
        return True

    ok = True
    for name, result in _RESULTS.items():
        ok = _check_metrics(name, result) and ok

    baseline = _RESULTS['3d_np1']
    for name in ('3d_np2', '3d_np4_mb8'):
        if name not in _RESULTS:
            continue
        candidate = _RESULTS[name]
        for field in ('density', 'velocity', 'magnetic_fluctuation'):
            error = float(np.max(np.abs(candidate[field] - baseline[field])))
            logger.info('%s %s decomposition error=% .8e', name, field, error)
            ok = (error < 2.0e-12) and ok
        mean_error = float(np.max(np.abs(
            candidate['magnetic_mean'] - baseline['magnetic_mean'])))
        logger.info('%s magnetic-mean decomposition error=% .8e',
                    name, mean_error)
        ok = (mean_error < 2.0e-12) and ok

    return ok
