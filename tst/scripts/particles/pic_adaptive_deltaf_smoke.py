import glob
import logging
import math
import os
import re
import subprocess

import numpy as np
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUT_DECK = 'tests/pic_adaptive_deltaf_smoke.athinput'
_KAPPA = 2.0
_STATES = np.asarray([[2.0, 1.0, 0.0], [1.0, 0.0, 2.0]])
_RESULTS = {}


def _athena_exe_dir():
    return os.path.join(os.getcwd(), 'build', 'src')


def _athena_input_path():
    return '../../' + athena.athena_rel_path + 'inputs/' + _INPUT_DECK


def _remove_outputs(basename):
    exe_dir = _athena_exe_dir()
    for pattern in [
            os.path.join(exe_dir, 'pvtk', basename + '.*.part.vtk'),
            os.path.join(exe_dir, 'rst', basename + '.*.rst*')]:
        for fname in glob.glob(pattern):
            os.remove(fname)


def _run_athena(label, arguments, restart_file=None, expect_fail=None):
    command = ['./athena']
    if restart_file is None:
        command += ['-i', _athena_input_path()]
    else:
        command += ['-r', restart_file]
    command += list(arguments)
    logger.info('Executing %s: %s', label, ' '.join(command))
    proc = subprocess.run(command, cwd=_athena_exe_dir(),
                          capture_output=True, text=True)
    output = (proc.stdout or '') + (proc.stderr or '')
    if expect_fail is None:
        if proc.returncode != 0:
            raise RuntimeError('Command failed for ' + label + '\n' + output)
        return output
    if proc.returncode == 0:
        raise RuntimeError('Expected failure for ' + label + ', but command passed')
    if expect_fail not in output:
        raise RuntimeError('Unexpected failure reason for ' + label + '\n' + output)
    return output


def _latest_file(pattern, label):
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No ' + label + ' files found for pattern: ' + pattern)
    return matches[-1]


def _latest_pvtk_file(basename):
    return _latest_file(os.path.join(_athena_exe_dir(), 'pvtk',
                                     basename + '.prtcl_all.*.part.vtk'),
                        'particle VTK')


def _latest_restart_file(basename):
    return _latest_file(os.path.join(_athena_exe_dir(), 'rst',
                                     basename + '.*.rst'), 'restart')


def _read_big_endian_floats(contents, offset, count, label):
    payload = contents[offset:offset + 4*count]
    if len(payload) != 4*count:
        raise RuntimeError('Truncated ' + label + ' payload in particle VTK output')
    return np.frombuffer(payload, dtype='>f4').astype(np.float64)


def _find_marker(contents, pattern, offset, label):
    match = re.search(pattern, contents[offset:])
    if match is None:
        raise RuntimeError('Could not find ' + label + ' marker in particle VTK output')
    return offset + match.start(), offset + match.end(), match


def _read_pvtk_snapshot(basename):
    path = _latest_pvtk_file(basename)
    with open(path, 'rb') as fp:
        contents = fp.read()
    _, offset, match = _find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS')
    npoint = int(match.group(1))
    offset += 4*3*npoint
    for name in ['gid', 'ptag', 'species']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' int\nLOOKUP_TABLE default\n')
        _, offset, _ = _find_marker(contents, pattern, offset, 'SCALARS ' + name)
        offset += 4*npoint
    scalars = {}
    for name in ['deltaf_f0', 'deltaf_weight']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' float\nLOOKUP_TABLE default\n')
        _, offset, _ = _find_marker(contents, pattern, offset, 'SCALARS ' + name)
        scalars[name] = _read_big_endian_floats(contents, offset, npoint, name)
        offset += 4*npoint
    _, offset, _ = _find_marker(contents, rb'\nVECTORS vel float\n',
                                offset, 'VECTORS vel')
    velocity = _read_big_endian_floats(contents, offset, 3*npoint, 'vel')
    return {
        'path': path,
        'npoint': npoint,
        'velocity': velocity.reshape(npoint, 3),
        'deltaf_f0': scalars['deltaf_f0'],
        'deltaf_weight': scalars['deltaf_weight'],
    }


def _expected_fit():
    parallel = np.sum(np.abs(_STATES[:, 0]))
    perpendicular = np.sum(np.linalg.norm(_STATES[:, 1:], axis=1))
    xi = (2.0/np.pi)*perpendicular/parallel
    shape = np.sqrt(xi**4*_STATES[:, 0]**2
                    + xi**2*np.sum(_STATES[:, 1:]**2, axis=1))
    prefactor = (np.sqrt(np.pi*_KAPPA)*(_KAPPA - 1.0)
                 * math.gamma(_KAPPA - 0.5)
                 / (2.0*math.gamma(_KAPPA + 1.0)))
    return xi, prefactor*np.mean(shape)


def _parse_fit(output):
    matches = re.findall(
        r'PIC adaptive delta-f fit: time=([^ ]+) bucket=([^ ]+) '
        r'xi=([^ ]+) p0=([^\n]+)', output)
    if len(matches) != 1:
        raise RuntimeError('Expected exactly one initial adaptive fit:\n' + output)
    time, bucket, xi, p0 = matches[0]
    return {
        'time': float(time),
        'bucket': int(bucket),
        'xi': float(xi),
        'p0': float(p0),
    }


def _parse_adaptive_timer_calls(output):
    match = re.search(
        r'q017\.telemetry\.timer\.particle\.adaptive_deltaf\.calls_rank_max=([^\n]+)',
        output)
    if match is None:
        raise RuntimeError('Missing adaptive delta-f Q-017 timer telemetry:\n' + output)
    return float(match.group(1))


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    full = 'pic_adaptive_deltaf_full'
    segment = 'pic_adaptive_deltaf_segment'
    restarted = 'pic_adaptive_deltaf_restart'
    for basename in [full, segment, restarted]:
        _remove_outputs(basename)

    full_output = _run_athena('full', ['job/basename=' + full, 'time/nlim=2'])
    _RESULTS['fit'] = _parse_fit(full_output)
    _RESULTS['adaptive_timer_calls'] = _parse_adaptive_timer_calls(full_output)
    _RESULTS['full'] = _read_pvtk_snapshot(full)

    _run_athena('segment', ['job/basename=' + segment, 'time/nlim=1'])
    restart_file = os.path.relpath(_latest_restart_file(segment),
                                   _athena_exe_dir())
    restart_output = _run_athena(
        'restart',
        ['job/basename=' + restarted, 'time/nlim=2',
         'output1/file_number=0', 'output2/dcycle=0'],
        restart_file=restart_file)
    _RESULTS['restart'] = _read_pvtk_snapshot(restarted)
    _RESULTS['restart_refit'] = 'PIC adaptive delta-f fit:' in restart_output

    _run_athena(
        'guard_restart_adapt_interval',
        ['particles/pic_deltaf_adapt_interval=11.0',
         'time/nlim=2', 'output1/dcycle=0', 'output2/dcycle=0'],
        restart_file=restart_file,
        expect_fail='Particle restart physical-model metadata mismatch')
    _run_athena(
        'guard_adapt_requires_extension_mode',
        ['particles/pic_physical_mode=paper_test_particle', 'time/nlim=0',
         'output1/dcycle=0', 'output2/dcycle=0'],
        expect_fail='requires <particles>/pic_physical_mode=extended_mhd_pic')
    _run_athena(
        'guard_adapt_requires_bounded_config',
        ['particles/pic_deltaf_mode=off', 'time/nlim=0',
         'output1/dcycle=0', 'output2/dcycle=0'],
        expect_fail='requires physical kappa_aniso delta-f')


def analyze():
    logger.debug('Analyzing test ' + __name__)
    fit = _RESULTS['fit']
    expected_xi, expected_p0 = _expected_fit()
    fit_errors = {
        'xi': abs(fit['xi'] - expected_xi),
        'p0': abs(fit['p0'] - expected_p0),
    }
    full = _RESULTS['full']
    restarted = _RESULTS['restart']
    restart_errors = {
        'velocity': float(np.max(np.abs(full['velocity'] - restarted['velocity']))),
        'deltaf_f0': float(np.max(np.abs(full['deltaf_f0']
                                          - restarted['deltaf_f0']))),
        'deltaf_weight': float(np.max(np.abs(full['deltaf_weight']
                                              - restarted['deltaf_weight']))),
    }
    measured = {
        'fit': fit,
        'expected_fit': {'xi': expected_xi, 'p0': expected_p0},
        'fit_errors': fit_errors,
        'restart_errors': restart_errors,
        'restart_refit': _RESULTS['restart_refit'],
        'adaptive_timer_calls': _RESULTS['adaptive_timer_calls'],
    }
    logger.info('adaptive delta-f metrics: %s', measured)
    _RESULTS['metrics'] = measured
    return (fit['time'] == 0.0 and fit['bucket'] == 0
            and max(fit_errors.values()) <= 1.0e-12
            and max(restart_errors.values()) <= 1.0e-6
            and not _RESULTS['restart_refit']
            and _RESULTS['adaptive_timer_calls'] > 0.0
            and full['npoint'] == restarted['npoint'] == 64
            and np.all(np.isfinite(full['deltaf_weight']))
            and np.max(np.abs(full['deltaf_weight'])) > 1.0e-6)
