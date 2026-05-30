import glob
import logging
import math
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from scripts.particles import pic_adaptive_deltaf_smoke as smoke

logger = logging.getLogger('athena' + __name__[7:])

_NLIM = 8
_FAST_INTERVAL = 0.2
_SLOW_INTERVAL = 10.0
_RESULTS = {}


def _particle_output_files(basename):
    pattern = os.path.join(
        smoke._athena_exe_dir(), 'pvtk', basename + '.prtcl_all.*.part.vtk')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No particle VTK files found for pattern: ' + pattern)
    return matches


def _read_deltaf_weights(path):
    with open(path, 'rb') as fp:
        contents = fp.read()
    _, offset, match = smoke._find_marker(
        contents, rb'\nPOINTS\s+([0-9]+)\s+float\n', 0, 'POINTS')
    npoint = int(match.group(1))
    offset += 4*3*npoint
    for name in ['gid', 'ptag', 'species']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' int\nLOOKUP_TABLE default\n')
        _, offset, _ = smoke._find_marker(contents, pattern, offset,
                                          'SCALARS ' + name)
        offset += 4*npoint
    for name in ['deltaf_f0', 'deltaf_weight']:
        pattern = (rb'\nSCALARS ' + name.encode('ascii') +
                   rb' float\nLOOKUP_TABLE default\n')
        _, offset, _ = smoke._find_marker(contents, pattern, offset,
                                          'SCALARS ' + name)
        values = smoke._read_big_endian_floats(contents, offset, npoint, name)
        offset += 4*npoint
    return values


def _parse_fits(output):
    matches = re.findall(
        r'PIC adaptive delta-f fit: time=([^ ]+) bucket=([^ ]+) '
        r'xi=([^ ]+) p0=([^\n]+)', output)
    return [
        {
            'time': float(time),
            'bucket': int(bucket),
            'xi': float(xi),
            'p0': float(p0),
        }
        for time, bucket, xi, p0 in matches
    ]


def _run_case(label, interval):
    basename = 'pic_adaptive_deltaf_cadence_variance_' + label
    smoke._remove_outputs(basename)
    output = smoke._run_athena(
        label,
        [
            'job/basename=' + basename,
            'time/nlim=' + str(_NLIM),
            'particles/pic_deltaf_adapt_interval=' + str(interval),
            'output2/dcycle=0',
        ])
    weights = [_read_deltaf_weights(path)
               for path in _particle_output_files(basename)]
    return {
        'fits': _parse_fits(output),
        'weights': weights,
        'variances': np.asarray([np.var(values) for values in weights]),
        'max_abs_weights': np.asarray([np.max(np.abs(values)) for values in weights]),
    }


def _expected_bucket(time, interval):
    ratio = time/interval
    roundoff = 64.0*np.finfo(float).eps*max(1.0, abs(ratio))
    return math.floor(ratio + roundoff)


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    _RESULTS.clear()
    _RESULTS['slow'] = _run_case('slow', _SLOW_INTERVAL)
    _RESULTS['fast'] = _run_case('fast', _FAST_INTERVAL)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    slow = _RESULTS['slow']
    fast = _RESULTS['fast']
    fast_times = np.asarray([row['time'] for row in fast['fits']])
    fast_buckets = np.asarray([row['bucket'] for row in fast['fits']])
    fast_xi = np.asarray([row['xi'] for row in fast['fits']])
    fast_p0 = np.asarray([row['p0'] for row in fast['fits']])
    expected_buckets = np.asarray(
        [_expected_bucket(time, _FAST_INTERVAL) for time in fast_times])
    final_variances = {
        'slow': float(slow['variances'][-1]),
        'fast': float(fast['variances'][-1]),
    }
    variance_separation = final_variances['fast'] - final_variances['slow']
    measured = {
        'slow_fit_count': len(slow['fits']),
        'fast_fit_times': fast_times.tolist(),
        'fast_fit_buckets': fast_buckets.tolist(),
        'fast_fit_xi': fast_xi.tolist(),
        'fast_fit_p0': fast_p0.tolist(),
        'expected_fast_buckets': expected_buckets.tolist(),
        'slow_variances': slow['variances'].tolist(),
        'fast_variances': fast['variances'].tolist(),
        'final_variances': final_variances,
        'variance_separation': variance_separation,
    }
    logger.info('adaptive delta-f cadence/variance metrics: %s', measured)
    _RESULTS['metrics'] = measured
    all_weights = slow['weights'] + fast['weights']
    return (
        len(slow['fits']) == 1
        and slow['fits'][0]['time'] == 0.0
        and slow['fits'][0]['bucket'] == 0
        and len(fast['fits']) == _NLIM
        and np.all(np.diff(fast_times) > 0.0)
        and np.all(np.diff(fast_buckets) > 0)
        and np.all(np.isfinite(fast_xi))
        and np.all(np.isfinite(fast_p0))
        and np.all(np.diff(fast_xi) < 0.0)
        and np.all(np.diff(fast_p0) < 0.0)
        and np.array_equal(fast_buckets, expected_buckets)
        and all(len(values) == 64 for values in all_weights)
        and all(np.all(np.isfinite(values)) for values in all_weights)
        and slow['variances'][0] == 0.0
        and fast['variances'][0] == 0.0
        and np.all(np.isfinite(slow['variances']))
        and np.all(np.isfinite(fast['variances']))
        and 0.20 <= final_variances['slow'] <= 0.30
        and 0.33 <= final_variances['fast'] <= 0.45
        and variance_separation >= 0.10
        and np.max(slow['max_abs_weights']) <= 1.30
        and np.max(fast['max_abs_weights']) <= 1.30
    )


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    run()
    if not analyze():
        raise SystemExit('pic_adaptive_deltaf_cadence_variance: FAIL')
    print('pic_adaptive_deltaf_cadence_variance: PASS')
