# Regression test for scalar_mixing_test blob advection-diffusion verification.
#
# This test is only active when AthenaK is configured with
#   --cmake=-DPROBLEM=scalar_mixing_test
#
# It exercises the calibrated 2D and 3D diffusive reference problems and also runs
# pure-advection smoke checks with scalar_diffusivity=0 at a conservative CFL.

import glob
import logging
import os
import sys

import scripts.utils.athena as athena

sys.path.insert(0, '../vis/python')
import scalar_mixing_test_verify as verify  # noqa

logger = logging.getLogger('athena' + __name__[7:])

_enabled = False

_DIFFUSIVE_CASES = [
    {
        'name': 'diffusive_2d',
        'input': 'tests/scalar_mixing_blob_diffusion_2d.athinput',
        'basename': 'scalar_blob_reg_diff_2d',
    },
    {
        'name': 'diffusive_3d',
        'input': 'tests/scalar_mixing_blob_diffusion_3d.athinput',
        'basename': 'scalar_blob_reg_diff_3d',
    },
]

_PURE_ADVECTION_CASES = [
    {
        'name': 'advection_2d',
        'input': 'tests/scalar_mixing_blob_diffusion_2d.athinput',
        'basename': 'scalar_blob_reg_adv_2d',
        'cfl': 0.2,
        'mass_rel_max': 1.0e-8,
        'centroid_max': 5.0e-4,
    },
    {
        'name': 'advection_3d',
        'input': 'tests/scalar_mixing_blob_diffusion_3d.athinput',
        'basename': 'scalar_blob_reg_adv_3d',
        'cfl': 0.2,
        'mass_rel_max': 1.0e-8,
        'centroid_max': 2.0e-3,
    },
]


def _compiled_problem():
    cfg_path = os.path.join('build', 'config.hpp')
    if not os.path.exists(cfg_path):
        return None
    with open(cfg_path, 'r', encoding='utf-8') as fobj:
        for line in fobj:
            if line.startswith('#define PROBLEM_GENERATOR'):
                return line.split('"')[1]
    return None


def _clear_old_outputs(basename):
    pattern = os.path.join('build', 'src', 'bin', basename + '.hydro_w.*.bin')
    for path in glob.glob(pattern):
        os.remove(path)


def _output_path(basename):
    pattern = os.path.join('build', 'src', 'bin', basename + '.hydro_w.*.bin')
    matches = sorted(glob.glob(pattern))
    if not matches:
        raise RuntimeError('No hydro_w output found for ' + basename)
    return matches[-1]


def _input_path(case):
    return os.path.join('..', 'inputs', case['input'])


def run(**kwargs):
    del kwargs
    global _enabled
    problem = _compiled_problem()
    _enabled = (problem == 'scalar_mixing_test')
    if not _enabled:
        logger.info('Skipping %s because build/config.hpp reports PROBLEM=%s',
                    __name__, problem)
        return

    for case in _DIFFUSIVE_CASES:
        _clear_old_outputs(case['basename'])
        athena.run(case['input'],
                   [f'job/basename={case["basename"]}',
                    'output1/dt=1.0'])

    for case in _PURE_ADVECTION_CASES:
        _clear_old_outputs(case['basename'])
        athena.run(case['input'],
                   [f'job/basename={case["basename"]}',
                    'hydro/scalar_diffusivity=0.0',
                    f'time/cfl_number={case["cfl"]}',
                    'output1/dt=1.0'])


def analyze():
    if not _enabled:
        return True

    analyze_status = True

    for case in _DIFFUSIVE_CASES:
        metrics = verify.analyze_dump(_input_path(case), _output_path(case['basename']))
        logger.info('%s: l1=%.3e l2=%.3e linf=%.3e mass=%.3e centroid=%.3e width=%.3e',
                    case['name'], metrics['field_errors']['l1'],
                    metrics['field_errors']['l2'], metrics['field_errors']['linf'],
                    metrics['mass']['rel_error'], metrics['blob']['centroid_error'],
                    metrics['blob']['width_rel_error'])
        if metrics['passes_thresholds'] is not True:
            logger.warning('%s failed calibrated thresholds', case['name'])
            analyze_status = False

    for case in _PURE_ADVECTION_CASES:
        metrics = verify.analyze_dump(_input_path(case), _output_path(case['basename']))
        logger.info('%s: mass=%.3e centroid=%.3e',
                    case['name'], metrics['mass']['rel_error'],
                    metrics['blob']['centroid_error'])
        if metrics['mass']['rel_error'] > case['mass_rel_max']:
            logger.warning('%s mass drift too large: %.3e',
                           case['name'], metrics['mass']['rel_error'])
            analyze_status = False
        if metrics['blob']['centroid_error'] > case['centroid_max']:
            logger.warning('%s centroid drift too large: %.3e',
                           case['name'], metrics['blob']['centroid_error'])
            analyze_status = False
        for comp in ('velx', 'vely', 'velz'):
            if metrics['velocity'][comp]['max_dev'] > 1.0e-12:
                logger.warning('%s velocity not uniform in %s: %.3e',
                               case['name'], comp, metrics['velocity'][comp]['max_dev'])
                analyze_status = False

    return analyze_status
