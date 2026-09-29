# Regression coverage for the CGL Landau-fluid STS closure.
#
# The built-in problem generator evaluates its quantitative reference comparison at
# finalization and terminates nonzero on failure.

import logging

import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])

_INPUTS = [
    'unit_tests/cgl_lf_sts_parallel.athinput',
    'unit_tests/cgl_lf_sts_perp.athinput',
    'unit_tests/cgl_lf_sts_grad_b.athinput',
    'unit_tests/cgl_lf_flux_limiter.athinput',
    'unit_tests/cgl_lf_limiter_heat_flux_suppression.athinput',
    'unit_tests/cgl_lf_limiter_mirror.athinput',
    'unit_tests/cgl_lf_limiter_firehose.athinput',
    'unit_tests/cgl_lf_field_wave.athinput',
]


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    for input_file in _INPUTS:
        athena.run(input_file, [])


def analyze():
    logger.debug('Analyzing test ' + __name__)
    return True
