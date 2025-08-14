# Smoke test for MHD with AMR refinement
#
# Runs a short MHD linear wave case with AMR enabled to exercise
# refine + prolongation paths. Analysis checks that the run completes
# and that a history file is produced with at least one data row.

import os
import logging
import scripts.utils.athena as athena

logger = logging.getLogger('athena' + __name__[7:])


def run(**kwargs):
    logger.debug('Running test ' + __name__)
    # Keep runtime modest; trigger refinement via input thresholds
    args = [
        'job/basename=LinWave',
        'time/tlim=0.2',
        'time/nlim=200',
        'mesh/nx1=128',
        'mesh/nx2=64',
        'mesh/nx3=1',
        'meshblock/nx1=16',
        'meshblock/nx2=16',
        'meshblock/nx3=1',
        # speed up I/O; ensure history writes at least once
        'output1/dt=-1.0',
        'output2/dt=-1.0',
        'output3/dt=-1.0',
        'output4/dt=0.05',
    ]
    athena.run('tests/linear_wave_mhd_amr.athinput', args)


def analyze():
    logger.debug('Analyzing test ' + __name__)
    # History file is written by rank 0 under build/src
    hst_path = os.path.join('build', 'src', 'LinWave.mhd.hst')
    if not os.path.exists(hst_path):
        logger.warning('History file not found: {}'.format(hst_path))
        return False

    # Require at least one non-comment row of data
    nrows = 0
    with open(hst_path, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            nrows += 1

    if nrows == 0:
        logger.warning('History file has no data rows: {}'.format(hst_path))
        return False

    return True

